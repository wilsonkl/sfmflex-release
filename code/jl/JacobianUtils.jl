module JacobianUtils

using SparseArrays
using LinearAlgebra
using Statistics
using Rotations
using StaticArrays
using HDF5
using Arpack
using DelimitedFiles

using LinAlgUtils

export BAProblemInstance
export readproblem, splitjacobian, estimate_parameter_scales
export estimate_camera_center_scale
export compute_measurement_covariance
export form_nullspace, build_schur_complement_matrix
export compute_covariance_evecs
export compute_evecs_deflation, compute_evecs_projection, compute_evecs_dropping


mutable struct BAProblemInstance
    num_cameras::Int
    num_points::Int
    num_observations::Int
    camera_parameters::Array{Float64}
    point_parameters::Array{Float64}
    residuals::Array{Float64}
    jacobian::SparseMatrixCSC{Float64,Int}
    scale_matrix::Diagonal{Float64}
end

"""
    readproblem(filename::String)

Unpack a HDF5 file containing a bundle adjustment problem.
"""
function readproblem(filename::String) :: BAProblemInstance
    num_cameras::Int                  = h5read(filename, "num_cameras")
    num_points::Int                   = h5read(filename, "num_points")
    num_observations::Int             = h5read(filename, "num_observations")

    camera_parameters::Array{Float64} = h5read(filename, "camera_parameters")
    point_parameters::Array{Float64}  = h5read(filename, "point_parameters")
    residuals::Array{Float64}         = h5read(filename, "residuals")

    rows::Array{Int}                  = h5read(filename, "jacobian_sparse_row")
    cols::Array{Int}                  = h5read(filename, "jacobian_sparse_column")
    vals::Array{Float64}              = h5read(filename, "jacobian_sparse_value")
    jacobian = sparse(rows .+ 1, cols .+ 1, vals)

    scale_matrix = Diagonal(ones(Float64, 6*num_cameras))

    return BAProblemInstance(
        num_cameras, num_points, num_observations,
        camera_parameters, point_parameters,
        residuals, jacobian,
        scale_matrix
    )
end

"""
Estimate the size of a SfM problem. Compute the distances from camera centers
to the mean camera location, and return the median of those distances.
"""
function estimate_camera_center_scale(problem::BAProblemInstance) :: Float64
    camera_centers = zeros(Float64, problem.num_cameras, 3)
    for k = 0:problem.num_cameras-1
        rot = problem.camera_parameters[9*k+1:9*k+3]
        θ = norm(rot)
        rot_aa = AngleAxis(θ, rot[1], rot[2], rot[3])
        t = problem.camera_parameters[9*k+4:9*k+6]
        camera_centers[k+1,:] = -rot_aa * t
    end

    average_cam_center = mean(camera_centers, dims=1)
    return median(sqrt.(sum((camera_centers .- average_cam_center).^2, dims=2)))
end

"""
Estimate the size of a SfM problem. Compute mean camera orientation, and
compute distances from that. Return the median orientation distance
"""
function estimate_camera_rot_scale(problem::BAProblemInstance) :: Float64
    # unpack the camera rotations from the problem params
    camera_rots = Vector{RotMatrix{3, Float64}}(undef, problem.num_cameras)
    for k = 0:problem.num_cameras-1
        rot = problem.camera_parameters[9*k+1:9*k+3]
        θ = norm(rot)
        camera_rots[k+1] = RotMatrix(AngleAxis(θ, rot[1], rot[2], rot[3]))
    end

    # compute the mean orientation
    Ravg = average_rotmat_l2(camera_rots)

    # compute distances to the mean rotation
    distances = Vector{Float64}(undef, problem.num_cameras)
    for k = 1:problem.num_cameras
        distances[k] = rotation_angle(Ravg' * camera_rots[k])
    end
    return median(distances)
end

"""
Split the jacobian into two. Jc are the columns representing camera params
and Jp are the columns representing point params. Discard all columns for
camera intrinsics.
"""
function splitjacobian(problem::BAProblemInstance) :: Tuple{SparseMatrixCSC{Float64,Int},SparseMatrixCSC{Float64,Int}}
    mask_cameras = falses(problem.jacobian.n)
    mask_cameras[1:9*problem.num_cameras] .= true

    mask_intrinsics = falses(problem.jacobian.n)
    for i in range(7, stop=9*problem.num_cameras, step=9)
        mask_intrinsics[i]   = true
        mask_intrinsics[i+1] = true
        mask_intrinsics[i+2] = true
    end

    Jc = problem.jacobian[:, mask_cameras .& .~mask_intrinsics]
    Jp = problem.jacobian[:, .~mask_cameras]
    return Jc, Jp
end

"""
Estimate the relative scales of rotations and translations, and save it as
a Diagonal matrix.
"""
function estimate_parameter_scales(problem::BAProblemInstance) :: Diagonal{Float64}
    scale_r = estimate_camera_rot_scale(problem)
    scale_x = estimate_camera_center_scale(problem)
    scale_factors = Vector{Float64}(undef, 6*problem.num_cameras)
    for i = 1:problem.num_cameras
        scale_factors[6*i-5:6*i-3] .= scale_r
        scale_factors[6*i-2:6*i-0] .= scale_x
    end
    camera_scale_matrix = Diagonal(scale_factors)
    problem.scale_matrix = camera_scale_matrix
    return camera_scale_matrix
end

"""
Estimate the covariance in the residual noise from the values of the
residuals. (This is very coarse.) Return a sparse diagonal matrix.

See Michal Polic's thesis (pg 10) for an overview of these methods.
method=
  :bishop: use diag(abs(residuals)) as an estimate of measurement cov.
  :lhuillier: Use a dof-adjusted norm of the residuals. Assume iid.
"""
function compute_measurement_covariance(
        problem::BAProblemInstance;
        method=:bishop) :: Diagonal{Float64}

    if method == :bishop
        return Diagonal(abs.(problem.residuals))

    elseif method == :lhuillier
        num_camera_params = 9
        gauge_rank = 7
        k = problem.num_observations
        n = problem.num_cameras
        σ = norm(problem.residuals) / (2*k - num_camera_params * n - gauge_rank)

        return Diagonal(σ .* ones(2*k))
        # return σ .* sparse(I, 2*k, 2*k)
    else
        throw(ArgumentError("method=$method is not supported"))
    end
end

"""
Create a num_params x 7 matrix that represents the kernel of a similarity
transformation. Vectors orthogonal to the columns of this matrix are
normal to the gauge subspace.
"""
function form_nullspace(problem::BAProblemInstance)
    num_params = 6 * problem.num_cameras
    nullspace = Array{Float64, 2}(undef, 6*problem.num_cameras, 7)
    for k in range(0, stop=problem.num_cameras-1)
        nullspace[6*k+1:6*k+6, 1:6] = Matrix{Float64}(I, 6, 6)
        # old code: this is a bug!
        # nullspace[6*k+1:6*k+3, 7] = problem.camera_parameters[9*k+1:9*k+3]
        # nullspace[6*k+4:6*k+6, 7] = zeros(Float64, 3, 1)
        rot = problem.camera_parameters[9*k+1:9*k+3]
        θ = norm(rot)
        rot_aa = AngleAxis(θ, rot[1], rot[2], rot[3])
        t = problem.camera_parameters[9*k+4:9*k+6]

        nullspace[6*k+1:6*k+3, 7] = zeros(Float64, 3, 1)
        nullspace[6*k+4:6*k+6, 7] = -rot_aa * t
    end

    # orthonormalize it
    fact = qr(problem.scale_matrix * nullspace)
    return Matrix(fact.Q)
end

"""
The Schur complement matrix is defined as D - (B' * A^-1 * B) where
    [A  B ] = J' Σ J
    [B' D ]
and the Jacobian J is partitioned: J = [Jp, Jc].

This function explicitly builds the Schur complement matrix as a dense matrix.
"""
function build_schur_complement_matrix(
        Jc :: SparseMatrixCSC{Float64,Int64},
        Jp :: SparseMatrixCSC{Float64,Int64},
        Sc :: Diagonal{Float64},
        Σ_input :: Diagonal{Float64}) :: Matrix{Float64}

    Σ_input_inv = inv(Σ_input)

    N  :: Int64 = size(Jc, 1)    # num rows
    Mp :: Int64 = size(Jp, 2)  # num point cols
    Mc :: Int64 = size(Jc, 2)  # num cam cols

    size(Jp,          1) == N || throw(DimensionMismatch())
    size(Σ_input_inv, 1) == N || throw(DimensionMismatch())
    size(Σ_input_inv, 2) == N || throw(DimensionMismatch())

    JpT = transpose(Jp)
    JcT = transpose(Jc)
    A = JpT * Σ_input_inv * Jp             # Mp x Mp
    B = JpT * Σ_input_inv * Jc             # Mp x Mc
    D = JcT * Σ_input_inv * Jc             # Mc x Mc

    # The A matrix is block-diagonal with 3x3 blocks
    # Create Ainv using a specialized 3x3 inversion function.
    # No nice way to reform a sparse matrix...
    num_pts = convert(Int64, Mp / 3)
    row = Vector{Int64}(undef, 9 * num_pts)
    col = Vector{Int64}(undef, 9 * num_pts)
    val = Vector{Float64}(undef, 9 * num_pts)
    for i = 1:num_pts
        # get the block and invert it
        Aidx = 3 * (i-1) + 1
        block = invertSym3x3(Array(A[Aidx:Aidx+2, Aidx:Aidx+2]))

        # stuff into the sparse storage (will this compile better?)
        Lidx = 9 * (i-1) + 1
        row[Lidx + 0] = Aidx + 0
        row[Lidx + 1] = Aidx + 0
        row[Lidx + 2] = Aidx + 0
        row[Lidx + 3] = Aidx + 1
        row[Lidx + 4] = Aidx + 1
        row[Lidx + 5] = Aidx + 1
        row[Lidx + 6] = Aidx + 2
        row[Lidx + 7] = Aidx + 2
        row[Lidx + 8] = Aidx + 2

        col[Lidx + 0] = Aidx + 0
        col[Lidx + 1] = Aidx + 1
        col[Lidx + 2] = Aidx + 2
        col[Lidx + 3] = Aidx + 0
        col[Lidx + 4] = Aidx + 1
        col[Lidx + 5] = Aidx + 2
        col[Lidx + 6] = Aidx + 0
        col[Lidx + 7] = Aidx + 1
        col[Lidx + 8] = Aidx + 2

        val[Lidx + 0] = block[1,1]
        val[Lidx + 1] = block[1,2]
        val[Lidx + 2] = block[1,3]
        val[Lidx + 3] = block[2,1]
        val[Lidx + 4] = block[2,2]
        val[Lidx + 5] = block[2,3]
        val[Lidx + 6] = block[3,1]
        val[Lidx + 7] = block[3,2]
        val[Lidx + 8] = block[3,3]
    end
    Ainv = sparse(row, col, val, Mp, Mp)

    # Compose the final operator: D - B'A^{-1}B
    Z = D - B' * Ainv * B

    return Matrix(Sc * Z * Sc)
end

"""
    compute_covariance_evecs

Compute the eigenvectors corresponding to the largest eigenvalues of the
covariance matrix.

Params:
  problem::BAProblemInstance - the problem instance
  num_vecs::Int - the number of eigenvectors to compute
Keywords:
  method:
    :projection - project the gauge space to null. Run ARPACK and discard
        the first few vectors.
    :deflation (default) - use Golub deflation to collapse the gauge space.
        Run ARPACK and then re-inflate the results.
    :dropping - just drop the lowest 7 eigenvectors
  debug::bool - display debugging output or not
"""
function compute_covariance_evecs(
    problem::BAProblemInstance,
    num_vecs::Int;
    method=:deflation::Symbol,
    debug=false::Bool
    )

    if method == :projection
        return compute_evecs_projection(problem, num_vecs, debug=debug)
    elseif method == :deflation
        return compute_evecs_deflation(problem, num_vecs, debug=debug)
    else
        return compute_evecs_dropping(problem, num_vecs, debug=debug)
    end
end

function compute_evecs_projection(
    problem::BAProblemInstance,
    num_vecs::Int;
    debug=false::Bool
    ) :: Tuple{Vector{Float64}, Matrix{Float64}, Vector{Float64}}

    Jc, Jp = splitjacobian(problem)
    Sc = estimate_parameter_scales(problem)
    Σ = compute_measurement_covariance(problem, method=:bishop)
    Z = build_schur_complement_matrix(Jc, Jp, Sc, Σ)
    N = form_nullspace(problem)

    # ensure symmetry
    Z = (Z + Z') / 2

    # project Z off of N
    P = Matrix{Float64}(I, size(Z, 1), size(Z, 1)) - N * N'
    PZP = Hermitian(P * Z * P)

    # now run ARPACK to get the eigensystem
    λ, V, nconv, niter, nmult, resid = eigs(PZP, nev=num_vecs+size(N,2), which=:SM, ritzvec=true, maxiter=1000)
    i = sortperm(λ)
    λ = λ[i]
    V = V[:,i]

    if debug
        println("Report: arpack + projection")
        println("  null eigs: $(λ[1:size(N, 2)])")
        println("  real eigs: $(λ[size(N,2)+1:end])")
        println("  nconv: $nconv")
        println("  niter: $niter")
        println("  nmult: $nmult")
        evec_accuracy_report(PZP, V, λ, constraints=N)
    end
    residuals = eigensystemresiduals(P*Z*P, V[:,size(N,2)+1:end], λ[size(N,2)+1:end])

    # Reverse the parameter scaling transformation
    V = problem.scale_matrix * V[:,size(N,2)+1:end]
    V = V ./ sqrt.(sum(abs2.(V), dims=1))

    return λ[size(N,2)+1:end].^-1, V, residuals
end

function compute_evecs_dropping(
    problem::BAProblemInstance,
    num_vecs::Int;
    debug=false::Bool
    ) :: Tuple{Vector{Float64}, Matrix{Float64}, Vector{Float64}}

    Jc, Jp = splitjacobian(problem)
    Sc = estimate_parameter_scales(problem)
    Σ = compute_measurement_covariance(problem, method=:bishop)
    Z = build_schur_complement_matrix(Jc, Jp, Sc, Σ)
    N = form_nullspace(problem) # only really using this for diagnostics

    # ensure symmetry
    Z = (Z + Z') / 2

    # now run ARPACK to get the eigensystem
    λ, V, nconv, niter, nmult, resid = eigs(Z, nev=num_vecs+size(N,2), which=:SM, ritzvec=true, maxiter=1000)
    i = sortperm(λ)
    λ = λ[i]
    V = V[:,i]

    if debug
        println("Report: arpack + projection")
        println("  null eigs: $(λ[1:size(N, 2)])")
        println("  real eigs: $(λ[size(N,2)+1:end])")
        println("  nconv: $nconv")
        println("  niter: $niter")
        println("  nmult: $nmult")
        evec_accuracy_report(Z, V, λ, constraints=N)
        err = grassmanian_metric_chordal(V[:,1:7], N)
        println("null space grassmanian chordal err:", err)
    end
    residuals = eigensystemresiduals(Z, V[:,size(N,2)+1:end], λ[size(N,2)+1:end])

    # Reverse the parameter scaling transformation
    V = problem.scale_matrix * V[:,size(N,2)+1:end]
    V = V ./ sqrt.(sum(abs2.(V), dims=1))

    return λ[size(N,2)+1:end].^-1, V, residuals
end

"""
Use ARPACK to solve a constrained eigenvalue problem, this time using Golub
deflation to do hard constraint enforcement.
"""
function compute_evecs_deflation(
    problem::BAProblemInstance,
    num_vecs::Int;
    debug=false::Bool
    ) :: Tuple{Vector{Float64}, Matrix{Float64}, Vector{Float64}}

    Jc, Jp = splitjacobian(problem)
    Sc = estimate_parameter_scales(problem)
    Σ = compute_measurement_covariance(problem, method=:bishop)
    Z = build_schur_complement_matrix(Jc, Jp, Sc, Σ)
    N = form_nullspace(problem)

    # ensure symmetry
    Z = (Z + Z') / 2

    # deflate the schur complement wrt the constraint space
    m = size(N, 1)
    p = size(N, 2)
    fact = qr(N)
    Q = fact.Q * Matrix(I, m, m) # get the full Q from the factorization
    Q2 = Q[:,p+1:end]
    G = Symmetric(Q2' * Z * Q2)

    # now run ARPACK to get the eigensystem
    λ, V, nconv, niter, nmult, resid = eigs(G, nev=num_vecs, which=:SM, ritzvec=true, maxiter=1000)

    # reinflate the answer
    V = Q2 * V

    P = Matrix(I, m, m) - N * N'
    if debug
        println("Report: arpack + golub deflation")
        println("  eigs: $λ")
        println("  nconv: $nconv")
        println("  niter: $niter")
        println("  nmult: $nmult")
        evec_accuracy_report(P*Z*P, V, λ, constraints=N)
    end
    residuals = eigensystemresiduals(P*Z*P, V, λ)

    # check convergences?
    e = eps(real(eltype(G)))/2
    f = cholesky(G)
    v = V[:,1]
    x = f / v
    println("|Ainv*v - λinv*v|:" , norm(x - v./λ[1]))

    # Reverse the parameter scaling transformation
    V = problem.scale_matrix * V
    V = V ./ sqrt.(sum(abs2.(V), dims=1))

    return λ.^-1, V, residuals
end



end # module
