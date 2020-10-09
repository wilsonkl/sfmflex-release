module LinAlgUtils

using SparseArrays
using LinearAlgebra
using StaticArrays
using Rotations

export normalizecolumns
export invertSym3x3
export average_rotmat_l2
export evec_accuracy_report, check_linops, iseigenpairaccurate, eigensystemresiduals
export grassmanian_metric_oper, grassmanian_metric_angle, grassmanian_metric_chordal


# divide all columns of a sparse matrix by their l2 norm.
# Return the new matrix and the vector of column norms
function normalizecolumns(A)
    # get the column sums of A
    S = sqrt.(sum(abs2.(A), dims=1))
    S = reshape(S, (:))

    # get the nonzero entries in A.
    # ei is row index, ej is col index, ev is the value in A
    ei, ej, ev = findnz(A)

    # get the number or rows and columns in A
    m,n = size(A)

    # create a new normalized matrix. For each nonzero index (ei,ej), its new value will be
    # the old value divided by the sum of that column, which can be obtained by S[ej]
    A = sparse(ei,ej,ev./S[ej],m,n)
    return A, S
end


"""
Compute and return the inverse of a Symmetric 3x3 matrix
(this has been carefully profiled / tweaked for performance)
"""
function invertSym3x3(A::Array{Float64,2}) :: SMatrix{3,3, Float64}

    A = SMatrix{3,3}(A) # 6x speedup from this
    detinv = 1 / det(A)

    a = (A[2,2]*A[3,3] - A[2,3]*A[3,2]) * detinv
    b = (A[2,3]*A[3,1] - A[2,1]*A[3,3]) * detinv
    c = (A[2,1]*A[3,2] - A[2,2]*A[3,1]) * detinv

    #d = (A[1,3]*A[3,2] - A[1,2]*A[3,3]) * detinv   # it's symmetric
    e = (A[1,1]*A[3,3] - A[1,3]*A[3,1]) * detinv
    f = (A[1,2]*A[3,1] - A[1,1]*A[3,2]) * detinv

    #g = (A[1,2]*A[2,3] - A[1,3]*A[2,2]) * detinv
    #h = (A[1,3]*A[2,1] - A[1,1]*A[2,3]) * detinv
    i = (A[1,1]*A[2,2] - A[1,2]*A[2,1]) * detinv

    # building array at end: 10% speedup with StaticArrays, slowdown with
    # dynamic memory Arrays.
    return @SMatrix [a b c;
                     b e f;
                     c f i]
end

"""
Compute the average rotation matrix under the l2 geodesic distance. This
algorithm only is guaranteed to work in the rotations lie in a π/2 ball. (This
is generally the case because a π/2 ball covers all but a measure 0 subset of
SO(3).)
"""
function average_rotmat_l2(rots::Vector{RotMatrix{3, Float64}}) :: RotMatrix{3, Float64}
    # Add up all of the rotation matrices using an R^3x3 sum.
    Rsum = SMatrix{3,3}(rots[1])
    for k = 2:length(rots)
        Rsum = Rsum + rots[k]
    end

    # Do an SVD factorization
    fact = svd(Rsum, full=true)

    # form the average matrix
    Ravg = fact.U * fact.Vt

    # correct for cheirality
    if det(Ravg) < 0
        Ravg = fact.U * Diagonal([1.0, 1.0, -1.0]) * fact.Vt
    end

    return Ravg
end

function evec_accuracy_report(operator, evecs, evals; constraints=nothing)

    println("Eigen accuracy report.")
    println("    evecs: ", typeof(evecs), size(evecs))
    N = size(evecs, 2)
    ortho = norm(evecs'*evecs - Matrix{Float64}(I, N, N))
    println("    evecs ortho: |V'V-I|: ", ortho)

    if constraints !== nothing
        Nc = size(constraints, 2)
        println("    constraints: ", typeof(constraints), size(constraints))
        ortho = norm(constraints'*constraints - Matrix{Float64}(I, Nc, Nc))
        println("    constraints ortho: |C'C-I|: ", ortho)
    end

    for (k, (v, λ)) in enumerate(zip(eachcol(evecs), evals))
        v ./= norm(v)
        println("$k. --->")
        # If v is a large eval of covariance
        # it's a small eval of J'J.
        product = operator * v
        λ_actual = norm(product)

        product ./= λ_actual


        println("    stat1: |AX - λX| : ", norm(λ*v - λ_actual * product))
        println("    stat2: |AX| / λ  : ", λ_actual / λ )

        θ = rad2deg(acos(clamp(dot(v, product), -1, 1)))
        println("    stat3: θ deg     : ", θ)
        println("    stat4: dot       : ", dot(v, product))

        println("    supplied λ:    ", λ)
        println("    actual λ=|Av|: ", λ_actual)
        println("    actual λ=v'Av: ", v'*operator*v)

        if constraints !== nothing
            println("    nullspace mass |v'C|: $(norm(v'*constraints))")
        end
    end
end

"""
Perform just a small part of the evec_accuracy_report computation: check the
accuracy of a single eigenvector under the |Av-λv| score. Return True if
the accuracy is satisfactory.
"""
function iseigenpairaccurate(operator, v, λ; ϵ=1e-6) :: Bool
    v ./= norm(v)
    return norm(λ*v - operator*v) < ϵ
end

"""
Return the residuals of each given eigenpair under the |Av-λv| score.
Normalize by a reference scale (Fro norm / num elements) to get a relative measure.
"""
function eigensystemresiduals(operator, evecs, evals) :: Vector{Float64}
    residuals = []
    for (k, (v, λ)) in enumerate(zip(eachcol(evecs), evals))
        v ./= norm(v)
        push!(residuals,  norm(λ*v - operator*v))
    end
    reference = norm(operator, 2) / prod(size(operator))
    return residuals / reference
end

"""
Check if two linear operators are the same.

Allow for them to be in different forms (i.e., dense matrix, SparseMatrixCSC,
    LinearMap, ...). Just assume that they support size() and mult!()

Rather than going element by element (which would take forever, given that we
    don't even have element-wise access) do a randomized search. Generate
    random vectors, apply the operator, and take the norm in the range-space.

Returns the maximum norm-difference across all trials.

Throws if operator dimensions are not identical.
"""
function check_linops(operatorA, operatorB; numtrials=1000)
    # I guess size() doesn't exist for some LinearMaps?
    N = size(operatorA, 1)
    M = size(operatorA, 2)
    # N == size(operatorB, 1) || throw("The operator dimensions are not the same.")
    # M == size(operatorB, 2) || throw("The operator dimensions are not the same.")


    diffs = SharedArray{Float64}(numtrials)
    Threads.@threads for i = 1:numtrials
        x = rand(M)
        diffs[i] = norm(operatorA*x - operatorB*x)
    end
    return maximum(diffs)
end

"""
Measure the distance between two linear subspaces. V and W are both n-by-k
matrices representing elements of the Grassmanian Gr(k, mathbb{R}^n).

Orthogonality is *not* assumed: V represents the subspaces spanned by the
columns of V, and likewise for W.

Return the metric d(V, W) = || P_V - P_W ||
    where P_V is the canonical orthogonal projection onto the column space of V
    where P_W is the canonical orthogonal projection onto the column space of W
    where || ⋅ || is the operator norm
Because the orthogonal projectors are real symmetric, the operator norm is the
largest eigenvalue of the operator.
"""
function grassmanian_metric_oper(V, W)
    size(V) == size(W) || throw("Dimensions don't match!")

    # orthnormalize V and W using a skinny QR
    function orthonormalize(A)
        fact = qr(A)
        return Matrix(fact.Q)
    end
    QV = LinearMap(orthonormalize(V))
    QW = LinearMap(orthonormalize(W))

    # create a matrix-free linear operator for the projection operation
    # PV - PW = (QV*QV') - (QW*QW')
    D = QV * transpose(QV) - QW * transpose(QW)

    # run a matrix-free eigensolver to find the largest eigenvalue
    λ, ϕ, nconv, niter, nmult, resid = eigs(D, nev=1, which=:LM, maxiter=100)
    println("λ: ", λ)
    println("nconv: ", nconv)
    println("niter: ", niter)
    println("nmult: ", nmult)

    return abs(λ[1])
end

"""
Measure the distance between two linear subspaces. V and W are both n-by-k
matrices representing elements of the Grassmanian Gr(k, mathbb{R}^n).

Orthogonality is *not* assumed: V represents the subspace spanned by the
columns of V, and likewise for W.

Return the cosine of the largest principle angle:
        cos θ = max(v in V, w in W, |v|=1, |w|=1)  v'w

Procedure: orthonormalize V and W, then find largest singular value of V'W

Note: this version is a lot faster than the operator norm. I'm seeing about .4ms
    for a 2000x5 problem.
"""
function grassmanian_metric_angle(V, W)
    size(V) == size(W) || throw("Dimensions don't match!")

    # orthnormalize V and W using a skinny QR
    function orthonormalize(A)
        fact = qr(A)
        return Matrix(fact.Q)
    end
    QV = orthonormalize(V)
    QW = orthonormalize(W)

    # assume dims are small: find all singular values with a dense solver.
    # (if #columns is not small, it would be better to use an iterative solver.)
    σ = svdvals(QV' * QW) # returned in descending order
    return σ[1]
end

"""
Measure the distance between two linear subspaces. V and W are both n-by-k
matrices representing elements of the Grassmanian Gr(k, mathbb{R}^n).

Orthogonality is *not* assumed: V represents the subspace spanned by the
columns of V, and likewise for W.

Return the sum of the sins of the principal angles.

Procedure: orthonormalize V and W, then find the singular values of V'W

Citation: Huang, Qui, Calderbank: The Role of Principal Angles in
Subspace Classification
"""
function grassmanian_metric_chordal(V, W)
    size(V) == size(W) || throw("Dimensions don't match!")

    # orthnormalize V and W using a skinny QR
    function orthonormalize(A)
        fact = qr(A)
        return Matrix(fact.Q)
    end
    QV = orthonormalize(V)
    QW = orthonormalize(W)

    # assume dims are small: find all singular values with a dense solver.
    # (if #columns is not small, it would be better to use an iterative solver.)
    σ = svdvals(QV' * QW) # returned in descending order
    # σ's are cos(θ). We want to return ∑sin^2(θ)
    return sum(1 .- σ.^2)
end

end # module
