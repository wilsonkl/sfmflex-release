#=
Run the method on all 1dsfm datasets! Save a .json file for viz, and also save the
eigenvalues, eigenvectors, and timing data for future analysis.
=#

include("../set_loadpath.jl")
using LinAlgUtils
using JacobianUtils

# JLD requires that packages whose types are being serialized be included, even
# if they aren't explicitly used otherwise in the file.
using SparseArrays
using LinearAlgebra
using Rotations

import Glob
using JSON
using JLD

########
# MAIN #
########
num_vecs = 20


dataset_root = "dataset/1dsfm"

# output location
outdir_root = "output/1dsfm"

# bin file location
convert_to_bal = normpath(joinpath(@__DIR__, "bundle2bal.jl"))
compute_jacobian = normpath(joinpath(@__DIR__, "..", "build", "bin", "compute_jacobian"))


problems = Glob.glob("$dataset_root/*")
for problem_file in problems

    dataset = splitext(basename(problem_file))[1]
    # find all of the individual problems in this dataset
    # only get post-bundle adjustment problems.
    outdir = "$outdir_root/$dataset"

    # check that the outdir exists. If not, make it.
    if !isdir(outdir)
        mkpath(outdir)
    end
    results_file = "$outdir/results.jld"
    if ~isfile(results_file)
        save(results_file, "all_runs", [])
    end

    println("Beginning computation on dataset: $dataset")



    problem_shortname = dataset
    bal_file = "$outdir/$problem_shortname.txt"
    json_file = "$outdir/$problem_shortname.json"
    jacobian_file = "$outdir/$problem_shortname.jac.h5"

    # check: does the bal file already exist? If so, skip this problem.
    if isfile(bal_file)
        println("Saved results found for $dataset/$problem_shortname. Skipping.")
        continue
    end

    # convert the bundle problem to a bal problem.
    println("Converting .out to BAL for problem $dataset/$problem_shortname")
    t_conv = @elapsed run(`julia $convert_to_bal "$problem_file" "$bal_file"`)
    println("  ---> that took $t_conv seconds.")

    # extract the jacobian from a problem instance
    println("Computing jacobian for problem $dataset/$problem_shortname")
    t_jac = @elapsed run(`$compute_jacobian "$bal_file" "$jacobian_file"`)
    println("  ---> that took $t_jac seconds.")

    # get the eigenvectors from this jacobian
    println("Computing eigenvectors for problem $dataset/$problem_shortname")
    problem = readproblem(jacobian_file)
    t_eig = @elapsed λ, X, residuals = compute_covariance_evecs(problem, num_vecs, method=:dropping)
    rm(jacobian_file)
    println("  ---> that took $t_eig seconds.")

    println("Saving output for problem $dataset/$problem_shortname")

    # TEMP: camera center calculations so we don't have to figure out how to do
    # this in Javascript
    camera_centers = zeros(Float64, 3*problem.num_cameras)
    for k = 0:problem.num_cameras-1
        rot = problem.camera_parameters[9*k+1:9*k+3]
        θ = norm(rot)
        rot_aa = AngleAxis(θ, rot[1], rot[2], rot[3])
        t = problem.camera_parameters[9*k+4:9*k+6]
        camera_centers[3*k+1 : 3*k+3] = -inv(rot_aa) * t
    end

    # save the vis information to a json file
    # stuff this into a json file
    json_data = Dict()
    json_data["num_cameras"] = problem.num_cameras
    json_data["num_points"] = problem.num_points
    json_data["camera_parameters"] = problem.camera_parameters
    json_data["point_parameters"] = problem.point_parameters
    json_data["camera_centers"] = camera_centers
    evecs = []
    for i=1:size(X,2)
        push!(evecs, X[:,i])
    end
    json_data["eigenvectors"] = evecs
    json_data["eigenvalues"] = λ
    open(json_file, "w") do f
        write(f, "sample_scene = ")
        JSON.print(f, json_data)
    end

    # prepare a structure with all of the results data
    data = Dict(
        "dataset_name"         => dataset,
        "problem_name"         => problem_shortname,
        "num_cameras"          => problem.num_cameras,
        "num_points"           => problem.num_points,
        "num_observations"     => problem.num_observations,
        "eigenvectors"         => evecs,
        "eigenvalues"          => λ,
        "computation_time"     => t_eig,
        "residuals"            => residuals,
        "source_problem_file"  => problem_file
    )
    # read the records file, append the new data, and push back to file
    all_records = load(results_file, "all_runs")
    rm(results_file)
    push!(all_records, data)
    save(results_file, "all_runs", all_records)

    println("")

end
