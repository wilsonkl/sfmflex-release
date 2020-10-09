include("../set_loadpath.jl")
using LinAlgUtils
using JacobianUtils
using JSON
using Rotations
using LinearAlgebra

###########
# OPTIONS #
###########
# bal_problem = "ladybug/problem-372-47423-post.txt"
# output_json = "vis/sample_scene.json"


bal_problem = "venice/problem-1490-935273-post.txt"
output_json = "vis/big_sample_scene.json"

# bal_problem = "trafalgar/problem-257-65132-post.txt"
# bal_problem = "dubrovnik/problem-253-163691-post.txt"
# output_json = "vis/other_sample_scene.json"
num_vecs = 20

##########
# SCRIPT #
##########

# extract the jacobian from a problem instance
problem_file = normpath(joinpath(@__DIR__, "..", "dataset", "bal", bal_problem))
compute_jacobian = normpath(joinpath(@__DIR__, "..", "build", "bin", "compute_jacobian"))
jacobian_file = tempname() * ".h5" # this is string concat in julia??? !!!
cmd = `$compute_jacobian "$problem_file" "$jacobian_file"`
run(cmd)

# get the eigenvectors from this jacobian
problem = readproblem(jacobian_file)
rm(jacobian_file)
λ, X = compute_covariance_evecs(problem, num_vecs, method=:dropping)

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

# write the json structure to a file
open(output_json, "w") do f
    write(f, "sample_scene = ")
    JSON.print(f, json_data)
end
