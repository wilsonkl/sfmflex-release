using Rotations
using StaticArrays
using Printf

function readBundleFile(bundle_file)
    fin = open(bundle_file)

    # read the header
    header = readline(fin)
    nCams, nPts = parse.(Int, split(readline(fin)))
    println("found $nCams cameras and $nPts points")

    # read the camera blocks
    camera_params = []
    for i=1:nCams
        f,   k1,  k2 = parse.(Float64, split(readline(fin)))
        r11, r12, r13 = parse.(Float64, split(readline(fin)))
        r21, r22, r23 = parse.(Float64, split(readline(fin)))
        r31, r32, r33 = parse.(Float64, split(readline(fin)))
        t1,  t2,  t3 = parse.(Float64, split(readline(fin)))
        R = @SMatrix [r11 r12 r13; r21 r22 r23; r31 r32 r33]
        aa = AngleAxis(RotMatrix{3}(R))
        v = rotation_angle(aa) .* rotation_axis(aa)
        push!(camera_params, [v[1], v[2], v[3], t1, t2, t3, f, k1, k2])
    end

    # read the point blocks
    observations = []
    point_params = []
    for i=1:nPts
        X1, X2, X3 = parse.(Float64, split(readline(fin))) # position
        readline(fin) # color (ignore)
        tokens = split(readline(fin))

        num_observations = parse(Int, tokens[1])
        for k=1:num_observations
            cam_id = parse(Int, tokens[(k-1)*4 + 2])
            key_id = parse(Int, tokens[(k-1)*4 + 3])
            xcoord = parse(Float64, tokens[(k-1)*4 + 4])
            ycoord = parse(Float64, tokens[(k-1)*4 + 5])

            push!(observations, [cam_id, i-1, xcoord, ycoord])
        end
        # Note to self: 0,0 is the center of the image.
        # (-w/2, -h/2) is the lower-left corner.
        push!(point_params, [X1, X2, X3])
    end

    close(fin)
    return camera_params, point_params, observations
end

function writeBALFile(bal_file, camera_params, point_params, observations)
    fout = open(bal_file, "w")
    @printf(fout, "%d %d %d\n", length(camera_params), length(point_params), length(observations))
    for obs in observations
        @printf(fout, "%d %d %lf %lf\n", obs...)
    end
    for cam in camera_params
        @printf(fout, "%lf %lf %lf %lf %lf %lf %lf %lf %lf\n", cam...)
    end
    for pt in point_params
        @printf(fout, "%lf %lf %lf\n", pt...)
    end
    close(fout)
end

function reindexProblem(camera_params, point_params, observations)

    changed_something = true
    while changed_something
        println("Start of pruning round.")
        changed_something = false

        # look for cameras that see fewer than 5 points
        # count how many obs per camera
        camera_counts = zeros(length(camera_params))
        for obs in observations
            cam_index = convert(Int, obs[1])
            camera_counts[cam_index+1] += 1
        end
        mask = camera_counts .< 5
        if any(mask)
            changed_something = true

            # compute the reindexing scheme
            old2new = Dict()
            count = 0
            for i=1:length(camera_params)
                if camera_counts[i] >= 5
                    old2new[i-1] = count
                    count += 1
                end
            end

            # apply the reindexing:
            new_observations = []
            for obs in observations
                cam_index = convert(Int, obs[1])
                if camera_counts[cam_index+1] >= 5
                    push!(new_observations, [old2new[cam_index], obs[2], obs[3], obs[4]])
                end
            end
            println("got rid of $(length(observations) - length(new_observations)) cameras.")
            observations = new_observations
            camera_params = camera_params[.!mask]
        end


        # look for points seen by fewer than two cameras
        # count how many obs per point
        point_counts = zeros(length(point_params))
        for obs in observations
            pt_index = convert(Int, obs[2])
            point_counts[pt_index+1] += 1
        end
        mask = point_counts .< 2
        if any(mask)
            changed_something = true

            # compute the reindexing scheme
            old2new = Dict()
            count = 0
            for i=1:length(point_params)
                if point_counts[i] >= 2
                    old2new[i-1] = count
                    count += 1
                end
            end

            # apply the reindexing:
            new_observations = []
            for obs in observations
                pt_index = convert(Int, obs[2])
                if point_counts[pt_index+1] >= 2
                    push!(new_observations, [obs[1], old2new[pt_index], obs[3], obs[4]])
                end
            end
            println("got rid of $(length(observations) - length(new_observations)) points.")
            observations = new_observations
            point_params = point_params[.!mask]
        end

    end
    println("Done pruning.")
    return camera_params, point_params, observations
end


bundle_file = ARGS[1]
println("reading in bundle file: $bundle_file")
bal_file = ARGS[2]
println("writing to bal file: $bal_file")
camera_params, point_params, observations = readBundleFile(bundle_file)
camera_params, point_params, observations = reindexProblem(camera_params, point_params, observations)
writeBALFile(bal_file, camera_params, point_params, observations)
