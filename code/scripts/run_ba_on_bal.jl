import Glob

dataset_dir = normpath(joinpath(@__DIR__, "..", "dataset", "bal"))
problems = Glob.glob("*/problem*-pre.txt", dataset_dir)

bundle_adjuster = normpath(joinpath(@__DIR__, "..", "build", "bin", "bundle_adjuster"))

for pre_name in problems
    post_name = replace(pre_name, "-pre.txt" => "-post.txt")
    log_name = replace(pre_name, "-pre.txt" => "-post.log")
    if !isfile(post_name)
        println("Bundle adjusting: ", joinpath(splitpath(pre_name)[end-1:end]...))
        cmd = `$bundle_adjuster --input="$pre_name" --final_txt="$post_name" --num_iterations=200 --num_threads=4 "2>&1"`
        write(log_name, read(cmd, String))
    end
end
