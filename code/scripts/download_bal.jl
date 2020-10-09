import HTTP
import Gumbo
import Cascadia

#=
Download the BAL datasets from https://grail.cs.washington.edu/projects/bal/.
Unzip all of the files.
=#

function scrape_downloads_for_category(category)
    downloads = String[]

    url = "https://grail.cs.washington.edu/projects/bal/$category.html"
    response = HTTP.get(url)

    html = Gumbo.parsehtml(String(response.body))
    links = Cascadia.eachmatch(
        Cascadia.Selector("a"),
        html.root
    )
    for link in links
        href = Gumbo.getattr(link, "href")

        # exclude the "return to index" link
        if occursin(category, href)
            push!(downloads, href)
        end
    end
    return downloads
end

function download_bal_problem(web_link, category)

    short_name = basename(web_link)

    if startswith(web_link, "http")
        long_name = web_link
    else
        long_name = "https://grail.cs.washington.edu/projects/bal/$web_link"
    end

    local_name = joinpath(joinpath(dataset_dir, category), short_name)
    local_unzipped_name = replace(local_name, ".txt.bz2" => ".txt")

    if !isfile(local_unzipped_name)

        if !isfile(local_name)
            println("Downloading ", category, "/", short_name)
            download(long_name, local_name)
        end

        println("Decompressing ", category, "/", short_name)
        run(`bzip2 -d $local_name`)
    end
end

function main()
    # put the BAL files here
    dataset_dir = normpath(joinpath(@__DIR__, "..", "dataset", "bal"))
    println("Downloading BAL dataset to ", dataset_dir)

    # make sure the dataset directory exists
    mkpath(dataset_dir)

    bal_categories = [
        "ladybug",
        "trafalgar",
        "dubrovnik",
        "venice",
        "final"
    ]
    for category in bal_categories

        category_dir = joinpath(dataset_dir, category)
        if !isdir(category_dir)
            mkdir(category_dir)
        end

        web_links = scrape_downloads_for_category(category)
        for link in web_links
            download_bal_problem(link, category)
        end
    end
end
