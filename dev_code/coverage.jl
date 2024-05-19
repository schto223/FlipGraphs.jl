using Coverage
using Pkg


clean_folder("local_code")
clean_folder("src")
clean_folder("test")

Pkg.test("FlipGraphs", coverage=true)

# process '*.cov' files
coverage = process_folder() # defaults to src/; alternatively, supply the folder name as argument
#coverage = append!(coverage, process_folder("deps"))  # useful if you want to analyze more than just src/
# process '*.info' files, if you collected them
coverage = merge_coverage_counts(coverage, filter!(
    let prefixes = (joinpath(pwd(), "src", "")#,
                    )#joinpath(pwd(), "deps", ""))
        c -> any(p -> startswith(c.filename, p), prefixes)
    end,
    LCOV.readfolder("test")))
# Get total coverage for all Julia files
covered_lines, total_lines = get_summary(coverage)
# Or process a single file
@show get_summary(process_file(joinpath("src", "FlipGraphs.jl")))

LCOV.writefile("local_code/lcov.info", coverage)

clean_folder("src")
clean_folder("test")