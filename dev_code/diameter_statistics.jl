using FlipGraphs
import XLSX
using StatsPlots
gr(size=(1000,800))
using DataFrames
import RDatasets
iris = RDatasets.dataset("datasets", "iris")
using Interact
using Blink

function collect_diameters(genus::Integer, points::Integer, n_diams::Integer, n_flips::Integer)
    return collect_diameters(deltacomplex(genus, points),n_diams,n_flips)
end
    
function collect_diameters(D::DeltaComplex, n_diams::Integer, n_flips::Integer)
    diams = Vector{Int}(undef, n_diams)
    for i in 1:n_diams
        random_flips!(D, n_flips)
        diams[i] = diameter(D)
    end
    return diams
end
    
function export_excel(genus::Integer, points::Integer, n_diams::Integer, n_flips::Integer, n_runs::Integer)
    XLSX.openxlsx("datapoints.xlsx", mode="w") do xf
        sheet = xf[1]
        XLSX.rename!(sheet, "new_sheet")
        for i in 1:n_runs
            sheet["A"*string(i)] = "Run "*string(i)
            sheet["B"*string(i)] = collect_diameters(genus,points, n_diams, n_flips)
        end
    end
end

function plot_average(genus::Integer, points::Integer, n_diams::Integer, n_flips::Integer, n_runs::Integer)
    mind = ones(Int, n_diams)*(4*genus + 2*points)
    maxd = zeros(Int, n_diams)
    cum = zeros(Int, n_diams)
    for i in 1:n_runs
        diams = collect_diameters(genus, points, n_diams, n_flips)
        for j in eachindex(diams)
            if diams[j] < mind[j]
                mind[j] = diams[j]
            end
            if diams[j] > maxd[j]
                maxd[j] = diams[j]
            end
            cum[j] += diams[j]
        end
    end
    avrg = cum./n_runs
    df = DataFrame(MIN = mind, AVERAGE = avrg, MAX = maxd)
    @df df plot(:MIN, :AVERAGE, :MAX, colour = [:blue :green :red])
    #w = Window()
    #body!(w, dataviewer(df))
end

function plot_averages(genus, points, n_diams::Integer, n_flips::Integer, n_runs::Integer)
    t = length(genus)*length(points)
    df = DataFrame()
    for g in eachindex(genus)
        for p in eachindex(points)
            cum = zeros(Int, n_diams)
            for i in 1:n_runs
                diams = collect_diameters(genus[g], points[p], n_diams, n_flips)
                for j in eachindex(diams)
                    cum[j] += diams[j]
                end
            end
            df[!,string("g",genus[g],"_p",points[p])] = cum./n_runs
        end
    end
    @df df plot(cols(1:t), palette = :prism)
    w = Window()
    body!(w, dataviewer(df))
end