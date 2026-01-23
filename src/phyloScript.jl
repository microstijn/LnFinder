#= 
Preamble & environment management
=# 

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

#= 
Packages
=# 

using Revise
using LnFinder
using CSV
using DataFrames

paths = (
    data_in  = "D:/hmmerStuff/out/final_deduplicated_data.csv",
    fasta_in = raw"D:\hmmerStuff\data\firstTry.fa",
    fasta_out = "D:/hmmerStuff/out/tree_input_known.fasta",
    itol_dir = "D:/hmmerStuff/out/itol_assets/"
)


mkpath(paths.itol_dir)


println("📂 Loading data...")
df = CSV.read(paths.data_in, DataFrame)


df_known = filter(row -> !ismissing(row.kingdom) && 
                         row.kingdom ∈ ["Bacteria", "Archaea", "Eukaryota"], df)

println("📉 Reduced dataset to $(nrow(df_known)) sequences with known taxonomy.")


println("🧬 Exporting FASTA...")
export_tree_fasta(df_known, paths.fasta_in, paths.fasta_out)

println("🎨 Generating iTOL annotation files...")

function generate_itol_habitat_strip(df, out_path)
    colors = Dict(
        "Aquatic" => "#0077be", "Terrestrial" => "#8b4513",
        "Host-associated" => "#228b22", "Anthropogenic" => "#ff4500",
        "Other" => "#808080", "_Not Mapped_" => "#d3d3d3"
    )
    open(out_path, "w") do io
        println(io, "DATASET_COLORSTRIP\nSEPARATOR TAB\nDATASET_LABEL\tHabitat\nCOLOR\t#ff0000\nDATA")
        for row in eachrow(df)
            header = "$(row.id)|$(row.kingdom)|$(replace(coalesce(row.report_group, "Other"), " " => "_"))|$(row.pos2_aa)"
            color = get(colors, row.report_group, "#808080")
            println(io, "$header\t$color\t$(row.report_group)")
        end
    end
end

generate_itol_habitat_strip(df_known, joinpath(paths.itol_dir, "habitat_strip.txt"))

function generate_itol_kingdom_strip(df, out_path)
    # Professional color palette for Kingdoms
    colors = Dict(
        "Bacteria"  => "#3498db", # Bright Blue
        "Archaea"   => "#e74c3c", # Vibrant Red
        "Eukaryota" => "#2ecc71", # Nature Green
        "Viruses"   => "#9b59b6"  # Purple (just in case)
    )

    open(out_path, "w") do io
        println(io, "DATASET_COLORSTRIP")
        println(io, "SEPARATOR TAB")
        println(io, "DATASET_LABEL\tKingdom")
        println(io, "COLOR\t#000000") # Legend color
        println(io, "DATA")
        
        for row in eachrow(df)
            # IMPORTANT: This must match your tree headers exactly
            # If your tree uses ID|Kingdom|Habitat|Pos2, use that here.
            # If it's just the ID, use: header = row.id
            header = "$(row.id)|$(row.kingdom)|$(replace(coalesce(row.report_group, "Other"), " " => "_"))|$(row.pos2_aa)"
            
            kingdom = row.kingdom
            color = get(colors, kingdom, "#808080")
            println(io, "$header\t$color\t$kingdom")
        end
    end
end

generate_itol_kingdom_strip(df_known, joinpath(paths.itol_dir, "itol_kingdoms.txt"))

using CSV, DataFrames


function generate_itol_function_strip_v2(df, out_path)
    open(out_path, "w") do io
        println(io, "DATASET_COLORSTRIP")
        println(io, "SEPARATOR TAB")
        println(io, "DATASET_LABEL\tActive_Site_Chemistry")
        println(io, "COLOR\t#ff00ff")
        println(io, "DATA")
        
        for row in eachrow(df)
            header = "$(row.id)|$(row.kingdom)|$(replace(coalesce(row.report_group, "Other"), " " => "_"))|$(row.pos2_aa)"
            
            # Robust classification check
            # D-x-D = Lanthanide (Magenta)
            # D-x-[ATS] = Calcium (Yellow)
            p3 = string(coalesce(row.pos3_aa, ""))
            
            color = if p3 == "D"
                "#ff00ff" # Magenta (Lanthanide)
            elseif p3 ∈ ["A", "T", "S"]
                "#f1c40f" # Yellow (Calcium)
            else
                "#7f8c8d" # Dark Grey (Other/Unknown)
            end
            
            label = p3 == "D" ? "Lanthanide" : (p3 ∈ ["A", "T", "S"] ? "Calcium" : "Other")
            println(io, "$header\t$color\t$label")
        end
    end
end

generate_itol_function_strip(df_known, joinpath(paths.itol_dir, "itol_mptifs.txt"))
println("✅ iTOL Function file generated at D:/hmmerStuff/out/itol_motifs.txt")

ln_users = filter(row -> row.classification == "Calcium", df_known)

