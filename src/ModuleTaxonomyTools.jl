module TaxonomyTools

export load_taxonomy_db, get_lineage, append_taxonomy!

using DataFrames
using ProgressMeter

"""
    TaxonomyDB
    Holds the parent-child relationships and names in memory.
"""
struct TaxonomyDB
    parents::Dict{Int, Int}
    ranks::Dict{Int, String}
    names::Dict{Int, String}
end

"""
    load_taxonomy_db(dump_path::String)
    Parses NCBI nodes.dmp and names.dmp.
"""
function load_taxonomy_db(dump_path::String)
    nodes_file = joinpath(dump_path, "nodes.dmp")
    names_file = joinpath(dump_path, "names.dmp")

    if !isfile(nodes_file) || !isfile(names_file)
        error("NCBI dump files not found in $dump_path. Ensure nodes.dmp and names.dmp are there.")
    end

    parents = Dict{Int, Int}()
    ranks = Dict{Int, String}()
    names_map = Dict{Int, String}()

    println("📂 Loading NCBI Taxonomy Nodes...")
    # nodes.dmp format: tax_id | parent_id | rank | ...
    open(nodes_file, "r") do io
        for line in eachline(io)
            parts = split(line, "|")
            tid = parse(Int, strip(parts[1]))
            pid = parse(Int, strip(parts[2]))
            rnk = strip(parts[3])
            parents[tid] = pid
            ranks[tid] = rnk
        end
    end

    println("📂 Loading NCBI Taxonomy Names...")
    # names.dmp format: tax_id | name_txt | unique_name | name_class |
    open(names_file, "r") do io
        for line in eachline(io)
            # Only keep "scientific name" to avoid synonyms
            if occursin("scientific name", line)
                parts = split(line, "|")
                tid = parse(Int, strip(parts[1]))
                nam = strip(parts[2])
                names_map[tid] = nam
            end
        end
    end

    println("✅ Taxonomy DB loaded ($(length(parents)) nodes).")
    return TaxonomyDB(parents, ranks, names_map)
end

"""
    get_lineage(taxid::Int, db::TaxonomyDB)
    Reconstructs the full lineage (Kingdom to Species) for a given TaxID.
"""
function get_lineage(taxid::Int, db::TaxonomyDB)
    # Ranks we care about for the final table
    target_ranks = ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]
    lineage = Dict{Symbol, String}()
    
    current_id = taxid
    # Safety counter to prevent infinite loops on corrupted dmp files
    max_depth = 50 
    
    while current_id != 1 && max_depth > 0
        rank = get(db.ranks, current_id, "no rank")
        
        if rank in target_ranks
            # Rename superkingdom to kingdom for common convention
            key = rank == "superkingdom" ? :kingdom : Symbol(rank)
            lineage[key] = get(db.names, current_id, "Unknown")
        end
        
        parent_id = get(db.parents, current_id, 1)
        if parent_id == current_id; break; end # Reached root
        current_id = parent_id
        max_depth -= 1
    end
    
    return isempty(lineage) ? missing : lineage
end

"""
    append_taxonomy!(df::DataFrame, db::TaxonomyDB; taxid_col=:TaxID)
    Modifies the DataFrame in-place to add taxonomy columns.
"""
function append_taxonomy!(df::DataFrame, db::TaxonomyDB; taxid_col=:TaxID)
    println("🧬 Appending Taxonomic Lineages...")
    
    # Initialize columns as Union{String, Missing}
    ranks = [:kingdom, :phylum, :class, :order, :family, :genus, :species]
    for r in ranks
        df[!, r] = Vector{Union{String, Missing}}(missing, nrow(df))
    end

    # Progress bar for long datasets
    @showprogress for i in 1:nrow(df)
        tid = df[i, taxid_col]
        
        # skipmissing handles our 'missing' TaxIDs from earlier steps
        if !ismissing(tid)
            lin = get_lineage(Int(tid), db)
            if !ismissing(lin)
                for r in ranks
                    if haskey(lin, r)
                        df[i, r] = lin[r]
                    end
                end
            end
        end
    end
    return df
end

end # module