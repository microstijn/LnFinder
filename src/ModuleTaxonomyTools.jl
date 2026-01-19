module TaxonomyTools

export load_taxonomy_db, get_lineage, append_taxonomy!, merge_taxonomy_with_hits, filter_centroids

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

using DataFrames

"""
    merge_taxonomy_with_hits(hits_df, tax_df; on_col=:TaxID)

Joins the HMM results with the OmniMicrobe/Taxonomy metadata. 
We use a left join to ensure no protein hits are lost during the merge.
"""
function merge_taxonomy_with_hits(hits_df::DataFrame, tax_df::DataFrame; on_col=:TaxID)
    # 1. Ensure the join column is the same type (Int64)
    # This prevents 'No matching method' errors during join
    hits_df[!, on_col] = Int64.(hits_df[!, on_col])
    tax_df[!, on_col] = Int64.(tax_dSf[!, on_col])

    println("🔗 Merging $(nrow(hits_df)) hits with taxonomic metadata...")

    # 2. Perform Left Join
    # hits_df is the 'Left' table (the primary data)
    # tax_df is the 'Right' table (the metadata)
    merged_df = leftjoin(hits_df, tax_df, on=on_col)

    # 3. Post-merge Cleanup
    # Sometimes tax_df has redundant columns; we ensure 'habitat_category' exists
    if !("habitat_category" in names(merged_df))
        @warn "Metadata columns not found in merged result. Check tax_df headers."
    end

    total_mapped = count(.!ismissing.(merged_df.habitat_specific))
    println("✅ Merge complete. $total_mapped / $(nrow(merged_df)) hits now have habitat data.")

    return merged_df
end

"""
    filter_centroids(df, uc_path)
Filters the main DataFrame to keep only the cluster representatives from USEARCH.
"""
function filter_centroids(df, uc_path)
    centroids = String[]
    for line in readlines(uc_path)
        parts = split(line, '\t')
        if parts[1] == "S"  # S indicates a Centroid (Seed)
            # Column 9 is the sequence label in UC format
            push!(centroids, parts[9])
        end
    end
    
    # Filter the original dataframe
    dedup_df = filter(row -> row.id in centroids, df)
    
    println("✂️ Deduplication Summary:")
    println("  Total sequences: $(nrow(df))")
    println("  Cluster Centroids: $(nrow(dedup_df))")
    println("  Reduction: $(round(100 * (1 - nrow(dedup_df)/nrow(df)), digits=1))%")
    
    return dedup_df
end

end # module