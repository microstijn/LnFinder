module EnvironmentTools

export fetch_omnicrobe_env!, agg_habitat, calculate_enrichment_stats, compare_archaea_bacteria

using HTTP, JSON3, DataFrames, ProgressMeter, HypothesisTests

# We cache OBT paths to avoid redundant API calls for the same habitat
const OBT_CACHE = Dict{String, String}()

function get_obt_path(obtid::String)
    if haskey(OBT_CACHE, obtid) return OBT_CACHE[obtid] end
    
    url = "https://omnicrobe.migale.inrae.fr/api/get/obt/$obtid"
    try
        r = HTTP.get(url, require_ssl_verification=false, status_exception=false)
        if r.status == 200
            json = JSON3.read(r.body)
            # The path is usually an array, e.g., ["OBT:000001/..."]
            path_str = !isempty(json.path) ? String(json.path[1]) : ""
            OBT_CACHE[obtid] = path_str
            return path_str
        end
    catch
    end
    return ""
end

function query_habitat_with_hierarchy(taxid::Int)
    base_url = "https://omnicrobe.migale.inrae.fr/api/search/relations?taxid=ncbi%3A$taxid&type=habitat"
    try
        r = HTTP.get(base_url, require_ssl_verification=false, status_exception=false)
        if r.status == 200
            json = JSON3.read(r.body)
            if !isempty(json)
                item = json[1]
                obtid = String(item.obtid)
                name = haskey(item, :obt_objects) ? String(item.obt_objects.name) : String(item.obt_forms[1])
                
                # CLIMB THE TOTEM POLE
                full_path = get_obt_path(obtid)
                
                # Standardized Mapping based on the path
                category = "Other"
                if occursin("OBT:000013", full_path) category = "Aquatic"
                elseif occursin("OBT:000006", full_path) category = "Terrestrial"
                elseif occursin("OBT:001660", full_path) || occursin("OBT:002027", full_path) category = "Host-associated"
                elseif occursin("OBT:000987", full_path) || occursin("OBT:000763", full_path) category = "Anthropogenic"
                end
                
                return (taxid, name, category)
            end
        end
    catch
    end
    return (taxid, missing, missing)
end

function fetch_omnicrobe_env!(df::DataFrame; taxid_col=:TaxID)
    unique_tids = unique(collect(skipmissing(df[!, taxid_col])))
    println("🌐 Fetching habitats and OBT hierarchies for $(length(unique_tids)) TaxIDs...")

    results = @showprogress asyncmap(query_habitat_with_hierarchy, unique_tids)
    
    tid_map = Dict(tid => (spec, cat) for (tid, spec, cat) in results)
    
    df.habitat_specific = [ismissing(t) ? missing : get(tid_map, Int(t), (missing, missing))[1] for t in df[!, taxid_col]]
    df.habitat_category = [ismissing(t) ? missing : get(tid_map, Int(t), (missing, missing))[2] for t in df[!, taxid_col]]
    
    return df
end

"""
    agg_habitat(name::Union{String, Missing})
Heuristic to roll up specific OBT leaf nodes into five primary ecological pillars.
"""
function agg_habitat(name::Union{String, Missing})
    if ismissing(name) || name == "Unknown"
        return "_Not Mapped_"
    end
    
    n = lowercase(name)
    # Priority 1: Anthropogenic
    if any(occursin.(["contaminated", "food", "broth", "agar", "industrial", "waste", "beer", "wine", "storage", "medium", "sewage", "polluted", "mine", "arsenic", "oil-", "hch-", "experimental"], n))
        return "Anthropogenic"
    # Priority 2: Host-associated
    elseif any(occursin.(["human", "plant", "root", "leaf", "host", "animal", "feces", "gut", "blood", "rhizosphere", "insect", "patient", "sputum", "lung", "skin", "nodule", "seed", "fruit", "tree", "biota", "organism", "cell", "intestine", "respiratory"], n))
        return "Host-associated"
    # Priority 3: Aquatic
    elseif any(occursin.(["water", "marine", "lake", "sea", "spring", "ocean", "ice", "glacier", "lagoon", "freshwater", "brine", "hydrothermal", "aquatic", "hotspring", "trench", "coastal"], n))
        return "Aquatic"
    # Priority 4: Terrestrial
    elseif any(occursin.(["soil", "paddy", "field", "permafrost", "sediment", "mud", "rhizoplane", "terrestrial"], n))
        return "Terrestrial"
    else
        return "Other"
    end
end

function calculate_enrichment_stats(df)
    # Filter for hits that actually have a habitat and a valid Pos 2 AA
    # We focus on the Top 5 AAs to keep the table readable
    top_aas = ["Y", "F", "W", "A", "G"]
    sub = filter(row -> row.report_group ∉ ["_Not Mapped_", "Other"] && 
                        row.pos2_aa ∈ top_aas, df)
    
    # Create contingency table
    ct = combine(groupby(sub, [:report_group, :pos2_aa]), nrow => :Count)
    matrix_df = unstack(ct, :pos2_aa, :report_group, :Count, fill=0)
    
    # Math for Z-scores
    obs = Array(matrix_df[:, 2:end])
    row_sums = sum(obs, dims=2)
    col_sums = sum(obs, dims=1)
    total = sum(obs)

    # Calculate Z-score matrix
    expected = (row_sums * col_sums) ./ total
    z_scores = (obs .- expected) ./ sqrt.(expected)
    
    # Pearson Chi-Squared Test
    pval = pvalue(ChisqTest(obs))
    
    return matrix_df, obs, z_scores, pval
end

"""
    compare_archaea_bacteria(df::DataFrame)
Vergelijkt de actieve site architectuur tussen Domains[cite: 3, 11, 12].
"""
function compare_archaea_bacteria(df::DataFrame)
    # Filter voor Lanthanide centroids in de twee hoofdgroepen [cite: 3, 22]
    sub = filter(row -> row.classification == :Lanthanide && 
                        row.kingdom in ["Bacteria", "Archaea"], df)
    
    # Maak contingentietabel voor Positie 2 [cite: 37, 38]
    ct = combine(groupby(sub, [:kingdom, :pos2_aa]), nrow => :Count)
    matrix_df = unstack(ct, :pos2_aa, :kingdom, :Count, fill=0)
    
    # Voer Chi-Squared test uit op de counts
    obs = Matrix{Int}(matrix_df[:, [:Bacteria, :Archaea]])
    p_val = pvalue(ChisqTest(obs))
    
    return matrix_df, p_val
end

end