module EnvironmentTools

export fetch_omnicrobe_env!

using HTTP
using JSON3
using DataFrames
using ProgressMeter

# --- ONTOLOGY HUB MAPPING ---
# These OBT IDs are the "Major Roots" in the OntoBiotope hierarchy.
# We use these to group specific terms (like "creek sediment") into broad categories.
const HABITAT_MAP = Dict(
    "OBT:000006" => "Soil",                 # Terrestrial
    "OBT:001773" => "Soil",                 # Soil variant
    "OBT:000013" => "Water",                # Aquatic
    "OBT:001645" => "Water",                # Aquatic variant
    "OBT:000008" => "Host-associated",      # Living organisms (Animal/Plant)
    "OBT:001660" => "Host-associated",      # Animal variant
    "OBT:000987" => "Artificial/Industrial",# Constructed habitats
    "OBT:001254" => "Food/Agriculture",     # Processed/Cultivated
    "OBT:000763" => "Polluted/Waste"        # Contaminated sites/Sludge
)

"""
    query_omnicrobe(taxid::Int, base_url::String)
    
Internal helper to fetch data for a single TaxID.
"""
function query_omnicrobe(taxid::Int, base_url::String)
    # We filter specifically for 'habitat' to avoid phenotypic noise
    full_url = "$base_url$taxid&type=habitat"
    
    try
        # require_ssl_verification=false is required for migale.inrae.fr
        r = HTTP.get(full_url, 
                     require_ssl_verification = false, 
                     status_exception = false, 
                     readtimeout = 10)

        if r.status == 200
            json_obj = JSON3.read(r.body)
            if isempty(json_obj)
                return (taxid, missing, missing)
            end
            
            # Use the first (most relevant) result
            item = json_obj[1]
            
            # 1. Get Specific Habitat (e.g., "rhizosphere")
            # Prefer canonical name, fallback to raw surface form
            specific = missing
            if haskey(item, :obt_objects) && haskey(item.obt_objects, :name)
                specific = String(item.obt_objects.name)
            elseif haskey(item, :obt_forms) && !isempty(item.obt_forms)
                specific = String(item.obt_forms[1])
            end

            # 2. Map to Broad Category using the OBT Hubs
            category = "Other"
            # The 'obt_root' usually contains the full OBT path
            path = haskey(item, :obt_root) ? String(item.obt_root) : ""
            
            for (hub_id, label) in HABITAT_MAP
                if occursin(hub_id, path) || (haskey(item, :obtid) && item.obtid == hub_id)
                    category = label
                    break
                end
            end
            
            return (taxid, specific, category)
        end
        return (taxid, missing, missing)
    catch e
        # Silently fail for individual taxids to allow the pipeline to continue
        return (taxid, missing, missing)
    end
end

"""
    fetch_omnicrobe_env!(df::DataFrame; taxid_col=:TaxID)
    
Queries OmniMicrobe for unique TaxIDs in the DataFrame and appends
standardized environment columns.
"""
function fetch_omnicrobe_env!(df::DataFrame; taxid_col=:TaxID)
    if string(taxid_col) ∉ names(df)
        @error "TaxID column '$taxid_col' not found."
        return df
    end

    base_url = "https://omnicrobe.migale.inrae.fr/api/search/relations?taxid=ncbi%3A"
    
    # Extract unique, non-missing TaxIDs
    unique_tids = unique(collect(skipmissing(df[!, taxid_col])))
    
    println("🌐 Querying OmniMicrobe for $(length(unique_tids)) unique habitats...")
    
    # Process concurrently using asyncmap
    results = @showprogress asyncmap(tid -> query_omnicrobe(Int(tid), base_url), unique_tids)
    
    # Create lookup dictionaries
    tid_to_spec = Dict{Int, Union{String, Missing}}()
    tid_to_cat = Dict{Int, Union{String, Missing}}()
    
    for (tid, spec, cat) in results
        tid_to_spec[tid] = spec
        tid_to_cat[tid] = cat
    end

    # Initialize new columns with 'missing'
    df.habitat_specific = Vector{Union{String, Missing}}(missing, nrow(df))
    df.habitat_category = Vector{Union{String, Missing}}(missing, nrow(df))

    # Assign results back to the main DataFrame
    for i in 1:nrow(df)
        tid = df[i, taxid_col]
        if !ismissing(tid) && haskey(tid_to_spec, Int(tid))
            df[i, :habitat_specific] = tid_to_spec[Int(tid)]
            df[i, :habitat_category] = tid_to_cat[Int(tid)]
        end
    end

    # Cleanup: If specific habitat is missing, category should be missing too
    # (prevents "Other" appearing for items with no data at all)
    for i in 1:nrow(df)
        if ismissing(df[i, :habitat_specific])
            df[i, :habitat_category] = missing
        end
    end

    println("✅ Environmental metadata appended.")
    return df
end

end # module