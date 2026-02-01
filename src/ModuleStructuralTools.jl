module StructuralTools

using HTTP, FASTX, DataFrames, JSON3, ProgressMeter, Random

export download_alphafold_structures, prepare_foldmason_subsample

"""
    get_af_metadata(accession::AbstractString)
Queries the AlphaFold API to get metadata for a specific UniProt accession.
"""
function get_af_metadata(accession::AbstractString)
    # Clean ID: A0A8T4GX75_9EURY -> A0A8T4GX75
    clean_id = split(split(accession, '_')[1], '-')[1]
    url = "https://alphafold.ebi.ac.uk/api/prediction/$(clean_id)"
    
    try
        response = HTTP.get(url, status_exception=false)
        if response.status == 200
            return JSON3.read(response.body)
        end
    catch
    end
    return nothing
end

"""
    download_alphafold_structures(df::DataFrame, out_dir::AbstractString; acc_col=:accession)
Uses the AlphaFold API to retrieve the correct PDB URLs and download them.
"""
function download_alphafold_structures(df::DataFrame, out_dir::AbstractString; acc_col=:accession)
    if !isdir(out_dir) mkpath(out_dir) end

    unique_ids = unique(collect(skipmissing(df[!, acc_col])))
    results = Dict{String, Union{String, Nothing}}()
    
    println("📡 Querying AlphaFold API for $(length(unique_ids)) structures...")

    @showprogress for id in unique_ids
        full_id_str = String(id)
        
        # 1. Get Metadata from API
        meta = get_af_metadata(full_id_str)
        
        if isnothing(meta) || isempty(meta)
            results[full_id_str] = nothing
            continue
        end

        # Extract the PDB URL from the first model result
        # The API returns an array of objects; we want pdbUrl
        pdb_url = meta[1][:pdbUrl]
        clean_acc = split(split(full_id_str, '_')[1], '-')[1]
        dest_path = joinpath(out_dir, "$(clean_acc).pdb")

        # Download if not already present
        if isfile(dest_path)
            results[full_id_str] = dest_path
            continue
        end

        try
            dl_res = HTTP.get(pdb_url, status_exception=false)
            if dl_res.status == 200
                open(dest_path, "w") do io
                    write(io, dl_res.body)
                end
                results[full_id_str] = dest_path
            else
                results[full_id_str] = nothing
            end
        catch
            results[full_id_str] = nothing
        end
    end

    success_count = count(!isnothing, values(results))
    println("✅ Found and downloaded $success_count / $(length(unique_ids)) structures.")
    return results
end


function prepare_foldmason_subsample(input_file, output_file; n_bact=60, n_euk=20, seed=42)
    # 1. Read all records
    records = []
    open(FASTA.Reader, input_file) do reader
        for record in reader
            push!(records, record)
        end
    end

    # 2. Filter by Kingdom
    archaea = filter(r -> occursin("|Archaea|", FASTA.identifier(r)), records)
    bacteria = filter(r -> occursin("|Bacteria|", FASTA.identifier(r)), records)
    eukarya  = filter(r -> occursin("|Eukaryota|", FASTA.identifier(r)), records)

    println("Found: $(length(archaea)) Archaea, $(length(bacteria)) Bacteria, $(length(eukarya)) Eukarya.")

    # 3. Create reproducible selection
    Random.seed!(seed) 
    
    final_selection = []
    append!(final_selection, archaea) # All 761 Archaeal hits for maximum coverage
    append!(final_selection, shuffle(bacteria)[1:min(n_bact, end)])
    append!(final_selection, shuffle(eukarya)[1:min(n_euk, end)])

    # 4. Corrected Writer Block (Removed the "w")
    open(FASTA.Writer, output_file) do writer
        for record in final_selection
            write(writer, record)
        end
    end

    println("Success! $(length(final_selection)) sequences saved to $output_file")
end

end # module