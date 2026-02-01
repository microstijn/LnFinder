module AccessionMapper

export parse_fasta_headers, fetch_taxids_async, fetch_uniprot_annotations

using HTTP
using JSON3
using DataFrames
using CodecZlib

# ------------------------------------------------------------------
# 1. FASTA Parsing
# ------------------------------------------------------------------
function parse_fasta_headers(fasta_path::String)
    headers = String[]
    accessions = String[]
    
    println("📖 Reading FASTA file: $fasta_path")
    if !isfile(fasta_path)
        error("File not found: $fasta_path")
    end

    open(fasta_path, "r") do io
        for line in eachline(io)
            if startswith(line, ">")
                raw_header = strip(line[2:end])
                acc = ""
                
                # Handle >sp|ACC|ID, >UniRef_ACC, and >ACC_ID formats
                if occursin('|', raw_header)
                    parts = split(raw_header, '|')
                    acc = length(parts) >= 2 ? String(parts[2]) : String(parts[1])
                elseif startswith(raw_header, "UniRef")
                    parts = split(raw_header, '_')
                    acc = length(parts) >= 2 ? String(parts[2]) : raw_header
                else
                    m = match(r"^([A-Z0-9]+)", raw_header)
                    acc = isnothing(m) ? raw_header : String(m.captures[1])
                end

                acc = split(acc, ' ')[1] # Remove trailing descriptions
                push!(headers, raw_header)
                push!(accessions, acc)
            end
        end
    end
    
    u_len = length(unique(accessions))
    println("✅ Found $(length(accessions)) sequences ($u_len unique).")
    return DataFrame(original_header = headers, accession = accessions)
end

# ------------------------------------------------------------------
# 2. Async Job Handling
# ------------------------------------------------------------------
function submit_job(ids::Vector{String}, from_db::String, to_db::String)
    url = "https://rest.uniprot.org/idmapping/run"
    params = Dict("from" => from_db, "to" => to_db, "ids" => join(ids, ","))
    
    resp = HTTP.post(url, body=HTTP.Form(params))
    return JSON3.read(resp.body).jobId
end

function wait_for_job(job_id::String)
    status_url = "https://rest.uniprot.org/idmapping/status/$job_id"
    print("⏳ Waiting for UniProt ($job_id)...")
    
    while true
        resp = HTTP.get(status_url, redirect=false, status_exception=false)
        if resp.status == 303
            println(" Done!")
            # Convert SubString to String to avoid type errors
            return String(HTTP.header(resp, "Location"))
        elseif resp.status == 200
            status = JSON3.read(resp.body).jobStatus
            if status == "FINISHED"
                println(" Done!")
                return "https://rest.uniprot.org/idmapping/results/$job_id"
            elseif status == "FAILED"
                error("UniProt Job Failed.")
            else
                sleep(2); print(".")
            end
        else
            error("HTTP Error: $(resp.status)")
        end
    end
end

function download_json_stream(result_url::AbstractString)
    # Ensure we get JSON and Compression
    stream_url = replace(result_url, "/results/" => "/results/stream/")
    final_url = "$stream_url?format=json&compressed=true"
    
    println("⬇️  Downloading response...")
    resp = HTTP.get(final_url, status_exception=false)
    
    if resp.status != 200
        error("Download failed: $(resp.status)")
    end
    
    return JSON3.read(transcode(GzipDecompressor, resp.body))
end




function fetch_taxids_async(accessions::Vector{String})
    unique_accs = unique(accessions)
    if isempty(unique_accs); return Dict{String, Int}(); end

    println("\n🚀 PHASE 1: Mapping $(length(unique_accs)) Accessions to UniProtKB...")
    
    # --- PHASE 1: Map to UniProtKB ---
    job1 = submit_job(unique_accs, "UniProtKB_AC-ID", "UniProtKB")
    url1 = wait_for_job(job1)
    json1 = download_json_stream(url1)
    
    acc_to_taxid = Dict{String, Int}()
    rescue_map = Dict{String, String}() # UPI -> Original_Accession
    
    for item in json1.results
        acc = String(item.from)
        
        if haskey(item, :to)
            target = item.to
            
            # Case A: Active Entry
            if haskey(target, :organism) && haskey(target.organism, :taxonId)
                acc_to_taxid[acc] = target.organism.taxonId
            
            # Case B: Inactive/Deleted Entry (Save UPI for Rescue)
            elseif get(target, :entryType, "") == "Inactive" && haskey(target, :extraAttributes)
                if haskey(target.extraAttributes, :uniParcId)
                    upi = String(target.extraAttributes.uniParcId)
                    rescue_map[upi] = acc
                end
            end
        end
    end

    # --- PHASE 2: Rescue via UniParc ---
    if !isempty(rescue_map)
        upis_to_fetch = collect(keys(rescue_map))
        println("\n⚠️  Found $(length(upis_to_fetch)) deleted entries. Attempting rescue via UniParc...")
        
        # FIX: We map "UniParc" -> "UniParc" (because mapping to TaxID directly is forbidden)
        job2 = submit_job(upis_to_fetch, "UniParc", "UniParc")
        url2 = wait_for_job(job2)
        json2 = download_json_stream(url2)
        
        rescued_count = 0
        for item in json2.results
            upi = String(item.from)
            
            # The result is the full UniParc entry
            if haskey(item, :to)
                uni_entry = item.to
                
                # UniParc entries have an 'organisms' array
                if haskey(uni_entry, :organisms) && !isempty(uni_entry.organisms)
                    # We take the first organism as the "best guess" for the deleted sequence
                    taxon_obj = uni_entry.organisms[1]
                    
                    if haskey(taxon_obj, :taxonId)
                        taxid = taxon_obj.taxonId
                        
                        # Map back to original accession
                        if haskey(rescue_map, upi)
                            original_acc = rescue_map[upi]
                            acc_to_taxid[original_acc] = taxid
                            rescued_count += 1
                        end
                    end
                end
            end
        end
        println("🩹 Rescued $rescued_count IDs via UniParc!")
    end
    
    mapped = length(acc_to_taxid)
    total = length(unique_accs)
    println("\n✅ Final Report: Mapped $mapped / $total IDs.")
    
    return acc_to_taxid
end


function fetch_uniprot_annotations(fasta_path::String)
    # 1. Parse IDs from the specific FASTA header format (>ID|Kingdom|Habitat|Pos2)
    accessions = String[]
    open(fasta_path, "r") do io
        for line in eachline(io)
            if startswith(line, ">")
                push!(accessions, String(split(strip(line[2:end]), "|")[1]))
            end
        end
    end
    
    # 2. Submit Mapping Job 
    job_id = AccessionMapper.submit_job(accessions, "UniProtKB_AC-ID", "UniProtKB")
    results_url = AccessionMapper.wait_for_job(job_id) 
    json_data = AccessionMapper.download_json_stream(results_url) 

    records = []
    for item in json_data.results
        entry = get(item, :to, nothing)
        isnothing(entry) && continue
        
        # --- NAME EXTRACTION ---
        p_desc = get(entry, :proteinDescription, nothing)
        name = "Unknown Protein"
        if !isnothing(p_desc)
            if haskey(p_desc, :recommendedName)
                name = p_desc.recommendedName.fullName.value
            elseif haskey(p_desc, :submissionNames) && !isempty(p_desc.submissionNames)
                # Unreviewed entries store the name here
                name = p_desc.submissionNames[1].fullName.value
            end
        end

        # --- GENE / ORF EXTRACTION ---
        gene_id = "N/A"
        if haskey(entry, :genes) && !isempty(entry.genes)
            g_obj = entry.genes[1]
            if haskey(g_obj, :geneName)
                gene_id = g_obj.geneName.value
            elseif haskey(g_obj, :orfNames) && !isempty(g_obj.orfNames)
                # Common in Archaea/Bacteria samples provided
                gene_id = g_obj.orfNames[1].value
            end
        end

        # --- FUNCTIONAL DESCRIPTION (Cross-references) ---
        # We look for InterPro or Pfam EntryNames to act as a functional description
        functional_desc = "PQQ-dependent dehydrogenase" # Default based on your study
        if haskey(entry, :uniProtKBCrossReferences)
            for xref in entry.uniProtKBCrossReferences
                if xref.database in ["InterPro", "Pfam", "SUPFAM"]
                    if haskey(xref, :properties)
                        # Find the EntryName property in the array
                        prop = findfirst(p -> p.key == "EntryName", xref.properties)
                        if !isnothing(prop)
                            functional_desc = xref.properties[prop].value
                            break # Take the first good match
                        end
                    end
                end
            end
        end

        # --- METADATA ---
        org = get(get(entry, :organism, Dict()), :scientificName, "Unknown")
        
        push!(records, (
            accession = item.from,
            protein_name = name,
            functional_description = functional_desc,
            gene_id = gene_id,
            organism = org
        ))
    end
    
    return DataFrame(records)
end

end