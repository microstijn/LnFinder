using FASTX
using DataFrames

using FASTX
using DataFrames
using Statistics

using FASTX
using DataFrames

"""
    analyze_alignment(fasta_path::AbstractString; sample_size=500)

Optimized for large files (350MB+).
1. reads a small sample to find the active site column.
2. Streams the rest of the file to classify sequences without loading them all into RAM.
"""
function analyze_alignment(fasta_path::AbstractString; sample_size=500)
    # Find Consensus Column via Sampling
    best_idx = 0
    max_score = 0
    
    # We only read the first `sample_size` sequences to find the column
    # This turns the O(N*L) search into O(sample_size*L), which is instant.
    FASTA.Reader(open(fasta_path)) do reader
        sample_seqs = String[]
        count = 0
        
        for record in reader
            count += 1
            push!(sample_seqs, FASTA.sequence(String, record))
            if count >= sample_size
                break 
            end
        end
        
        if isempty(sample_seqs)
            error("File is empty or could not be read.")
        end

        # Calculate consensus on the sample
        seq_len = length(sample_seqs[1])
        
        # Scan columns (heuristic: skip first/last 10% to save time)
        start_search = div(seq_len, 10)
        end_search = seq_len - div(seq_len, 10)

        for i in start_search:end_search
            score = 0
            for seq in sample_seqs
                # Bounds check
                if i+2 > length(seq); continue; end
                
                # Check for D-x-[D/A/T/S]
                if seq[i] == 'D'
                    switch = seq[i+2]
                    if switch == 'D' || switch in ['A', 'T', 'S']
                        score += 1
                    end
                end
            end
            
            if score > max_score
                max_score = score
                best_idx = i
            end
        end
    end

    # Stream & Classify 
    ids = String[]
    motifs = MotifType[]
    contexts = String[] # Store the 3-letter motif (e.g. "DLD")

    # Re-open file to stream from the start
    FASTA.Reader(open(fasta_path)) do reader
        for record in reader
            seq = FASTA.sequence(String, record)
            
            # Fast fail if sequence is too short or weird
            if best_idx == 0 || best_idx+2 > length(seq)
                push!(ids, FASTA.identifier(record))
                push!(motifs, LnFinder.Unknown)
                push!(contexts, "---")
                continue
            end

            # Direct O(1) Access - No searching!
            anchor = seq[best_idx]
            switch = seq[best_idx+2]
            
            # Classification
            if anchor == 'D' && switch == 'D'
                push!(motifs, LnFinder.Lanthanide)
            elseif anchor == 'D' && (switch == 'A' || switch == 'T' || switch == 'S')
                push!(motifs, LnFinder.Calcium)
            else
                push!(motifs, LnFinder.Unknown)
            end
            
            # Store ID and Motif context
            push!(ids, FASTA.identifier(record))
            push!(contexts, "$anchor$(seq[best_idx+1])$switch")
        end
    end

    return DataFrame(id=ids, motif=motifs, motif_seq=contexts), best_idx
end

