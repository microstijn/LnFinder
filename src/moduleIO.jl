using FASTX
using DataFrames

"""
    classify_hits_from_fasta(fasta_path::AbstractString)

Reads a FASTA file (downloaded from HMMER web server) containing PQQ-ADH hits.
Classifies each sequence as Lanthanide, Calcium, or Unknown.

Returns a DataFrame with columns: `id`, `motif`, `sequence`.
"""
function classify_hits_from_fasta(fasta_path::AbstractString)
    # Initialize output vectors
    ids = String[]
    motifs = MotifType[]
    sequences = String[]

    # Open the FASTA file using FASTX.FASTA.Reader
    reader = FASTA.Reader(open(fasta_path))

    for record in reader
        # Extract data
        seq_id = FASTA.identifier(record)
        seq_data = FASTA.sequence(String, record) # Get sequence as String
        
        # Run our classifier (The "Funnel" Step 2)
        # Note: We convert to LongAA inside the classifier wrapper we wrote earlier
        motif = classify_pqq_adh(seq_data)

        # Store results
        push!(ids, seq_id)
        push!(motifs, motif)
        push!(sequences, seq_data)
    end
    
    close(reader)

    return DataFrame(id=ids, motif=motifs, sequence=sequences)
end