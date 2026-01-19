
module SequenceTools

using BioSequences
using DataFrames
using FASTX

export count_residues_per_position, export_tree_fasta

"""
    count_residues_per_position(motif_column::Vector{String})
    
Returns a DataFrame where rows are Amino Acids and columns are 
the probability/count at each of the 3 positions in the D-x-D motif.
"""
function count_residues_per_position(motif_column)
    # 1. Filter out any remaining gaps or truncated entries
    valid_motifs = [m for m in motif_column if length(m) == 3 && !occursin("-", m)]
    total = length(valid_motifs)
    
    # 2. Get the standard amino acid alphabet
    amino_acids = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", 
                   "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]
    
    # 3. Initialize counts
    counts = Dict(aa => [0, 0, 0] for aa in amino_acids)
    
    # 4. Populate counts
    for motif in valid_motifs
        for pos in 1:3
            aa = string(motif[pos])
            if haskey(counts, aa)
                counts[aa][pos] += 1
            end
        end
    end
    
    # 5. Convert to DataFrame
    df_counts = DataFrame(AminoAcid = amino_acids)
    df_counts.Pos1 = [counts[aa][1] for aa in amino_acids]
    df_counts.Pos2 = [counts[aa][2] for aa in amino_acids]
    df_counts.Pos3 = [counts[aa][3] for aa in amino_acids]
    
    # 6. Add percentages for the presentation
    df_counts.Pct1 = round.((df_counts.Pos1 ./ total) .* 100, digits=2)
    df_counts.Pct2 = round.((df_counts.Pos2 ./ total) .* 100, digits=2)
    df_counts.Pct3 = round.((df_counts.Pos3 ./ total) .* 100, digits=2)
    
    return sort(df_counts, :Pos2, rev=true)
end


"""
    export_tree_fasta(df_dedup::DataFrame, input_fasta::String, output_fasta::String)
Streams sequences van schijf en schrijft een FASTA met metadata-headers. 
Zorgt voor correcte regeleinden tussen records.
"""
function export_tree_fasta(df_dedup::DataFrame, input_fasta::String, output_fasta::String)
    # Maak een set van de ID's die we willen behouden (centroids) [cite: 22]
    target_ids = Set(df_dedup.id)
    
    # Map ID's naar metadata voor de header-verrijking [cite: 10, 15]
    metadata = Dict(row.id => (
        kingdom = coalesce(row.kingdom, "Unknown"),
        habitat = replace(coalesce(row.report_group, "Other"), " " => "_"),
        pos2 = row.pos2_aa
    ) for row in eachrow(df_dedup))

    println("🧵 Streaming sequences van $input_fasta...")
    
    count = 0
    # Gebruik de FASTX Reader om de originele fasta te lezen 
    FASTA.Reader(open(input_fasta)) do reader
        # Open de output file in write-modus
        open(output_fasta, "w") do writer
            for record in reader
                id = FASTA.identifier(record)
                
                if id in target_ids
                    meta = metadata[id]
                    # Formatteer de nieuwe header: ID|Kingdom|Habitat|Pos2 [cite: 59, 63]
                    new_header = "$(id)|$(meta.kingdom)|$(meta.habitat)|$(meta.pos2)"
                    
                    # Haal de sequentie op uit het huidige record 
                    seq = FASTA.sequence(record)
                    
                    # Schrijf het record handmatig om volledige controle over regeleinden te hebben
                    write(writer, ">", new_header, "\n")
                    write(writer, seq, "\n")
                    
                    count += 1
                end
            end
        end
    end
    println("✅ Export voltooid. $count centroids weggeschreven naar $output_fasta")
end

end