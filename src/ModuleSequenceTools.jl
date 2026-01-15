
module SequenceTools

using BioSequences
using DataFrames

export count_residues_per_position

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

end