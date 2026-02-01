module StructuralUtils

using FASTX
using DataFrames

export organize_pdb_files

"""
    organize_pdb_files(source_dir, output_dir, df)

Sorts AlphaFold PDB files based on the 'accession', 'kingdom', and 
'classification' columns in the provided DataFrame.
"""
function organize_pdb_files(source_dir, output_dir, df::DataFrame)
    # Ensure source exists
    if !isdir(source_dir)
        error("Source directory not found: $source_dir")
    end

    # AlphaFold naming convention: AF-[ACCESSION]-F1-model_v4.pdb
    id_regex = r"AF-(.*?)-F\d"

    for file_name in readdir(source_dir)
        if endswith(file_name, ".pdb")
            
            # 1. Extract Accession from filename
            m = match(id_regex, file_name)
            accession = isnothing(m) ? split(file_name, ".")[1] : m.captures[1]
            
            # 2. Lookup in DataFrame
            # We use the 'accession' column for the most reliable match
            row = filter(r -> r.accession == accession, df)
            
            if !isempty(row)
                # Handle "missing" or actual missing values
                k = row[1, :kingdom]
                kingdom = (ismissing(k) || k == "missing") ? "Unknown" : k
                
                classification = row[1, :classification]
                
                # 3. Create target directory: root/Kingdom/Classification
                target_dir = joinpath(output_dir, string(kingdom), string(classification))
                mkpath(target_dir)
                
                # 4. Copy file
                cp(joinpath(source_dir, file_name), joinpath(target_dir, file_name), force=true)
            else
                # Optional: log files not found in the DF
                # println("Accession $accession not found in DataFrame.")
            end
        end
    end
    println("✅ PDB files organized into $output_dir based on DataFrame classification.")
end
end # module