using DataFrames
using CSV

"""
    run_hmmer(hmm_profile, target_sequences, output_table; cpu=1)

Runs `hmmsearch` via the command line.
Arguments are ordered: Input -> Output.
"""
function run_hmmer(hmm_profile::AbstractString, target_sequences::AbstractString, output_table::AbstractString; cpu::Integer=1)
    if isnothing(Sys.which("hmmsearch"))
        error("hmmsearch not found! Please install HMMER.")
    end

    # Use `cmd` backticks for safe shell execution
    cmd = `hmmsearch --tblout $output_table --cpu $cpu -E 1e-10 $hmm_profile $target_sequences`
    run(cmd)
    
    return output_table
end

"""
    parse_hmmer_tblout(file_path)

Parses the whitespace-aligned HMMER output into a DataFrame.
"""
function parse_hmmer_tblout(file_path::AbstractString)
    # Initialize typed vectors instead of a Vector{Any} to avoid "elaborate container types"
    targets = String[]
    queries = String[]
    evalues = Float64[]
    scores = Float64[]
    
    for line in eachline(file_path)
        # Skip comments
        startswith(line, "#") && continue
        
        # Split by whitespace
        parts = split(line)
        
        # Ensure the line has enough parts to avoid bounds errors (basic safety)
        if length(parts) >= 6
            push!(targets, String(parts[1]))
            push!(queries, String(parts[3]))
            push!(evalues, parse(Float64, parts[5]))
            push!(scores, parse(Float64, parts[6]))
        end
    end
    
    # Return a clean DataFrame
    return DataFrame(
        target = targets, 
        query = queries, 
        evalue = evalues, 
        score = scores
    )
end