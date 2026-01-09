module LnFinder

using BioSequences
using FASTX
using DataFrames
using CSV

# Include sub-modules
#include("ModuleHmmer.jl")
include("ModuleMotifSearch.jl")
include("ModuleIO.jl")

# Export public interface methods
#export run_hmmer, parse_hmmer_tblout
export classify_pqq_adh, MotifType
export analyze_alignment
#export classify_hits_from_fasta

end # module LnFinder
