#= 
Preamble & environment management
=# 

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

#= 
Packages
=# 
using LnFinder

#=
prelim scripting
=# 

f = raw"C:\Users\peete074\OneDrive - Wageningen University & Research\programming\LnFinder\dat\bc0f8812-ae10-4feb-ae88-69b3e5628280-5f14c9c6783ff9c0.fa"

df = classify_hits_from_fasta(f)

# 3. Filter for Lanthanide users (The "Where To Find Them" part)
ln_users = filter(row -> row.motif == LnFinder.Lanthanide, df)