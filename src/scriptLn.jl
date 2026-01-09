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

f = raw"C:\Users\peete074\OneDrive - Wageningen University & Research\programming\LnFinder\dat\hmmResults\6768f934-391a-4ffa-a8f2-f34d41af7082-5f14c9c6783ff9c0.afa"

df, col_idx = analyze_alignment(f, sample_size=500)

# Filter for Lanthanide
ln_users = filter(row -> row.motif == LnFinder.Lanthanide, df)



