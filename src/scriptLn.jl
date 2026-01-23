#= 
Preamble & environment management
=# 

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

#= 
Packages
=# 

using Revise
using LnFinder
using CSV
using DataFrames 

#=
prelim scripting
=# 

#=
setup
=# 

afa = raw"D:\hmmerStuff\data\aligned_to_model.afa"
fa = raw"D:\hmmerStuff\data\firstTry.fa"
taxdump_dir = raw"D:\ncbi_downloads\taxdump"
output_csv = "D:/hmmerStuff/out/lanthipeptide_lineages.csv"

#=
active site analysis
=# 

idx = calibrate_active_site(afa)
df = parse_hmmalign_results(afa, idx)
ln_users = filter(row -> row.classification == :Lanthanide, df)
stats_df = count_residues_per_position(df.active_site)

#=
taxonomy
=# 

df_hits = parse_fasta_headers(fa)

# fetch. Real slow. 
acc_map = fetch_taxids_async(df_hits.accession)

# Map the TaxIDs back to the DataFrame
# If an accession wasn't found (e.g., deleted), we fill 0 or missing
df_hits.TaxID = [get(acc_map, acc, missing) for acc in df_hits.accession]

# get taxonomy from database
db = load_taxonomy_db(taxdump_dir)

# apply local db to taxonomy 
append_taxonomy!(df_hits, db; taxid_col=:TaxID)


#=
 merge results
=# 

final_df = leftjoin(df, df_hits, on = :id => :original_header)
# filter lanthanide specifically. 
ln_users = filter(row -> row.classification == :Lanthanide, final_df)

sum(ismissing.(ln_users.TaxID))

# OUTPUT
CSV.write(output_csv, final_df)


df = CSV.File(output_csv) |> DataFrame

df.classification = Symbol.(df.classification )

#df = df[1:1000, :]

fetch_omnicrobe_env!(df; taxid_col=:TaxID)

CSV.write(output_csv, df)

df = CSV.File(output_csv) |> DataFrame
df.classification = Symbol.(df.classification )
df.report_group = strip.(agg_habitat.(df.habitat_specific))





# Prepare a sub-dataframe for statistics
sub = copy(df)
# Ensure Pos 2 is a string to match our target list
sub.pos2_aa = [string(s[2]) for s in sub.active_site]

target_aas = Set(["Y", "F", "W", "A", "G"])
target_groups = Set(["Aquatic", "Terrestrial", "Host-associated", "Anthropogenic"])

# Filter for valid ecological pillars and target residues
sub = filter(row -> row.report_group in target_groups && 
                        row.pos2_aa in target_aas, sub)





unique(df.habitat_specific)
unique(df.habitat_category)

dropmissing!(df)
ln_users = filter(row -> row.classification == "Lanthanide", df)

ln_users[:, :habitat_specific]



df.TaxID
#figure gene

# 2. Get your results
results_df = parse_hmmalign_results(afa, idx)

# 3. Extract the 3-residue motifs for BioSequences visualization
# We only take the :Lanthanide ones for the logo
lanthanide_hits = filter(row -> row.classification == :Lanthanide, results_df)

# Create LongAminoAcidSeqs for Logo generation
# We replace the "-x-" string with the actual 3 characters



function rerun_niche_stats(df)
    # 1. Clean data for stats
    # Focus only on target Aromatics and valid habitat pillars
    target_aas = Set(["Y", "F", "W", "A", "G"])
    target_groups = Set(["Aquatic", "Terrestrial", "Host-associated", "Anthropogenic"])
    
    # Extract Pos 2 AA from the active_site string
    df.pos2_aa = [length(s) >= 2 ? string(s[2]) : "-" for s in df.active_site]
    
    sub = filter(row -> row.report_group in target_groups && 
                        row.pos2_aa in target_aas, df)

    # 2. Build Contingency Table
    ct = combine(groupby(sub, [:report_group, :pos2_aa]), nrow => :Count)
    matrix_df = unstack(ct, :pos2_aa, :report_group, :Count, fill=0)
    
    # 3. Perform Chi-Squared Test
    obs = Matrix{Int}(matrix_df[:, 2:end])
    chisq = ChiSquaredTest(obs)
    pval = pvalue(chisq)
    
    # 4. Calculate Z-Scores
    row_sums = sum(obs, dims=2)
    col_sums = sum(obs, dims=1)
    total = sum(obs)
    
    println("\n📈 RECALCULATED STATISTICS (Deduplicated)")
    println("Chi-Squared p-value: $pval")
    println("-"^40)
    
    # Print a quick text-based heatmap
    headers = names(matrix_df)[2:end]
    print(rpad("AA", 6))
    for h in headers; print(rpad(h, 15)); end
    println()
    
    for i in 1:size(obs, 1)
        print(rpad("**$(matrix_df[i, 1])**", 6))
        for j in 1:size(obs, 2)
            expected = (row_sums[i] * col_sums[j]) / total
            z = (obs[i,j] - expected) / sqrt(expected)
            
            # Format output with indicators
            indicator = z > 2.0 ? "🔥" : (z < -2.0 ? "❄️" : "  ")
            print(rpad("$indicator $(round(z, digits=1))", 15))
        end
        println()
    end
    
    return pval, matrix_df
end

# Usage:
pval, results_matrix = rerun_niche_stats(df_dedup)


df_dedup2 = dropmissing(df_dedup)
println(df_dedup2[df_dedup2.kingdom .== "Bacteria", :species])

"""
    summarize_genera(df)
Maakt een top-10 lijst van de meest voorkomende bacteriële genera.
"""
function summarize_genera(df)
    # Filter alleen bacteriën uit de Lanthanide groep
    bac_df = filter(row -> row.kingdom == "Bacteria" && row.classification == :Lanthanide, df)
    
    # Extraheer het eerste woord van de species naam (het Genus)  
    genus_summary = combine(groupby(bac_df, [:genus, :kingdom]), nrow => :Count)
    sort!(genus_summary, :Count, rev=true)
    
    println("🏆 Top 10 Lanthanide-Dependent Genera:")
    display(first(genus_summary, 40))
    
    return genus_summary
end

summarize_genera(df_dedup2)

#= 
Preamble & environment management
=# 

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

#=
    Packages 
=#

using LnFinder
using CSV
using DataFrames
using Serialization

#=
    Config 
=#

paths = (
    afa = raw"D:\hmmerStuff\data\aligned_to_model.afa",
    fa = raw"D:\hmmerStuff\data\firstTry.fa",
    taxdump = raw"D:\ncbi_downloads\taxdump",
    uc_file = "D:/usearch/clusters.uc",
    uc50_file = "D:/usearch/clusters_50.uc",
    out_dir = "D:/hmmerStuff/out/"
) 

#=
    Active Site Calibration and Parsing 
=#

# Identifies the model index for the D-x-D motif 
idx = calibrate_active_site(paths.afa)
results_df = parse_hmmalign_results(paths.afa, idx)


#=
    taxnomy mapping
    Extract accessions and fetch TaxIDs via UniProt 
=#

df_headers = parse_fasta_headers(paths.fa)
acc_map = fetch_taxids_async(df_headers.accession)
df_headers.TaxID = [get(acc_map, acc, missing) for acc in df_headers.accession]

# Load local NCBI database and append lineages 
db = load_taxonomy_db(paths.taxdump)
append_taxonomy!(df_headers, db; taxid_col=:TaxID)
results_df
final_df = leftjoin(results_df, df_headers[:, Not(:original_header)], on=:accession)

rename!(final_df, :original_header => :id)
#=
    Merge with usearch results to prevent overrepresenting ubiquitous sequences (Lineage bias, cutoff 0.9)
=#

#df_dedup = filter_centroids(final_df, paths.uc_file)
df_dedup = filter_centroids(final_df, paths.uc50_file)

#=
    Fetch OmniCrobe env ontologies and merge 
=#

fetch_omnicrobe_env!(df_dedup; taxid_col=:TaxID)

# Heuristic rollup into ecological pillars 
df_dedup.report_group = agg_habitat.(df_dedup.habitat_specific)

#=
    Statistics ecological / tax vs. motif second AA.  
    and a small export for a tree
=#

# Calculate residue frequencies 

stats_df = count_residues_per_position(df_dedup.active_site)

# Add column for only pos 2
df_dedup.pos2_aa = [string(s[2]) for s in df_dedup.active_site]

# tree fa export
export_tree_fasta(df_dedup, paths.fa, joinpath(paths.out_dir, "tree_input_deduped.fasta"))

# drop missing for furtyer analysis. 
dropmissing!(df_dedup, :kingdom)

# Arc | Bacteria vs AA
ln_users = filter(row -> row.classification == :Lanthanide, df_dedup)

arch_vs_bac_df, p_val_arch = compare_archaea_bacteria(ln_users)

# Env vs AA
matrix_df, obs, z_scores, pval = calculate_enrichment_stats(ln_users)

# convert to z-score df
z_df = copy(matrix_df)
for (i, colname) in enumerate(names(z_df)[2:end])
    z_df[!, colname] = round.(z_scores[:, i], digits=2)
end


CSV.write(joinpath(paths.out_dir, "archaea_comparison.csv"), arch_vs_bac_df)
CSV.write(joinpath(paths.out_dir, "niche_enrichment.csv"), matrix_df)
CSV.write(joinpath(paths.out_dir, "niche_z_scores.csv"), z_df)
# Final Data Output
CSV.write(joinpath(paths.out_dir, "final_deduplicated_data.csv"), df_dedup)
CSV.write(joinpath(paths.out_dir, "motif_conservation_stats.csv"), stats_df)
CSV.write(joinpath(paths.out_dir, "niche_enrichment_matrix.csv"), matrix_df)

# Serialize the raw statistical objects for later visualization in Julia/R
serialize(joinpath(paths.out_dir, "statistical_results.jls"), 
    Dict(:p_value => p_val, :obs => obs_matrix, :total => grand_total))


#=
    Create a FastTree 
=#

f = readdir(raw"D:\teamsDownoad\Monthly_values")

println(f)