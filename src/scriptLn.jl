
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
    uc_fasta = raw"D:\usearch\centroids.fa",
    out_dir = "D:/hmmerStuff/out/",
    itol_dir = "D:/hmmerStuff/out/itol_assets/"
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
append_taxonomy!(df_headers, db; taxid_col = :TaxID)

final_df = leftjoin(results_df, df_headers[:, Not(:original_header)], on=:accession)

rename!(final_df, :original_header => :id)

CSV.write(joinpath(paths.out_dir, "final_df.csv"), final_df)

#=
    Merge with usearch results to prevent overrepresenting ubiquitous sequences (Lineage bias, cutoff 0.9)
=#

#df_dedup = filter_centroids(final_df, paths.uc_file)
df_dedup = filter_centroids(final_df, paths.uc_file)

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
    get og uniprot annotation
=#

# get earlier data
final_df = CSV.File(joinpath(paths.out_dir, "final_df.csv")) |> DataFrame
final_df.classification = Symbol.(final_df.classification)

# filter for ln archaea
ln = filter(row -> row.classification == :Lanthanide, final_df)
dropmissing!(ln, :kingdom)
ln = filter(row -> row.kingdom == "Archaea", ln)
CSV.write(joinpath(paths.out_dir, "lnArch_uniprot_annotations.csv"), ln)

# annotate with omnicrobe and uniprot data
fetch_omnicrobe_env!(ln; taxid_col=:TaxID)
ln.report_group = agg_habitat.(ln.habitat_specific)
ln.pos2_aa = [string(s[2]) for s in ln.active_site]

export_tree_fasta(ln, paths.fa, joinpath(paths.out_dir, "lnArch.fasta"))

anno_df = fetch_uniprot_annotations(joinpath(paths.out_dir, "lnArch.fasta"))

#=
    get pdf files
=#

# this works now we will apply it to all hits. 
download_alphafold_structures(
    final_df, 
    joinpath(paths.out_dir, "pdb")
)

# Organize AlphaFold models
organize_pdb_files(
    joinpath(paths.out_dir, "pdb"),
    joinpath(paths.out_dir, "organized_pdbs"), 
    final_df
)

#=
    prep foldmason subsample
=#
final_df = CSV.File(joinpath(paths.out_dir, "final_df.csv")) |> DataFrame
final_df.classification = Symbol.(final_df.classification)

input_file = raw"D:\hmmerStuff\out\tree_input.fasta"
output_file = joinpath(paths.out_dir, "foldmason_subsample.fasta")




prepare_foldmason_subsample(
    input_file,
    output_file;
    n_bact = 60,
    n_euk  = 20,
    seed   = 42
)

using JSON

j = JSON.parsefile(raw"D:\hmmerStuff\out\foldmason\foldmason.json")

masonTree = j.tree

writefile = joinpath(paths.out_dir, "foldmason_subsample_tree.nwk")
open(writefile, "w") do io
    write(io, masonTree)
end 