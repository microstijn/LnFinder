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

final_df = merge_taxonomy_with_hits(df, df_hits)

# filter lanthanide specifically. 
ln_users = filter(row -> row.classification == :Lanthanide, final_df)

sum(ismissing.(ln_users.TaxID))

# OUTPUT
CSV.write(output_csv, final_df)





df = CSV.File(output_csv) |> DataFrame


t = df[1:100, :]

fetch_omnicrobe_env!(df; taxid_col=:TaxID)

unique(df.habitat_specific)
unique(df.habitat_category)

dropmissing!(df)
ln_users = filter(row -> row.classification == "Lanthanide", df)

ln_users[:, :habitat_specific]




#figure gene

# 2. Get your results
results_df = parse_hmmalign_results(afa, idx)

# 3. Extract the 3-residue motifs for BioSequences visualization
# We only take the :Lanthanide ones for the logo
lanthanide_hits = filter(row -> row.classification == :Lanthanide, results_df)

# Create LongAminoAcidSeqs for Logo generation
# We replace the "-x-" string with the actual 3 characters
