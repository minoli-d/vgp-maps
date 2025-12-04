species="$1"
shapefile="data/shapes/mammals_combined.gpkg"
metadata="data/vgp_iucn_loc_metadata.csv"
sp_sane=$(echo "$species" | tr ' ' '_')
outdir="results/species/${sp_sane}"
goode_template="data/goode_template.tif"

mkdir -p "$outdir"

Rscript scripts/per_species_data.R "$species" "$shapefile" "$metadata" "$outdir" "$goode_template"

# RUN WITH THIS COMMAND :
# cat data/sample_mammals.txt | parallel -j 4 scripts/run_per_species_data.sh