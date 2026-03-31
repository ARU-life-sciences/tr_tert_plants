/software/team301/cdhit/cd-hit -i plant_TERT_regions.faa -o plant_TERT_cdhit.faa -c 0.9 -n 5

# requires clstr
clstr tofasta plant_TERT_cdhit.faa.clstr ./plant_TERT_regions.faa

mkdir clusters
mv plant_TERT_cdhit.* clusters
