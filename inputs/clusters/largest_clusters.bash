mmft num *.fasta \
  | sort -k2,2nr \
  | tail -n+4 \
  | cut -f1 \
  | xargs --no-run-if-empty rm

for file in ./plant_TERT_cdhit.faa.*.fasta; do
  nm=$(basename $file .fasta)
  echo $nm
  aa2nucaln $file > "${nm}_as_nucl.fasta"
done

for file in ./*_as_nucl.fasta; do
  echo $file
  /software/team301/mafft-7.525-with-extensions/core/mafft $file > "${file%.fasta}_aligned.fasta"
done
