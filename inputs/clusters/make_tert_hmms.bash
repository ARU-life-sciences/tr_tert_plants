for file in ./*as_nucl_aligned.fasta; do
  base=$(basename "$file" .fasta)
  /software/team301/hmmer-3.4/src/hmmbuild --dna "${base}.hmm" "$file"
  echo "Built HMM: ${base}.hmm"
done

 mv plant_TERT_cdhit.faa.telomerase_reverse_transcriptase_[Oryza_sativa_Japonica_Group]_as_nucl_aligned.hmm tert1.hmm
 mv plant_TERT_cdhit.faa.telomerase_reverse_transcriptase_catalytic_subunit_[Zea_mays]_as_nucl_aligned.hmm catalytic_subunit_tert.hmm
 mv plant_TERT_cdhit.faa.telomerase_reverse_transcriptase_isoform_X1_[Nicotiana_sylvestris]_as_nucl_aligned.hmm isoform_X1_tert.hmm
