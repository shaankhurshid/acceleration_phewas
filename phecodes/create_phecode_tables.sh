for file in /mnt/ml4cvd/projects/skhurshid/accel_phewas/new_phecode_tabs3/*
do
  ./ukbb2disease.linux --project ukbb-analyses --usegp=true -database "ukbb-analyses.ukbb7089_202106" -tabfile "$file" > "$file".tsv
done