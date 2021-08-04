# Depends
library(data.table)

# Load sample files
genetic <- fread(file='/Volumes/medpop_afib/skhurshid/acceleration_phewas/cox_data_mvpa_sqrt_cutoff_prs.csv')
accel <- fread(file='/Volumes/medpop_afib/skhurshid/acceleration_phewas/cox_data_mvpa_rate_bmi.csv')
self <- fread(file='/Volumes/medpop_afib/skhurshid/acceleration_phewas/cox_data_self_bmi.csv')

# Genetic + accel
ga <- genetic[sample_id %in% accel$sample_id]

# Genetic + self
gs <- genetic[sample_id %in% self$sample_id]

# Accel + self
as <- accel[sample_id %in% self$sample_id]