# Script to create seed files for phecode cox models

# Depends
library(data.table)
library(stringr)

# Load censor data
censor_data <- fread(file='/Volumes/medpop_afib/skhurshid/phenotypes/2020_06/censor_202006.csv')

# Load PRS data
prs <- fread(file='/Volumes/medpop_afib/skhurshid/acceleration_phewas/prs/scores/prs_no_accel_white/exposure_prs_no_accel_white.csv')

# Load PC data
pcs <- fread(file='/Volumes/medpop_afib/skhurshid/phenotypes/genetic/pcs15.csv')

# Merges
setkey(prs,userId); setkey(censor_data,sample_id); setkey(pcs,sample_id)
censor_data[prs,prs_std := i.exposure]
censor_data <- censor_data[pcs,nomatch=0]

#### Fix censor data in censor file
# Load center categories
center <- fread(file='/Volumes/medpop_afib/skhurshid/phenotypes/center0.csv')
center_lookup <- fread(file='/Volumes/medpop_afib/skhurshid/phenotypes/enrollment_correspondences.csv')

# Add center value to dataset
setkey(censor_data,sample_id); setkey(center,sample_id)
censor_data[center,':='(center_code = i.value)]

setkey(censor_data,center_code); setkey(center_lookup,Code)
censor_data[center_lookup,':='(center_location = i.Region)]

# Now correct censor dates based on location
censor_data[,':='(phenotype_censor_date = as.Date(ifelse(center_location=='England',phenotype_censor_date,
                                                         ifelse(center_location=='Scotland',pmin(phenotype_censor_date,as.Date('2016-10-31',format='%Y-%m-%d')),
                                                                pmin(phenotype_censor_date,as.Date('2016-02-29',format='%Y-%m-%d')))),origin='1970-01-01'))]

# And set censor date to date of death for those who died
censor_data[,':='(phenotype_censor_date = as.Date(ifelse(!is.na(death_date),pmin(death_date,phenotype_censor_date),phenotype_censor_date),origin='1970-01-01'))]

# Remove people with no PRS
prs <- censor_data[!is.na(prs_std)]

# Scope columns for PRS analysis
prs <- prs[,c('sample_id','enroll_age','enroll_date','sex','prs_std','phenotype_censor_date',paste0('pc',1:5))]

# Write out
write.csv(prs,file='/Volumes/medpop_afib/skhurshid/acceleration_phewas/cox_data_prs_white.csv',row.names = F)
