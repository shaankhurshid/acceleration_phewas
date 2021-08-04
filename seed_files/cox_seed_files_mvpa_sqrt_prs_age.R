# Script to create seed files for phecode cox models

# Depends
library(data.table)
library(stringr)

# Load censor data
censor_data <- fread(file='/Volumes/medpop_afib/skhurshid/phenotypes/2021_06/censor_202106.csv')

# Load pre-seed data
prs <- fread(file='/Volumes/medpop_afib/skhurshid/acceleration_phewas/no_accel_mvpa_sqrt_prs_phenotype.csv')

# Merges
setkey(prs,sample_id); setkey(censor_data,sample_id)
censor_data[prs,':='(inferred_accel = i.inferred_accel_std, PC1 = i.PC1, PC2 = i.PC2, PC3 = i.PC3, PC4 = i.PC4, PC5 = i.PC5)]

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
                                                         ifelse(center_location=='Scotland',pmin(phenotype_censor_date,as.Date('2021-03-31',format='%Y-%m-%d')),
                                                                pmin(phenotype_censor_date,as.Date('2018-02-28',format='%Y-%m-%d')))),origin='1970-01-01'))]

# And set censor date to date of death for those who died
censor_data[,':='(phenotype_censor_date = as.Date(ifelse(!is.na(death_date),pmin(death_date,phenotype_censor_date),phenotype_censor_date),origin='1970-01-01'))]

# Remove people with no PRS
prs <- censor_data[!is.na(inferred_accel)]

### Age subgroups
young <- prs[enroll_age < 55]
young[,inferred_accel := (inferred_accel - mean(inferred_accel,na.rm=T))/sd(inferred_accel,na.rm=T)]

medium <- prs[enroll_age >= 55 & enroll_age < 65]
medium[,inferred_accel := (inferred_accel - mean(inferred_accel,na.rm=T))/sd(inferred_accel,na.rm=T)]

old <- prs[enroll_age >= 65]
old[,inferred_accel := (inferred_accel - mean(inferred_accel,na.rm=T))/sd(inferred_accel,na.rm=T)]

# Scope columns for cutoff BMI
prs_young <- young[,c('sample_id','enroll_age','enroll_date','sex','inferred_accel','phenotype_censor_date',paste0('PC',1:5))]
prs_medium <- medium[,c('sample_id','enroll_age','enroll_date','sex','inferred_accel','phenotype_censor_date',paste0('PC',1:5))]
prs_old <- old[,c('sample_id','enroll_age','enroll_date','sex','inferred_accel','phenotype_censor_date',paste0('PC',1:5))]

# Write out
write.csv(prs_young,file='/Volumes/medpop_afib/skhurshid/acceleration_phewas/cox_data_sqrt_prs_young.csv',row.names = F)
write.csv(prs_medium,file='/Volumes/medpop_afib/skhurshid/acceleration_phewas/cox_data_sqrt_prs_medium.csv',row.names = F)
write.csv(prs_old,file='/Volumes/medpop_afib/skhurshid/acceleration_phewas/cox_data_sqrt_prs_old.csv',row.names = F)