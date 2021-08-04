# Script to create seed files for phecode cox models

# Depends
library(data.table)
library(stringr)
library(sktools)

# Load censor data
censor_data <- fread(file='/Volumes/medpop_afib/skhurshid/phenotypes/2020_06/censor_202006.csv')

# Load pre-seed data
prs <- fread(file='/Volumes/medpop_afib/skhurshid/acceleration_phewas/no_accel_prs_phenotype.csv')

# Merges
setkey(prs,sample_id); setkey(censor_data,sample_id)
censor_data[prs,':='(inferred_accel = i.inferred_accel_std, PC1 = i.PC1, PC2 = i.PC2, PC3 = i.PC3, PC4 = i.PC4, PC5 = i.PC5)]

# Fix dates
for (j in (c('enroll_date','phenotype_censor_date'))){set(acceleration,j=j,value=as.Date(acceleration[[j]],format='%Y-%m-%d'))}

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

# Remove people with no PRS
prs <- censor_data[!is.na(inferred_accel)]

# Decile transformation
prs[,inferred_accel_decile := quantilize(inferred_accel,ncuts=10)]

# Scope columns for PRS analysis
prs <- prs[,c('sample_id','enroll_age','enroll_date','sex','inferred_accel_decile','phenotype_censor_date',paste0('PC',1:5))]

# Write out
write.csv(prs,file='/Volumes/medpop_afib/skhurshid/acceleration_phewas/cox_data_prs_decile.csv',row.names = F)
