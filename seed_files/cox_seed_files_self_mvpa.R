# Script to create seed files for phecode cox models

# Depends
library(data.table)
library(stringr)

# Load censor data
censor_data <- fread(file='/Volumes/medpop_afib/skhurshid/phenotypes/2021_06/censor_202106.csv')

# Load self data
self <- fread(file='/Volumes/medpop_afib/skhurshid/acceleration_phewas/self_phenotype.csv')

# Load BMI data
bmi <- fread(file='/Volumes/medpop_afib/skhurshid/acceleration_phewas/self_bmi_confounder.csv')

# Load withdrawals
withdrawals <- fread(file='/Volumes/medpop_afib/skhurshid/phenotypes/withdrawals/w7089_20210201.csv')

# Merges
setkey(self,sample_id); setkey(censor_data,sample_id); setkey(bmi,userId)
censor_data[self,':='(self_mvpa = i.self_mvpa, self_guidelines = i.self_guidelines)]
censor_data[bmi,bmi := i.x21001_0_0]

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

# Remove withdrawals
censor_data <- censor_data[!(sample_id %in% withdrawals$V1)]

# Scope columns for PRS analysis
self <- censor_data[c(!is.na(self_mvpa) & !is.na(self_guidelines)),
                    c('sample_id','enroll_age','enroll_date','sex','self_mvpa','self_guidelines',
                       'bmi','phenotype_censor_date')] 

# Complete covariates
self <- self[!is.na(self_mvpa) & !is.na(self_guidelines) & !is.na(enroll_age) & !is.na(bmi)]
# Final N=457633

# Standardize self MVPA
self[,':='(self_mvpa = (self_mvpa - mean(self_mvpa,na.rm=T))/sd(self_mvpa,na.rm=T))]

# Write out
write.csv(self,file='/Volumes/medpop_afib/skhurshid/acceleration_phewas/cox_data_self_bmi.csv',row.names = F)
