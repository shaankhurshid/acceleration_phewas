# Script to create seed files for phecode cox models

# Depends
library(data.table)
library(stringr)

# Load censor data
censor_data <- fread(file='/Volumes/medpop_afib/skhurshid/phenotypes/2021_06/censor_202106.csv')

# Load self data
self <- fread(file='/Volumes/medpop_afib/skhurshid/acceleration_phewas/self_phenotype.csv')

### BMI
bmi <- fread(file='/Volumes/medpop_afib/skhurshid/phenotypes/homegrown/bmi_instance0_202109.csv')

## Smoking
tob <- fread('/Volumes/medpop_afib/skhurshid/phenotypes/tobacco_instance_0.csv')
tob[,tobacco_accel_selfreport := ifelse(value==2,"Current",
                                        ifelse(value==1,"Former","Never"))]

## TDI
tdi <- fread(file='/Volumes/medpop_afib/skhurshid/phenotypes/2020_06/townsend_0.csv')

## EtOH
etoh <- fread(file='/Volumes/medpop_afib/skhurshid/phenotypes/homegrown/etoh_0_012022.csv')

## BP
sbp <- fread('/Volumes/medpop_afib/skhurshid/phenotypes/homegrown/sbp_combined_instance0_202006.csv')
dbp <- fread('/Volumes/medpop_afib/skhurshid/phenotypes/homegrown/dbp_combined_instance0_202006.csv')

## anti-HTN
bpmed <- fread('/Volumes/medpop_afib/skhurshid/phenotypes/homegrown/bpmed_combined_instance0.csv')

## Education
edu <- fread('/Volumes/medpop_afib/skhurshid/phenotypes/homegrown/education_0.csv')

## Diet
diet <- fread('/Volumes/medpop_afib/skhurshid/phenotypes/homegrown/diet_instance0.csv')

# Load withdrawals
withdrawals <- fread(file='/Volumes/medpop_afib/skhurshid/phenotypes/withdrawals/w7089_20220222.csv')

# Merges
setkey(self,sample_id); setkey(censor_data,sample_id); setkey(bmi,sample_id)
censor_data[self,':='(self_mvpa = i.self_mvpa, self_guidelines = i.self_guidelines)]
censor_data[bmi,bmi := i.bmi]

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

# Remove missing exposure data
censor_data <- censor_data[!is.na(self_mvpa) & !is.na(self_guidelines)] #502459 - 42597 = 459862

# Merges
setkey(censor_data, sample_id); setkey(bpmed,sample_id); setkey(tob,sample_id)
setkey(tdi,sample_id); setkey(dbp,sample_id); setkey(sbp,sample_id); setkey(etoh,sample_id)
setkey(diet,sample_id); setkey(edu,sample_id)
censor_data[bmi, bmi := i.bmi]
censor_data[bpmed, bpmed := i.prev_bpmed]
censor_data[tob, tob := i.tobacco_accel_selfreport]
censor_data[tdi, tdi := i.value]
censor_data[dbp, dbp := i.dbp0]
censor_data[sbp, sbp := i.sbp0]
censor_data[diet, diet := i.diet_quality]
censor_data[edu, qual_ea := i.qual_ea]
censor_data[etoh, ':='(etoh_grams = i.etoh_grams, etoh_status = i.f_1558)]
censor_data[is.na(bpmed)]$bpmed <- 0
censor_data[is.na(qual_ea)]$qual_ea <- median(censor_data$qual_ea,na.rm=T)
censor_data[is.na(diet)]$diet <- 'intermediate' # 105 or 0.1%

# Missing data removal
censor_data <- censor_data[!is.na(bmi)] # 459862 - 2228 = 457634
censor_data <- censor_data[!is.na(tdi)] # 457634 - 564 = 457070
censor_data <- censor_data[!is.na(sbp)] # 457070 - 383 = 456687
censor_data <- censor_data[!is.na(dbp)] # 456687 - 0 = 456687

# Impute for alcohol 
censor_data <- censor_data[!is.na(etoh_status)] # 456687 - 244 = 456443

# 61871 + 19 to impute
## Fix 1 outlier
censor_data <- censor_data[,etoh_grams := ifelse(etoh_grams > 2000, NA, etoh_grams)]

monthly_mean <- round(median(censor_data[c(!is.na(etoh_grams) & etoh_status=='Monthly')]$etoh_grams),0)
weekly_mean <- round(median(censor_data[c(!is.na(etoh_grams) & etoh_status=='Weekly')]$etoh_grams),0)

censor_data[c(etoh_status=='Weekly' & is.na(etoh_grams)),
             etoh_grams := (weekly_mean)]
censor_data[c(etoh_status=='Monthly' & is.na(etoh_grams)),
             etoh_grams := (monthly_mean)]

# Remove withdrawals
censor_data <- censor_data[!(sample_id %in% withdrawals$V1)] #456443 - 69 = 456374

# Scope columns for PRS analysis
self <- censor_data[,c('sample_id','enroll_age','enroll_date','sex','self_mvpa','self_guidelines',
                       'bmi','phenotype_censor_date',
                      'sbp','dbp','bpmed','tob','tdi','etoh_grams','diet','qual_ea')]

# Complete covariates
# Final N=456374

# Write out raw
write.csv(self,file='/Volumes/medpop_afib/skhurshid/acceleration_phewas/cox_data_self_bmi_raw_covar.csv',row.names = F)

# Standardize self MVPA
#self[,':='(self_mvpa = (self_mvpa - mean(self_mvpa,na.rm=T))/sd(self_mvpa,na.rm=T))]
self[,':='(self_mvpa = self_mvpa/150)]

# Write out
write.csv(self,file='/Volumes/medpop_afib/skhurshid/acceleration_phewas/cox_data_self_bmi_covar.csv',row.names = F)
