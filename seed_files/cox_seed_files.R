# Script to create seed files for phecode cox models

# Depends
library(data.table)
library(stringr)
library(sktools)

# Load BMI confounder data
bmi <- fread(file='/Volumes/medpop_afib/skhurshid/acceleration_phewas/bmi_confounder.csv')
setnames(bmi,'x21001_0_0','bmi')

# Load censor data
censor_data <- fread(file='/Volumes/medpop_afib/skhurshid/phenotypes/2021_06/censor_202106.csv')
old_censor_data <- fread(file='/Volumes/medpop_afib/skhurshid/phenotypes/2020_06/censor_202006.csv')

# Load accel data
load(file='/Volumes/medpop_afib/skhurshid/acceleration_phewas/accel_phewas_110820.RData')

# Load fixed MVPA and convert to bouts
bouts <- list()
for (i in 1:21){
  load(file=paste0('/Volumes/medpop_afib/skhurshid/accel/bouts_v2/out_bouts',i,'.RData'))
  bouts[[i]] <- mvpa_5
}
file_list <- fread('/Volumes/medpop_afib/skhurshid/accel/accel_file_list.csv')
ids <- as.numeric(str_extract(file_list$x,'\\d+'))
bouts <- data.table(do.call(rbind,bouts))
bouts <- cbind(bouts,ids)
bouts <- bouts[,id := NULL]
setnames(bouts,'ids','sample_id')
bouts[,sample_id := as.numeric(sample_id)]

# Load the key to translate into sek ID
linker <- fread('/Volumes/medpop_afib/skhurshid/accel/7089_17488_linker.csv')

# Joins
setkey(linker,app17488); setkey(bouts,sample_id)
bouts[linker,sek_id := i.app7089]

setkey(bouts,sek_id); setkey(acceleration,sample_id)
acceleration[bouts,mvpa_bouts := (i.mvpa)*60*5]

# Create guideline-based cutoffs
acceleration[,':='(cutoff300 = ifelse((mvpa_bouts/60)/(total/86400) >= 300/7,1,0),
                  cutoff75 = ifelse((mvpa_bouts/60)/(total/86400) >= 75/7,1,0),
                  who_acc_rate_bouted = ifelse((mvpa_bouts/60)/(total/86400) >= 21.42857,1,0),
                  mvpa_rate_bouted = (mvpa_bouts/60)/(total/86400))]

#### Fix censor data in censor file
# Load center categories
center <- fread(file='/Volumes/medpop_afib/skhurshid/phenotypes/center0.csv')
center_lookup <- fread(file='/Volumes/medpop_afib/skhurshid/phenotypes/enrollment_correspondences.csv')

# Add center value to dataset
setkey(censor_data,sample_id); setkey(old_censor_data,sample_id); setkey(center,sample_id)
censor_data[center,':='(center_code = i.value)]
old_censor_data[center,':='(center_code = i.value)]

setkey(censor_data,center_code); setkey(old_censor_data,center_code); setkey(center_lookup,Code)
censor_data[center_lookup,':='(center_location = i.Region)]
old_censor_data[center_lookup,':='(center_location = i.Region)]

# Now correct censor dates based on location
censor_data[,':='(phenotype_censor_date = as.Date(ifelse(center_location=='England',phenotype_censor_date,
                                                         ifelse(center_location=='Scotland',pmin(phenotype_censor_date,as.Date('2021-03-31',format='%Y-%m-%d')),
                                                                pmin(phenotype_censor_date,as.Date('2018-02-28',format='%Y-%m-%d')))),origin='1970-01-01'))]
old_censor_data[,':='(phenotype_censor_date = as.Date(ifelse(center_location=='England',phenotype_censor_date,
                                                         ifelse(center_location=='Scotland',pmin(phenotype_censor_date,as.Date('2016-10-31',format='%Y-%m-%d')),
                                                                pmin(phenotype_censor_date,as.Date('2016-02-29',format='%Y-%m-%d')))),origin='1970-01-01'))]

# And set censor date to date of death for those who died
censor_data[,':='(phenotype_censor_date = as.Date(ifelse(!is.na(death_date),pmin(death_date,phenotype_censor_date),phenotype_censor_date),origin='1970-01-01'))]
old_censor_data[,':='(phenotype_censor_date = as.Date(ifelse(!is.na(death_date),pmin(death_date,phenotype_censor_date),phenotype_censor_date),origin='1970-01-01'))]

# Merges
setkey(bmi,userId); setkey(censor_data,sample_id)
acceleration[bmi, bmi := i.bmi]
acceleration[censor_data,':='(enroll_date = i.enroll_date, phenotype_censor_date = i.phenotype_censor_date)]
setnames(acceleration,'end_date','accel_date')
acceleration[,value_std := (value - mean(value))/sd(value)]
acceleration[,mvpa_std := (mvpa_rate_bouted - mean(mvpa_rate_bouted,na.rm=T))/sd(mvpa_rate_bouted,na.rm=T)]

# Fix dates
for (j in (c('enroll_date','phenotype_censor_date'))){set(acceleration,j=j,value=as.Date(acceleration[[j]],format='%Y-%m-%d'))}

# Remove individuals with missing BMI
acceleration <- acceleration[!is.na(bmi)]

# Remove individuals in the latest withdrawal file
withdrawals <- fread('/Volumes/medpop_afib/skhurshid/phenotypes/withdrawals/w7089_20210201.csv') # UKBB withdrawals
acceleration <- acceleration[!(sample_id %in% withdrawals$V1)]

# Save out phenotype file
write.csv(acceleration,file='/Volumes/medpop_afib/skhurshid/acceleration_phewas/acceleration_phenotype.csv',row.names = F)

### OUTPUT 1: VALUE BMI

# Scope columns for value BMI
value_bmi <- acceleration[,c('sample_id','accel_date','age_accel','sex','bmi','value_std','phenotype_censor_date')]

# Write out
write.csv(value_bmi,file='/Volumes/medpop_afib/skhurshid/acceleration_phewas/cox_data_value_bmi.csv',row.names = F)

### OUTPUT 2: INVNORM VALUE BMI

# Scope columns for value BMI invnorm
value_bmi_invnorm <- acceleration[,c('sample_id','accel_date','age_accel','sex','bmi','value_std','phenotype_censor_date')]

# Inverse normal transformation
value_bmi_invnorm[,value_invnorm := qnorm((rank(value_std,na.last="keep")-0.5)/sum(!is.na(value_std)))]

# Write out
write.csv(value_bmi_invnorm,file='/Volumes/medpop_afib/skhurshid/acceleration_phewas/cox_data_value_bmi_invnorm.csv',row.names = F)

### OUTPUT 3: DECILE VALUE BMI
# Scope columns for value BMI invnorm
value_bmi_decile <- acceleration[,c('sample_id','accel_date','age_accel','sex','bmi','value_std','phenotype_censor_date')]

# Inverse normal transformation
value_bmi_decile[,value_decile := quantilize(value_std,ncuts=10)]

# Write out
write.csv(value_bmi_decile,file='/Volumes/medpop_afib/skhurshid/acceleration_phewas/cox_data_value_bmi_decile.csv',row.names = F)

### OUTPUT 4: MVPA RATE BMI

# Scope columns for value BMI
mvpa_bmi <- acceleration[,c('sample_id','accel_date','age_accel','sex','bmi','mvpa_std','phenotype_censor_date')]

# Write out
write.csv(mvpa_bmi,file='/Volumes/medpop_afib/skhurshid/acceleration_phewas/cox_data_mvpa_rate_bmi.csv',row.names = F)

### OUTPUT 5: MVPA RATE CUTOFF BMI
# Scope columns for cutoff BMI
mvpa_bmi <- acceleration[,c('sample_id','accel_date','age_accel','sex','bmi','who_acc_rate_bouted','phenotype_censor_date')]

# Write out
write.csv(mvpa_bmi,file='/Volumes/medpop_afib/skhurshid/acceleration_phewas/cox_data_mvpa_cutoff_bmi.csv',row.names = F)

### OUTPUT 6: LOW CUTOFF BMI
# Scope columns for cutoff BMI
mvpa_bmi <- acceleration[,c('sample_id','accel_date','age_accel','sex','bmi','cutoff75','phenotype_censor_date')]

# Write out
write.csv(mvpa_bmi,file='/Volumes/medpop_afib/skhurshid/acceleration_phewas/cox_data_mvpa_cutoff75_bmi.csv',row.names = F)

### OUTPUT 7: HIGH CUTOFF BMI
# Scope columns for cutoff BMI
mvpa_bmi <- acceleration[,c('sample_id','accel_date','age_accel','sex','bmi','cutoff300','phenotype_censor_date')]

# Write out
write.csv(mvpa_bmi,file='/Volumes/medpop_afib/skhurshid/acceleration_phewas/cox_data_mvpa_cutoff300_bmi.csv',row.names = F)

### OUTPUT 8: AGE SUBGROUPS
young <- acceleration[age_accel < 55]
young[,mvpa_std := (mvpa_rate_bouted - mean(mvpa_rate_bouted,na.rm=T))/sd(mvpa_rate_bouted,na.rm=T)]

medium <- acceleration[age_accel >= 55 & age_accel < 65]
medium[,mvpa_std := (mvpa_rate_bouted - mean(mvpa_rate_bouted,na.rm=T))/sd(mvpa_rate_bouted,na.rm=T)]

old <- acceleration[age_accel >= 65]
old[,mvpa_std := (mvpa_rate_bouted - mean(mvpa_rate_bouted,na.rm=T))/sd(mvpa_rate_bouted,na.rm=T)]

# Scope columns for cutoff BMI
mvpa_bmi_young <- young[,c('sample_id','accel_date','age_accel','sex','bmi','mvpa_std','phenotype_censor_date')]
mvpa_bmi_medium <- medium[,c('sample_id','accel_date','age_accel','sex','bmi','mvpa_std','phenotype_censor_date')]
mvpa_bmi_old <- old[,c('sample_id','accel_date','age_accel','sex','bmi','mvpa_std','phenotype_censor_date')]

# Write out
write.csv(mvpa_bmi_young,file='/Volumes/medpop_afib/skhurshid/acceleration_phewas/cox_data_mvpa_bmi_young.csv',row.names = F)
write.csv(mvpa_bmi_medium,file='/Volumes/medpop_afib/skhurshid/acceleration_phewas/cox_data_mvpa_bmi_medium.csv',row.names = F)
write.csv(mvpa_bmi_old,file='/Volumes/medpop_afib/skhurshid/acceleration_phewas/cox_data_mvpa_bmi_old.csv',row.names = F)


