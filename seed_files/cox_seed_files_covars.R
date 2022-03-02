# Script to create seed files for phecode cox models

# Depends
library(data.table)
library(stringr)
library(sktools)

# Load confounder data
### BMI
bmi <- fread(file='bmi_instance01_202109.csv')

## Smoking
tob0 <- fread('tobacco_instance_0.csv')
tob1 <- fread('tobacco_instance1.csv')
setkey(tob0,sample_id); setkey(tob1,sample_id)
tob <- tob0[tob1,':='(value1 = i.value)]
tob[,value_unified := ifelse(!is.na(value1),value1,value)]
tob[,tobacco_accel_selfreport := ifelse(value_unified==2,"Current",
                                        ifelse(value_unified==1,"Former","Never"))]

## TDI
tdi <- fread(file='2020_06/townsend_0.csv')

## EtOH
etoh <- fread(file='etoh_01_012022.csv')

## BP
sbp <- fread('sbp_combined_instance01_202006.csv')
dbp <- fread('dbp_combined_instance01_202006.csv')

## anti-HTN
bpmed <- fread('bpmed_combined_instance01.csv')

# Load censor data
censor_data <- fread(file='censor_202106.csv')

# Load accel data
load(file='accel_phewas_110820.RData')

# Load fixed MVPA and convert to bouts
bouts <- list()
for (i in 1:21){
  load(file=paste0('/out_bouts',i,'.RData'))
  bouts[[i]] <- mvpa_5
}
file_list <- fread('accel_file_list.csv')
ids <- as.numeric(str_extract(file_list$x,'\\d+'))
bouts <- data.table(do.call(rbind,bouts))
bouts <- cbind(bouts,ids)
bouts <- bouts[,id := NULL]
setnames(bouts,'ids','sample_id')
bouts[,sample_id := as.numeric(sample_id)]

# Load the key to translate into sek ID
linker <- fread('7089_17488_linker.csv')

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
center <- fread(file='center0.csv')
center_lookup <- fread(file='enrollment_correspondences.csv')

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

# Merges
setkey(bmi,sample_id); setkey(censor_data,sample_id); setkey(bpmed,sample_id); setkey(tob,sample_id)
setkey(tdi,sample_id); setkey(dbp,sample_id); setkey(sbp,sample_id); setkey(etoh,sample_id)
setnames(acceleration,'end_date','accel_date')
acceleration[bmi, bmi := i.bmi]
acceleration[bpmed, bpmed := i.prev_bpmed]
acceleration[tob, tob := i.tobacco_accel_selfreport]
acceleration[tdi, tdi := i.value]
acceleration[dbp, dbp := i.dbp_combined]
acceleration[sbp, sbp := i.sbp_combined]
acceleration[etoh, ':='(etoh_grams = i.etoh_grams, etoh_status = i.f_1558)]
acceleration[censor_data,':='(enroll_date = i.enroll_date, phenotype_censor_date = i.phenotype_censor_date)]
acceleration[,value_std := (value - mean(value))/sd(value)]
acceleration[,mvpa_std := (mvpa_rate_bouted - mean(mvpa_rate_bouted,na.rm=T))/sd(mvpa_rate_bouted,na.rm=T)]
acceleration[is.na(bpmed)]$bpmed <- 0

# Fix dates
for (j in (c('enroll_date','phenotype_censor_date'))){set(acceleration,j=j,value=as.Date(acceleration[[j]],format='%Y-%m-%d'))}

# Remove individuals with missing BMI
acceleration <- acceleration[!is.na(bmi)] # 96685 - 205 = 96480
acceleration <- acceleration[!is.na(tdi)] # 96480 - 107 = 96373
acceleration <- acceleration[!is.na(sbp)] # 96373 - 55 = 96318
acceleration <- acceleration[!is.na(dbp)] # 96318 - 0 = 96318
acceleration <- acceleration[!is.na(etoh_status)] # 96318 - 65 = 96253
acceleration[is.na(tob)]$tob <- "Never"

# Remove individuals in the latest withdrawal file
withdrawals <- fread('withdrawals/w7089_20210809.csv') # 96253 - 4 = 96249
acceleration <- acceleration[!(sample_id %in% withdrawals$V1)]

# Impute for alcohol (N=10454)
## Fix 1 outlier
acceleration <- acceleration[,etoh_grams := ifelse(etoh_grams > 2000, NA, etoh_grams)]

monthly_mean <- round(median(acceleration[c(!is.na(etoh_grams) & etoh_status=='Monthly')]$etoh_grams),0)
weekly_mean <- round(median(acceleration[c(!is.na(etoh_grams) & etoh_status=='Weekly')]$etoh_grams),0)

acceleration[c(etoh_status=='Weekly' & is.na(etoh_grams)),
             etoh_grams := (weekly_mean)]
acceleration[c(etoh_status=='Monthly' & is.na(etoh_grams)),
             etoh_grams := (monthly_mean)]

# Save out phenotype file
write.csv(acceleration,file='acceleration_phenotype_covar.csv',row.names = F)

### OUTPUT 1: VALUE BMI

# Scope columns for value BMI
value_bmi <- acceleration[,c('sample_id','accel_date','age_accel','sex','bmi','value_std','phenotype_censor_date',
                             'sbp','dbp','bpmed','tob','tdi','etoh_grams')]

# Write out
write.csv(value_bmi,file='cox_data_value_bmi_covar.csv',row.names = F)

### OUTPUT 2: INVNORM VALUE BMI

# Scope columns for value BMI invnorm
value_bmi_invnorm <- acceleration[,c('sample_id','accel_date','age_accel','sex','bmi','value_std','phenotype_censor_date',
                                     'sbp','dbp','bpmed','tob','tdi','etoh_grams')]

# Inverse normal transformation
value_bmi_invnorm[,value_invnorm := qnorm((rank(value_std,na.last="keep")-0.5)/sum(!is.na(value_std)))]

# Write out
write.csv(value_bmi_invnorm,file='cox_data_value_bmi_invnorm_covar.csv',row.names = F)

### OUTPUT 3: DECILE VALUE BMI
# Scope columns for value BMI invnorm
value_bmi_decile <- acceleration[,c('sample_id','accel_date','age_accel','sex','bmi','value_std','phenotype_censor_date',
                                    'sbp','dbp','bpmed','tob','tdi','etoh_grams')]

# Inverse normal transformation
value_bmi_decile[,value_decile := quantilize(value_std,ncuts=10)]

# Write out
write.csv(value_bmi_decile,file='cox_data_value_bmi_decile_covar.csv',row.names = F)

### OUTPUT 4: MVPA RATE BMI

# Scope columns for value BMI
mvpa_bmi <- acceleration[,c('sample_id','accel_date','age_accel','sex','bmi','mvpa_std','phenotype_censor_date',
                            'sbp','dbp','bpmed','tob','tdi','etoh_grams')]

# Write out
write.csv(mvpa_bmi,file='cox_data_mvpa_rate_bmi_covar.csv',row.names = F)

### OUTPUT 5: MVPA RATE CUTOFF BMI
# Scope columns for cutoff BMI
mvpa_bmi <- acceleration[,c('sample_id','accel_date','age_accel','sex','bmi','who_acc_rate_bouted','phenotype_censor_date',
                            'sbp','dbp','bpmed','tob','tdi','etoh_grams')]

# Write out
write.csv(mvpa_bmi,file='cox_data_mvpa_cutoff_bmi_covar.csv',row.names = F)

### OUTPUT 6: LOW CUTOFF BMI
# Scope columns for cutoff BMI
mvpa_bmi <- acceleration[,c('sample_id','accel_date','age_accel','sex','bmi','cutoff75','phenotype_censor_date',
                            'sbp','dbp','bpmed','tob','tdi','etoh_grams')]

# Write out
write.csv(mvpa_bmi,file='cox_data_mvpa_cutoff75_bmi_covar.csv',row.names = F)

### OUTPUT 7: HIGH CUTOFF BMI
# Scope columns for cutoff BMI
mvpa_bmi <- acceleration[,c('sample_id','accel_date','age_accel','sex','bmi','cutoff300','phenotype_censor_date',
                            'sbp','dbp','bpmed','tob','tdi','etoh_grams')]

# Write out
write.csv(mvpa_bmi,file='cox_data_mvpa_cutoff300_bmi_covar.csv',row.names = F)

### OUTPUT 8: AGE SUBGROUPS
young <- acceleration[age_accel < 55]
young[,mvpa_std := (mvpa_rate_bouted - mean(mvpa_rate_bouted,na.rm=T))/sd(mvpa_rate_bouted,na.rm=T)]

medium <- acceleration[age_accel >= 55 & age_accel < 65]
medium[,mvpa_std := (mvpa_rate_bouted - mean(mvpa_rate_bouted,na.rm=T))/sd(mvpa_rate_bouted,na.rm=T)]

old <- acceleration[age_accel >= 65]
old[,mvpa_std := (mvpa_rate_bouted - mean(mvpa_rate_bouted,na.rm=T))/sd(mvpa_rate_bouted,na.rm=T)]

# Scope columns for cutoff BMI
mvpa_bmi_young <- young[,c('sample_id','accel_date','age_accel','sex','bmi','mvpa_std','phenotype_censor_date',
                           'sbp','dbp','bpmed','tob','tdi','etoh_grams')]
mvpa_bmi_medium <- medium[,c('sample_id','accel_date','age_accel','sex','bmi','mvpa_std','phenotype_censor_date',
                             'sbp','dbp','bpmed','tob','tdi','etoh_grams')]
mvpa_bmi_old <- old[,c('sample_id','accel_date','age_accel','sex','bmi','mvpa_std','phenotype_censor_date',
                       'sbp','dbp','bpmed','tob','tdi','etoh_grams')]

# Write out
write.csv(mvpa_bmi_young,file='cox_data_mvpa_bmi_young_covar.csv',row.names = F)
write.csv(mvpa_bmi_medium,file='cox_data_mvpa_bmi_medium_covar.csv',row.names = F)
write.csv(mvpa_bmi_old,file='cox_data_mvpa_bmi_old_covar.csv',row.names = F)

### OUTPUT 9: BLANKING PERIOD SENSITIVITY ANALYSIS
# Create blanked date (2 years after accelerometer)
acceleration[,blanked_date := accel_date + 365.25*2]

# Scope columns for value BMI
mvpa_bmi <- acceleration[,c('sample_id','accel_date','blanked_date','age_accel','sex','bmi','mvpa_std','phenotype_censor_date',
                            'sbp','dbp','bpmed','tob','tdi','etoh_grams')]

# Write out
write.csv(mvpa_bmi,file='cox_data_mvpa_rate_bmi_blanked_covar.csv',row.names = F)

### OUTPUT 10: DECILES FOR DECILE FINDER
# Scope columns for value BMI
acceleration[,mvpa_decile := factor(quantilize(mvpa_rate_bouted,10),levels = as.character(1:10))]

# Scope columns
mvpa_bmi <- acceleration[,c('sample_id','accel_date','age_accel','sex','bmi','mvpa_decile','phenotype_censor_date',
                            'sbp','dbp','bpmed','tob','tdi','etoh_grams')]

# Write out
write.csv(mvpa_bmi,file='cox_data_mvpa_rate_bmi_decile_covar.csv',row.names = F)

### OUTPUT 11: QUINTILES 
# Scope columns for value BMI
acceleration[,mvpa_decile := factor(quantilize(mvpa_rate_bouted,5),levels = as.character(1:5))]

# Scope columns
mvpa_bmi <- acceleration[,c('sample_id','accel_date','age_accel','sex','bmi','mvpa_decile','phenotype_censor_date',
                            'sbp','dbp','bpmed','tob','tdi','etoh_grams')]

# Write out
write.csv(mvpa_bmi,file='cox_data_mvpa_rate_bmi_quintile_covar.csv',row.names = F)

### OUTPUT 11: RAW MVPA RATE BOUTED
# Scope columns
mvpa_bmi <- acceleration[,c('sample_id','accel_date','age_accel','sex','bmi','mvpa_rate_bouted','phenotype_censor_date',
                            'sbp','dbp','bpmed','tob','tdi','etoh_grams')]

# Write out
write.csv(mvpa_bmi,file='cox_data_mvpa_rate_bmi_nonstd_covar.csv',row.names = F)

### OUTPUT 12: VPA RATE BOUTED
# Read in vig files
vig_files <- list.files('/accel/')[str_detect(list.files('/accel/'),'summary_vig')]
vig <- list()
for (i in vig_files){
  load(paste0('/accel/',i))
  vig[[i]] <- output
}
vig_all <- data.table(do.call(rbind,vig))
setnames(vig_all,c('V1','V2','V3'),c('sample_id','vpa','total'))

# Join to main datasets
setkey(linker,app17488); setkey(vig_all,sample_id)
vig_all[linker,sek_id := i.app7089]

setkey(vig_all,sek_id); setkey(acceleration,sample_id)
acceleration[vig_all,vpa := i.vpa]

# VPA variables
# Rate in minutes per week
acceleration[,":="(vpa_rate = (vpa/total)/(total/604800)*604800/60)]
acceleration[,vpa_std := (vpa_rate - mean(vpa_rate,na.rm=T))/sd(vpa_rate,na.rm=T)]

# Deciles
acceleration[,vpa_decile := factor(quantilize(vpa_std,5),levels = as.character(1:5))]

# Scope columns
vpa_bmi <- acceleration[,c('sample_id','accel_date','age_accel','sex','bmi','vpa_std','vpa_decile','phenotype_censor_date',
                           'sbp','dbp','bpmed','tob','tdi','etoh_grams')]

# Write out
write.csv(vpa_bmi,file='cox_data_vpa_rate_bmi_covar.csv',row.names = F)

# Quintiles
acceleration[,vpa_decile := factor(quantilize(vpa_std,5),levels = as.character(1:5))]

# Scope columns
vpa_bmi <- acceleration[,c('sample_id','accel_date','age_accel','sex','bmi','vpa_std','vpa_decile','phenotype_censor_date',
                           'sbp','dbp','bpmed','tob','tdi','etoh_grams')]

# Write out
write.csv(vpa_bmi,file='cox_data_vpa_rate_bmi_quintile_covar.csv',row.names = F)





