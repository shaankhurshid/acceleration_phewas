# Depends
library(data.table)
library(plyr)

########### Set 1: Acceleration/BMI
# Load phenotype files
value_bmi <- fread(file='/Volumes/medpop_afib/skhurshid/acceleration_phewas/cox_data_value_bmi.csv')

# Load phenotypes
ht0 <- fread('/Volumes/medpop_afib/skhurshid/phenotypes/2020_06/height0.csv')
wt0 <- fread('/Volumes/medpop_afib/skhurshid/phenotypes/2020_06/weight0.csv')
ht1 <- fread('/Volumes/medpop_afib/skhurshid/accel/height_instance1.csv')
wt1 <- fread('/Volumes/medpop_afib/skhurshid/accel/weight_instance1.csv')
sex <- fread('/Volumes/medpop_afib/skhurshid/phenotypes/2020_06/sex0.csv')
censor <- fread('/Volumes/medpop_afib/skhurshid/phenotypes/2021_06/censor_202106.csv')
htn <- fread('/Volumes/medpop_afib/skhurshid/phenotypes/2021_06/htn_202106.csv')
dm <- fread('/Volumes/medpop_afib/skhurshid/phenotypes/2021_06/any_dm_202106.csv')
race <- fread('/Volumes/medpop_afib/skhurshid/phenotypes/2020_06/race_0.csv')
tob0 <- fread('/Volumes/medpop_afib/skhurshid/phenotypes/tobacco_instance_0.csv')
tob1 <- fread('/Volumes/medpop_afib/skhurshid/accel/tobacco_instance1.csv')

# Create combo fields
setkey(ht0,sample_id); setkey(ht1,sample_id); setkey(wt0,sample_id); setkey(wt1,sample_id)
setkey(tob0,sample_id); setkey(tob1,sample_id)

ht0[ht1,ht1 := i.value]
wt0[wt1,wt1 := i.value]

ht0[,ht_unified := ifelse(is.na(ht1),height,ht1)]
wt0[,wt_unified := ifelse(is.na(wt1),weight,wt1)]

tob0[,tobacco_accel_selfreport := ifelse(value==2,"Current",
                                        ifelse(value==1,"Former","Never"))]

tob <- tob0[tob1,':='(value1 = i.value)]
tob[,value_unified := ifelse(!is.na(value1),value1,value)]
tob[,tobacco_accel_selfreport := ifelse(value_unified==2,"Current",
                                        ifelse(value_unified==1,"Former","Never"))]

# Joins
setkey(censor,sample_id); setkey(value_bmi,sample_id); setkey(sex,sample_id); setkey(acceleration,sample_id)
setkey(htn,sample_id); setkey(dm,sample_id)
value_bmi[ht0,ht := i.ht_unified]
value_bmi[wt0,wt := i.wt_unified]
value_bmi[tob,tob := i.tobacco_accel_selfreport]
value_bmi[censor,':='(birthdate = i.birthdate,enroll_date = i.enroll_date)]
value_bmi[acceleration,':='(age_accel = i.age_accel, sex = i.sex)]

# Prevalent disease
value_bmi[htn,':='(has_htn = i.has_disease, htn_censor_age = i.censor_age)]
value_bmi[dm,':='(has_dm = i.has_dm, dm_censor_age = i.any_dm_age)]
value_bmi[,':='(prev_dm = ifelse(c(has_dm == 1 & dm_censor_age <= age_accel),1,0),
                prev_htn = ifelse(c(has_htn == 1 & htn_censor_age <= age_accel),1,0))]

# Granular race
race[,race_categories := ifelse(value %in% c(1,1001,1002,1003),'white',
                                ifelse(value %in% c(2,2001,2002,2003,2004),'mixed',
                                       ifelse(value %in% c(3,3001,3002,3003,3004,5),'asian',
                                              ifelse(value %in% c(4,4001,4002,4003),'black',
                                                     ifelse(value %in% 6,'other',NA)))))]
setkey(race,sample_id)
value_bmi[race,race_category := i.race_categories]

########### Set 2: PRS NO ACCEL
prs_no_accel <- fread(file='/Volumes/medpop_afib/skhurshid/acceleration_phewas/cox_data_mvpa_sqrt_prs.csv')
setkey(prs_no_accel,sample_id)
prs_no_accel[ht0,ht := i.height]
prs_no_accel[wt0,wt := i.weight]
prs_no_accel[sex,sex := i.height]
prs_no_accel[censor,':='(birthdate = i.birthdate,enroll_date = i.enroll_date,enroll_age = i.enroll_age)]
prs_no_accel[race,race_category := i.race_categories]
prs_no_accel[tob0,tob := i.tobacco_accel_selfreport]
prs_no_accel[,bmi := wt/((ht/100)**2)]

# Prevalent disease
prs_no_accel[htn,':='(has_htn = i.has_disease, htn_censor_date = i.censor_date)]
prs_no_accel[dm,':='(has_dm = i.has_disease, dm_censor_date = i.censor_date)]
prs_no_accel[,':='(prev_dm = ifelse(c(has_dm == 1 & dm_censor_date <= enroll_date),1,0),
                prev_htn = ifelse(c(has_htn == 1 & htn_censor_date <= enroll_date),1,0))]

########### Set 3: Self
self <- fread(file='/Volumes/medpop_afib/skhurshid/acceleration_phewas/cox_data_self_bmi.csv')
setkey(self,sample_id)
self[ht0,ht := i.height]
self[wt0,wt := i.weight]
self[sex,sex := i.height]
self[tob0,tob := i.tobacco_accel_selfreport]
self[censor,':='(birthdate = i.birthdate,enroll_date = i.enroll_date,enroll_age = i.enroll_age)]
self[race,race_category := i.race_categories]
self[,bmi := wt/((ht/100)**2)]

# Prevalent disease
self[htn,':='(has_htn = i.has_disease, htn_censor_date = i.censor_date)]
self[dm,':='(has_dm = i.has_disease, dm_censor_date = i.censor_date)]
self[,':='(prev_dm = ifelse(c(has_dm == 1 & dm_censor_date <= enroll_date),1,0),
                   prev_htn = ifelse(c(has_htn == 1 & htn_censor_date <= enroll_date),1,0))]


