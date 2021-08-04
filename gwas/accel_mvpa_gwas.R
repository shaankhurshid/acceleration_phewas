# Script for formatting acceleration data for GWAS

# Depends
library(data.table)

# Load raw acceleration data
acceleration <- fread('/Volumes/medpop_afib/skhurshid/accel/ukbb_accel2.csv')
load(file='/Volumes/medpop_afib/skhurshid/accel/accel_raw_cat.RData')

# Load withdrawals
withdrawals <- fread('/Volumes/medpop_afib/skhurshid/phenotypes/withdrawals/w7089_20210201.csv') # UKBB withdrawals

# Load censor master file (for exclusions)
censor <- fread(file='/Volumes/medpop_afib/skhurshid/phenotypes/2020_06/censor_202006.csv')

# Recommended QC for accelerometer data
## Could not calibrate
no_calibrate <- acceleration[calibration==1] # 6996
## Inadequate wear time
acceleration <- no_calibrate[wear_time==1] # 4

## Remove revoked consent
acceleration <- acceleration[!(sample_id %in% withdrawals$V1)] # 7

# Remove in plink but not imputed
not_imputed <- fread('/Volumes/medpop_afib/skhurshid/lvm_gwas/qc/bolt.in_plink_but_not_imputed.FID_IID.968.txt') # UKBB withdrawals
acceleration <- acceleration[!(sample_id %in% not_imputed$V1)]

# Keep only in PLINK
gt <- fread(file='/Volumes/medpop_afib/skhurshid/central_qc/in_v2_gt.csv')
acceleration <- acceleration[(sample_id %in% gt$FID)]

# Remove individual missing GT (>0.1)
missing_indiv <- fread(file='/Volumes/medpop_afib/skhurshid/central_qc/v2_gt_missing_indiv10.csv')
acceleration <- acceleration[!(sample_id %in% missing_indiv$IID)]

# Join total time
setkey(dt_all,sek_id); setkey(acceleration,sample_id)
acceleration[dt_all,total := i.total]

# Load bouts files
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

# Load the ID key
linker <- fread('/Volumes/medpop_afib/skhurshid/accel/7089_17488_linker.csv')

# Joins
setkey(linker,app17488); setkey(bouts,sample_id)
bouts[linker,sek_id := i.app7089]

setkey(bouts,sek_id); setkey(acceleration,sample_id)
acceleration[bouts,mvpa_bouts := (i.mvpa)*60*5]

# Create guideline-based cutoffs
acceleration[,':='(mvpa_rate_bouted = (mvpa_bouts/60)/(total/86400))]

## Add in age at accel
setkey(censor,sample_id); setkey(acceleration,sample_id)
acceleration[censor,':='(birthdate = i.birthdate)]

# Fix dates
for (j in c('start_date','end_date','birthdate')){
  set(acceleration,j=j,value=as.Date(substr(acceleration[[j]],1,10),format='%Y-%m-%d'))}

# Age at accel
acceleration[,age_accel := as.numeric(end_date - birthdate)/365.25]

## PLINK output function
create<-function(trait,exclude_all_both=NULL,exclude_all_cases=NULL,exclude_all_controls=NULL,
                 exclude_incident_cases=NULL,exclude_flexible=NULL,data){
  ##phenotype file
  a <- data
  a<-a[!is.na(get(trait))]
  print(paste0('Total N:',nrow(a)))
  
  ##sample qc file
  b<-fread("/Volumes/medpop_esp2/pradeep/UKBiobank/v2data/ukb_sqc_v2_7089.tsv",header=T) 
  setkey(a,'sample_id'); setkey(b,'eid')
  ab <- a[b,nomatch=0]
  ab[,':='(array_UKBB = ifelse(genotyping_array=='UKBB',1,0))]
  print(paste0('N after merge with sample QC file:',nrow(ab)))
  
  ##remove poor quality
  ab[,':='(ex_poor = ifelse(het_missing_outliers==1 | putative_sex_chromosome_aneuploidy==1,1,0),
           ex_sex = ifelse(Submitted_Gender==Inferred_Gender,0,1),
           ex_misKin = ifelse(ab$excluded_from_kinship_inference==1,1,0))]
  
  #high quality data
  ab <- ab[ab$ex_sex==0]
  print(paste0('N after removal of sex mismatch:',nrow(ab)))
  ab <- ab[ab$ex_poor==0]
  print(paste0('N after removal of poor:',nrow(ab)))
  
  #Renames
  ab[,':='(male = ifelse(Inferred_Gender=='M',1,0))]
  ab1<-ab
  pcs <- paste0('PC',1:5)  #First 5 PCs
  
  #######
  ##create phenotype file
  #######
  ## Choose columns
  pheno<-ab1[,c("sample_id",trait,"age_accel",pcs,"array_UKBB","male"),with=F]
  ## Format for PLINK
  setnames(pheno,"sample_id","FID")
  pheno[,':='(IID = FID)]
  setcolorder(pheno,c('FID','IID'))
  print(paste0('Final phenotype N: ',nrow(pheno)))
  write.table(pheno,file=paste0('/Volumes/medpop_afib/skhurshid/acceleration_phewas/',trait,".tsv"),sep="\t",col.names =T,row.names = F,quote = F)
}

create(trait="mvpa_rate_bouted",data=acceleration)

# Exclusions list for BOLT
## Read processed phenotype
acceleration <- fread(file='/Volumes/medpop_afib/skhurshid/acceleration_phewas/mvpa_rate_bouted.tsv')
## Load list of all individuals in UKBB (includes people in plink GT set that are not in censor files)
ukbb_all <- fread(file='/Volumes/medpop_afib/skhurshid/central_qc/ukbb_all.csv')
exclusions <- ukbb_all[!(FID %in% acceleration$FID)]
## Add individuals that are in the plink dataset but not the censor file
write.table(exclusions,file='/Volumes/medpop_afib/skhurshid/acceleration_phewas/mvpa_rate_bouted_exclusions.tsv',sep="\t",col.names =T,row.names = F,quote = F)