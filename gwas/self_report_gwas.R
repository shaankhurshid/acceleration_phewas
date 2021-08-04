# Script to prep a self-report activity gwas

# Load activity data
freq_walk <- fread(file='/Volumes/medpop_afib/skhurshid/accel/self_report/freq_walk.csv')
freq_walk1 <- fread(file='/Volumes/medpop_afib/skhurshid/accel/self_report/freq_walk1.csv')
freq_mod <- fread(file='/Volumes/medpop_afib/skhurshid/accel/self_report/freq_mod.csv')
freq_mod1 <- fread(file='/Volumes/medpop_afib/skhurshid/accel/self_report/freq_mod1.csv')
freq_vigorous <- fread(file='/Volumes/medpop_afib/skhurshid/accel/self_report/freq_vigorous.csv')
freq_vigorous1 <- fread(file='/Volumes/medpop_afib/skhurshid/accel/self_report/freq_vigorous1.csv')
duration_walk <- fread(file='/Volumes/medpop_afib/skhurshid/accel/self_report/duration_walks.csv')
duration_walk1 <- fread(file='/Volumes/medpop_afib/skhurshid/accel/self_report/duration_walks1.csv')
duration_mod <- fread(file='/Volumes/medpop_afib/skhurshid/accel/self_report/duration_mod.csv')
duration_mod1 <- fread(file='/Volumes/medpop_afib/skhurshid/accel/self_report/duration_mod1.csv')
duration_vigorous <- fread(file='/Volumes/medpop_afib/skhurshid/accel/self_report/duration_vigorous.csv')
duration_vigorous1 <- fread(file='/Volumes/medpop_afib/skhurshid/accel/self_report/duration_vigorous1.csv')

# Join into single table
## Set keys
setkey(freq_walk,sample_id); setkey(freq_walk1,sample_id); setkey(freq_mod,sample_id); setkey(freq_mod1,sample_id)
setkey(freq_vigorous,sample_id); setkey(freq_vigorous1,sample_id); setkey(duration_walk,sample_id); setkey(duration_walk1,sample_id)
setkey(duration_mod,sample_id); setkey(duration_mod1,sample_id); setkey(duration_vigorous,sample_id); setkey(duration_vigorous1,sample_id)

## Join
freq_walk <- freq_walk[freq_walk1,freq_walk1 := i.freq_walk]
freq_mod <- freq_mod[freq_mod1,freq_mod1 := i.freq_mod]
freq_vigorous <- freq_vigorous[freq_vigorous1,freq_vigorous1 := i.freq_vigorous]
duration_walk <- duration_walk[duration_walk1,duration_walk1 := i.duration_walk]
duration_mod <- duration_mod[duration_mod1,duration_mod1 := i.duration_mod]
duration_vigorous <- duration_vigorous[duration_vigorous1,duration_vigorous1 := i.duration_vigorous]

activity <- freq_walk[freq_mod,':='(freq_mod = i.freq_mod, freq_mod1 = i.freq_mod1)]
activity <- activity[freq_vigorous,':='(freq_vigorous = i.freq_vigorous, freq_vigorous1 = i.freq_vigorous1)]
activity <- activity[duration_walk,':='(duration_walk = i.duration_walks, duration_walk1 = i.duration_walk1)]
activity <- activity[duration_mod,':='(duration_mod = i.duration_moderate, duration_mod1 = i.duration_mod1)]
activity <- activity[duration_vigorous,':='(duration_vigorous = i.duration_vigorous, duration_vigorous1 = i.duration_vigorous1)]

## Unify variables
activity[,':='(freq_walk_unified = ifelse(is.na(freq_walk1),freq_walk,freq_walk1),
               freq_mod_unified = ifelse(is.na(freq_mod1),freq_mod,freq_mod1),
               freq_vigorous_unified = ifelse(is.na(freq_vigorous1),freq_vigorous,freq_vigorous1),
               duration_walk_unified = ifelse(is.na(duration_walk1),duration_walk,duration_walk1),
               duration_mod_unified = ifelse(is.na(duration_mod1),duration_mod,duration_mod1),
               duration_vigorous_unified = ifelse(is.na(duration_vigorous1),duration_vigorous,duration_vigorous1))]

## Set non-answers to NA (-1 = Don't know, -2 = Cant' walk, and -3 = Prefer not to answer)
names <- c('freq_walk_unified','freq_mod_unified','freq_vigorous_unified','duration_walk_unified','duration_mod_unified','duration_vigorous_unified')
for (j in (names)){set(activity,i=which(activity[[j]] %in% c(-1,-2,-3)),j=j,value=NA)}

names <- c('freq_walk1','freq_mod1','freq_vigorous1','duration_walk1','duration_mod1','duration_vigorous1')
for (j in (names)){set(activity,i=which(activity[[j]] %in% c(-1,-2,-3)),j=j,value=NA)}

## Generate MET summary columns
activity[,':='(met_walk = duration_walk_unified*freq_walk_unified*3.3,
               met_walk0 = duration_walk*freq_walk*3.3,
               met_walk1 = duration_walk1*freq_walk1*3.3,
               met_mod = duration_mod_unified*freq_mod_unified*4,
               met_mod0 = duration_mod*freq_mod*4,
               met_mod1 = duration_mod1*freq_mod1*4,
               met_vigorous = duration_vigorous_unified*freq_vigorous_unified*8,
               met_vigorous0 = duration_vigorous*freq_vigorous*8,
               met_vigorous1 = duration_vigorous1*freq_vigorous1*8)]

activity[!c(is.na(met_walk) & is.na(met_mod) & is.na(met_vigorous)),
         total_met := apply(.SD,1,sum,na.rm=T),.SDcols=c('met_walk','met_mod','met_vigorous')]

activity[!c(is.na(met_walk1) & is.na(met_mod1) & is.na(met_vigorous1)),
         total_met1 := apply(.SD,1,sum,na.rm=T),.SDcols=c('met_walk1','met_mod1','met_vigorous1')]

activity[!c(is.na(met_walk0) & is.na(met_mod0) & is.na(met_vigorous0)),
         total_met0 := apply(.SD,1,sum,na.rm=T),.SDcols=c('met_walk0','met_mod0','met_vigorous0')]

# Load withdrawals
withdrawals <- fread('/Volumes/medpop_afib/skhurshid/phenotypes/withdrawals/w7089_20210201.csv') # UKBB withdrawals

# Reduce dataset
self <- activity[!is.na(total_met),c('sample_id','total_met')]
self <- self[!(sample_id %in% withdrawals$V1)]

# Load censor master file (for exclusions)
censor <- fread(file='/Volumes/medpop_afib/skhurshid/phenotypes/2020_06/censor_202006.csv')

## Load phenotype file to get age
load(file='/Volumes/medpop_afib/skhurshid/accel/phenos_060520.RData')
setkey(phenos,ID); setkey(self,sample_id)
self[phenos,':='(agevisit0 = i.agevisit0)]

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
  #ab <- ab[ab$ex_misKin==0]
  #print(paste0('N after removal of missing kinship inference:',nrow(ab)))
  
  # Loop over "exclude all both" phenotypes - all individuals with exclusion phenotype at any time removed for both cases/controls
  if (length(exclude_all_both)!=0){
    for (i in exclude_all_both){
      exclude <- fread(paste0("/Volumes/medpop_afib/skhurshid/svt_brady_gwas/phenotypes/",i,'.tab.tsv'),header=T)
      setkey(ab,sample_id); setkey(exclude,sample_id)
      ab[exclude,':='(exclude_prev=i.prevalent_disease,exclude_incd=i.incident_disease,exclude_censor = i.censor_age)]
      ab[,exclude := ifelse(c(c(!is.na(exclude_prev) & exclude_prev==1) | c(!is.na(exclude_incd) & exclude_incd == 1)),1,0)]
      print(paste0('I am going to exclude ',sum(ab$exclude),' individuals for diagnosis: ',i,' occurring at any time for cases/controls'))
      ab <- ab[exclude==0]
      ab <- ab[,!(c('exclude_prev','exclude_incd','exclude_censor','exclude'))]
    }}
  
  # Loop over "exclude all cases" phenotypes - all individuals with exclusion phenotype at any time removed for cases
  if (length(exclude_all_cases)!=0){
    for (i in exclude_all_cases){
      exclude <- fread(paste0("/Volumes/medpop_afib/skhurshid/svt_brady_gwas/phenotypes/",i,'.tab.tsv'),header=T)
      setkey(ab,sample_id); setkey(exclude,sample_id)
      ab[exclude,':='(exclude_prev=i.prevalent_disease,exclude_incd=i.incident_disease,exclude_censor = i.censor_age)]
      ab[,exclude := ifelse(c(c(!is.na(get(trait)) & get(trait)==1) & 
                                c(c(!is.na(exclude_prev) & exclude_prev==1) | c(!is.na(exclude_incd) & exclude_incd==1))),1,0)]
      print(paste0('I am going to exclude ',sum(ab$exclude),' individuals for diagnosis: ',i,' occurring at any time for cases only'))
      ab <- ab[exclude==0]
      ab <- ab[,!(c('exclude_prev','exclude_incd','exclude_censor','exclude'))]
    }}
  
  # Loop over "exclude all controls" phenotypes - all individuals with exclusion phenotype at any time removed for controls
  if (length(exclude_all_controls)!=0){
    for (i in exclude_all_controls){
      exclude <- fread(paste0("/Volumes/medpop_afib/skhurshid/svt_brady_gwas/phenotypes/",i,'.tab.tsv'),header=T)
      setkey(ab,sample_id); setkey(exclude,sample_id)
      ab[exclude,':='(exclude_prev=i.prevalent_disease,exclude_incd=i.incident_disease,exclude_censor = i.censor_age)]
      ab[,exclude := ifelse(c(c(get(trait)==0 | is.na(get(trait))) & 
                                c(c(!is.na(exclude_prev) & exclude_prev==1) | c(!is.na(exclude_incd) & exclude_incd==1))),1,0)]
      print(paste0('I am going to exclude ',sum(ab$exclude),' individuals for diagnosis: ',i,' occurring at any time for controls only'))
      ab <- ab[exclude==0]
      ab <- ab[,!(c('exclude_prev','exclude_incd','exclude_censor','exclude'))]
    }}
  
  # Loop over "exclude incident" - only cases with exclusion phenotype before disease removed
  if (length(exclude_incident_cases)!=0){
    for (i in exclude_incident_cases){
      exclude <- fread(paste0("/Volumes/medpop_afib/skhurshid/svt_brady_gwas/phenotypes/",i,'.tab.tsv'),header=T)
      setkey(ab,sample_id); setkey(exclude,sample_id)
      ab[exclude,':='(exclude_disease = i.has_disease, exclude_prev = i.prevalent_dsease, exclude_censor = i.censor_date)]
      ab[,exclude := ifelse(c(c(!is.na(get(trait)) & get(trait)==1) & 
                                c(c(!is.na(exclude_disease) & (exclude_censor <= censor_date)) |
                                    c(!is.na(exclude_disease) & exclude_prev==1))),1,0)]
      print(paste0('I am going to exclude ',sum(ab$exclude),' individuals for diagnosis: ',i,' occurring before case diagnosis'))
      ab <- ab[exclude==0]
      ab <- ab[,!(c('exclude_prev','exclude_censor','exclude_disease','exclude'))]
    }}
  
  # Loop over "exclude flexible" - excludes any instance of exclusion phenotype among controls, and only exclusion phenotype prior to disease for cases
  if (length(exclude_flexible)!=0){
    for (i in exclude_flexible){
      exclude <- fread(paste0("/Volumes/medpop_afib/skhurshid/svt_brady_gwas/phenotypes/",i,'.tab.tsv'),header=T)
      setkey(ab,sample_id); setkey(exclude,sample_id)
      ab[exclude,':='(exclude_incd = i.incident_disease, exclude_disease = i.has_disease, exclude_prev = i.prevalent_disease, exclude_censor = i.censor_date)]
      ab[,exclude := ifelse(c(!is.na(get(trait)) & get(trait)==1),
                            ifelse(c(!is.na(exclude_prev) & exclude_prev==1),1,
                                   ifelse(c(!is.na(exclude_incd) & (exclude_incd==1) & (exclude_censor <= censor_date)),1,0)),
                            ifelse(c(!is.na(exclude_disease) & exclude_disease==1),1,0))]
      print(paste0('I am going to exclude ',sum(ab$exclude),' individuals for diagnosis: ',i,' occurring before case diagnosis or at any time for controls'))
      ab <- ab[exclude==0]
      ab <- ab[,!(c('exclude_disease','exclude_prev','exclude_incd','exclude_censor','exclude'))]
    }}
  
  #This file lists the pairs of individuals related up to the third degree in the data set. It is a plaintext file with space separated columns.
  rel<-fread("/Volumes/medpop_esp2/pradeep/UKBiobank/v2data/ukb708_rel_chr1_s488374.dat",header=T)
  pheno<-ab[,.SD,.SDcols=c('sample_id',trait)]
  names(pheno)<-c("ID","pheno")
  
  # Merge phenotypes with relative index
  setkey(rel,ID1); setkey(pheno,ID)
  rel <- rel[pheno,pheno.ID1:=i.pheno]
  setkey(rel,ID2)
  rel <- rel[pheno,pheno.ID2:=i.pheno]
  
  # Now create a list of people to remove for relatedness
  rel[,rel_id := ifelse(is.na(pheno.ID1) & !is.na(pheno.ID2),ID1,
                        ifelse(is.na(pheno.ID2) & !is.na(pheno.ID1),ID2,
                               ifelse(is.na(pheno.ID2) & is.na(pheno.ID2),ID2,
                                      ifelse(!is.na(pheno.ID1) & !is.na(pheno.ID2),ID2,ID2))))]
  
  rel_id<-unique(rel$rel_id) #81838 ->81847 
  
  ab[,':='(male = ifelse(Inferred_Gender=='M',1,0),
           ex_rel = ifelse(sample_id %in% rel_id,1,0))]
  
  print(paste0('N related:',sum(ab$ex_rel)))
  
  #######
  ###test AF related PCs in each cleaned dataset
  #######
  
  #dim(subset(ab, used_in_pca_calculation==1 & ex_sex==0 & ex_poor==0 & ex_misKin==0 & white==1))
  #######
  ##all white, but relatives OK for BOLT-LMM
  #######
  ab1<-ab
  #print(paste0('N white:',sum(ab$in_white_British_ancestry_subset)))
  
  #######
  ###test AF related PCs in each cleaned dataset
  #######
  
  #dim(subset(ab1, used_in_pca_calculation==1 & ex_sex==0 & ex_poor==0 & ex_misKin==0 & white==1))
  #######
  ##all white, no relatives
  #######
  form1<-formula(paste0(trait,"~agevisit0 + ",paste0("PC",1:40,collapse="+"),"+ array_UKBB + male",collapse="+"))
  s1<-summary(lm(form1,data=ab1))$coefficients
  s1<-s1[substring(rownames(s1),1,2)=="PC",]
  ab1.1<-ab1[,.SD,.SDcols=c("agevisit0",trait)]
  allN_1<-nrow(ab1.1)
  male_1<-nrow(ab1[ab1[["male"]]==1 & !is.na(ab1.1[[trait]])])
  pcs1<-paste(rownames(s1)[1:5],collapse=",") # First 5 PCs
  
  #######
  ##create summary file
  #######
  t1<-c("all",allN_1,"mean_age",round(mean(ab1[!is.na(trait)]$agevisit0),2),"sd_age",round(sd(ab1[!is.na(trait)]$agevisit0),2),"male_N",male_1,"male%",round(mean(ab1[!is.na(trait)]$male)*100,2),"related-PCs",ifelse((pcs1!=""),pcs1,"None"))
  write.table(t1,file=paste0("/Volumes/medpop_afib/skhurshid/accel/summary_",trait,".txt"),row.names=F,quote=F,sep="\t")
  
  #######
  ##create phenotype file
  #######
  ## Choose columns
  pheno<-ab1[,c("sample_id",trait,"agevisit0",rownames(s1)[1:5],"array_UKBB","male"),with=F]
  ## Format for PLINK
  setnames(pheno,"sample_id","FID")
  setnames(pheno,"agevisit0","enroll_age")
  pheno[,':='(IID = FID)]
  setcolorder(pheno,c('FID','IID'))
  print(paste0('Final phenotype N: ',nrow(pheno)))
  write.table(pheno,file=paste0('/Volumes/medpop_afib/skhurshid/acceleration_phewas/',trait,".tsv"),sep="\t",col.names =T,row.names = F,quote = F)
}

#sqc<-ab14[,c(1,29,31:70,74,75,80)]#486553
#write.tab1le(sqc,file="/medpop/afib/lcweng/UKBB_all/QC/sqc.tsv",sep="\t",col.names =T,row.names = F,quote = F)

#exclusion phenotypes
# SND - exclude valve disease, cardiac surgery, MI at or prior to SND (cases)
# SND - exclude valve disease, cardiac surgery, MI, SND-inclusive, DCD-inclusive, PM
# DCD - exclude valve disease, cardiac surgery, MI at or prior to DCD cases)
# DCD - exclude valve disease, cardiac surgery, MI, SND-inclusive, DCD-inclusive, PM (controls)
# PM - exclude valve disease, cardiac surgery, MI at or prior to PM (cases)
# PM - exclude valve disease, cardiac surgery, MI, SND-inclusive, DCD-inclusive, PM (controls)
# WPW - exclude HCM, Ebstein (cases)
# WPW - exclude HCM, Ebstein (controls)
# SVT - exclude HCM, Ebstein (cases)
# SVT - exclude HCM, Ebstein (controls)

create(trait="total_met",data=self)

# # Inverse normal version
# total_met <- fread(file='/Volumes/medpop_afib/skhurshid/acceleration_phewas/total_met.tsv')
# total_met[,total_met_invnorm := qnorm((rank(total_met,na.last="keep")-0.5)/sum(!is.na(total_met)))]
# total_met[,total_met := NULL]
# write.table(total_met,file='/Volumes/medpop_afib/skhurshid/acceleration_phewas/total_met_invnorm.tsv',sep="\t",col.names =T,row.names = F,quote = F)

# Exclusions list for BOLT
exclusions <- fread(file='/Volumes/medpop_afib/skhurshid/acceleration_phewas/total_met.tsv')
exclusions <- censor[!(sample_id %in% exclusions$FID)]$sample_id
exclusions <- data.table(IID = exclusions, FID = exclusions)
write.table(exclusions,file='/Volumes/medpop_afib/skhurshid/acceleration_phewas/self_exclusions.tsv',sep="\t",col.names =T,row.names = F,quote = F)


