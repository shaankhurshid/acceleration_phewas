# Depends
library(data.table)
library(stringr)

# Load master phecode library
all <- fread(file='~/Documents/MGH Research/accel_phewas/phecode_definitions/icd9_10cm_unrolled.txt')

# Load unrolled ICD-10 (not CM)
icd10 <- fread(file='~/Documents/MGH Research/accel_phewas/phecode_definitions/icd10_phecode_map_unrolled_sk.csv')

# Isolate to codes not already in the library
icd10_extra <- icd10[!(ICD10 %in% all[flag==10]$ICD)]
icd10_extra <- icd10_extra[!is.na(PheCode) & !is.na(ICD10) & !is.na(`ICD10 String`)]
icd10_extra <- data.table(ICD = icd10_extra$ICD10, flag=10, phecode = icd10_extra$PheCode)
all <- rbind(all,icd10_extra)

# Test
n <- 1
for (i in unique(all$phecode)){
  display_name <- str_replace_all(as.character(i),'\\.','\\_')
  tens <- all[phecode==i & flag==10]$ICD
  nines <- all[phecode==i & flag==9]$ICD
  if (c(length(tens)==0 & length(nines)==0)){next}
  else if (length(nines)==0){
    output <- data.table(Field=c(41202,41204,40001,40002),
                         Coding=rep(paste0(tens,collapse=','),4),
                         exclude=rep(0,4))
    fwrite(output,paste0('~/Documents/MGH Research/accel_phewas/phecode_tabs/phecode_',display_name,'.tab'),row.names=F,quote=F,sep='\t')
  }
  else if (length(tens)==0){
    output <- data.table(Field=c(41203,41205),
                         Coding=rep(paste0(nines,collapse=','),2),
                         exclude=c(0,0))
    fwrite(output,paste0('~/Documents/MGH Research/accel_phewas/phecode_tabs/phecode_',display_name,'.tab'),row.names=F,quote=F,sep='\t')
}
  else {
    output <- data.table(Field=c(41203,41205,41202,41204,40001,40002),
                         Coding=c(rep(paste0(nines,collapse=','),2),rep(paste0(tens,collapse=','),4)),
                         exclude=c(rep(0,6)))
    fwrite(output,paste0('~/Documents/MGH Research/accel_phewas/phecode_tabs/phecode_',display_name,'.tab'),row.names=F,quote=F,sep='\t')
  }
if (n %% 50 == 0){print(paste0("Just finished tabfile ",n," out of ",length(all),"!"))}
n <- n + 1
}