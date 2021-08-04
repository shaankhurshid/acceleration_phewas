# Script to generate unrolled ICD10 phecode map
# Maps obtained for Phecode Map v1.2 (https://phewascatalog.org/phecodes_icd10)
library(stringr)

# Load existing maps
icd10 = fread('~/Documents/MGH Research/accel_phewas/icd10_phecode_map.csv')
all = fread('~/Documents/MGH Research/accel_phewas/icd9_10_phecode_unrolled.txt')
icd10_cm <- all[flag==10]

# Identify parent categories
## Get number of entries per phecode
icd10_cm_ct <- icd10_cm[,.N,by='phecode']
## Isolate to N > 1
parents <- icd10_cm_ct[N > 1]

# Now loop through and gather ICD10 codes for each parent
output <- data.table()
for (i in parents$phecode){
  codes <- icd10[PheCode==i | str_detect(PheCode,paste0('^',i,'\\.\\d+'))]$ICD10
  result <- data.table(ICD10 = codes, `ICD10 String` = rep(NA,length(codes)), PheCode = rep(i,length(codes)),
                       Phenotype = rep(NA,length(codes)),
                       `Excl. Phecodes` = rep(NA,length(codes)), `Excl. Phenotypes` = rep(NA,length(codes)))
  output <- rbind(output,result)
}

# Add back to ICD 10 phecode map
icd10 <- rbind(icd10,output)

# Remove duplicates within PheCode
icd10 = unique(icd10,by=c('PheCode','ICD10'))

# Write out
write.csv(icd10,file='~/Documents/MGH Research/accel_phewas/icd10_phecode_map_unrolled_sk.csv',row.names=F)


