# Script to prepare PRS for PHB replication

# Depends
library(data.table)
library(stringr)

# Load PRS
prs <- fread(file='/data/arrhythmia_source/lcweng/Acceleration/score_white_ID.txt')

# Load PCs
pcs <- fread(file='/data/arrhythmia_source/data/phsBio/phsBio_201812/genotyped/sal44/Biobank_sal44_20181227_103631_x30717_VCF/preIMPUTE_QC/PHB_PCA_round2.txt')

# Load whites
white <- fread(file='/data/arrhythmia_source/data/phsBio/phsBio_201812/genotyped/sal44/Biobank_sal44_20181227_103631_x30717_VCF/preIMPUTE_QC/PHB_EUR_ID.txt')

# Load GWAS linker
gwas_linker <- fread(file='/data/arrhythmia_source/data/phsBio/phsBio_201812/imputed/ID/matchID.txt')

# Load MGB linker
linker <- fread('/data/arrhythmia_source/data/phsBio_201910/linker/patient_linker_file_phsbio.csv')

# Add individual ID to PCs file
setkey(pcs,IID); setkey(gwas_linker,V2)
pcs[gwas_linker,ID := i.V4]

# Separate IDs in PRS file
prs[,':='(FID = str_extract(ID,'FAM\\d+'),
          IID = str_extract(ID,'IND\\d+'))]

# Link PCs file and PRS file
setkey(pcs,ID); setkey(prs,IID)
prs[pcs,paste0('PC',1:5) := .(i.PC1,i.PC2,i.PC3,i.PC4,i.PC5)]
prs[pcs,lcw_id := i.IID]

# Isolate to whites
prs <- prs[lcw_id %in% white$V2]

# Link
setkey(linker,EMPI); setkey(prs,EMPI)
prs[linker,linker_id := i.linker_id]

# Standardize PRS variable
prs[,prs_std := (SCORE_SUM - mean(SCORE_SUM))/sd(SCORE_SUM)]

# Remove identifiers
prs <- prs[,c('linker_id','collect_age','prs_std',paste0('PC',1:5))]

# Save out
write.csv(prs,file='/data/genes/skhurshid/accel_prs_mgb_eur.csv',row.names=F,quote=F)