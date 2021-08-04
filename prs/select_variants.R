# Script to select variants per chromosome from clump file 

# Depends
library(data.table)

# Loop function
divide <- function(clump,gwas,output_path){
  for (i in 1:22){
  chr_list <- gwas[CHR==i] 
  select_list <- chr_list[SNP %in% clump$V1]$SNP
  write.table(select_list,file=paste0(output_path,'variants_chr',i,'.tsv'),sep='\t',quote=F,row.names=F)
  }
}

# Load clump file
clump <- fread(file='/Volumes/medpop_afib/skhurshid/acceleration_phewas/prs/clumps/clump4_rsid.tsv',header = F)

# Load gwas
gwas <- fread(file='/medpop/afib/skhurshid/accel/gwas/bolt_output/v2/value_bolt.imputed.unambiguous.tsv')

divide(clump=clump,gwas=gwas,output_path='/medpop/afib/skhurshid/acceleration_phewas/no_accel_prs_assoc/variant_list/')

