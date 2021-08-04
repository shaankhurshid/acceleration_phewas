library(stringr)
library(data.table)

post_clump <- function(results_dir,input_gwas,clump,name=''){
  
  # Name (optional)
  if (name != ''){header <- paste0('_',name)}
  
  # STEP 3: Generate list of rsIDs from clump data
  clump_number <- ifelse(is.na(str_extract(clump,'\\d+')),1,str_extract(clump,'\\d+'))
  clump <- fread(paste0(results_dir,clump))
  gwas <- fread(input_gwas)
  gwas[,':='(ID = paste0(CHR,':',BP,':',ALLELE0,':',ALLELE1))]
  setkey(gwas,ID); setkey(clump,SNP)
  clump[gwas,':='(rsID = i.SNP)]
  id <- clump$rsID
  write.table(id,file=paste0(results_dir,'clump',clump_number,header,'_rsid.tsv'),sep='\t',row.names=F,quote=F,col.names=F)
  
  # STEP 4: Format summary GWAS for score function
  bolt_clumped <- gwas[SNP %in% id]
  
  # Isolate duplicates
  duplicated <- bolt_clumped[SNP %in% bolt_clumped[duplicated(SNP)]$SNP]
  duplicated[,BETA := lapply(.SD,mean),by='SNP',.SDcols='BETA']
  special <- unique(duplicated,by='SNP')
  
  # Add back to master dataset
  bolt_clumped <- bolt_clumped[!(SNP %in% bolt_clumped[duplicated(SNP)]$SNP)]
  bolt_clumped <- rbind(bolt_clumped,special)
  
  # Reduce to necessary rows
  out <- bolt_clumped[,c('SNP','ALLELE1','BETA')]
  write.table(out,file=paste0(results_dir,'bolt',header,'_clumped',clump_number,'.tsv'),sep='\t',row.names=F,quote=F)
}

post_clump(results_dir='/medpop/afib/skhurshid/acceleration_phewas/prs/clumps/',
          input_gwas='/medpop/afib/skhurshid/acceleration_phewas/gwas/v24_mvpa_rate_bouted_sqrt_sensi/v24_mvpa_rate_bouted_sqrt_sensi_bolt.imputed.filtered.tsv',
          clump='out_mvpa_sqrt_sensi_prs4.clumped',name='mvpa_sqrt_sensi')