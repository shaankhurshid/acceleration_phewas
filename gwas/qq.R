# Script to generate QQ plots manually

# Depends
library(data.table)

# Plotting function
qq_plot <- function(input_gwas,output_path,pval_column='P_BOLT_LMM',
                    chr_column=NULL,pos_column=NULL,label_chr=NULL,label_pos=NULL){
  gwas <- fread(input_gwas)
  
  if (length(label_pos)==0){
    col <- 'black'
  } else {
    col <- rep('black',nrow(gwas))
    col[which(c(gwas[,get(chr_column)]==label_chr & gwas[,get(pos_column)] >= label_pos[1] 
                & gwas[,get(pos_column)] <= label_pos[2]))] <- 'red'
  }
  
  obv <- gwas[,get(pval_column)]
  exp <- (rank(gwas[,get(pval_column)],ties.method='first')+0.5)/(length(obv)+1)
  ylim <- (max(-log10(obv)) %/% 2) * 2 + 2
  xlim <- (max(-log10(exp)) %/% 2) * 2 + 2
  png(paste0(output_path),width=2000, height=2000, res=300)
  par(mar=c(3,3,1,1),oma=c(3,3,1,1))
  plot(-log10(exp),-log10(obv),col=col,xlim=c(0,xlim),ylim=c(0,ylim),xaxt='n',yaxt='n',xlab='',ylab='')
  abline(0,1,col='red')
  axis(1,at=seq(0,xlim,2),cex.axis=1.2)
  axis(2,at=seq(0,ylim,2),cex.axis=1.2,las=2)
  mtext("Expected -log(p)",1,cex=1.5,line=2.5)
  mtext("Observed -log(p)",2,cex=1.5,line=2.5)
  dev.off()
}

# Load data
qq_plot(input_gwas='/medpop/afib/skhurshid/acceleration_phewas/gwas/v23_mvpa_rate_bouted_sqrt/v23_mvpa_rate_bouted_sqrt_bolt.imputed.filtered.tsv',
        output_path='/medpop/afib/skhurshid/acceleration_phewas/qq_mvpa_sqrt.png')
