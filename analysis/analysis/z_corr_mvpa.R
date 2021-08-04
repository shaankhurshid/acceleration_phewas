library(fdrtool)

# Load summary results
prs_result <- fread(file='~/Documents/MGH Research/accel_phewas/phecode_outputs/cox_mvpa_sqrt_prs_processed_fdr.csv')
mvpa_result <- fread(file='~/Documents/MGH Research/accel_phewas/phecode_outputs/cox_mvpa_bmi_processed_fdr.csv')

# Beta column in mvpa set
mvpa_result[,beta := log(hr)]

# Loop over significance levels
corr_threshold <- function(data1,data2,statistic='p',cuts){
  out <- list(); n <- 1
  fdr <- fdrtool(x=data1[,get(statistic)],statistic='pvalue')
  for (i in cuts){
    subset <- data1[which(fdr$qval<=i)]
    setkey(subset,phecode); setkey(data2,phecode)
    compare_set <- subset[data2,nomatch=0]
    corr <- cor.test(compare_set$beta,compare_set$i.beta)
    out[[n]] <- c(i,nrow(compare_set),max(subset$p),corr$estimate,corr$conf.int[1],corr$conf.int[2])
    n <- n + 1
  }
output <- data.table(do.call(rbind,out))
names(output) <- c('fdr','n_diseases','max_p','corr','lower','upper')
return(output)
}

corr_threshold_sp <- function(data1,data2,statistic='p',cuts){
  out <- list(); n <- 1
  fdr <- fdrtool(x=data1[,get(statistic)],statistic='pvalue')
  for (i in cuts){
    subset <- data1[which(fdr$qval<=i)]
    setkey(subset,phecode); setkey(data2,phecode)
    compare_set <- subset[data2,nomatch=0]
    corr <- cor.test(subset$beta,compare_set$i.beta,method='spearman')
    out[[n]] <- c(i,nrow(compare_set),max(subset$p),corr$estimate)
    n <- n + 1
  }
  output <- data.table(do.call(rbind,out))
  names(output) <- c('fdr','n_diseases','max_p','corr')
  return(output)
}

# Loop over significance levels
corr_p_threshold <- function(data1,data2,statistic='p',cuts){
  out <- list(); n <- 1
  fdr <- fdrtool(x=data1[,get(statistic)],statistic='pvalue')
  for (i in cuts){
    subset <- data1[which(fdr$qval<=i)]
    setkey(subset,phecode); setkey(data2,phecode)
    compare_set <- subset[data2,nomatch=0]
    corr <- cor.test(subset$p,compare_set$i.p)
    out[[n]] <- c(i,nrow(compare_set),max(subset$p),corr$estimate,corr$conf.int[1],corr$conf.int[2])
    n <- n + 1
  }
  output <- data.table(do.call(rbind,out))
  names(output) <- c('fdr','n_diseases','max_p','corr','lower','upper')
  return(output)
}

corrs = corr_threshold(data1=mvpa_result,data2=prs_result,cuts=c(0.30,0.25,0.20,0.15,0.10,0.05,0.01,0.005,
                                                                    0.001,5e-4,1e-4,5e-5,
                                                                    1e-5,5e-6,1e-6,5e-7,1e-7,5e-8,1e-8,
                                                                    5e-9,1e-9,5e-10,1e-10))
corrs2 = corr_threshold(data1=mvpa_result,data2=prs_result,cuts=c(0.30,0.25,0.20,0.15,0.10,0.05,0.01,0.005,
                                                                  0.001,5e-4,1e-4,5e-5,
                                                                  1e-5,5e-6,1e-6,1e-7,1e-8,1e-9,1e-10,1e-11,1e-12))
corrs_p = corr_p_threshold(data1=mvpa_result,data2=prs_result,cuts=c(0.30,0.25,0.20,0.15,0.10,0.05,0.01,0.005,
                                                                  0.001,5e-4,1e-4,5e-5,
                                                                  1e-5,5e-6,1e-6,5e-7,1e-7,5e-8,1e-8,
                                                                  5e-9,1e-9,5e-10,1e-10))
corrs_p2 = corr_p_threshold(data1=mvpa_result,data2=prs_result,cuts=c(0.30,0.25,0.20,0.15,0.10,0.05,0.01,0.005,
                                                                   0.001,5e-4,1e-4,5e-5,
                                                                   1e-5,5e-6,1e-6))
corrs_sp = corr_threshold_sp(data1=mvpa_result,data2=prs_result,cuts=c(0.30,0.25,0.20,0.15,0.10,0.05,0.01,0.005,
                                                                  0.001,5e-4,1e-4,5e-5,
                                                                  1e-5,5e-6,1e-6,5e-7,1e-7,5e-8,1e-8,
                                                                  5e-9,1e-9,5e-10,1e-10))
# Plot
pdf(file='~/Documents/MGH Research/accel_phewas/fdr_corr_ct.pdf',height=4,width=6,pointsize=5)
par(mar=c(6,5,1,1),oma=c(1,1,1,1))
plot(x=0:(nrow(corrs)-1),y=corrs$corr,xlim=c(0,nrow(corrs)-1),ylim=c(-0.2,1),col='#000000',
     pch=19,xlab='',ylab='',xaxt='n',yaxt='n',bty='n',cex=1.5)
axis(1,at=seq(0,nrow(corrs)-1),cex.axis=1.5,labels=rep('',nrow(corrs)))
axis(2,at=seq(-0.2,1,0.1),cex.axis=1.5,las=2,pos=-0.5)
mtext("FDR threshold",1,cex=1.8,line=5.5)
mtext("Correlation of betas for mean versus inferred acceleration",2,cex=1.8,line=3.8)
text(x=0:(nrow(corrs)-1),y=corrs$upper+0.07,labels=as.character(corrs$n_disases),cex=1.2)
polygon(x=c(-0.2,1:(nrow(corrs)-2),(nrow(corrs)-0.5),(nrow(corrs)-0.5),(nrow(corrs)-2):1,-0.2),
        y=c(corrs$upper,rev(corrs$lower)),col="#b3e2cd4D",border=NA)
text(x = 0:(nrow(corrs)-1),
     y = par("usr")[3] - 0.05,
     labels = as.character(corrs$fdr),
     xpd = NA,
     srt = -45,
     adj = 0,
     cex = 1.5)

dev.off()

# Plot
pdf(file='~/Documents/MGH Research/accel_phewas/fdr_corr2_ct.pdf',height=4,width=9,pointsize=5)
par(mar=c(6,5,1.5,1),oma=c(1,1,1,1),xpd=TRUE)
plot(x=0:(nrow(corrs2)-1),y=corrs2$corr,xlim=c(0,nrow(corrs2)-1),ylim=c(0,1),col='#000000',
     pch=19,xlab='',ylab='',xaxt='n',yaxt='n',bty='n',cex=1.5)
axis(1,at=seq(0,nrow(corrs2)-1),cex.axis=1.5,labels=rep('',nrow(corrs2)))
axis(2,at=seq(0,1,0.1),cex.axis=1.8,las=2)
mtext("FDR threshold",1,cex=2,line=6)
mtext("Correlation of betas for measured vs inferred MVPA",2,cex=1.8,line=4.5)
text(x=0:(nrow(corrs2)-1),y=corrs2$upper+0.07,labels=as.character(corrs2$n_diseases),cex=1.6)
polygon(x=c(-0.2,1:(nrow(corrs2)-2),(nrow(corrs2)-0.5),(nrow(corrs2)-0.5),(nrow(corrs2)-2):1,-0.2),
        y=c(corrs2$upper,rev(corrs2$lower)),col="#fdcdac4D",border=NA)
text(x = 0:(nrow(corrs2)-1),
     y = par("usr")[3] - 0.04,
     labels = as.character(corrs2$fdr),
     xpd = NA,
     srt = -45,
     adj = 0,
     cex = 1.6)
dev.off()

# Plot
pdf(file='~/Documents/MGH Research/accel_phewas/fdr_corr_p.pdf',height=4,width=10,pointsize=5)
par(mar=c(5,5,1,1),oma=c(1,1,1,1))
plot(x=0:(nrow(corrs_p)-1),y=corrs_p$corr,xlim=c(0,nrow(corrs_p)-1),ylim=c(0,1),col='darkgray',
     pch=19,xlab='',ylab='',xaxt='n',yaxt='n',bty='n',cex=1.5)
axis(1,at=seq(0,nrow(corrs_p)-1),cex.axis=1.3,labels=as.character(corrs_p$fdr))
axis(2,at=seq(0,1,0.1),cex.axis=1.3,las=2)
mtext("FDR threshold",1,cex=1.5,line=3.2)
mtext("Correlation of betas for mean versus inferred acceleration",2,cex=1.5,line=4)
text(x=0:(nrow(corrs_p)-1),y=corrs_p$upper+0.05,labels=as.character(corrs_p$n_diseases),cex=1.2)
segments(0:(nrow(corrs_p)-1),corrs_p$lower,0:(nrow(corrs_p)-1),corrs_p$upper,col='darkgray')
dev.off()

# Plot
pdf(file='~/Documents/MGH Research/accel_phewas/fdr_corr_p2.pdf',height=4,width=10,pointsize=5)
par(mar=c(5,5,1,1),oma=c(1,1,1,1))
plot(x=0:(nrow(corrs_p2)-1),y=corrs_p2$corr,xlim=c(0,nrow(corrs_p2)-1),ylim=c(0,1),col='darkgray',
     pch=19,xlab='',ylab='',xaxt='n',yaxt='n',bty='n',cex=1.5)
axis(1,at=seq(0,nrow(corrs_p2)-1),cex.axis=1.3,labels=as.character(corrs_p2$fdr))
axis(2,at=seq(0,1,0.1),cex.axis=1.3,las=2)
mtext("FDR threshold",1,cex=1.5,line=3.2)
mtext("Correlation of betas for mean versus inferred acceleration",2,cex=1.5,line=4)
text(x=0:(nrow(corrs_p2)-1),y=corrs_p2$upper+0.05,labels=as.character(corrs_p2$n_diseases),cex=1.2)
segments(0:(nrow(corrs_p2)-1),corrs_p2$lower,0:(nrow(corrs_p2)-1),corrs_p2$upper,col='darkgray')
dev.off()