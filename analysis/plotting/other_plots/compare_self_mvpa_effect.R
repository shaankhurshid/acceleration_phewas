# Depends
library(Cairo)

# Load processed results
mvpa_cutoff <- fread(file='~/Documents/MGH Research/accel_phewas/phecode_outputs/cox_cutoff_bmi_processed_fdr_covar.csv')
mvpa <- fread(file='~/Documents/MGH Research/accel_phewas/phecode_outputs/cox_self_mvpa_cutoff_processed_fdr_covar.csv')

# Color list
col_corr <- fread(file='~/Documents/MGH Research/accel_phewas/col_corr.csv')

# Isolate to common diseases
setkey(mvpa_cutoff,phecode); setkey(mvpa,phecode)
both <- mvpa_cutoff[mvpa,nomatch=0]

# Isolate to significant for both
both_sig <- both[p<0.05 & i.p<0.05 & hr<1 & i.hr<1]

# Now obtain median HR per category
med_mvpa <- both_sig[,lapply(.SD,median),.SDcols='hr',by='category']
med_self <- both_sig[,lapply(.SD,median),.SDcols='i.hr',by='category']

# Combine
setkey(med_mvpa,category); setkey(med_self,category)
med_mvpa[med_self,self_hr := i.i.hr]

# Add color indicator
setkey(col_corr,all_categories)
med_mvpa[col_corr,col := i.cat_col]

# Set order
setorder(med_mvpa,order=-category)

# Plot
CairoPDF(file='~/Documents/MGH Research/accel_phewas/compare_self_mvpa_covar.pdf',height=4,width=5,pointsize=5)
par(oma=c(1,1,1,1),mar=c(5,13,1,1),xpd=TRUE)
plot(x=c(med_mvpa$hr,med_mvpa$self_hr),y=rep(seq(1,nrow(med_mvpa)*2,2),2),col=rep(med_mvpa$col,2),bty='n',
     yaxt='n',xaxt='n',ylab='',xlab='',pch=c(rep(19,nrow(med_mvpa)),rep(17,nrow(med_mvpa))),
     cex=c(rep(1.5,nrow(med_mvpa)),rep(2.2,nrow(med_mvpa))),xlim=c(0.5,1.1))
axis(1,at=seq(0.5,1.1,0.1),cex.axis=1.5)
axis(2,at=seq(1,nrow(med_mvpa)*2,2),labels=paste0(med_mvpa$category),las=2,cex.axis=1.5)
segments(1,0,1,nrow(med_mvpa)*2,lty=5)
segments(0.49,seq(3,31,4),
         1.085,seq(3,31,4),col=rgb(0,0,0,alpha=0.05),
         lty=1,lend=2,lwd=13)
mtext("Median within-category hazard ratio",1,cex=1.8,line=3.5)
legend(x=1.01,y=33,legend=c('Measured','Self-reported'),pch=c(19,17),
       pt.cex=c(1.3,2),cex=1.2,col='darkgray',y.intersp=1.3,
       box.lwd=0.1,box.col='gray',bg=rgb(1,1,1,alpha=1))
dev.off()


