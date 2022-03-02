# Script to make stacked barplots

# Depends
library(data.table)
library(RColorBrewer)
library(Cairo)

# Load processed phewas results
value_bmi <- fread(file='~/Documents/MGH Research/accel_phewas/phecode_outputs/cox_mvpa_bmi_processed_fdr_covar.csv')
self_bmi <- fread(file='~/Documents/MGH Research/accel_phewas/phecode_outputs/cox_self_bmi_processed_fdr_covar.csv')

# Number of significant hits per category over total diseases in that category
hits_cutoff <- value_bmi[sig==1 & hr < 1,.N,by='category']; total_cutoff <- value_bmi[,.N,by='category']
hits_self <- self_bmi[sig==1 & hr < 1,.N,by='category']; total_self <- self_bmi[,.N,by='category']

setkey(hits_cutoff,category); setkey(total_cutoff,category)
setkey(hits_self,category); setkey(total_self,category)

hits_cutoff[total_cutoff,total := i.N]
hits_self[total_self,total := i.N]

hits_cutoff[,fraction := N/total]
hits_self[,fraction := N/total]

setorder(hits_cutoff,-'fraction'); setorder(hits_self,-'fraction')

# Generate a parent frame
col_corr <- fread(file='~/Documents/MGH Research/accel_phewas/col_corr.csv')

setkey(col_corr,all_categories); setkey(hits_self,category); setkey(hits_cutoff,category)
col_corr[hits_cutoff,':='(c_hits_cutoff_fraction = i.fraction)]
col_corr[hits_self,':='(d_hits_self_fraction = i.fraction)]

# Category ordering
col_corr[,category_ordered := ifelse(all_categories=='other','a_other',all_categories)]

# Now melt
melted <- melt(col_corr,id.vars=c('all_categories','cat_col','category_ordered'),
               measure.vars = c('c_hits_cutoff_fraction',
                                'd_hits_self_fraction'))
setkey(melted,'category_ordered','variable')
melted[is.na(value)]$value <- 0

################################## Plot 1: Positive
CairoPDF(file='~/Documents/MGH Research/accel_phewas/pos_fraction_wide.pdf',height=6,width=18,
         pointsize=5)
par(oma=c(1,1,1,1),mar=c(22,7,2.5,10),xpd=TRUE)

coords <- barplot(melted$value,space=c(rep(0.05,2),1,rep(c(rep(0.05,1),1),16),rep(0.05,1)),border=NA,
                  col=rep(c('#fdb462','#b3de69'),18),ylim=c(0,0.6),
                  yaxt='n',xaxt='n')

axis(2,cex.axis=3,las=2,at=seq(0,0.6,0.1),pos=-1,labels=as.character(0:6))
mtext("% of Category Associated (Lower Risk)",2,cex=3.2,line=3)

legend(x=max(coords)-7,y=0.65,legend=c("Accelerometer-Derived",
                                       "Self-Reported"),
       col=c('#fdb462','#b3de69'),lwd=3,bty='n',cex=3.2)

offset=0
label_coord <- apply(rbind(coords[seq(1,35,2)],coords[seq(2,36,2)]),FUN=mean,MARGIN=2)
text(x = label_coord,
     y = par("usr")[3] - 0.03,
     labels = as.character(unique(melted$all_categories)),
     xpd = NA,
     srt = -45,
     adj = 0,
     cex = 3)

dev.off()

#######################################
# Number of significant hits per category over total diseases in that category
hits_cutoff <- cutoff_bmi[sig==1 & hr > 1,.N,by='category']; total_cutoff <- value_bmi[,.N,by='category']
hits_self <- self_bmi[sig==1 & hr > 1,.N,by='category']; total_self <- value_bmi[,.N,by='category']

setkey(hits_cutoff,category); setkey(total_cutoff,category)
setkey(hits_self,category); setkey(total_self,category)

hits_cutoff[total_cutoff,total := i.N]
hits_self[total_self,total := i.N]

hits_cutoff[,fraction := N/total]
hits_self[,fraction := N/total]

setorder(hits_cutoff,-'fraction'); setorder(hits_self,-'fraction')

# Generate a parent frame
col_corr <- fread(file='~/Documents/MGH Research/accel_phewas/col_corr.csv')

setkey(col_corr,all_categories); setkey(hits_self,category); setkey(hits_cutoff,category)
col_corr[hits_cutoff,':='(c_hits_cutoff_fraction = i.fraction)]
col_corr[hits_self,':='(d_hits_self_fraction = i.fraction)]

# Category ordering
col_corr[,category_ordered := ifelse(all_categories=='other','a_other',all_categories)]

# Now melt
melted <- melt(col_corr,id.vars=c('all_categories','cat_col','category_ordered'),
               measure.vars = c('c_hits_cutoff_fraction',
                                'd_hits_self_fraction'))
setkey(melted,'category_ordered','variable')
melted[is.na(value)]$value <- 0

################################## Plot 2: Negative
CairoPDF(file='~/Documents/MGH Research/accel_phewas/neg_fraction_wide.pdf',height=6,width=18,
         pointsize=5)
par(oma=c(1,1,1,1),mar=c(21,7,2.5,10),xpd=TRUE)

coords <- barplot(melted$value,space=c(rep(0.05,2),1,rep(c(rep(0.05,1),1),16),rep(0.05,1)),border=NA,
                  col=rep(c('#fdb462','#b3de69'),18),ylim=c(0.5,0),
                  yaxt='n',xaxt='n')

axis(2,cex.axis=3,las=2,at=seq(0.5,0,-0.1),pos=-1,labels=as.character(5:0))
mtext("% of Category Associated (Higher Risk)",2,cex=3.2,line=2.5)

# legend(x=max(coords)-7,y=0.42,legend=c("Accelerometer-Derived",
#                                        "Self-Reported","Genetically Inferred"),
#        col=c('#fdb462','#b3de69','#fccde5'),lwd=3,bty='n',cex=2)

offset=0
# text(x = coords[seq(2,53,3)],
#      y = par("usr")[3] - 0.03,
#      labels = as.character(unique(melted$all_categories)),
#      xpd = NA,
#      srt = -45,
#      adj = 0,
#      cex = 2.5)

dev.off()
