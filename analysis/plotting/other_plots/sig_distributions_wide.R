# Script to make stacked barplots

# Depends
library(data.table)
library(RColorBrewer)
library(Cairo)

# Load processed phewas results
cutoff_bmi <- fread(file='~/Documents/MGH Research/accel_phewas/phecode_outputs/cox_cutoff_bmi_processed_fdr.csv')
self_bmi <- fread(file='~/Documents/MGH Research/accel_phewas/phecode_outputs/cox_self_mvpa_cutoff_processed_fdr.csv')
prs <- fread(file='~/Documents/MGH Research/accel_phewas/phecode_outputs/cox_mvpa_sqrt_cutoff_prs_processed_fdr.csv')

# Number of significant hits per category over total diseases in that category
hits_cutoff <- cutoff_bmi[sig==1 & hr < 1,.N,by='category']; total_cutoff <- value_bmi[,.N,by='category']
hits_self <- self_bmi[sig==1 & hr < 1,.N,by='category']; total_self <- value_bmi[,.N,by='category']
hits_prs <- prs[sig==1 & hr < 1,.N,by='category']; total_prs <- value_bmi[,.N,by='category']

# # Chi-square null distribution test
# setkey(hits,category); setkey(total,category)
# hits[total,total := i.N]
# hits[,not_hit := total-N]
# chisq.test(matrix(rbind(hits$N,hits$not_hit),nrow=2))

# Proceed
setorder(hits_cutoff,-'N'); setorder(hits_self,-'N'); setorder(hits_prs,-'N')

# Generate a parent frame
col_corr <- fread(file='~/Documents/MGH Research/accel_phewas/col_corr.csv')

setkey(col_corr,all_categories); setkey(hits_self,category); setkey(hits_cutoff,category); setkey(hits_prs,category)
col_corr[hits_cutoff,':='(c_hits_cutoff = i.N)]
col_corr[hits_self,':='(d_hits_self = i.N)]
col_corr[hits_prs,':='(e_hits_prs = i.N)]

# Category ordering
col_corr[,category_ordered := ifelse(all_categories=='other','a_other',all_categories)]

# Now melt
melted <- melt(col_corr,id.vars=c('all_categories','cat_col','category_ordered'),
               measure.vars = c('c_hits_cutoff',
                                'd_hits_self','e_hits_prs'))
setkey(melted,'category_ordered','variable')
melted[is.na(value)]$value <- 0

################################## Plot 1: Positive
CairoPDF(file='~/Documents/MGH Research/accel_phewas/pos_wide.pdf',height=6,width=18,
         pointsize=5)
par(oma=c(1,1,1,1),mar=c(21,7,2.5,10),xpd=TRUE)

coords <- barplot(melted$value,space=c(rep(0.05,3),1,rep(c(rep(0.05,2),1),16),rep(0.05,2)),border=NA,
                  col=rep(c('#fdb462','#b3de69','#fccde5'),18),ylim=c(0,70),
                  yaxt='n',xaxt='n')

axis(2,cex.axis=2.5,las=2,at=seq(0,70,5),pos=-1)
mtext("# of Associations (Lower Risk)",side=2,cex=3,line=2.5)

legend(x=max(coords)-15,y=68,legend=c("Accelerometer-Derived",
                                    "Self-Reported","Genetically Inferred"),
       col=c('#fdb462','#b3de69','#fccde5'),lwd=3,bty='n',cex=2.5)

offset=0
text(x = coords[seq(2,53,3)],
     y = par("usr")[3] - 3,
     labels = as.character(unique(melted$all_categories)),
     xpd = NA,
     srt = -45,
     adj = 0,
     cex = 2.5)

dev.off()

#######################################
# Number of significant hits per category over total diseases in that category
hits_cutoff <- cutoff_bmi[sig==1 & hr > 1,.N,by='category']; total_cutoff <- value_bmi[,.N,by='category']
hits_self <- self_bmi[sig==1 & hr > 1,.N,by='category']; total_self <- value_bmi[,.N,by='category']
hits_prs <- prs[sig==1 & hr > 1,.N,by='category']; total_prs <- value_bmi[,.N,by='category']

setorder(hits_cutoff,-'N')
setorder(hits_self,-'N'); setorder(hits_prs,-'N')

# Generate a parent frame
col_corr <- fread(file='~/Documents/MGH Research/accel_phewas/col_corr.csv')

setkey(col_corr,all_categories); setkey(hits_self,category); setkey(hits_cutoff,category); setkey(hits_prs,category)
col_corr[hits_cutoff,':='(c_hits_cutoff = i.N)]
col_corr[hits_self,':='(d_hits_self = i.N)]
col_corr[hits_prs,':='(e_hits_prs = i.N)]

# Category ordering
col_corr[,category_ordered := ifelse(all_categories=='other','a_other',all_categories)]

# Now melt
melted <- melt(col_corr,id.vars=c('all_categories','cat_col','category_ordered'),
               measure.vars = c('c_hits_cutoff',
                                'd_hits_self','e_hits_prs'))
setkey(melted,'category_ordered','variable')
melted[is.na(value)]$value <- 0

################################## Plot 2: Negative
CairoPDF(file='~/Documents/MGH Research/accel_phewas/neg_wide.pdf',height=6,width=18,
         pointsize=5)
par(oma=c(1,1,1,1),mar=c(21,7,2.5,10),xpd=TRUE)

coords <- barplot(melted$value,space=c(rep(0.05,3),1,rep(c(rep(0.05,2),1),16),rep(0.05,2)),border=NA,
                  col=rep(c('#fdb462','#b3de69','#fccde5'),18),ylim=c(20,0),
                  yaxt='n',xaxt='n')

axis(2,cex.axis=2.5,las=2,at=seq(20,0,-5),pos=-1)
mtext("# of Associations (Higher Risk)",side=2,cex=3,line=2.5)

# offset=0
# text(x = coords[seq(2,53,3)],
#      y = par("usr")[3] - 1.5,
#      labels = as.character(unique(melted$all_categories)),
#      xpd = NA,
#      srt = -45,
#      adj = 0,
#      cex = 2.5)

dev.off()
