# Script to generate visualizations for within-group effect sizes

# Dependences
library(data.table)
library(viridis)

# Load processed phewas results
mvpa_bmi <- fread('~/Documents/MGH Research/accel_phewas/phecode_outputs/cox_mvpa_bmi_processed_fdr_medium.csv')

# Get median HR for each group
mvpa_bmi[,group_median_hr := lapply(.SD,median,na.rm=T),.SDcols='hr',by='category']
setkey(mvpa_bmi,group_median_hr)

### BOX AND WHISKER PLOTS
## ALL
# Plot
pdf('~/Documents/MGH Research/accel_phewas/phecode_plots/effect_size_box_mvpa_bmi_medium.pdf',pointsize=6,
    height=6,width=9)
par(oma=c(1,1,1,1))
par(mar=c(15,4,2,8.5))

# Plot
boxplot(formula = mvpa_bmi$hr ~ mvpa_bmi$group_median_hr,
        col=unique(mvpa_bmi[order(mvpa_bmi$group_median_hr)]$col),xlab='',ylab='',xaxt='n',yaxt='n',frame=F,
        cex=1.4,outline=FALSE,ylim=c(0.2,1.4),xlim=c(0,length(unique(mvpa_bmi$category))),
        at=1:length(unique(mvpa_bmi$category)))

# Axes
axis(1,cex.axis=2,at=1:length(unique(mvpa_bmi$category)),labels = NA)
axis(2,cex.axis=2,las=2,pos=0,at=seq(0.2,1.4,0.2))
mtext(side=2,'Hazard ratio',cex=2,line=2.8)

text(x = 1:length(unique(mvpa_bmi$category)),
     y = par("usr")[3] - 0.04,
     labels = unique(mvpa_bmi$category),
     xpd = NA,
     srt = -45,
     adj = 0,
     cex = 1.8)

segments(-0.2,1,18.5,1,lty=5,col='#bd0026')

dev.off()