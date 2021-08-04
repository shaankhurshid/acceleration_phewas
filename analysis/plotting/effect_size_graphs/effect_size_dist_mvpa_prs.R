# Script to generate visualizations for within-group effect sizes

# Dependences
library(data.table)
library(viridis)

# Load processed phewas results
value_bmi <- fread('~/Documents/MGH Research/accel_phewas/phecode_outputs/cox_mvpa_prs_processed_fdr.csv')

# Get median HR for each group
value_bmi[,group_median_hr := lapply(.SD,median,na.rm=T),.SDcols='hr',by='category']
setkey(value_bmi,group_median_hr)

# Color
color  <- viridis(length(unique(value_bmi$category)))

### BOX AND WHISKER PLOTS
## ALL
# Plot
pdf('~/Documents/MGH Research/accel_phewas/phecode_plots/effect_size_box_mvpa_prs.pdf',pointsize=6,
    height=5,width=6.8)
par(oma=c(1,1,1,1))
par(mar=c(12,4,2,8))

# Plot
boxplot(formula = value_bmi$hr ~ value_bmi$group_median_hr,
        col=color,xlab='',ylab='',xaxt='n',yaxt='n',frame=F,
        cex=1.4,outline=FALSE,ylim=c(0.8,1.2),xlim=c(0,length(unique(value_bmi$category))),
        at=1:length(unique(value_bmi$category)))

# Axes
axis(1,cex.axis=1.8,at=1:length(value_bmi$category),labels = NA)
axis(2,cex.axis=1.8,las=2,pos=-0.2,at=seq(0.8,1.2,0.2))
mtext(side=2,'Hazard ratio',cex=2,line=2.8)

text(x = 1:length(unique(value_bmi$category)),
     y = par("usr")[3] - 0.02,
     labels = unique(value_bmi$category),
     xpd = NA,
     srt = -45,
     adj = 0,
     cex = 1.3)

dev.off()