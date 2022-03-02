# # Script to create incidence waterfalls using pre-processed data

# Depends
library(data.table)
library(viridis)
library(stringr)
library(fdrtool)

# Load data
mvpa_bmi <- fread('~/Documents/MGH Research/accel_phewas/phecode_outputs/cox_mvpa_bmi_processed_fdr_covar.csv')
self <- fread('~/Documents/MGH Research/accel_phewas/phecode_outputs/cox_self_bmi_processed_fdr_covar.csv')

# Top 3 conditions per category that are significant for measured
setkey(mvpa_bmi,p)
results <- mvpa_bmi[sig==1,.SD[1:min(3,length(.SD))],by='category']
results <- results[!is.na(n_events)]

# Sort
setkeyv(results,c('category_order','hr'))

# Generate x-axis coordinates
coord <- results[,.N,by='category']

x_coord <- -1
for (j in 1:length(coord$N)){
  group <- coord$N[j]
  for (i in 1:group){
    if (i == 1){start <- max(x_coord) + 1; new <- start + 1; x_coord <- c(x_coord,new)}
    else {start <- max(x_coord); new <- start+1; x_coord <- c(x_coord,new)}
  }
}
x_coord <- x_coord[-1]

# Friendly renames
results[str_detect(phenotype,'Adverse effects of sedatives')]$phenotype <- "Adverse effects of sedatives or other CNS depressants"
results[str_detect(phenotype,'Dizziness and giddiness')]$phenotype <- "Light-headedness and vertigo"
results[str_detect(phenotype,'Myalgia and myositis')]$phenotype <- "Myalgia and myositis"
results[str_detect(phenotype,'Vertiginous')]$phenotype <- "Vertiginous syndromes"
results[str_detect(phenotype,'Congenital anomalies of muscle')]$phenotype <- "Congenital anomalies of connective tissue"
results[str_detect(phenotype,'Otitis media')]$phenotype <- "Otitis media and Eustachian tube disorders"
results[str_detect(phenotype,'Symptoms involving nervous')]$phenotype <- "Nervous and musculoskeletal symptoms"

# Plot of measured
pdf(file='~/Documents/MGH Research/accel_phewas/phecode_plots/cox_mvpa_bmi_effect_zoomed.pdf',height=9,width=12,pointsize=5)
par(mar=c(30,4,1,14.5),oma=c(1,1,1,1),xpd=TRUE)

plot(x=x_coord,y=results$hr,col=results$col,
     bty='n',xaxt='n',yaxt='n',pch=19,xlim=c(0,max(x_coord)),
     ylim=c(0,1),xlab='',ylab='',cex=2.5)

axis(1,at=x_coord,cex.axis=2,labels=rep('',nrow(results)))
axis(2,cex.axis=2.2,at=seq(0,1,0.2),las=2,pos=0)

mtext("Hazard ratio",2,line=1.5,cex=2.8)

segments(0.6,1,max(x_coord)+4,1,col='#bd0026',lty=2)

segments(x_coord,results$lower,x_coord,results$upper,col=results$col,
         lwd=2)

legend(x=30,y=0.2,legend=unique(results$category),col=unique(results$col),ncol=4,
       bty='n',pch=19,cex=1.6)

text(x = x_coord,
     y = par("usr")[3] - 0.02,
     labels = results$phenotype,
     xpd = NA,
     srt = -45,
     adj = 0,
     cex = 1.6)

dev.off()

# Top 3 conditions per category that are significant for measured
setkey(self,p)
results <- self[sig==1,.SD[1:min(3,length(.SD))],by='category']
results <- results[!is.na(n_events)]

# Sort
setkeyv(results,c('category_order','hr'))

# Generate x-axis coordinates
coord <- results[,.N,by='category']

x_coord <- -1
for (j in 1:length(coord$N)){
  group <- coord$N[j]
  for (i in 1:group){
    if (i == 1){start <- max(x_coord) + 1; new <- start + 1; x_coord <- c(x_coord,new)}
    else {start <- max(x_coord); new <- start+1; x_coord <- c(x_coord,new)}
  }
}
x_coord <- x_coord[-1]

# Friendly renames
results[str_detect(phenotype,'Adverse effects of sedatives')]$phenotype <- "Adverse effects of sedatives or other CNS depressants"
results[str_detect(phenotype,'Dizziness and giddiness')]$phenotype <- "Light-headedness and vertigo"
results[str_detect(phenotype,'Myalgia and myositis')]$phenotype <- "Myalgia and myositis"
results[str_detect(phenotype,'Vertiginous')]$phenotype <- "Vertiginous syndromes"
results[str_detect(phenotype,'Congenital anomalies of muscle')]$phenotype <- "Congenital anomalies of connective tissue"
results[str_detect(phenotype,'Otitis media')]$phenotype <- "Otitis media and Eustachian tube disorders"
results[str_detect(phenotype,'Thoracic or lumbosacral neuritis or radiculitis')]$phenotype <- "Thoracic/lumbosacral neuritis, unspecified"

# Plot of self
pdf(file='~/Documents/MGH Research/accel_phewas/phecode_plots/cox_mvpa_prs_effect_zoomed_same.pdf',height=9,width=12,pointsize=5)
par(mar=c(30,4,1,14.5),oma=c(1,1,1,1),xpd=TRUE)

plot(x=x_coord,y=results$hr,col=results$col,
     bty='n',xaxt='n',yaxt='n',pch=19,xlim=c(0,max(x_coord)),
     ylim=c(0.2,1.4),xlab='',ylab='',cex=2.5)

axis(1,at=x_coord,cex.axis=2,labels=rep('',nrow(results)))
axis(2,cex.axis=2.2,at=seq(0.2,1.4,0.2),las=2,pos=0)

mtext("Hazard ratio",2,line=1.5,cex=2.8)

segments(0,1,max(x_coord)+3,1,col='#bd0026',lty=2)

segments(x_coord,results$lower,x_coord,results$upper,col=results$col,
         lwd=2)

text(x = x_coord,
     y = par("usr")[3] - 0.025,
     labels = results$phenotype,
     xpd = NA,
     srt = -45,
     adj = 0,
     cex = 1.6)

legend(x=26,y=0.43,legend=unique(results$category),col=unique(results$col),ncol=4,
       bty='n',pch=19,cex=1.6)

dev.off()

