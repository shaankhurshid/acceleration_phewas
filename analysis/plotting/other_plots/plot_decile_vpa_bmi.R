# # Script to create incidence waterfalls using pre-processed data

# Depends
library(data.table)
library(viridis)
library(stringr)
library(fdrtool)

# Loads
dict <- fread(file='~/Documents/MGH Research/accel_phewas/phecode_definitions/phecode_definitions1.2.csv')
results <- fread(file='~/Documents/MGH Research/accel_phewas/phecode_outputs/cox_vpa_bmi_decile.csv')

# Load color correspondences from Fig 2
col_corr <- fread(file='~/Documents/MGH Research/accel_phewas/col_corr.csv')

# Apply post-hoc event filter
results <- results[n_events>=20]

# Convert output format to phecode
results[,phecode := as.numeric(str_remove(str_replace_all(disease,'\\_','\\.'),'phecode\\.'))]

# Join on phecode to get names
setkey(dict,phecode); setkey(results,phecode)
results[dict,':='(phenotype = i.phenotype)]

# Limit to significant results on primary test
primary_results <- fread(file='~/Documents/MGH Research/accel_phewas/phecode_outputs/cox_vpa_bmi_processed_fdr.csv')
sigs <- primary_results[sig==1]
results <- results[phecode %in% sigs$phecode]
setkey(sigs,phecode)
results[sigs,':='(category = i.category,
                  category_order = i.category_order,
                  col = i.col)]

# Sort by category then randomly within category for better spread
setkeyv(results,c('category_order'))
results <- results[,.SD[sample(.N)],by='category']

# Create vector of midpoints per category for x-axis labels
x_counts <- results[,.N,by='category_order']

midpoint <- function(x){
  if(x %% 2==0){result <- x %/% 2
  } else {result <- (x %/% 2) + 1}
}

x_locs <- data.table(category=x_counts$category_order,x_count=sapply(x_counts$N,midpoint))

out <- rep(NA,length(x_locs$x_count))
for (i in 1:length(x_locs$x_count))
  if (i == 1){out[i] <- x_locs$x_count[i]
  } else {
    out[i] <- x_locs$x_count[i] + cumsum(x_counts$N)[i-1]
  }

x_locs$x_coord <- out

setkey(x_locs,category)

# Plot p-values
pdf(file='~/Documents/MGH Research/accel_phewas/phecode_plots/cox_vpa_bmi_min_decile.pdf',height=5,width=12.5,pointsize=5)
par(mar=c(15.6,10,1,2),oma=c(1,1,1,1))

plot(x=1:nrow(results),y=results$min_decile,col=results$col,
     bty='n',xaxt='n',yaxt='n',xlim=c(0,nrow(results)),
     ylim=c(0.5,10.5),xlab='',ylab='',cex=2.2,pch=19)

axis(1,at=x_locs$x_coord,cex.axis=2.4,labels=rep('',length(x_locs$x_coord)))
axis(2,cex.axis=2.4,at=seq(1,10,1),las=2,pos=-5,labels=rep('',10))

mvpa_labels <- c('<20','20-36','37-54','55-76','77-103',
                 '104-138','139-190','191-278','279-463','>463')

mtext("Weekly VPA With Lowest Risk (min)",2,line=8,cex=2.4)

text(x = par("usr")[1] - 0.5,
     y = 1:10,
     labels = mvpa_labels,
     xpd = NA,
     srt = 30,
     adj = 1,
     cex = 2)

text(x = x_locs$x_coord,
     y = par("usr")[3] - 0.4,
     labels = str_remove_all(x_locs$category,'[A-z]\\_'),
     xpd = NA,
     srt = -45,
     adj = 0,
     cex = 2.2)

dev.off()

# Plot p-values
pdf(file='~/Documents/MGH Research/accel_phewas/phecode_plots/cox_vpa_bmi_max_decile.pdf',height=5,width=12.5,pointsize=5)
par(mar=c(15.6,10,1,2),oma=c(1,1,1,1))

plot(x=1:nrow(results),y=results$max_decile,col=results$col,
     bty='n',xaxt='n',yaxt='n',xlim=c(0,nrow(results)),
     ylim=c(0.5,10.5),xlab='',ylab='',cex=2.2,pch=19)

axis(1,at=x_locs$x_coord,cex.axis=2.4,labels=rep('',length(x_locs$x_coord)))
axis(2,cex.axis=2.4,at=seq(1,10,1),las=2,pos=-5,labels=rep('',10))

mtext("Weekly VPA With Highest Risk (min)",2,line=8,cex=2.4)

mvpa_labels <- c('<20','20-36','37-54','55-76','77-103',
                 '104-138','139-190','191-278','279-463','>463')
  
text(x = par("usr")[1] - 0.5,
     y = 1:10,
     labels = mvpa_labels,
     xpd = NA,
     srt = 30,
     adj = 1,
     cex = 2)

text(x = x_locs$x_coord,
     y = par("usr")[3] - 0.4,
     labels = str_remove_all(x_locs$category,'[A-z]\\_'),
     xpd = NA,
     srt = -45,
     adj = 0,
     cex = 2.2)

dev.off()

write.csv(results,file='~/Documents/MGH Research/accel_phewas/phecode_outputs/cox_vpa_bmi_decile_processed.csv')
