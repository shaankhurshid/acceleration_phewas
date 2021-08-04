# # Script to create incidence waterfalls using pre-processed data

# Depends
library(data.table)
library(viridis)
library(stringr)
library(fdrtool)
library(RColorBrewer)

# Loads
dict <- fread(file='~/Documents/MGH Research/accel_phewas/phecode_definitions/phecode_definitions1.2.csv')
results <- fread(file='~/Documents/MGH Research/accel_phewas/phecode_outputs/cox_value_bmi.csv')

# Load color correspondences from Fig 2
col_corr <- fread(file='~/Documents/MGH Research/accel_phewas/fig2_col_list.csv')

# Add additional colors where needed
## Choose palette
n <- 6
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
set.seed(42)
col=sample(col_vector, n)
extra <- data.table(all_categories=c('congenital anomalies','hematopoietic','mental disorders',
                                     'infectious diseases','neurological','pregnancy complications'),
                    col=col)
col_corr <- rbind(col_corr,extra)
setnames(col_corr,'col','cat_col')
setkey(col_corr,all_categories)

# Save color palette
#write.csv(col_corr,file='~/Documents/MGH Research/accel_phewas/col_corr.csv',row.names=F)

# Apply post-hoc event filter
results <- results[n_events>=20]

# Convert output format to phecode
results[,phecode := as.numeric(str_remove(str_replace_all(disease,'\\_','\\.'),'phecode\\.'))]

# Join on phecode to get names
setkey(dict,phecode); setkey(results,phecode)
results[dict,':='(phenotype = i.phenotype,
                  category = i.category)]

# Remove NAs
results <- results[!is.na(hr)]

# Have to replace zeros with a really small number
results[p==0]$p <- 1*10^-10

fdr_out <- fdrtool(x=results$p,statistic='pvalue')

results[which(fdr_out$qval<=0.01),sig := 1]
results[which(fdr_out$qval<=0.05),sig_lax := 1]

# Set up plotting variables
setkey(results,p)
correction <- nrow(results)

# Unite miscellaneous category
results[is.na(category) | category=='NULL']$category <- 'other'

# Create dummy variable to control category order
results[,category_order := ifelse(category=='other','a_other',category)]

# color
setkey(results,category); setkey(col_corr,all_categories)
results[col_corr,col := i.cat_col]

# shape
results[,shape := ifelse(hr < 1,6,2)]

# Sort by category then disease
setkeyv(results,c('category_order','phenotype'))

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
pdf(file='~/Documents/MGH Research/accel_phewas/phecode_plots/cox_value_bmi_fdr.pdf',height=5,width=12,pointsize=5)
par(mar=c(15,3,1,6),oma=c(1,1,1,1))

plot(x=1:nrow(results),y=-log10(results$p),col=ifelse(!is.na(results$sig),results$col,paste0(results$col,'4D')),
     bty='n',xaxt='n',yaxt='n',xlim=c(0,nrow(results)),
     ylim=c(0,20),xlab='',ylab='',cex=2,pch=results$shape)

axis(1,at=x_locs$x_coord,cex.axis=2,labels=rep('',length(x_locs$x_coord)))
axis(2,cex.axis=2,at=seq(0,20,5),las=2,pos=-12)

mtext("-log(p)",2,line=1.5,cex=2)

segments(0,-log10(max(results[sig==1]$p)),nrow(results),-log10(max(results[sig==1]$p)),col='#bd0026',lty=5)

text(x = x_locs$x_coord,
     y = par("usr")[3] - 0.8,
     labels = str_remove_all(x_locs$category,'[A-z]\\_'),
     xpd = NA,
     srt = -45,
     adj = 0,
     cex = 2)

dev.off()

# Plot effect sizes
setkey(results,hr)

# color
sig_results <- results[sig==1]
sig_results$color <- rep('darkgray',nrow(sig_results))
color <- viridis(nrow(sig_results[sig==1]))
for (i in 1:length(color)){
  sig_results$color[i] <- color[i]
}

pdf(file='~/Documents/MGH Research/accel_phewas/phecode_plots/cox_value_bmi_effect_fdr.pdf',height=5,width=20,pointsize=5)
par(mar=c(18,1,1,0),oma=c(1,1,1,1))

plot(x=1:nrow(sig_results),y=sig_results$hr,col=sig_results$color,
     bty='n',xaxt='n',yaxt='n',pch=19,xlim=c(0,nrow(sig_results)),
     ylim=c(0,1.5),xlab='',ylab='',cex=1.2)

axis(1,at=1:nrow(sig_results),cex.axis=2,labels=rep('',nrow(sig_results)))
axis(2,cex.axis=2,at=seq(0,1.5,0.3),las=2,pos=0.6)

mtext("Hazard ratio",2,line=-1,cex=2.5)

segments(0.6,1,nrow(sig_results),1,col='black',lty=2)

segments(1:nrow(sig_results),sig_results$lower,1:nrow(sig_results),sig_results$upper,col=sig_results$color,
         lwd=1.5)

text(x = 1:nrow(sig_results),
     y = par("usr")[3] - 0.05,
     labels = sig_results$phenotype,
     xpd = NA,
     srt = -45,
     adj = 0,
     cex = 0.6)

dev.off()

write.csv(results,'~/Documents/MGH Research/accel_phewas/phecode_outputs/cox_value_bmi_processed_fdr.csv')

