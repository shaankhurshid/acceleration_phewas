# Script to create incidence waterfalls using pre-processed data

# Depends
library(data.table)
library(viridis)
library(stringr)
library(fdrtool)

# Loads
dict <- fread(file='~/Documents/MGH Research/accel_phewas/phecode_definitions/phecode_definitions1.2.csv')
results <- fread(file='~/Documents/MGH Research/accel_phewas/phecode_outputs/cox_mvpa_bmi_quintile_full.csv')

# Load color correspondences from Fig 2
col_corr <- fread(file='~/Documents/MGH Research/accel_phewas/col_corr.csv')

# Apply post-hoc event filter
results <- results[n_events>=50]

# Convert output format to phecode
results[,phecode := as.numeric(str_remove(str_replace_all(disease,'\\_','\\.'),'phecode\\.'))]

# Join on phecode to get names
setkey(dict,phecode); setkey(results,phecode)
results[dict,':='(phenotype = i.phenotype)]

# Limit to significant results on primary test
primary_results <- fread(file='~/Documents/MGH Research/accel_phewas/phecode_outputs/cox_mvpa_bmi_processed_fdr_covar.csv')
sigs <- primary_results[sig==1]
results <- results[phecode %in% sigs$phecode]
setkey(sigs,phecode)
results[sigs,':='(category = i.category,
                  category_order = i.category_order,
                  col = i.col)]

# Eliminate non-sensical hazard ratios (beta did not converge)
for (j in c('d2','d3','d4','d5')){
  set(results,i=which(results[[j]] < 0.01 | results[[j]] > 10),j=j,value=999)
}

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

gap = 10
int1 <- 1:x_counts$N[1]
int2 <- (gap+max(int1)):(gap+max(int1)+x_counts$N[2]-1)
int3 <- (gap+max(int2)):(gap+max(int2)+x_counts$N[3]-1)
int4 <- (gap+max(int3)):(gap+max(int3)+x_counts$N[4]-1)
int5 <- (gap+max(int4)):(gap+max(int4)+x_counts$N[5]-1)
int6 <- (gap+max(int5)):(gap+max(int5)+x_counts$N[6]-1)
int7 <- (gap+max(int6)):(gap+max(int6)+x_counts$N[7]-1)
int8 <- (gap+max(int7)):(gap+max(int7)+x_counts$N[8]-1)
int9 <- (gap+max(int8)):(gap+max(int8)+x_counts$N[9]-1)
int10 <- (gap+max(int9)):(gap+max(int9)+x_counts$N[10]-1)
int11 <- (gap+max(int10)):(gap+max(int10)+x_counts$N[11]-1)
int12 <- (gap+max(int11)):(gap+max(int11)+x_counts$N[12]-1)
int13 <- (gap+max(int12)):(gap+max(int12)+x_counts$N[13]-1)
int14 <- (gap+max(int13)):(gap+max(int13)+x_counts$N[14]-1)
int15 <- (gap+max(int14)):(gap+max(int14)+x_counts$N[15]-1)
int16 <- (gap+max(int15)):(gap+max(int15)+x_counts$N[16]-1)
int17 <- (gap+max(int16)):(gap+max(int16)+x_counts$N[17]-1)

x_locs$x_coord_new <- c(mean(int1),mean(int2),mean(int3),mean(int4),mean(int5),mean(int6),mean(int7),
                        mean(int8),mean(int9),mean(int10),mean(int11),mean(int12),mean(int13),mean(int14),
                        mean(int15),mean(int16),mean(int17))

# Plot p-values
pdf(file='~/Documents/MGH Research/accel_phewas/phecode_plots/cox_mvpa_bmi_quintile_full.pdf',height=5,width=12,pointsize=5)
par(mar=c(15.6,10,1,10),oma=c(1,1,1,1))

y_vals <- c(results$d2,results$d3,results$d4,results$d5)
col_vals <- c(rep('#d7191c8C',355),rep('#fdae618C',355),rep('#a6d96a8C',355),rep('#1a96418C',355))
x_vals <- c(int1,int2,int3,int4,int5,int6,int7,int8,int9,int10,
            int11,int12,int13,int14,int15,int16,int17)
pch_vals <- ifelse(y_vals<1,25,24)

plot(x=rep(x_vals,4),y=y_vals,col=col_vals,bg=col_vals,
     bty='n',xaxt='n',yaxt='n',xlim=c(0,max(int17)),
     ylim=c(0,1.5),xlab='',ylab='',cex=2.2,pch=19)

axis(1,at=x_locs$x_coord_new,cex.axis=2.4,labels=rep('',length(x_locs$x_coord_new)))
axis(2,cex.axis=2.4,at=c(0,0.5,1,1.5),las=2,pos=-8)

segments(-8,1,max(int17),1,col='black',lty=5)

mtext("Hazard ratio versus lowest MVPA quintile",2,line=5,cex=2.4)

text(x = x_locs$x_coord_new,
     y = par("usr")[3] - 0.07,
     labels = str_remove_all(x_locs$category,'[A-z]\\_'),
     xpd = NA,
     srt = -45,
     adj = 0,
     cex = 2.2)

par(xpd=TRUE)
legend(x=max(int17)+10,y=1.3,legend=c('','','',''),pch=19,cex=2.4,y.intersp = 2,
       col=c('#d7191c8C','#fdae618C','#a6d96a8C','#1a96418C'),bty='n')
dev.off()

write.csv(results,file='~/Documents/MGH Research/accel_phewas/phecode_outputs/cox_mvpa_bmi_quintile_full_covar_processed.csv',row.names=F)