# Depends
library(data.table)
library(prodlim)
library(Cairo)

#### STEP 1: Identify candidate diseases
# Load results files
cutoff_result <- fread(file='cox_self_mvpa_cutoff_processed_fdr_covar.csv')

# Will choose: 496.2 (Chronic bronchitis), 250 (Type 2 Diabetes), 574 (Cholelithiasis), 428.1 (Heart failure)

# Load seed files
value_bmi <- fread(file='cox_data_self_bmi_covar.csv')

# Load color file
col_corr <- fread(file='col_corr.csv')

# Create cutoff variables
setkey(value_bmi,sample_id)

# Load phecode tables
copd <- read.table(file='phecode_496_2.tab.tsv',sep='\t',header = TRUE); setDT(copd)
gastritis <- read.table(file='phecode_574_1.tab.tsv',sep='\t',header = TRUE); setDT(gastritis)
pad <- read.table(file='phecode_428_1.tab.tsv',sep='\t',header = TRUE); setDT(pad)
dm <- read.table(file='phecode_250.tab.tsv',sep='\t',header = TRUE); setDT(dm)

# Format dates
for (j in (c('enroll_date','phenotype_censor_date'))){set(inferred,j=j,value=as.Date(inferred[[j]],format='%Y-%m-%d'))}
for (j in (c('enroll_date','phenotype_censor_date'))){set(value_bmi,j=j,value=as.Date(value_bmi[[j]],format='%Y-%m-%d'))}

################### Phenotype 1: COPD
setkey(copd,sample_id)
copd_value <- value_bmi[copd,nomatch=0]

# Format variables
copd_value[,censor_date := as.Date(censor_date,format='%Y-%m-%d')]

# Create analysis variables
copd_value[,time_to_copd := ifelse(c(has_disease == 0 | is.na(has_disease)),pmin(as.numeric(censor_date - enroll_date)/365.25,as.numeric(phenotype_censor_date - enroll_date)/365.25),
                                      as.numeric(censor_date - enroll_date)/365.25)]

# Remove prevalent disease or no follow-up
copd_analysis_value <- copd_value[!is.na(time_to_copd) & time_to_copd > 0]

# Graphical variables
copd_analysis_value[,guidelines := factor(ifelse(self_guidelines==1,'Meeting','Not Meeting'),levels=c('Not Meeting','Meeting'))]

# Fit prodlim models
# Models
mod_value <- prodlim(Hist(time_to_copd,has_disease)~guidelines,data=copd_analysis_value)

# Plot
CairoPDF(file='km_copd_cutoff_self.pdf',height=3,width=3.7,
         pointsize=5)
par(oma=c(3,1,1,1),mar=c(4,2,1,1))
plot(mod_value,"cuminc",ylim=c(0,0.01),xlim=c(0,5), # Remove cuminc if you want survival
     lwd=1.5, # width of lines
     background=F, # background horizontal lines, generally set to false
     axis2.at=seq(0,0.01,0.002),axis2.las=2,axis2.cex.axis=1.8, #y-axis labeling parameters
     axis1.at=seq(0,5,1),axis1.labels=as.character(seq(0,5,1)),axis1.padj=0.5,axis1.cex.axis=1.8, #x-axis labeling parameters
     col=c(col_corr[all_categories=='respiratory']$cat_col,'#8dd3c7'), # color of curves
     atrisk.col='black',
     confint=FALSE, # whether you want CI on the curves
     legend.x=0,legend.y=0.01*1.05, # legend parameters
     legend.cex=1.8,legend.legend=c("Not Meeting Guidelines","Meeting Guidelines"), # more legend parameters
     atrisk.title='                  ',atrisk.pos=-0.5,atrisk.line=c(2,3.5), # position of the N at risk rows
     atrisk.cex=1.6,atrisk.interspace=1.4, # more N at risk parameters
     atrisk.times=c(0,2.5,5), # x-axis points where you want the N at risk to show
     xlab='',ylab='') # Empty axis labels, using mtext instead
mtext("Cumulative risk (%)",side=2,line=-1,at=0.005,cex=2) # y-axis label
mtext("Years",side=1, line=-0.6,cex=2) # x-axis label
mtext('Stratum',side=1,line=0.5,cex=1.6,at=-1.05) # descriptor for N at risk
dev.off()

################### Phenotype 2: PAD
setkey(pad,sample_id)
pad_value <- value_bmi[pad,nomatch=0]

# Format variables
pad_value[,censor_date := as.Date(censor_date,format='%Y-%m-%d')]

# Create analysis variables
pad_value[,time_to_pad := ifelse(c(has_disease == 0 | is.na(has_disease)),pmin(as.numeric(censor_date - enroll_date)/365.25,as.numeric(phenotype_censor_date - enroll_date)/365.25),
                                   as.numeric(censor_date - enroll_date)/365.25)]

# Remove prevalent disease or no follow-up
pad_analysis_value <- pad_value[!is.na(time_to_pad) & time_to_pad > 0]

# Graphical variables
pad_analysis_value[,guidelines := factor(ifelse(self_guidelines==1,'Meeting','Not Meeting'),levels=c('Not Meeting','Meeting'))]

# Fit prodlim models
# Models
mod_value <- prodlim(Hist(time_to_pad,has_disease)~guidelines,data=pad_analysis_value)

# Plot
CairoPDF(file='km_pad_cutoff_self.pdf',height=3,width=3.7,
         pointsize=5)
par(oma=c(3,1,1,1),mar=c(4,2,1,1))
plot(mod_value,"cuminc",ylim=c(0,0.012),xlim=c(0,5), # Remove cuminc if you want survival
     lwd=1.5, # width of lines
     background=F, # background horizontal lines, generally set to false
     axis2.at=seq(0,0.012,0.002),axis2.las=2,axis2.cex.axis=1.8, #y-axis labeling parameters
     axis1.at=seq(0,5,1),axis1.labels=as.character(seq(0,5,1)),axis1.padj=0.5,axis1.cex.axis=1.8, #x-axis labeling parameters
     col=c(col_corr[all_categories=='circulatory system']$cat_col,'#8dd3c7'), # color of curves
     atrisk.col='black',
     confint=FALSE, # whether you want CI on the curves
     legend.x=0,legend.y=0.012*1.05, # legend parameters
     legend.cex=1.8,legend.legend=c("Not Meeting Guidelines","Meeting Guidelines"), # more legend parameters
     atrisk.title='                  ',atrisk.pos=-0.5,atrisk.line=c(2,3.5), # position of the N at risk rows
     atrisk.cex=1.6,atrisk.interspace=1.4, # more N at risk parameters
     atrisk.times=c(0,2.5,5), # x-axis points where you want the N at risk to show
     xlab='',ylab='') # Empty axis labels, using mtext instead
mtext("Cumulative risk (%)",side=2,line=-1,at=0.006,cex=2) # y-axis label
mtext("Years",side=1, line=-0.6,cex=2) # x-axis label
mtext('Stratum',side=1,line=0.5,cex=1.6,at=-1.05) # descriptor for N at risk
dev.off()

################### Phenotype 3: DM
setkey(dm,sample_id)
dm_value <- value_bmi[dm,nomatch=0]

# Format variables
dm_value[,censor_date := as.Date(censor_date,format='%Y-%m-%d')]

# Create analysis variables
dm_value[,time_to_dm := ifelse(c(has_disease == 0 | is.na(has_disease)),pmin(as.numeric(censor_date - enroll_date)/365.25,as.numeric(phenotype_censor_date - enroll_date)/365.25),
                                 as.numeric(censor_date - enroll_date)/365.25)]

# Remove prevalent disease or no follow-up
dm_analysis_value <- dm_value[!is.na(time_to_dm) & time_to_dm > 0]

# Graphical variables
dm_analysis_value[,guidelines := factor(ifelse(self_guidelines==1,'Meeting','Not Meeting'),levels=c('Not Meeting','Meeting'))]

# Fit prodlim models
# Models
mod_value <- prodlim(Hist(time_to_dm,has_disease)~guidelines,data=dm_analysis_value)

# Plot
CairoPDF(file='km_dm_cutoff_self.pdf',height=3,width=3.7,
         pointsize=5)
par(oma=c(3,1,1,1),mar=c(4,2,1,1))
plot(mod_value,"cuminc",ylim=c(0,0.04),xlim=c(0,5), # Remove cuminc if you want survival
     lwd=1.5, # width of lines
     background=F, # background horizontal lines, generally set to false
     axis2.at=seq(0,0.04,0.01),axis2.las=2,axis2.cex.axis=1.8, #y-axis labeling parameters
     axis1.at=seq(0,5,1),axis1.labels=as.character(seq(0,5,1)),axis1.padj=0.5,axis1.cex.axis=1.8, #x-axis labeling parameters
     col=c(col_corr[all_categories=='endocrine/metabolic']$cat_col,'#8dd3c7'), # color of curves
     atrisk.col='black',
     confint=FALSE, # whether you want CI on the curves
     legend.x=0,legend.y=0.04*1.05, # legend parameters
     legend.cex=1.8,legend.legend=c("Not Meeting Guidelines","Meeting Guidelines"), # more legend parameters
     atrisk.title='                  ',atrisk.pos=-0.5,atrisk.line=c(2,3.5), # position of the N at risk rows
     atrisk.cex=1.6,atrisk.interspace=1.4, # more N at risk parameters
     atrisk.times=c(0,2.5,5), # x-axis points where you want the N at risk to show
     xlab='',ylab='') # Empty axis labels, using mtext instead
mtext("Cumulative risk (%)",side=2,line=-1,at=0.02,cex=2) # y-axis label
mtext("Years",side=1, line=-0.6,cex=2) # x-axis label
mtext('Stratum',side=1,line=0.5,cex=1.6,at=-1.05) # descriptor for N at risk
dev.off()

################### Phenotype 4: Gastritis
setkey(gastritis,sample_id)
gastritis_value <- value_bmi[gastritis,nomatch=0]

# Format variables
gastritis_value[,censor_date := as.Date(censor_date,format='%Y-%m-%d')]

# Create analysis variables
gastritis_value[,time_to_gastritis := ifelse(c(has_disease == 0 | is.na(has_disease)),pmin(as.numeric(censor_date - enroll_date)/365.25,as.numeric(phenotype_censor_date - enroll_date)/365.25),
                               as.numeric(censor_date - enroll_date)/365.25)]

# Remove prevalent disease or no follow-up
gastritis_analysis_value <- gastritis_value[!is.na(time_to_gastritis) & time_to_gastritis > 0]

# Graphical variables
gastritis_analysis_value[,guidelines := factor(ifelse(self_guidelines==1,'Meeting','Not Meeting'),levels=c('Not Meeting','Meeting'))]

# Fit prodlim models
# Models
mod_value <- prodlim(Hist(time_to_gastritis,has_disease)~guidelines,data=gastritis_analysis_value)

# Plot
CairoPDF(file='km_gastritis_cutoff_self.pdf',height=3,width=3.7,
         pointsize=5)
par(oma=c(3,1,1,1),mar=c(4,2,1,1))
plot(mod_value,"cuminc",ylim=c(0,0.02),xlim=c(0,5), # Remove cuminc if you want survival
     lwd=1.5, # width of lines
     background=F, # background horizontal lines, generally set to false
     axis2.at=seq(0,0.02,0.005),axis2.las=2,axis2.cex.axis=1.8, #y-axis labeling parameters
     axis1.at=seq(0,5,1),axis1.labels=as.character(seq(0,5,1)),axis1.padj=0.5,axis1.cex.axis=1.8, #x-axis labeling parameters
     col=c(col_corr[all_categories=='digestive']$cat_col,'#8dd3c7'), # color of curves
     atrisk.col='black',
     confint=FALSE, # whether you want CI on the curves
     legend.x=0,legend.y=0.02*1.05, # legend parameters
     legend.cex=1.8,legend.legend=c("Not Meeting Guidelines","Meeting Guidelines"), # more legend parameters
     atrisk.title='                  ',atrisk.pos=-0.5,atrisk.line=c(2,3.5), # position of the N at risk rows
     atrisk.cex=1.6,atrisk.interspace=1.4, # more N at risk parameters
     atrisk.times=c(0,2.5,5), # x-axis points where you want the N at risk to show
     xlab='',ylab='') # Empty axis labels, using mtext instead
mtext("Cumulative risk (%)",side=2,line=-1,at=0.01,cex=2) # y-axis label
mtext("Years",side=1, line=-0.6,cex=2) # x-axis label
mtext('Stratum',side=1,line=0.5,cex=1.6,at=-1.05) # descriptor for N at risk
dev.off()