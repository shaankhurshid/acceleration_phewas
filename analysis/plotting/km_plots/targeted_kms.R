# Depends
library(data.table)
library(prodlim)

#### STEP 1: Identify candidate diseases
# Load results files
cutoff_prs_result <- fread(file='~/Documents/MGH Research/accel_phewas/phecode_outputs/cox_mvpa_sqrt_cutoff_prs_processed_fdr.csv')
cutoff_result <- fread(file='~/Documents/MGH Research/accel_phewas/phecode_outputs/cox_cutoff_bmi_processed_fdr.csv')

# Find common ground among sigs
prs_sig <- cutoff_prs_result[sig==1]; mvpa_sig <- cutoff_result[sig==1]
setkey(prs_sig,phecode); setkey(mvpa_sig,phecode)
sig_both <- prs_sig[mvpa_sig,nomatch=0]
setkey(sig_both,hr) # Sort by genetic hr
dz <- sig_both[,.SD[1:10],by='category'] # Lowest HRs for each category

# Will choose: 496.21 (Obstructive chronic bronchitis), 250 (Type 2 Diabetes), 571.5 (Non-Alcholic Liver Diseease), 443 (Peripheral Vascular Disease)

# Load seed files
value_bmi <- fread(file='/Volumes/medpop_afib/skhurshid/acceleration_phewas/cox_data_mvpa_cutoff_bmi.csv')
inferred <- fread(file='/Volumes/medpop_afib/skhurshid/acceleration_phewas/cox_data_mvpa_sqrt_cutoff_prs.csv')

# Load color file
col_corr <- fread(file='~/Documents/MGH Research/accel_phewas/col_corr.csv')

# Create cutoff variables
setkey(value_bmi,sample_id); setkey(inferred,sample_id)

# Load phecode tables
copd <- read.table(file='/Volumes/medpop_afib/skhurshid/phenotypes/2020_06/phecode_tables/phecode_496_21.tab.tsv',sep='\t',header = TRUE); setDT(copd)
liver <- read.table(file='/Volumes/medpop_afib/skhurshid/phenotypes/2020_06/phecode_tables/phecode_571_5.tab.tsv',sep='\t',header = TRUE); setDT(liver)
pad <- read.table(file='/Volumes/medpop_afib/skhurshid/phenotypes/2020_06/phecode_tables/phecode_443.tab.tsv',sep='\t',header = TRUE); setDT(pad)
dm <- read.table(file='/Volumes/medpop_afib/skhurshid/phenotypes/2020_06/phecode_tables/phecode_250.tab.tsv',sep='\t',header = TRUE); setDT(dm)

# Format dates
for (j in (c('enroll_date','phenotype_censor_date'))){set(inferred,j=j,value=as.Date(inferred[[j]],format='%Y-%m-%d'))}
for (j in (c('accel_date','phenotype_censor_date'))){set(value_bmi,j=j,value=as.Date(value_bmi[[j]],format='%Y-%m-%d'))}

################### Phenotype 1: COPD
setkey(copd,sample_id)
copd_value <- value_bmi[copd,nomatch=0]
copd_prs <- inferred[copd,nomatch=0]

# Format variables
copd_value[,censor_date := as.Date(censor_date,format='%Y-%m-%d')]
copd_prs[,censor_date := as.Date(censor_date,format='%Y-%m-%d')]

# Create analysis variables
copd_value[,time_to_copd := ifelse(c(has_disease == 0 | is.na(has_disease)),pmin(as.numeric(censor_date - accel_date)/365.25,as.numeric(phenotype_censor_date - accel_date)/365.25),
                                      as.numeric(censor_date - accel_date)/365.25)]
copd_prs[,time_to_copd := ifelse(c(has_disease == 0 | is.na(has_disease)),pmin(as.numeric(censor_date - enroll_date)/365.25,as.numeric(phenotype_censor_date - enroll_date)/365.25),
                                    as.numeric(censor_date - enroll_date)/365.25)]

# Remove prevalent disease or no follow-up
copd_analysis_value <- copd_value[!is.na(time_to_copd) & time_to_copd > 0]
copd_analysis_prs <- copd_prs[!is.na(time_to_copd) & time_to_copd > 0]

# Graphical variables
copd_analysis_value[,guidelines := factor(ifelse(who_acc_rate_bouted==1,'Meeting','Not Meeting'),levels=c('Not Meeting','Meeting'))]
copd_analysis_prs[,guidelines := factor(ifelse(inferred_accel==1,'Meeting','Not Meeting'),levels=c('Not Meeting','Meeting'))]

# Fit prodlim models
# Models
mod_value <- prodlim(Hist(time_to_copd,has_disease)~guidelines,data=copd_analysis_value)
mod_prs <- prodlim(Hist(time_to_copd,has_disease)~guidelines,data=copd_analysis_prs)

# Plot
CairoPDF(file='~/Documents/MGH Research/accel_phewas/km_copd_cutoff.pdf',height=3,width=3.7,
         pointsize=5)
par(oma=c(3,1,1,1),mar=c(4,2,1,1))
plot(mod_value,"cuminc",ylim=c(0,0.02),xlim=c(0,5), # Remove cuminc if you want survival
     lwd=1.5, # width of lines
     background=F, # background horizontal lines, generally set to false
     axis2.at=seq(0,0.02,0.005),axis2.las=2,axis2.cex.axis=1.8, #y-axis labeling parameters
     axis1.at=seq(0,5,1),axis1.labels=as.character(seq(0,5,1)),axis1.padj=0.5,axis1.cex.axis=1.8, #x-axis labeling parameters
     col=c(col_corr[all_categories=='respiratory']$cat_col,'#8dd3c7'), # color of curves
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

# Plot
CairoPDF(file='~/Documents/MGH Research/accel_phewas/km_copd_prs.pdf',height=3,width=3.7,
         pointsize=5)
par(oma=c(3,1,1,1),mar=c(4,2,1,1))
plot(mod_prs,"cuminc",ylim=c(0,0.02),xlim=c(0,10), # Remove cuminc if you want survival
     lwd=1.5, # width of lines
     background=F, # background horizontal lines, generally set to false
     axis2.at=seq(0,0.02,0.005),axis2.las=2,axis2.cex.axis=1.8, #y-axis labeling parameters
     axis1.at=seq(0,10,2),axis1.labels=as.character(seq(0,10,2)),axis1.padj=0.5,axis1.cex.axis=1.8, #x-axis labeling parameters
     col=c(col_corr[all_categories=='respiratory']$cat_col,'#8dd3c7'), # color of curves
     atrisk.col='black',
     confint=FALSE, # whether you want CI on the curves
     legend.x=0,legend.y=0.02*1.05, # legend parameters
     legend.cex=1.8,legend.legend=c("Not Meeting Guidelines","Meeting Guidelines"), # more legend parameters
     atrisk.title='                  ',atrisk.pos=-1,atrisk.line=c(2,3.5), # position of the N at risk rows
     atrisk.cex=1.6,atrisk.interspace=1.4, # more N at risk parameters
     atrisk.times=c(0,2.5,5,7.5,10), # x-axis points where you want the N at risk to show
     xlab='',ylab='') # Empty axis labels, using mtext instead
mtext("Cumulative risk (%)",side=2,line=-1,at=0.01,cex=2) # y-axis label
mtext("Years",side=1, line=-0.5,cex=2) # x-axis label
mtext('Stratum',side=1, line=0.5,cex=1.6,at=-2.1) # descriptor for N at risk
dev.off()

################### Phenotype 2: PAD
setkey(pad,sample_id)
pad_value <- value_bmi[pad,nomatch=0]
pad_prs <- inferred[pad,nomatch=0]

# Format variables
pad_value[,censor_date := as.Date(censor_date,format='%Y-%m-%d')]
pad_prs[,censor_date := as.Date(censor_date,format='%Y-%m-%d')]

# Create analysis variables
pad_value[,time_to_pad := ifelse(c(has_disease == 0 | is.na(has_disease)),pmin(as.numeric(censor_date - accel_date)/365.25,as.numeric(phenotype_censor_date - accel_date)/365.25),
                                   as.numeric(censor_date - accel_date)/365.25)]
pad_prs[,time_to_pad := ifelse(c(has_disease == 0 | is.na(has_disease)),pmin(as.numeric(censor_date - enroll_date)/365.25,as.numeric(phenotype_censor_date - enroll_date)/365.25),
                                 as.numeric(censor_date - enroll_date)/365.25)]
# Remove prevalent disease or no follow-up
pad_analysis_value <- pad_value[!is.na(time_to_pad) & time_to_pad > 0]
pad_analysis_prs <- pad_prs[!is.na(time_to_pad) & time_to_pad > 0]

# Graphical variables
pad_analysis_value[,guidelines := factor(ifelse(who_acc_rate_bouted==1,'Meeting','Not Meeting'),levels=c('Not Meeting','Meeting'))]
pad_analysis_prs[,guidelines := factor(ifelse(inferred_accel==1,'Meeting','Not Meeting'),levels=c('Not Meeting','Meeting'))]

# Fit prodlim models
# Models
mod_value <- prodlim(Hist(time_to_pad,has_disease)~guidelines,data=pad_analysis_value)
mod_prs <- prodlim(Hist(time_to_pad,has_disease)~guidelines,data=pad_analysis_prs)

# Plot
CairoPDF(file='~/Documents/MGH Research/accel_phewas/km_pad_cutoff.pdf',height=3,width=3.7,
         pointsize=5)
par(oma=c(3,1,1,1),mar=c(4,2,1,1))
plot(mod_value,"cuminc",ylim=c(0,0.025),xlim=c(0,5), # Remove cuminc if you want survival
     lwd=1.5, # width of lines
     background=F, # background horizontal lines, generally set to false
     axis2.at=seq(0,0.025,0.005),axis2.las=2,axis2.cex.axis=1.8, #y-axis labeling parameters
     axis1.at=seq(0,5,1),axis1.labels=as.character(seq(0,5,1)),axis1.padj=0.5,axis1.cex.axis=1.8, #x-axis labeling parameters
     col=c(col_corr[all_categories=='circulatory system']$cat_col,'#8dd3c7'), # color of curves
     atrisk.col='black',
     confint=FALSE, # whether you want CI on the curves
     legend.x=0,legend.y=0.025*1.05, # legend parameters
     legend.cex=1.8,legend.legend=c("Not Meeting Guidelines","Meeting Guidelines"), # more legend parameters
     atrisk.title='                  ',atrisk.pos=-0.5,atrisk.line=c(2,3.5), # position of the N at risk rows
     atrisk.cex=1.6,atrisk.interspace=1.4, # more N at risk parameters
     atrisk.times=c(0,2.5,5), # x-axis points where you want the N at risk to show
     xlab='',ylab='') # Empty axis labels, using mtext instead
mtext("Cumulative risk (%)",side=2,line=-1,at=0.0125,cex=2) # y-axis label
mtext("Years",side=1, line=-0.6,cex=2) # x-axis label
mtext('Stratum',side=1,line=0.5,cex=1.6,at=-1.05) # descriptor for N at risk
dev.off()

# Plot
CairoPDF(file='~/Documents/MGH Research/accel_phewas/km_pad_prs.pdf',height=3,width=3.7,
         pointsize=5)
par(oma=c(3,1,1,1),mar=c(4,2,1,1))
plot(mod_prs,"cuminc",ylim=c(0,0.025),xlim=c(0,10), # Remove cuminc if you want survival
     lwd=1.5, # width of lines
     background=F, # background horizontal lines, generally set to false
     axis2.at=seq(0,0.025,0.005),axis2.las=2,axis2.cex.axis=1.8, #y-axis labeling parameters
     axis1.at=seq(0,10,2),axis1.labels=as.character(seq(0,10,2)),axis1.padj=0.5,axis1.cex.axis=1.8, #x-axis labeling parameters
     col=c(col_corr[all_categories=='circulatory system']$cat_col,'#8dd3c7'), # color of curves
     atrisk.col='black',
     confint=FALSE, # whether you want CI on the curves
     legend.x=0,legend.y=0.025*1.05, # legend parameters
     legend.cex=1.8,legend.legend=c("Not Meeting Guidelines","Meeting Guidelines"), # more legend parameters
     atrisk.title='                  ',atrisk.pos=-1,atrisk.line=c(2,3.5), # position of the N at risk rows
     atrisk.cex=1.6,atrisk.interspace=1.4, # more N at risk parameters
     atrisk.times=c(0,2.5,5,7.5,10), # x-axis points where you want the N at risk to show
     xlab='',ylab='') # Empty axis labels, using mtext instead
mtext("Cumulative risk (%)",side=2,line=-1,at=0.0125,cex=2) # y-axis label
mtext("Years",side=1, line=-0.5,cex=2) # x-axis label
mtext('Stratum',side=1, line=0.5,cex=1.6,at=-2.1) # descriptor for N at risk
dev.off()

################### Phenotype 3: DM
setkey(dm,sample_id)
dm_value <- value_bmi[dm,nomatch=0]
dm_prs <- inferred[dm,nomatch=0]

# Format variables
dm_value[,censor_date := as.Date(censor_date,format='%Y-%m-%d')]
dm_prs[,censor_date := as.Date(censor_date,format='%Y-%m-%d')]

# Create analysis variables
dm_value[,time_to_dm := ifelse(c(has_disease == 0 | is.na(has_disease)),pmin(as.numeric(censor_date - accel_date)/365.25,as.numeric(phenotype_censor_date - accel_date)/365.25),
                                 as.numeric(censor_date - accel_date)/365.25)]
dm_prs[,time_to_dm := ifelse(c(has_disease == 0 | is.na(has_disease)),pmin(as.numeric(censor_date - enroll_date)/365.25,as.numeric(phenotype_censor_date - enroll_date)/365.25),
                               as.numeric(censor_date - enroll_date)/365.25)]

# Remove prevalent disease or no follow-up
dm_analysis_value <- dm_value[!is.na(time_to_dm) & time_to_dm > 0]
dm_analysis_prs <- dm_prs[!is.na(time_to_dm) & time_to_dm > 0]

# Graphical variables
dm_analysis_value[,guidelines := factor(ifelse(who_acc_rate_bouted==1,'Meeting','Not Meeting'),levels=c('Not Meeting','Meeting'))]
dm_analysis_prs[,guidelines := factor(ifelse(inferred_accel==1,'Meeting','Not Meeting'),levels=c('Not Meeting','Meeting'))]

# Fit prodlim models
# Models
mod_value <- prodlim(Hist(time_to_dm,has_disease)~guidelines,data=dm_analysis_value)
mod_prs <- prodlim(Hist(time_to_dm,has_disease)~guidelines,data=dm_analysis_prs)

# Plot
CairoPDF(file='~/Documents/MGH Research/accel_phewas/km_dm_cutoff.pdf',height=3,width=3.7,
         pointsize=5)
par(oma=c(3,1,1,1),mar=c(4,2,1,1))
plot(mod_value,"cuminc",ylim=c(0,0.08),xlim=c(0,5), # Remove cuminc if you want survival
     lwd=1.5, # width of lines
     background=F, # background horizontal lines, generally set to false
     axis2.at=seq(0,0.08,0.02),axis2.las=2,axis2.cex.axis=1.8, #y-axis labeling parameters
     axis1.at=seq(0,5,1),axis1.labels=as.character(seq(0,5,1)),axis1.padj=0.5,axis1.cex.axis=1.8, #x-axis labeling parameters
     col=c(col_corr[all_categories=='endocrine/metabolic']$cat_col,'#8dd3c7'), # color of curves
     atrisk.col='black',
     confint=FALSE, # whether you want CI on the curves
     legend.x=0,legend.y=0.08*1.05, # legend parameters
     legend.cex=1.8,legend.legend=c("Not Meeting Guidelines","Meeting Guidelines"), # more legend parameters
     atrisk.title='                  ',atrisk.pos=-0.5,atrisk.line=c(2,3.5), # position of the N at risk rows
     atrisk.cex=1.6,atrisk.interspace=1.4, # more N at risk parameters
     atrisk.times=c(0,2.5,5), # x-axis points where you want the N at risk to show
     xlab='',ylab='') # Empty axis labels, using mtext instead
mtext("Cumulative risk (%)",side=2,line=-1,at=0.04,cex=2) # y-axis label
mtext("Years",side=1, line=-0.6,cex=2) # x-axis label
mtext('Stratum',side=1,line=0.5,cex=1.6,at=-1.05) # descriptor for N at risk
dev.off()

# Plot
CairoPDF(file='~/Documents/MGH Research/accel_phewas/km_dm_prs.pdf',height=3,width=3.7,
         pointsize=5)
par(oma=c(3,1,1,1),mar=c(4,2,1,1))
plot(mod_prs,"cuminc",ylim=c(0,0.08),xlim=c(0,10), # Remove cuminc if you want survival
     lwd=1.5, # width of lines
     background=F, # background horizontal lines, generally set to false
     axis2.at=seq(0,0.08,0.02),axis2.las=2,axis2.cex.axis=1.8, #y-axis labeling parameters
     axis1.at=seq(0,10,2),axis1.labels=as.character(seq(0,10,2)),axis1.padj=0.5,axis1.cex.axis=1.8, #x-axis labeling parameters
     col=c(col_corr[all_categories=='endocrine/metabolic']$cat_col,'#8dd3c7'), # color of curves
     atrisk.col='black',
     confint=FALSE, # whether you want CI on the curves
     legend.x=0,legend.y=0.08*1.05, # legend parameters
     legend.cex=1.8,legend.legend=c("Not Meeting Guidelines","Meeting Guidelines"), # more legend parameters
     atrisk.title='                  ',atrisk.pos=-1,atrisk.line=c(2,3.5), # position of the N at risk rows
     atrisk.cex=1.6,atrisk.interspace=1.4, # more N at risk parameters
     atrisk.times=c(0,2.5,5,7.5,10), # x-axis points where you want the N at risk to show
     xlab='',ylab='') # Empty axis labels, using mtext instead
mtext("Cumulative risk (%)",side=2,line=-1,at=0.04,cex=2) # y-axis label
mtext("Years",side=1, line=-0.5,cex=2) # x-axis label
mtext('Stratum',side=1, line=0.5,cex=1.6,at=-2.1) # descriptor for N at risk
dev.off()

################### Phenotype 4: Liver disease
setkey(liver,sample_id)
liver_value <- value_bmi[liver,nomatch=0]
liver_prs <- inferred[liver,nomatch=0]

# Format variables
liver_value[,censor_date := as.Date(censor_date,format='%Y-%m-%d')]
liver_prs[,censor_date := as.Date(censor_date,format='%Y-%m-%d')]

# Create analysis variables
liver_value[,time_to_liver := ifelse(c(has_disease == 0 | is.na(has_disease)),pmin(as.numeric(censor_date - accel_date)/365.25,as.numeric(phenotype_censor_date - accel_date)/365.25),
                               as.numeric(censor_date - accel_date)/365.25)]
liver_prs[,time_to_liver := ifelse(c(has_disease == 0 | is.na(has_disease)),pmin(as.numeric(censor_date - enroll_date)/365.25,as.numeric(phenotype_censor_date - enroll_date)/365.25),
                             as.numeric(censor_date - enroll_date)/365.25)]
# Remove prevalent disease or no follow-up
liver_analysis_value <- liver_value[!is.na(time_to_liver) & time_to_liver > 0]
liver_analysis_prs <- liver_prs[!is.na(time_to_liver) & time_to_liver > 0]

# Graphical variables
liver_analysis_value[,guidelines := factor(ifelse(who_acc_rate_bouted==1,'Meeting','Not Meeting'),levels=c('Not Meeting','Meeting'))]
liver_analysis_prs[,guidelines := factor(ifelse(inferred_accel==1,'Meeting','Not Meeting'),levels=c('Not Meeting','Meeting'))]

# Fit prodlim models
# Models
mod_value <- prodlim(Hist(time_to_liver,has_disease)~guidelines,data=liver_analysis_value)
mod_prs <- prodlim(Hist(time_to_liver,has_disease)~guidelines,data=liver_analysis_prs)

# Plot
CairoPDF(file='~/Documents/MGH Research/accel_phewas/km_liver_cutoff.pdf',height=3,width=3.7,
         pointsize=5)
par(oma=c(3,1,1,1),mar=c(4,2,1,1))
plot(mod_value,"cuminc",ylim=c(0,0.015),xlim=c(0,5), # Remove cuminc if you want survival
     lwd=1.5, # width of lines
     background=F, # background horizontal lines, generally set to false
     axis2.at=seq(0,0.015,0.005),axis2.las=2,axis2.cex.axis=1.8, #y-axis labeling parameters
     axis1.at=seq(0,5,1),axis1.labels=as.character(seq(0,5,1)),axis1.padj=0.5,axis1.cex.axis=1.8, #x-axis labeling parameters
     col=c(col_corr[all_categories=='digestive']$cat_col,'#8dd3c7'), # color of curves
     atrisk.col='black',
     confint=FALSE, # whether you want CI on the curves
     legend.x=0,legend.y=0.015*1.05, # legend parameters
     legend.cex=1.8,legend.legend=c("Not Meeting Guidelines","Meeting Guidelines"), # more legend parameters
     atrisk.title='                  ',atrisk.pos=-0.5,atrisk.line=c(2,3.5), # position of the N at risk rows
     atrisk.cex=1.6,atrisk.interspace=1.4, # more N at risk parameters
     atrisk.times=c(0,2.5,5), # x-axis points where you want the N at risk to show
     xlab='',ylab='') # Empty axis labels, using mtext instead
mtext("Cumulative risk (%)",side=2,line=-1,at=0.0075,cex=2) # y-axis label
mtext("Years",side=1, line=-0.6,cex=2) # x-axis label
mtext('Stratum',side=1,line=0.5,cex=1.6,at=-1.05) # descriptor for N at risk
dev.off()

# Plot
CairoPDF(file='~/Documents/MGH Research/accel_phewas/km_liver_prs.pdf',height=3,width=3.7,
         pointsize=5)
par(oma=c(3,1,1,1),mar=c(4,2,1,1))
plot(mod_prs,"cuminc",ylim=c(0,0.015),xlim=c(0,10), # Remove cuminc if you want survival
     lwd=1.5, # width of lines
     background=F, # background horizontal lines, generally set to false
     axis2.at=seq(0,0.015,0.005),axis2.las=2,axis2.cex.axis=1.8, #y-axis labeling parameters
     axis1.at=seq(0,10,2),axis1.labels=as.character(seq(0,10,2)),axis1.padj=0.5,axis1.cex.axis=1.8, #x-axis labeling parameters
     col=c(col_corr[all_categories=='digestive']$cat_col,'#8dd3c7'), # color of curves
     atrisk.col='black',
     confint=FALSE, # whether you want CI on the curves
     legend.x=0,legend.y=0.015*1.05, # legend parameters
     legend.cex=1.8,legend.legend=c("Not Meeting Guidelines","Meeting Guidelines"), # more legend parameters
     atrisk.title='                  ',atrisk.pos=-1,atrisk.line=c(2,3.5), # position of the N at risk rows
     atrisk.cex=1.6,atrisk.interspace=1.4, # more N at risk parameters
     atrisk.times=c(0,2.5,5,7.5,10), # x-axis points where you want the N at risk to show
     xlab='',ylab='') # Empty axis labels, using mtext instead
mtext("Cumulative risk (%)",side=2,line=-1,at=0.0075,cex=2) # y-axis label
mtext("Years",side=1, line=-0.5,cex=2) # x-axis label
mtext('Stratum',side=1, line=0.5,cex=1.6,at=-2.1) # descriptor for N at risk
dev.off()