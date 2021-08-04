# Depends
library(data.table)
library(prodlim)
library(survival)

#### STEP 1: Identify candidate diseases
# Load results files
cutoff_prs_result <- fread(file='~/Documents/MGH Research/accel_phewas/phecode_outputs/cox_mvpa_sqrt_cutoff_prs_processed_fdr.csv')
cutoff_result <- fread(file='~/Documents/MGH Research/accel_phewas/phecode_outputs/cox_cutoff_bmi_processed_fdr.csv')

# Find common ground among sigs
prs_sig <- cutoff_prs_result[sig==1]; mvpa_sig <- cutoff_result[sig==1]
setkey(prs_sig,phecode); setkey(mvpa_sig,phecode)
sig_both <- prs_sig[mvpa_sig,nomatch=0]
setkey(sig_both,p) # Sort by genetic p-value
dz <- sig_both[,.SD[1:10],by='category'] # Lowest p-value for each category

# Will choose: 496.21 (Obstructive chronic bronchitis), 250 (Type 2 Diabetes), 530 (Diseases of esophagus), 401.1 (Essential hypertension)

# Load seed files
value_bmi <- fread(file='/Volumes/medpop_afib/skhurshid/acceleration_phewas/cox_data_mvpa_cutoff_bmi.csv')
inferred <- fread(file='/Volumes/medpop_afib/skhurshid/acceleration_phewas/cox_data_mvpa_sqrt_cutoff_prs.csv')

# Load color file
col_corr <- fread(file='~/Documents/MGH Research/accel_phewas/col_corr.csv')

# Create cutoff variables
setkey(value_bmi,sample_id); setkey(inferred,sample_id)

# Load phecode tables
copd <- read.table(file='/Volumes/medpop_afib/skhurshid/phenotypes/2021_06/phecode_tables/phecode_496_21.tab.tsv',sep='\t',header = TRUE); setDT(copd)
gastritis <- read.table(file='/Volumes/medpop_afib/skhurshid/phenotypes/2021_06/phecode_tables/phecode_535.tab.tsv',sep='\t',header = TRUE); setDT(gastritis)
pad <- read.table(file='/Volumes/medpop_afib/skhurshid/phenotypes/2021_06/phecode_tables/phecode_443.tab.tsv',sep='\t',header = TRUE); setDT(pad)
dm <- read.table(file='/Volumes/medpop_afib/skhurshid/phenotypes/2021_06/phecode_tables/phecode_250.tab.tsv',sep='\t',header = TRUE); setDT(dm)

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

############ A. COPD
# Fit models for adjustment
# Relevel
mod_value <- coxph(Surv(time_to_copd,has_disease) ~ strata(guidelines) + age_accel + sex,data=copd_analysis_value)
mod_prs <- coxph(Surv(time_to_copd,has_disease) ~ strata(guidelines) + enroll_age + sex,data=copd_analysis_prs)

weights <- data.frame(guidelines = levels(copd_analysis_value$guidelines),
                      age_accel = rep(mean(copd_analysis_value$age_accel),2),
                      sex = rep('Female',2))
weights2 <- data.frame(guidelines = levels(copd_analysis_prs$guidelines),
                      enroll_age = rep(mean(copd_analysis_prs$enroll_age),2),
                      sex = rep('Female',2))

# Plot
fit <- survfit(mod_value,newdata=weights)

CairoPDF(file='~/Documents/MGH Research/accel_phewas/km_copd_cutoff_adjusted_f.pdf',height=3,width=3.5,
         pointsize=6)
par(oma=c(3,4,1,1),mar=c(3,4,1,1))
plot(fit,fun=function(x) 1-x,ylim=c(0,0.005),xlim=c(0,5),
     col=c(col_corr[all_categories=='respiratory']$cat_col,'#8dd3c7'),lwd=1.5,
     bty='n',xaxt='n',yaxt='n',xlab='',ylab='')

axis(1,at=0:5,cex.axis=1.5,labels=as.character(0:5))
axis(2,at=seq(0,0.005,0.001),las=2,cex.axis=1.5,labels=c("0 %","0.1 %",'0.2 %','0.3 %','0.4 %','0.5 %'))

mtext("Cumulative risk (%)",side=2,line=6,at=0.0025,cex=1.5)
mtext("Years",side=1, line=2.5,cex=1.5)
legend(-0.1,0.005*1.05,legend=c("Not Meeting Guidelines","Meeting Guidelines"),
       col=c(col_corr[all_categories=='respiratory']$cat_col,'#8dd3c7'),lty=1,lwd=1.5,bty='n',cex=1.5)
dev.off()

# Plot
fit <- survfit(mod_prs,newdata=weights2)

CairoPDF(file='~/Documents/MGH Research/accel_phewas/km_copd_prs_adjusted_f.pdf',height=3,width=3.5,
         pointsize=6)
par(oma=c(3,4,1,1),mar=c(3,4,1,1))
plot(fit,fun=function(x) 1-x,ylim=c(0,0.005),xlim=c(0,5),
     col=c(col_corr[all_categories=='respiratory']$cat_col,'#8dd3c7'),lwd=1.5,
     bty='n',xaxt='n',yaxt='n',xlab='',ylab='')

axis(1,at=0:5,cex.axis=1.5,labels=as.character(0:5))
axis(2,at=seq(0,0.005,0.001),las=2,cex.axis=1.5,labels=c("0 %","0.1 %",'0.2 %','0.3 %','0.4 %','0.5 %'))

mtext("Cumulative risk (%)",side=2,line=6,at=0.0025,cex=1.5)
mtext("Years",side=1, line=2.5,cex=1.5)
legend(-0.1,0.005*1.05,legend=c("Not Meeting Guidelines","Meeting Guidelines"),
       col=c(col_corr[all_categories=='respiratory']$cat_col,'#8dd3c7'),lty=1,lwd=1.5,bty='n',cex=1.5)
dev.off()

############ B. PAD
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

# Fit models for adjustment
# Relevel
mod_value <- coxph(Surv(time_to_pad,has_disease) ~ strata(guidelines) + age_accel + sex,data=pad_analysis_value)
mod_prs <- coxph(Surv(time_to_pad,has_disease) ~ strata(guidelines) + enroll_age + sex,data=pad_analysis_prs)

weights <- data.frame(guidelines = levels(pad_analysis_value$guidelines),
                      age_accel = rep(mean(pad_analysis_value$age_accel),2),
                      sex = rep('Female',2))
weights2 <- data.frame(guidelines = levels(pad_analysis_prs$guidelines),
                       enroll_age = rep(mean(pad_analysis_prs$enroll_age),2),
                       sex = rep('Female',2))

# Plot
fit <- survfit(mod_value,newdata=weights)

CairoPDF(file='~/Documents/MGH Research/accel_phewas/km_pad_cutoff_adjusted_f.pdf',height=3,width=3.5,
         pointsize=6)
par(oma=c(3,4,1,1),mar=c(3,4,1,1))
plot(fit,fun=function(x) 1-x,ylim=c(0,0.008),xlim=c(0,5),
     col=c(col_corr[all_categories=='circulatory system']$cat_col,'#8dd3c7'),lwd=1.5,
     bty='n',xaxt='n',yaxt='n',xlab='',ylab='')

axis(1,at=0:5,cex.axis=1.5,labels=as.character(0:5))
axis(2,at=seq(0,0.008,0.002),las=2,cex.axis=1.5,labels=c('0 %','0.2 %','0.4 %','0.6 %','0.8 %'))

mtext("Cumulative risk (%)",side=2,line=6,at=0.004,cex=1.5)
mtext("Years",side=1, line=2.5,cex=1.5)
legend(-0.1,0.008*1.05,legend=c("Not Meeting Guidelines","Meeting Guidelines"),
       col=c(col_corr[all_categories=='circulatory system']$cat_col,'#8dd3c7'),lty=1,lwd=1.5,bty='n',cex=1.5)
dev.off()

# Plot
fit <- survfit(mod_prs,newdata=weights2)

CairoPDF(file='~/Documents/MGH Research/accel_phewas/km_pad_prs_adjusted_f.pdf',height=3,width=3.5,
         pointsize=6)
par(oma=c(3,4,1,1),mar=c(3,4,1,1))
plot(fit,fun=function(x) 1-x,ylim=c(0,0.008),xlim=c(0,5),
     col=c(col_corr[all_categories=='circulatory system']$cat_col,'#8dd3c7'),lwd=1.5,
     bty='n',xaxt='n',yaxt='n',xlab='',ylab='')

axis(1,at=0:5,cex.axis=1.5,labels=as.character(0:5))
axis(2,at=seq(0,0.008,0.002),las=2,cex.axis=1.5,labels=c('0 %','0.2 %','0.4 %','0.6 %','0.8 %'))

mtext("Cumulative risk (%)",side=2,line=6,at=0.004,cex=1.5)
mtext("Years",side=1, line=2.5,cex=1.5)
legend(-0.1,0.008*1.05,legend=c("Not Meeting Guidelines","Meeting Guidelines"),
       col=c(col_corr[all_categories=='circulatory system']$cat_col,'#8dd3c7'),lty=1,lwd=1.5,bty='n',cex=1.5)
dev.off()

############ C. DM
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

# Fit models for adjustment
# Relevel
mod_value <- coxph(Surv(time_to_dm,has_disease) ~ strata(guidelines) + age_accel + sex,data=dm_analysis_value)
mod_prs <- coxph(Surv(time_to_dm,has_disease) ~ strata(guidelines) + enroll_age + sex,data=dm_analysis_prs)

weights <- data.frame(guidelines = levels(dm_analysis_value$guidelines),
                      age_accel = rep(mean(dm_analysis_value$age_accel),2),
                      sex = rep('Female',2))
weights2 <- data.frame(guidelines = levels(dm_analysis_prs$guidelines),
                       enroll_age = rep(mean(dm_analysis_prs$enroll_age),2),
                       sex = rep('Female',2))

# Plot
fit <- survfit(mod_value,newdata=weights)

CairoPDF(file='~/Documents/MGH Research/accel_phewas/km_dm_cutoff_adjusted_f.pdf',height=3,width=3.5,
         pointsize=6)
par(oma=c(3,4,1,1),mar=c(3,4,1,1))
plot(fit,fun=function(x) 1-x,ylim=c(0,0.04),xlim=c(0,5),
     col=c(col_corr[all_categories=='endocrine/metabolic']$cat_col,'#8dd3c7'),lwd=1.5,
     bty='n',xaxt='n',yaxt='n',xlab='',ylab='')

axis(1,at=0:5,cex.axis=1.5,labels=as.character(0:5))
axis(2,at=seq(0,0.04,0.01),las=2,cex.axis=1.5,labels=c('0 %','1 %','2 %','3 %','4 %'))

mtext("Cumulative risk (%)",side=2,line=6,at=0.02,cex=1.5)
mtext("Years",side=1, line=2.5,cex=1.5)
legend(-0.1,0.04*1.05,legend=c("Not Meeting Guidelines","Meeting Guidelines"),
       col=c(col_corr[all_categories=='endocrine/metabolic']$cat_col,'#8dd3c7'),lty=1,lwd=1.5,bty='n',cex=1.5)
dev.off()

# Plot
fit <- survfit(mod_prs,newdata=weights2)

CairoPDF(file='~/Documents/MGH Research/accel_phewas/km_dm_prs_adjusted_f.pdf',height=3,width=3.5,
         pointsize=6)
par(oma=c(3,4,1,1),mar=c(3,4,1,1))
plot(fit,fun=function(x) 1-x,ylim=c(0,0.04),xlim=c(0,5),
     col=c(col_corr[all_categories=='endocrine/metabolic']$cat_col,'#8dd3c7'),lwd=1.5,
     bty='n',xaxt='n',yaxt='n',xlab='',ylab='')

axis(1,at=0:5,cex.axis=1.5,labels=as.character(0:5))
axis(2,at=seq(0,0.04,0.01),las=2,cex.axis=1.5,labels=c('0 %','1 %','2 %','3 %','4 %'))

mtext("Cumulative risk (%)",side=2,line=6,at=0.02,cex=1.5)
mtext("Years",side=1, line=2.5,cex=1.5)
legend(-0.1,0.04*1.05,legend=c("Not Meeting Guidelines","Meeting Guidelines"),
       col=c(col_corr[all_categories=='endocrine/metabolic']$cat_col,'#8dd3c7'),lty=1,lwd=1.5,bty='n',cex=1.5)
dev.off()

############ D. gastritis
setkey(gastritis,sample_id)
gastritis_value <- value_bmi[gastritis,nomatch=0]
gastritis_prs <- inferred[gastritis,nomatch=0]

# Format variables
gastritis_value[,censor_date := as.Date(censor_date,format='%Y-%m-%d')]
gastritis_prs[,censor_date := as.Date(censor_date,format='%Y-%m-%d')]

# Create analysis variables
gastritis_value[,time_to_gastritis := ifelse(c(has_disease == 0 | is.na(has_disease)),pmin(as.numeric(censor_date - accel_date)/365.25,as.numeric(phenotype_censor_date - accel_date)/365.25),
                                   as.numeric(censor_date - accel_date)/365.25)]
gastritis_prs[,time_to_gastritis := ifelse(c(has_disease == 0 | is.na(has_disease)),pmin(as.numeric(censor_date - enroll_date)/365.25,as.numeric(phenotype_censor_date - enroll_date)/365.25),
                                 as.numeric(censor_date - enroll_date)/365.25)]
# Remove prevalent disease or no follow-up
gastritis_analysis_value <- gastritis_value[!is.na(time_to_gastritis) & time_to_gastritis > 0]
gastritis_analysis_prs <- gastritis_prs[!is.na(time_to_gastritis) & time_to_gastritis > 0]

# Graphical variables
gastritis_analysis_value[,guidelines := factor(ifelse(who_acc_rate_bouted==1,'Meeting','Not Meeting'),levels=c('Not Meeting','Meeting'))]
gastritis_analysis_prs[,guidelines := factor(ifelse(inferred_accel==1,'Meeting','Not Meeting'),levels=c('Not Meeting','Meeting'))]

# Fit models for adjustment
# Relevel
mod_value <- coxph(Surv(time_to_gastritis,has_disease) ~ strata(guidelines) + age_accel + sex,data=gastritis_analysis_value)
mod_prs <- coxph(Surv(time_to_gastritis,has_disease) ~ strata(guidelines) + enroll_age + sex,data=gastritis_analysis_prs)

weights <- data.frame(guidelines = levels(gastritis_analysis_value$guidelines),
                      age_accel = rep(mean(gastritis_analysis_value$age_accel),2),
                      sex = rep('Female',2))
weights2 <- data.frame(guidelines = levels(gastritis_analysis_prs$guidelines),
                       enroll_age = rep(mean(gastritis_analysis_prs$enroll_age),2),
                       sex = rep('Female',2))

# Plot
fit <- survfit(mod_value,newdata=weights)

CairoPDF(file='~/Documents/MGH Research/accel_phewas/km_gastritis_cutoff_adjusted_f.pdf',height=3,width=3.5,
         pointsize=6)
par(oma=c(3,4,1,1),mar=c(3,4,1,1))
plot(fit,fun=function(x) 1-x,ylim=c(0,0.04),xlim=c(0,5),
     col=c(col_corr[all_categories=='digestive']$cat_col,'#8dd3c7'),lwd=1.5,
     bty='n',xaxt='n',yaxt='n',xlab='',ylab='')

axis(1,at=0:5,cex.axis=1.5,labels=as.character(0:5))
axis(2,at=seq(0,0.04,0.01),las=2,cex.axis=1.5,labels=c('0 %','1 %','2 %','3 %','4 %'))

mtext("Cumulative risk (%)",side=2,line=6,at=0.02,cex=1.5)
mtext("Years",side=1, line=2.5,cex=1.5)
legend(-0.1,0.04*1.05,legend=c("Not Meeting Guidelines","Meeting Guidelines"),
       col=c(col_corr[all_categories=='digestive']$cat_col,'#8dd3c7'),lty=1,lwd=1.5,bty='n',cex=1.5)
dev.off()

# Plot
fit <- survfit(mod_prs,newdata=weights2)

CairoPDF(file='~/Documents/MGH Research/accel_phewas/km_gastritis_prs_adjusted_f.pdf',height=3,width=3.5,
         pointsize=6)
par(oma=c(3,4,1,1),mar=c(3,4,1,1))
plot(fit,fun=function(x) 1-x,ylim=c(0,0.04),xlim=c(0,5),
     col=c(col_corr[all_categories=='digestive']$cat_col,'#8dd3c7'),lwd=1.5,
     bty='n',xaxt='n',yaxt='n',xlab='',ylab='')

axis(1,at=0:5,cex.axis=1.5,labels=as.character(0:5))
axis(2,at=seq(0,0.04,0.01),las=2,cex.axis=1.5,labels=c('0 %','1 %','2 %','3 %','4 %'))

mtext("Cumulative risk (%)",side=2,line=6,at=0.02,cex=1.5)
mtext("Years",side=1, line=2.5,cex=1.5)
legend(-0.1,0.04*1.05,legend=c("Not Meeting Guidelines","Meeting Guidelines"),
       col=c(col_corr[all_categories=='digestive']$cat_col,'#8dd3c7'),lty=1,lwd=1.5,bty='n',cex=1.5)
dev.off()