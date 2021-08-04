# Script to analyze PRS-pheWAS

# Depends 
library(data.table)
library(ggplot2)

# Load PRS
prs <- fread(file='/Volumes/medpop_afib/skhurshid/acceleration_phewas/prs/scores/mvpa_prs4/mvpa_prs4_processed.csv')

# Load accel phenotype file
acceleration_pheno <- fread(file='/Volumes/medpop_afib/skhurshid/acceleration_phewas/acceleration_phenotype.csv')

# Load accel genetic file
acceleration <- fread(file='/Volumes/medpop_afib/skhurshid/acceleration_phewas/mvpa_rate_bouted.tsv')

# Join
setkey(prs,sample_id); setkey(acceleration,FID)
acceleration[prs,':='(prs_std = i.prs_std)]

# Remove no PRS (poor quality genetic data)
prs_set <- acceleration[!is.na(prs_std)]

# Check exposure ~ PRS
mod <- lm(mvpa_rate_bouted ~ prs_std,data=prs_set)

# Now fit genetic acceleration model
mod_fitted <- lm(mvpa_rate_bouted ~ prs_std + PC1 + PC2 + PC3 + PC4 + PC5 + array_UKBB + male + age_accel,data=prs_set)

# Now load the no accel sets
prs_na <- fread(file='/Volumes/medpop_afib/skhurshid/acceleration_phewas/prs/scores/mvpa_prs_no_accel/mvpa_prs_no_accel_processed.csv')
prs_covar <- fread(file='/Volumes/medpop_afib/skhurshid/acceleration_phewas/prs/no_accel_seed.tsv')

# Merge
setkey(prs_na,sample_id); setkey(prs_covar,sample_id)
prs_no_accel <- prs_covar[prs_na,nomatch=0]

# Reduce feature space
prs_no_accel <- prs_no_accel[,c('sample_id','sex','enroll_age','array_UKBB','prs_std',paste0('PC',1:5))]
setnames(prs_no_accel,c('enroll_age','sex'),c('age_accel','male'))
prs_no_accel[,male := ifelse(male=='Male',1,0)]

# Get inferred acceleration for analysis
prs_set[,inferred_accel := predict(mod_fitted,newdata=prs_set)]

# Compare distributions
# Actual vs Inferred
x <- list(v1=prs_set$mvpa_rate_bouted,v2=prs_set$inferred_accel)
data <- melt(x)

# Plot
ggplot() + geom_density(data=data,aes(x=value,fill=L1),alpha=0.55) + 
  scale_x_continuous(breaks=seq(0,150,50),expand=c(0,0),limits=c(0,150)) +
  scale_y_continuous(breaks=seq(0,0.06,0.01),expand=c(0,0),limits=c(0,0.06)) +
  scale_fill_manual(values=c("#2b8cbe","#f03b20"),name='',labels=c('MVPA','Inferred MVPA')) +
  theme(panel.background=element_blank(),axis.line=element_line(color='black'),legend.position=c(0.80,0.90),
        axis.text=element_text(size=20,color='black'),plot.margin=unit(c(0.5,1,0.5,0.5),'cm'),
        axis.title.y = element_text(size=20,margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(size=20),legend.text=element_text(size=20)) +
  geom_vline(xintercept=150/7,linetype=5,color='black') +
  labs(x='MVPA minutes/day',y='Density') 
ggsave('~/Documents/MGH Research/accel_phewas/accel_mvpa_density_compare.pdf',
       height=2,width=2.5,units='in',scale=4)

# Figure out expected counts for guidelines
prs_set[,':='(guideline = ifelse(mvpa_rate_bouted >= 150/7,1,0),
              guideline_inferred = ifelse(inferred_accel >= 150/7,1,0))]

# Standardize inferred acceleration for analysis
prs_no_accel[,inferred_accel_std := (inferred_accel - mean(inferred_accel))/sd(inferred_accel)]

# Save out for seed file
write.csv(prs_no_accel,file='/Volumes/medpop_afib/skhurshid/acceleration_phewas/no_accel_mvpa_prs_phenotype.csv',row.names=F)

# Create inferred cutoff variable
prs_no_accel[,inferred_cutoff := ifelse(inferred_accel >= (150/7),1,0)]
prs_no_accel[,guideline_log := ifelse(exp(inferred_accel_log) >= (150/7),1,0)]

# Save out for seed file
write.csv(prs_no_accel,file='/Volumes/medpop_afib/skhurshid/acceleration_phewas/no_accel_inferred_cutoff_phenotype.csv',row.names=F)

