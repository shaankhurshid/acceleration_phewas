# Script to analyze PRS-pheWAS

# Depends 
library(data.table)
library(ggplot2)

# Load PRS
prs <- fread(file='/Volumes/medpop_afib/skhurshid/acceleration_phewas/prs/scores/mvpa_sqrt_sensi_prs4/mvpa_sqrt_sensi_prs4_processed.csv')

# Load accel genetic file
acceleration <- fread(file='/Volumes/medpop_afib/skhurshid/acceleration_phewas/mvpa_rate_bouted_sqrt_sensi.tsv')

# Join
setkey(prs,sample_id); setkey(acceleration,FID)
acceleration[prs,':='(prs_std = i.prs_std)]

# Remove no PRS (poor quality genetic data)
prs_set <- acceleration[!is.na(prs_std)]

# Check exposure ~ PRS
mod <- lm(mvpa_rate_bouted_sqrt ~ prs_std,data=prs_set)

# Now fit genetic acceleration model
mod_fitted <- lm(mvpa_rate_bouted_sqrt ~ prs_std + PC1 + PC2 + PC3 + PC4 + PC5 + array_UKBB + male + age_accel,data=prs_set)

# Now load the no accel sets
prs_na <- fread(file='/Volumes/medpop_afib/skhurshid/acceleration_phewas/prs/scores/mvpa_sqrt_sensi_prs_no_accel/mvpa_sqrt_sensi_prs_no_accel_processed.csv')
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
prs_no_accel[,inferred_accel := predict(mod_fitted,newdata=prs_no_accel)]

# Compare distributions
# Actual vs Inferred
x <- list(v1=prs_set$mvpa_rate_bouted_sqrt**2,v2=prs_set$inferred_accel**2)
data <- melt(x)

# Plot
ggplot() + geom_density(data=data,aes(x=value,fill=L1),alpha=0.55) + 
  scale_x_continuous(breaks=seq(0,150,50),expand=c(0,0),limits=c(0,150)) +
  scale_y_continuous(breaks=seq(0,0.07,0.01),expand=c(0,0),limits=c(0,0.07)) +
  scale_fill_manual(values=c("#2b8cbe","#f03b20"),name='',labels=c('MVPA','Inferred MVPA')) +
  theme(panel.background=element_blank(),axis.line=element_line(color='black'),legend.position=c(0.80,0.90),
        axis.text=element_text(size=20,color='black'),plot.margin=unit(c(0.5,1,0.5,0.5),'cm'),
        axis.title.y = element_text(size=20,margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(size=20),legend.text=element_text(size=20)) +
  geom_vline(xintercept=150/7,linetype=5,color='black') +
  labs(x='MVPA minutes/day',y='Density') 
ggsave('~/Documents/MGH Research/accel_phewas/accel_mvpa_sqrt_density_compare.pdf',
       height=2,width=2.5,units='in',scale=4)

# Inferred
ggplot() + geom_density(data=prs_no_accel,aes(x=prs_no_accel$inferred_accel),fill="#4575b4",alpha=0.55) +
  scale_x_continuous(breaks=seq(0,60,10),expand=c(0.01,0),limits=c(0,60)) +
  scale_y_continuous(breaks=seq(0,0.15,0.05),expand=c(0,0),limits=c(0,0.15)) +
  theme(panel.background=element_blank(),axis.line=element_line(color='black'),legend.position='None',
        axis.text=element_text(size=20,color='black'),plot.margin=unit(c(0.5,0.6,0.5,0.5),'cm'),
        axis.title.y = element_text(size=20,margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(size=20),legend.text=element_text(size=20)) +
  labs(x=expression(paste("Mean acceleration")),y='Density')
ggsave(filename='~/Documents/MGH Research/accel_phewas/accel_inferred_value_density.pdf',height=1.8,width=2.53,
       scale=4,device='pdf')

# Square and standardize inferred acceleration for analysis
prs_no_accel[,inferred_accel := prs_no_accel$inferred_accel**2]
prs_no_accel[,inferred_accel_std := (inferred_accel - mean(inferred_accel))/sd(inferred_accel)]

# Figure out expected counts for guidelines
prs_no_accel[,':='(guideline_inferred = ifelse(inferred_accel >= 150/7,1,0))]
prs_no_accel[,':='(inferred_75 = ifelse(inferred_accel >= 75/7,1,0))]
prs_no_accel[,':='(inferred_300 = ifelse(inferred_accel >= 300/7,1,0))]

# Save out for seed file
write.csv(prs_no_accel,file='/Volumes/medpop_afib/skhurshid/acceleration_phewas/no_accel_mvpa_sqrt_sensi_prs_phenotype.csv',row.names=F)

