# Script to analyze PRS-pheWAS

# Depends 
library(data.table)
library(ggplot2)

# Load PRS
prs <- fread(file='/Volumes/medpop_afib/skhurshid/acceleration_phewas/prs/scores/prs4_white/prs4_white_processed.csv')

# Load accel genetic file
acceleration <- fread(file='/Volumes/medpop_afib/skhurshid/acceleration_phewas/value_white.tsv')

# Join
setkey(prs,sample_id); setkey(acceleration,FID)
acceleration[prs,':='(prs_std = i.prs_std)]

# Check exposure ~ PRS
mod <- lm(value ~ prs_std,data=acceleration)

# Now fit genetic acceleration model
mod_fitted <- lm(value ~ prs_std + PC1 + PC2 + PC3 + PC4 + PC5 + array_UKBB + male + age_accel,data=acceleration)

# Now load the no accel sets
prs_na <- fread(file='/Volumes/medpop_afib/skhurshid/acceleration_phewas/prs/scores/prs_no_accel_white/prs_no_accel_white_processed.csv')
prs_covar <- fread(file='/Volumes/medpop_afib/skhurshid/acceleration_phewas/prs/no_accel_white_seed.tsv')

# Merge
setkey(prs_na,sample_id); setkey(prs_covar,sample_id)
prs_no_accel <- prs_covar[prs_na,nomatch=0]

# Reduce feature space
prs_no_accel <- prs_no_accel[,c('sample_id','sex','enroll_age','array_UKBB','prs_std',paste0('PC',1:5))]
setnames(prs_no_accel,c('enroll_age','sex'),c('age_accel','male'))
prs_no_accel[,male := ifelse(male=='Male',1,0)]

# Get inferred acceleration for analysis
prs_no_accel[,inferred_accel := predict(mod_fitted,newdata=prs_no_accel)]

# Compare distributions
# Actual
ggplot() + geom_density(data=acceleration_pheno,aes(x=acceleration_pheno$value),fill="#f03b20",alpha=0.55) +
  scale_x_continuous(breaks=seq(0,60,10),expand=c(0.01,0),limits=c(0,60)) +
  scale_y_continuous(breaks=seq(0,0.07,0.01),expand=c(0,0),limits=c(0,0.07)) +
  theme(panel.background=element_blank(),axis.line=element_line(color='black'),legend.position='None',
        axis.text=element_text(size=20,color='black'),plot.margin=unit(c(0.5,0.6,0.5,0.5),'cm'),
        axis.title.y = element_text(size=20,margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(size=20),legend.text=element_text(size=20)) +
  labs(x=expression(paste("Mean acceleration")),y='Density')
ggsave(filename='~/Documents/MGH Research/accel_phewas/accel_value_density.pdf',height=1.8,width=2.53,
       scale=4,device='pdf')

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

# Save out for seed file
write.csv(prs_no_accel,file='/Volumes/medpop_afib/skhurshid/acceleration_phewas/no_accel_prs_phenotype.csv',row.names=F)

