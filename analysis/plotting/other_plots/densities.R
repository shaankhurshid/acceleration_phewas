# Density plots

# Depends 
library(data.table)
library(ggplot2)

# Load data
cutoff <- fread(file='cox_data_mvpa_cutoff_bmi_covar.csv')
self_mvpa <- fread(file='cox_data_self_bmi_raw_covar.csv')
mvpa <- fread(file='mvpa_rate_bouted.tsv')
inferred <- fread(file='no_accel_mvpa_sqrt_prs_phenotype.csv')

# Compare distributions
# Actual
ggplot() + geom_density(data=mvpa,aes(x=mvpa$mvpa_rate_bouted*7),fill="#4575b4",alpha=0.55) +
  scale_x_continuous(breaks=seq(0,1000,100),expand=c(0,0),limits=c(0,1000)) +
  scale_y_continuous(breaks=seq(0,0.005,0.001),expand=c(0,0),limits=c(0,0.005)) +
  theme(panel.background=element_blank(),axis.line=element_line(color='black'),legend.position=c(0.80,0.90),
        axis.text=element_text(size=20,color='black'),plot.margin=unit(c(0.5,1,0.5,0.5),'cm'),
        axis.title.y = element_text(size=20,margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(size=20),legend.text=element_text(size=20)) +
  labs(x='MVPA minutes/week',y='Density') 
ggsave('accel_mvpa_density.pdf',
       height=2,width=2.5,units='in',scale=4)

# Actual self
ggplot() + geom_density(data=self_mvpa,adjust=1.4,aes(x=self_mvpa$mvpa_raw),fill="#a8ddb5",alpha=0.55) +
  scale_x_continuous(breaks=seq(0,1600,200),expand=c(0,0),limits=c(0,1600)) +
  scale_y_continuous(breaks=seq(0,0.005,0.001),expand=c(0,0),limits=c(0,0.005)) +
  theme(panel.background=element_blank(),axis.line=element_line(color='black'),legend.position=c(0.80,0.90),
        axis.text=element_text(size=20,color='black'),plot.margin=unit(c(0.5,1,0.5,0.5),'cm'),
        axis.title.y = element_text(size=20,margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(size=20),legend.text=element_text(size=20)) +
  geom_vline(xintercept=150,linetype=5,color='black') +
  labs(x='MVPA minutes/week',y='Density') 
ggsave('accel_self_mvpa_density.pdf',
       height=2,width=2.5,units='in',scale=4)

# Actual cutoff
ggplot() + geom_density(data=mvpa,aes(x=mvpa$mvpa_rate_bouted*7),fill="#4575b4",alpha=0.55) +
  scale_x_continuous(breaks=seq(0,1000,100),expand=c(0,0),limits=c(0,1000)) +
  scale_y_continuous(breaks=seq(0,0.005,0.001),expand=c(0,0),limits=c(0,0.005)) +
  theme(panel.background=element_blank(),axis.line=element_line(color='black'),legend.position=c(0.80,0.90),
        axis.text=element_text(size=20,color='black'),plot.margin=unit(c(0.5,1,0.5,0.5),'cm'),
        axis.title.y = element_text(size=20,margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(size=20),legend.text=element_text(size=20)) +
  geom_vline(xintercept=150,linetype=5,color='black') +
  labs(x='MVPA minutes/week',y='Density') 
ggsave('accel_mvpa_cutoff_density.pdf',
       height=2,width=2.5,units='in',scale=4)

# Inferred cutoff
ggplot() + geom_density(data=inferred,aes(x=inferred$inferred_accel*7),fill="#a8ddb5",alpha=0.55) +
  scale_x_continuous(breaks=seq(0,500,100),expand=c(0,0),limits=c(0,500)) +
  scale_y_continuous(breaks=seq(0,0.008,0.002),expand=c(0,0),limits=c(0,0.008)) +
  theme(panel.background=element_blank(),axis.line=element_line(color='black'),legend.position=c(0.80,0.90),
        axis.text=element_text(size=20,color='black'),plot.margin=unit(c(0.5,1,0.5,0.5),'cm'),
        axis.title.y = element_text(size=20,margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(size=20),legend.text=element_text(size=20)) +
  geom_vline(xintercept=150,linetype=5,color='black') +
  labs(x='MVPA minutes/week',y='Density') 
ggsave('accel_mvpa_density_inferred_sqrt.pdf',
       height=2,width=2.5,units='in',scale=4)
