# Script to create Table 2

# Depends
library(stringr)

# Load data
sig_both <- fread(file='~/Documents/MGH Research/accel_phewas/sig_both_071121.csv')

# Start formatting
sig_both[,hr_formatted := paste0(round(i.hr,2),' (',round(i.lower,2),'-',round(i.upper,2),')')]
sig_both[,hr_genetic_formatted := paste0(round(hr,2),' (',round(lower,2),'-',round(upper,2),')')]

# Write out
write.csv(sig_both,file='~/Documents/MGH Research/accel_phewas/sig_both_formatted_071121.csv',row.names=F)