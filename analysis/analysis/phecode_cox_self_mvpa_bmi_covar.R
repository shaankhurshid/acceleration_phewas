# Script to process phecode tables

# Depends
library(data.table)
library(stringr)
library(survival)

# Index of file names
list <- paste0('/mnt/ml4cvd/projects/skhurshid/accel_phewas/new_phecode_tables/',list.files('/mnt/ml4cvd/projects/skhurshid/accel_phewas/new_phecode_tables'))

# Load exposure/covariate data
value_bmi <- fread('/mnt/ml4cvd/projects/skhurshid/accel_phewas/cox_data_self_bmi_covar.csv')
setkey(value_bmi,sample_id)

# Format dates
for (j in (c('enroll_date','phenotype_censor_date'))){set(value_bmi,j=j,value=as.Date(value_bmi[[j]],format='%Y-%m-%d'))}

# Init vars
out <- data.table(); n <- 1

# Looping cox model
for (i in list){
# Merge
  phecode <- NULL; analysis_set <- NULL
  phecode <- read.table(i,sep='\t',header = TRUE); setDT(phecode)
  setkey(phecode,sample_id)
  analysis_set <- value_bmi[phecode,nomatch=0]
# Format variables
  analysis_set[,censor_date := as.Date(censor_date,format='%Y-%m-%d')]
# Create analysis variables
  analysis_set[,time_to_event := ifelse(c(has_disease == 0 | is.na(has_disease)),pmin(as.numeric(censor_date - enroll_date)/365.25,as.numeric(phenotype_censor_date - enroll_date)/365.25),
                                        as.numeric(censor_date - enroll_date)/365.25)]
# Remove prevalent disease or no follow-up
  analysis_set <- analysis_set[!is.na(time_to_event) & time_to_event > 0]
# Define events and follow-up
  disease <- analysis_set$disease[1]
  n_events <- nrow(analysis_set[has_disease==1])
  fu_median <- quantile(analysis_set$time_to_event,0.50); fu_q1 <- quantile(analysis_set$time_to_event,0.25); fu_q3 <- quantile(analysis_set$time_to_event,0.75)
# If less than 10 cases, abort
  if (n_events < 10){
    hr <- NA; lower <- NA; upper <- NA; p <- NA; z <- NA
    result <- data.table(disease,n_events,fu_median,fu_q1,fu_q3,hr,lower,upper,z,p)
    out <- rbind(out,result)
    print(paste0("Skipping phenotype ",analysis_set$disease[1]," since < 10 cases"))
    if (n %% 50 == 0){print(paste0("Just finished model ",n," out of ",length(list),"!"))}
    n <- n+1; next}
# Fit cox model
  model <- coxph(Surv(time_to_event,has_disease) ~ self_mvpa + enroll_age + sex + bmi + sbp + dbp + bpmed + tob + tdi + etoh_grams + diet + qual_ea, data=analysis_set)
  hr <- summary(model)$coefficients[1,2]; lower <- summary(model)$conf.int[1,3]
  upper <- summary(model)$conf.int[1,4]; z <- summary(model)$coefficients[1,4]
  p <- 2*pnorm(abs(z),lower.tail=FALSE)
  result <- data.table(disease,n_events,fu_median,fu_q1,fu_q3,hr,lower,upper,z,p)
  out <- rbind(out,result)
  if (n %% 50 == 0){print(paste0("Just finished model ",n," out of ",length(list),"!"))}
  n <- n+1
}

# Save out
write.csv(out,file='/mnt/ml4cvd/projects/skhurshid/accel_phewas/phecode_processed/cox_self_mvpa_bmi_covar.csv',row.names=F)