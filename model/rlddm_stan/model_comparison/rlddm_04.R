.libPaths()

#assign(".lib.loc", "C:/Program Files/R/R-3.5.2/library", envir = environment(.libPaths))

# after successful installation, load required packages
library("StanHeaders")
library("rstan")
options(mc.cores = 4)
library("Rcpp")
library("hBayesDM")
library("boot")
library("readr")
library("tidyr")
library("dplyr")
library("loo")

library("bayesplot")
library("rstanarm")

##########################
# LOADING AND PREPARING  #
##########################

gather_data <- function(files){
  # summarize all data in 1 data frame
  datalist <- list()
  for (i in 1:length(files)){
    no_col <- max(count.fields(files[i], sep = "\t"))
    D <- read_delim(
      files[i],"\t", escape_double = FALSE, locale = locale(), trim_ws = TRUE)
    # change 1,12 for kids data!
    D <- cbind(rep(substr(files[i],1,12),dim(D)[1]),D)
    datalist[[i]] <- D
  }
  transformed <- data.table::rbindlist(datalist) # combine all data frames in on
  return(transformed)
}

## define paths
path <- "C:/Rscripts"
model_path <- paste0(path,"/model/rlddm_stan/rlddm_04.stan")
data_path <- paste0(path,"/data/piloting/pilots_biokurs/subjects")


### load data
setwd(data_path)
files <- dir(pattern=".txt", recursive=TRUE)
raw_data <- gather_data(files)

colnames(raw_data)[1] <- "subjID"
raw_data$rt <- raw_data$rt/1000
names(raw_data)[names(raw_data)=="rt"] <- "RT"

DT_trials <- raw_data[, .N, by = subjID]
subjs <- DT_trials$subjID
n_subj    <- length(subjs)

# automatically filter missed responses (since RT = 0)
raw_data <- raw_data[which(raw_data$RT > 0.15),]
raw_data$trial <- as.integer(raw_data$trial)
#raw_data$trial_subj <- rep("NA",nrow(raw_data))
# since we discarded some observations, we have to assign new trial numbers 
for (subj in subjs){
  sub <- which(raw_data$subjID==subj)
  raw_data[sub,]$trial <- seq.int(nrow(raw_data[sub,]))
}

DT_trials_per_block <- raw_data[, .N, by = list(subjID,block)]
raw_data$block <- as.factor(raw_data$block)

# rename blocks
for (subj in subjs){
  sub <- which(raw_data$subjID==subj)
  n1 <- DT_trials_per_block[which(DT_trials_per_block==subj),]$N[1]
  b1 <-rep("1",n1)
  n2 <- DT_trials_per_block[which(DT_trials_per_block==subj),]$N[2]
  b2 <- rep("2",n2)
  n3 <- DT_trials_per_block[which(DT_trials_per_block==subj),]$N[3]
  b3 <-rep("3",n3)
  raw_data[sub]$block <- c(b1,b2,b3)
}
raw_data$block <- as.integer(raw_data$block)


# raw data: fb = 0 incorrect, fb = 1 correct, (fb = 2 missed)
# encoding for simulation: lower (incorrect) response=1, upper (correct) response =2 
raw_data$response = raw_data$fb+1
# raw_data$nonresponse = abs(raw_data$fb-2) # not used atm

raw_data$aStim <- as.double(raw_data$aStim)
# split vstim columns
raw_data <- raw_data %>% separate(vStims, c("vStim1", "vStim2"),sep="\\_")
raw_data$vStim1 <- as.double(raw_data$vStim1)
raw_data$vStim2 <- as.double(raw_data$vStim2)

# assign every stimulus pair for each block a unique number
raw_data[which(raw_data$block==2),]$aStim = raw_data[which(raw_data$block==2),]$aStim + 8
raw_data[which(raw_data$block==2),]$vStim1 = raw_data[which(raw_data$block==2),]$vStim1 + 8
raw_data[which(raw_data$block==2),]$vStim2 = raw_data[which(raw_data$block==2),]$vStim2 + 8

raw_data[which(raw_data$block==3),]$aStim = raw_data[which(raw_data$block==3),]$aStim + 16
raw_data[which(raw_data$block==3),]$vStim1 = raw_data[which(raw_data$block==3),]$vStim1 + 16
raw_data[which(raw_data$block==3),]$vStim2 = raw_data[which(raw_data$block==3),]$vStim2 + 16

raw_data$vStimNassoc <- ifelse(raw_data$aStim==raw_data$vStim1,raw_data$vStim2,raw_data$vStim1)

DT_trials <- raw_data[, .N, by = subjID]

# get minRT
minRT <- with(raw_data, aggregate(RT, by = list(y = subjID), FUN = min)[["x"]])
ifelse(is.null(dim(minRT)),minRT<-as.array(minRT))


first <- which(raw_data$trial==1)
# if N=1 transform int to 1-d array
first<-as.array(first)
# last is a Sx1 matrix identifying all last trials of a subject for each choice
last <- (first + DT_trials$N - 1)
# if N=1 transform int to 1-d array
last<-as.array(last)
# define the values for the rewards: if upper resp, value = 1
value <- ifelse(raw_data$response==2, 1, 0)
n_trials <- nrow(raw_data)

#blocks <- tapply(raw_data$block,raw_data$subjID, max,simplify = TRUE)
blocks <- aggregate( raw_data$block ~ raw_data$subjID, FUN = max )
blocks <- blocks$`raw_data$block`
# if N=1 transform int to 1-d array
ifelse(is.null(dim(blocks)),blocks<-as.array(blocks))

stims_per_block <- 8
dat <- list("N" = n_subj, "T"=n_trials,"RTbound" = 0.15,"minRT" = minRT, "iter" = raw_data$trial, "response" = raw_data$response, 
            "stim_assoc" = raw_data$aStim, "stim_nassoc" = raw_data$vStimNassoc, "RT" = raw_data$RT, "first" = first, "last" = last, "value"=value, "n_stims"=stims_per_block*blocks)  # names list of numbers


###############
# LOAD MODEL  #
###############

stanmodel <- rstan::stan_model(model_path)

fit <- rstan::sampling(object  = stanmodel,
                       data    = dat,
                       pars    = c("alpha","a_mod","v_mod","tau","eta_pos","eta_neg","assoc_active_pair","assoc_inactive_pair","delta_hat","pe_hat","mu_alpha","mu_eta_neg","mu_eta_pos","mu_v_mod","mu_a_mod","mu_tau","log_lik"),
                       #init    = "random",
                       chains  = 3,
                       iter    = 3000,
                       warmup  = 1000,
                       thin    = 1,
                       init_r = 0.05,
                       save_warmup = FALSE,
                       control = list(adapt_delta   = 0.999,
                                      stepsize      = 0.05,
                                      max_treedepth = 30),
                       verbose =FALSE)

# save model
saveRDS(fit, "fit_03_notconverged.rds")

# fit <- readRDS("fit1.rds") 

parVals <- rstan::extract(fit, permuted = TRUE)
fit_summary <- rstan::summary(fit)

parVals$lp__

head(fit_summary$summary)
tail(fit_summary$summary)

lp_rlddm03 <- log_posterior(fit)
head(lp_rlddm03)

posterior_rlddm02 <- as.array(fit)
color_scheme_set("darkgray")
#mcmc_parcoord(posterior_rlddm02)


#####################################
# EXTRACT AND WRITE OUT PARAMETERS  #
#####################################


alpha <- rep("NA",dat$N)
for (i in 1:dat$N){
  index <- paste0("alpha[",i,"]")
  alpha[i] <- fit_summary$summary[index,1]
}
alpha <- as.double(alpha)

alpha_mod <- rep("NA",dat$N)
for (i in 1:dat$N){
  index <- paste0("a_mod[",i,"]")
  alpha_mod[i] <- fit_summary$summary[index,1]
}
alpha_mod <- as.double(alpha_mod)

drift_mod <- rep("NA",dat$N)
for (i in 1:dat$N){
  index <- paste0("v_mod[",i,"]")
  drift_mod[i] <- fit_summary$summary[index,1]
}
drift_mod <- as.double(drift_mod)

tau <- rep("NA",dat$N)
for (i in 1:dat$N){
  index <- paste0("tau[",i,"]")
  tau[i] <- fit_summary$summary[index,1]
}
tau <- as.double(tau)

eta_pos <- rep("NA",dat$N)
for (i in 1:dat$N){
  index <- paste0("eta_pos[",i,"]")
  eta_pos[i] <- fit_summary$summary[index,1]
}
eta_pos <- as.double(eta_pos)

eta_neg <- rep("NA",dat$N)
for (i in 1:dat$N){
  index <- paste0("eta_neg[",i,"]")
  eta_neg[i] <- fit_summary$summary[index,1]
}
eta_neg <- as.double(eta_neg)

# delta = average of upper and lower response boundary
deltas <- rep("NA",dat$T)
for (i in 1:dat$T){
  index <- paste0("delta_hat[",i,"]")
  deltas[i] <- fit_summary$summary[index,1]
}
deltas <- as.double(deltas)

# association strength of active (i.e. correct) stimulus pair
assoc_active_pair <- rep("NA",dat$T)
for (i in 1:dat$T){
  index <- paste0("assoc_active_pair[",i,"]")
  assoc_active_pair[i] <- fit_summary$summary[index,1]
}
assoc_active_pair <- as.double(assoc_active_pair)

# association strength of inactive (i.e. incorrect) stimulus pair
assoc_inactive_pair <- rep("NA",dat$T)
for (i in 1:dat$T){
  index <- paste0("assoc_inactive_pair[",i,"]")
  assoc_inactive_pair[i] <- fit_summary$summary[index,1]
}
assoc_inactive_pair <- as.double(assoc_inactive_pair)

# prediction errors
pe_hat <- rep("NA",dat$T)
for (i in 1:dat$T){
  index <- paste0("pe_hat[",i,"]")
  pe_hat[i] <- fit_summary$summary[index,1]
}
pe_hat <- as.double(pe_hat)


# add the trial-by-trial regressors to the dataset
tbtregs <- cbind(raw_data[,1:9],raw_data[,18],raw_data[,14:16],pe=as.array(pe_hat),delta=as.array(deltas),drift=rep(0,dat$T),
                 assoc_active=as.array(assoc_active_pair),assoc_inactive=as.array(assoc_inactive_pair))
for(i in 1:n_subj){
  subj = as.character(subjs[i])
  tbtregs[which(tbtregs$subjID==subj),]$drift = tbtregs[which(tbtregs$subjID==subj),]$delta * drift_mod[i]
}

# write out parameters for each subject separately
for(i in 1:n_subj){
  subj = as.character(subjs[i])
  write.table(tbtregs[which(tbtregs$subjID==subj),], paste0(subj,"_pars_rlddm03",".csv"),
              quote=FALSE, sep=",", row.names=FALSE, col.names=TRUE)
}

# compute mean drift for each subject
mean_drift <-  rep("NA",dat$N)
mean_drift <- aggregate(list(drift=tbtregs$drift),list(subjID=tbtregs$subjID), mean)
subj_params <- cbind(mean_drift,alpha,alpha_mod,tau,drift_mod,eta_pos,eta_neg)
# write out subject parameters
write.table(subj_params, "subj_params_rlddm03_notconverged.csv", sep=",", row.names=FALSE, col.names=TRUE,quote = FALSE)

##########################
##### COMPARE MODELS #####
##########################

saved_model_path <- paste0(path,"/modeling_data/")
setwd(saved_model_path)

fit01 <- readRDS("fit_01_converged.rds")
loo1 <- loo(fit01, pars = "log_lik")

fit02 <- readRDS("fit_02_converged.rds")
loo1 <- loo(fit02, pars = "log_lik")
fit03 <- readRDS("fit_03_notconverged.rds")