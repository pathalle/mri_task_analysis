.libPaths()

assign(".lib.loc", "C:/Program Files/R/R-3.5.3/library", envir = environment(.libPaths))

# after successful installation, load required packages
library("StanHeaders", lib.loc="C:/Program Files/R/R-3.5.3/library")
library("rstan", lib.loc="C:/Program Files/R/R-3.5.3/library")
#Sys.setenv(LOCAL_CPPFLAGS = '-march=native')
library("Rcpp", lib.loc="C:/Program Files/R/R-3.5.3/library")
library("hBayesDM", lib.loc="C:/Program Files/R/R-3.5.3/library")
library("boot", lib.loc="C:/Program Files/R/R-3.5.3/library")
library("readr", lib.loc="C:/Program Files/R/R-3.5.3/library")
library("tidyr", lib.loc="C:/Program Files/R/R-3.5.3/library")
library("dplyr", lib.loc="C:/Program Files/R/R-3.5.3/library")

library("bayesplot", lib.loc="C:/Program Files/R/R-3.5.3/library")
library("rstanarm", lib.loc="C:/Program Files/R/R-3.5.3/library")

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

# not in use 
# get_astim_trials <- function(data){
#   df_subj <- list()
#   #data<-data[data$choice!=0,] # remove 'too slow ' responses 
#   new_data <- data[FALSE,]
#   new_col_trial <- data.frame(trial_astim = integer(0))
#   new_data <- cbind(new_data,new_col_trial)
#   for(i in unique(data$subjID)){
#     df_subj[[i]] <- subset(data, subjID == i)
#     for(j in unique(df_subj[[i]]$block)){
#       df_subj_block <- subset(df_subj[[i]],block==j)
#       # compute cumulative sum for each auditory stimulus in a given block of a subject
#       for(k in unique(df_subj_block$aStim)){
#         df_subj_block_astim <- list()
#         df_subj_block_astim[[k]] <- subset(df_subj_block, aStim==k )
#         new_col_trial <- 1:nrow(df_subj_block_astim[[k]])
#         df_subj_block_astim[[k]]<- cbind(df_subj_block_astim[[k]],new_col_trial)
#         colnames(df_subj_block_astim[[k]])[ncol(df_subj_block_astim[[k]])] <- "trial_astim"
#         new_data <- rbind(new_data,df_subj_block_astim[[k]])
#         # reorder data
#         new_data <- new_data[
#           with(new_data, order(subjID, block,trial)),
#           ]
#       }
#     }
#   }
#   return(new_data)
# }

## define paths
path <- "N:/Users/phaller/mri_task_analysis"
model_path <- paste0(path,"/model/rlddm_stan/rlddm_blocks.stan")
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



#################
# SIMULATE DATA #
#################

stanmodel <- rstan::stan_model(model_path)



###############
# LOAD MODEL  #
###############

stanmodel <- rstan::stan_model(model_path)

fit <- rstan::sampling(object  = stanmodel,
                       data    = dat,
                       init    = "random",
                       chains  = 4,
                       iter    = 2000,
                       warmup  = 1000,
                       thin    = 1,
                       control = list(adapt_delta   = 0.95,
                                      stepsize      = 1,
                                      max_treedepth = 10),
                       verbose =TRUE)

# save model
saveRDS(fit, "fit_vp01-17_october.rds")
fit <- readRDS("N:/Users/phaller/modeling_data/fit1.rds")
# fit <- readRDS("fit1.rds") 

parVals <- rstan::extract(fit, permuted = TRUE)
fit_summary <- rstan::summary(fit)

head(fit_summary$summary)
tail(fit_summary$summary)


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
tbtregs <- cbind(raw_data[,1:9],raw_data[,18],raw_data[,14:16],pe=as.array(pe_hat),delta=as.array(deltas),drift=rep(0,dat$T),assoc_active=as.array(assoc_active_pair),assoc_inactive=as.array(assoc_inactive_pair))
for(i in 1:n_subj){
  subj = as.character(subjs[i])
  tbtregs[which(tbtregs$subjID==subj),]$drift = tbtregs[which(tbtregs$subjID==subj),]$delta * drift_mod[i]
}

# write out parameters for each subject separately
for(i in 1:n_subj){
  subj = as.character(subjs[i])
  write.table(tbtregs[which(tbtregs$subjID==subj),], paste0(subj,"_pars",".csv"),
              quote=FALSE, sep=",", row.names=FALSE, col.names=TRUE)
}

# compute mean drift for each subject
mean_drift <-  rep("NA",dat$N)
mean_drift <- aggregate(list(drift=tbtregs$drift),list(subjID=tbtregs$subjID), mean)
subj_params <- cbind(mean_drift,alpha,alpha_mod,tau,drift_mod)
# write out subject parameters
write.table(subj_params, "subj_params_vp01-17.csv", sep=",", row.names=FALSE, col.names=TRUE,quote = FALSE)


### Analysis
path <- "N:/Users/phaller/modeling_data"
setwd(path)
fit <- readRDS("fit_pilots_biokurs_vp01-17.rds") 

