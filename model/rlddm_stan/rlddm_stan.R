.libPaths()

assign(".lib.loc", "C:/Program Files/R/R-3.5.2/library", envir = environment(.libPaths))


#install.packages("rstan", lib="C:\\Program Files\\R\\R-3.5.2\\library")
#install.packages("StanHeaders", lib="C:\\Program Files\\R\\R-3.5.2\\library")
#install.packages("rstantools", lib="C:\\Program Files\\R\\R-3.5.2\\library")
#install.packages("hBayesDM", lib="C:\\Program Files\\R\\R-3.5.2\\library")
#install.packages("Rcpp", lib="C:\\Program Files\\R\\R-3.5.2\\library")

# after successful installation, load required packages
library("StanHeaders", lib.loc="C:/Program Files/R/R-3.5.2/library")
library("rstan", lib.loc="C:/Program Files/R/R-3.5.2/library")
library("Rcpp", lib.loc="C:/Program Files/R/R-3.5.2/library")
library("hBayesDM", lib.loc="C:/Program Files/R/R-3.5.2/library")
library("boot")
require("tidyr")

###
gather_data <- function(files){
  # summarize all data in 1 data frame
  datalist <- list()
  for (i in 1:length(files)){
    no_col <- max(count.fields(files[i], sep = "\t"))
    D <- read_delim(
      files[i],"\t", escape_double = FALSE, locale = locale(), trim_ws = TRUE)
    D <- cbind(rep(substr(files[i],1,12),dim(D)[1]),D)
    datalist[[i]] <- D
  }
  transformed <- data.table::rbindlist(datalist) # combine all data frames in on
  return(transformed)
}

get_astim_trials <- function(data){
  df_subj <- list()
  #data<-data[data$choice!=0,] # remove 'too slow ' responses 
  new_data <- data[FALSE,]
  new_col_trial <- data.frame(trial_astim = integer(0))
  new_data <- cbind(new_data,new_col_trial)
  for(i in unique(data$subjID)){
    df_subj[[i]] <- subset(data, subjID == i)
    for(j in unique(df_subj[[i]]$block)){
      df_subj_block <- subset(df_subj[[i]],block==j)
      # compute cumulative sum for each auditory stimulus in a given block of a subject
      for(k in unique(df_subj_block$aStim)){
        df_subj_block_astim <- list()
        df_subj_block_astim[[k]] <- subset(df_subj_block, aStim==k )
        new_col_trial <- 1:nrow(df_subj_block_astim[[k]])
        df_subj_block_astim[[k]]<- cbind(df_subj_block_astim[[k]],new_col_trial)
        colnames(df_subj_block_astim[[k]])[ncol(df_subj_block_astim[[k]])] <- "trial_astim"
        new_data <- rbind(new_data,df_subj_block_astim[[k]])
        # reorder data
        new_data <- new_data[
          with(new_data, order(subjID, block,trial)),
          ]
      }
    }
  }
  return(new_data)
}

## define paths
path <- "N:/Users/phaller/mri_task_analysis"
model_path <- paste0(path,"/model/rlddm_stan/rlddm_per_stimulus_v2.stan")
#data_path <- paste0(path,"/test_input.txt")
data_path <- paste0(path,"/data/piloting/pilots_biokurs/subjects")


### load data
setwd(data_path)
files <- dir(pattern=".txt", recursive=TRUE)
raw_data <- gather_data(files)

#data_path <- paste0(path,"/test_input.txt")
#setwd(path)
#raw_data <- data.table::fread(file = data_path, header = TRUE, sep = "\t", data.table = TRUE,
 #                             fill = TRUE, stringsAsFactors = TRUE, logical01 = FALSE)

#raw_data <- cbind(rep(substr(data_path,1,12),dim(raw_data)[1]),raw_data)
#D<-D[D$resp!=0,] # remove 'too slow ' responses 
### Rename and transform some columns
colnames(raw_data)[1] <- "subjID"

raw_data$rt <- raw_data$rt/1000
names(raw_data)[names(raw_data)=="rt"] <- "RT"
# automatically filter missed responses (since RT = 0)
raw_data <- raw_data[which(raw_data$RT > 0.15),]
#raw_data$subjID = rep('01',nrow(raw_data))

# raw data: fb = 0 incorrect, fb = 1 correct, (fb = 2 missed)
# encoding for simulation: lower (incorrect) response=1, upper (correct) response =2 
raw_data$response = raw_data$fb+1
# raw_data$nonresponse = abs(raw_data$fb-2) # not used atm

# split vstim columns
raw_data <- raw_data %>% separate(vStims, c("vStim1", "vStim2"),sep="\\_")
# get new column with non-associated stimulus
raw_data$vStimNassoc <- ifelse(raw_data$aStim==raw_data$vStim1,as.integer(raw_data$vStim2),as.integer(raw_data$vStim1))

raw_data <- get_astim_trials(raw_data)

## prepare data for jags
#raw_data$row <- seq.int(nrow(raw_data))
DT_trials <- raw_data[, .N, by = "subjID"]
subjs     <- DT_trials$subjID
n_subj    <- length(subjs)
# get minRT
minRT <- with(raw_data, aggregate(RT, by = list(y = subjID), FUN = min)[["x"]])
ifelse(is.null(dim(minRT)),minRT<-as.array(minRT))

# change blocks
for (subj in subjs){
  sub <- which(raw_data$subjID==subj)
  raw_data[sub,]$block <- as.factor(raw_data[sub,]$block)
  levels(raw_data[sub,]$block) <- c("1","2","3")
  raw_data[sub,]$block <- as.integer(raw_data[sub,]$block)
}

raw_data <- raw_data[which(raw_data$block==1),]
# since we discarded some observations, we have to assign new trial numbers 
for (subj in subjs){
  sub <- which(raw_data$subjID==subj)
  raw_data[sub,]$trial <- seq.int(nrow(raw_data[sub,]))
}

DT_trials <- raw_data[, .N, by = "subjID"]
DT_trials

# first is Sx1 matrix identifying all first trials of a subject for each choice
first <- which(raw_data$trial==1)
# if N=1 transform int to 1-d array
ifelse(is.null(dim(first)),first<-as.array(first))
# last is a Sx1 matrix identifying all last trials of a subject for each choice
last <- as.integer(first + DT_trials$N - 1)
ifelse(is.null(dim(last)),last<-as.array(last))
# incorrect is the inverse vector of choice and is needed to update the ev for the non-choices
raw_data$incorrect <- as.integer(ifelse(raw_data$correct==1, 0, 1))
# define the values for the rewards: if upper resp, value = 1
value <- ifelse(raw_data$response==2, 1, 0)
## all RT with negative choices -> -1
#new_RT <- ifelse(raw_data$correct==1, raw_data$RT*-1, raw_data$RT)
## # obs
n_trials <- nrow(raw_data)
## 


dat <- list("N" = n_subj, "T"=n_trials,"RTbound" = 0.15,"minRT" = minRT, "iter" = raw_data$trial, "response" = raw_data$response, "trial_astim" = raw_data$trial_astim,
            "stim_assoc" = raw_data$aStim, "stim_nassoc" = raw_data$vStimNassoc, "RT" = raw_data$RT, "first" = first, "last" = last, "value"=value, "n_stims"=8)  # names list of numbers

### with fixed learning rates ###

stanmodel_per_stimulus <- rstan::stan_model(model_path)

fit_invlog <- rstan::sampling(object  = stanmodel_per_stimulus,
                              data    = dat,
                              init    = "random",
                              chains  = 2,
                              iter    = 8000,
                              warmup  = 2000,
                              thin    = 1,
                              control = list(adapt_delta   = 0.95,
                                             stepsize      = 1,
                                             max_treedepth = 10),
                              verbose =TRUE)

parValsinvl <- rstan::extract(fit_invlog, permuted = TRUE)

fit_summary_invlog <- rstan::summary(fit_invlog)


### with fixed learning rates and updates for both associated and non-associated stimulus###

stanmodel_v2 <- rstan::stan_model(model_path)

fit_v2 <- rstan::sampling(object  = stanmodel_v2,
                          data    = dat,
                          init    = "random",
                          chains  = 3,
                          iter    = 10000,
                          warmup  = 3000,
                          thin    = 1,
                          control = list(adapt_delta   = 0.95,
                                         stepsize      = 1,
                                         max_treedepth = 10),
                          verbose =TRUE)

parVals_v2 <- rstan::extract(fit_v2, permuted = TRUE)

fit_summary_v2 <- rstan::summary(fit_v2)

head(fit_summary_v2$summary)
tail(fit_summary_v2$summary)

ev_mean <-  matrix(data= NA, nrow=dat$T, ncol=8)
ev_mean[1,] <- fit_summary_v2$summary[27:34,1]
for(i in 0:38){
  ev_mean[i,] <- fit_summary_v2$summary[(19+i*8):(26+i*8),1]
}


assoc_active_pair <- rep("NA",dat$T)
for (i in 1:dat$T){
  index <- paste0("assoc_active_pair[",i,"]")
  assoc_active_pair[i] <- fit_summary_v2$summary[index,1]
}

assoc_inactive_pair <- rep("NA",dat$T)
for (i in 1:dat$T){
  index <- paste0("assoc_inactive_pair[",i,"]")
  assoc_inactive_pair[i] <- fit_summary_v2$summary[index,1]
}

pe_pos_hat <- rep("NA",dat$T)
for (i in 1:dat$T){
  index <- paste0("pe_pos_hat[",i,"]")
  pe_pos_hat[i] <- fit_summary_v2$summary[index,1]
}

pe_neg_hat <- rep("NA",dat$T)
for (i in 1:dat$T){
  index <- paste0("pe_neg_hat[",i,"]")
  pe_neg_hat[i] <- fit_summary_v2$summary[index,1]
}

deltas <- rep("NA",dat$T)
for (i in 1:dat$T){
  index <- paste0("delta_hat[",i,"]")
  deltas[i] <- fit_summary_v2$summary[index,1]
}

v_mod <- rep("NA",dat$N)
for (i in 1:dat$N){
  index <- paste0("v_mod[",i,"]")
  v_mod[i] <- fit_summary_v2$summary[index,1]
}

write.csv(assoc_active_pair,paste(data_path,"/assoc_active_pair.csv",sep=""),quote=FALSE,row.names=FALSE)
write.csv(assoc_inactive_pair,paste(data_path,"/assoc_inactive_pair.csv",sep=""),quote=FALSE,row.names=FALSE)
write.csv(pe_hat,paste(data_path,"/pe_hat.csv",sep=""),quote=FALSE,row.names=FALSE)
write.csv(deltas,paste(data_path,"/deltas.csv",sep=""),quote=FALSE,row.names=FALSE)

v_mods <- rep("NA", 15)
for (i in 1:15){
  index <- paste0("v_mod[",i,"]")
  v_mods[i] <- fit_summary_v2$summary[index,1]
}
#########

rstan::stan_diag(fit_invlog, info = 'sample') # shows three plots together
rstan::stan_par(fit_invlog, par = "alpha")

parVals <- rstan::extract(fit, permuted = TRUE)

## now access the data of the fit 
print(names(parVals))

head(parVals$mu_pr)


fit_summary <- summary(fit)
# In fit_summary$summary all chains are merged whereas fit_summary$c_summary contains summaries for each chain individually. 
# Typically we want the summary for all chains merged
print(names(fit_summary))
print(fit_summary$summary)


mean(parVals$alpha)
mean(parVals$tau)
# eta is 3 dimensional: 1-d: iteration, 2-d: subject, 3-d: pos or neg eta
mean(parVals$eta_pos)
mean(parVals$eta_neg)
