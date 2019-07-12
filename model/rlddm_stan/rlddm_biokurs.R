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
library("readr")
library("tidyr")
library("dplyr")

###
gather_data <- function(files){
  # summarize all data in 1 data frame
  datalist <- list()
  for (i in 1:length(files)){
    no_col <- max(count.fields(files[i], sep = "\t"))
    D <- read_delim(
      files[i],"\t", escape_double = FALSE, locale = locale(), trim_ws = TRUE)
    #D <- cbind(rep(substr(files[i],1,12),dim(D)[1]),D) for adult pilots
    D <- cbind(rep(substr(files[i],17,22),dim(D)[1]),D)
    datalist[[i]] <- D
  }
  transformed <- data.table::rbindlist(datalist) # combine all data frames in on
  return(transformed)
}

#get_astim_trials <- function(data){
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
model_path <- paste0(path,"/model/rlddm_stan/rlddm_blocks.stan")
#data_path <- paste0(path,"/test_input.txt")
data_path <- paste0(path,"/data/piloting/piloting_kloten/2x4")
#data_path <- paste0(path,"/data/piloting/pilots_biokurs")


### load data
setwd(data_path)
files <- dir(pattern=".txt", recursive=TRUE)
raw_data <- gather_data(files)

###################### for test file
setwd(path)

colnames(raw_data)
raw_data <- data.table::fread(file = data_path, header = TRUE, sep = "\t", data.table = TRUE,
                              fill = TRUE, stringsAsFactors = TRUE, logical01 = FALSE)

raw_data <- cbind(rep(substr(data_path,1,12),dim(raw_data)[1]),raw_data)
#D<-D[D$resp!=0,] # remove 'too slow ' responses 
### Rename and transform some columns
#######################


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


# rename blocks
for (subj in subjs){
  sub <- which(raw_data$subjID==subj)
  raw_data[sub,]$block <- as.factor(raw_data[sub,]$block)
  #levels(raw_data[sub,]$block) <- c("1","2","3")
  levels(raw_data[sub,]$block) <- c("1","2")
  raw_data[sub,]$block <- as.integer(raw_data[sub,]$block)
}


# raw data: fb = 0 incorrect, fb = 1 correct, (fb = 2 missed)
# encoding for simulation: lower (incorrect) response=1, upper (correct) response =2 
raw_data$response = raw_data$fb+1
# raw_data$nonresponse = abs(raw_data$fb-2) # not used atm

raw_data$aStim <- as.double(raw_data$aStim)
# split vstim columns
raw_data <- raw_data %>% separate(vStims, c("vStim1", "vStim2"),sep="\\_")
raw_data$vStim1 <- as.double(raw_data$vStim1)
raw_data$vStim2 <- as.double(raw_data$vStim2)
# get new column with non-associated stimulus

raw_data[which(raw_data$block==2),]$aStim = raw_data[which(raw_data$block==2),]$aStim + 4
raw_data[which(raw_data$block==2),]$vStim1 = raw_data[which(raw_data$block==2),]$vStim1 + 4
raw_data[which(raw_data$block==2),]$vStim2 = raw_data[which(raw_data$block==2),]$vStim2 + 4

#raw_data[which(raw_data$block==2),]$aStim = raw_data[which(raw_data$block==2),]$aStim + 8
#raw_data[which(raw_data$block==2),]$vStim1 = raw_data[which(raw_data$block==2),]$vStim1 + 8
#raw_data[which(raw_data$block==2),]$vStim2 = raw_data[which(raw_data$block==2),]$vStim2 + 8

#raw_data[which(raw_data$block==3),]$aStim = raw_data[which(raw_data$block==3),]$aStim + 16
#raw_data[which(raw_data$block==3),]$vStim1 = raw_data[which(raw_data$block==3),]$vStim1 + 16
#raw_data[which(raw_data$block==3),]$vStim2 = raw_data[which(raw_data$block==3),]$vStim2 + 16

raw_data$vStimNassoc <- ifelse(raw_data$aStim==raw_data$vStim1,raw_data$vStim2,raw_data$vStim1)

#raw_data <- get_astim_trials(raw_data)

DT_trials <- raw_data[, .N, by = subjID]

# get minRT
minRT <- with(raw_data, aggregate(RT, by = list(y = subjID), FUN = min)[["x"]])
ifelse(is.null(dim(minRT)),minRT<-as.array(minRT))


first <- which(raw_data$trial==1)
# if N=1 transform int to 1-d array
first<-as.array(first)
# last is a Sx1 matrix identifying all last trials of a subject for each choice
last <- (first + DT_trials$N - 1)
last<-as.array(last)
# incorrect is the inverse vector of choice and is needed to update the ev for the non-choices
#raw_data$incorrect <- as.integer(ifelse(raw_data$correct==1, 0, 1))
# define the values for the rewards: if upper resp, value = 1
value <- ifelse(raw_data$response==2, 1, 0)


n_trials <- nrow(raw_data)

#blocks <- tapply(raw_data$block,raw_data$subjID, max,simplify = TRUE)
blocks <- aggregate( raw_data$block ~ raw_data$subjID, FUN = max )
blocks <- blocks$`raw_data$block`
ifelse(is.null(dim(blocks)),blocks<-as.array(blocks))


stims_per_block <- 4
dat <- list("N" = n_subj, "T"=n_trials,"RTbound" = 0.15,"minRT" = minRT, "iter" = raw_data$trial, "response" = raw_data$response, 
            "stim_assoc" = raw_data$aStim, "stim_nassoc" = raw_data$vStimNassoc, "RT" = raw_data$RT, "first" = first, "last" = last, "value"=value, "n_stims"=stims_per_block*blocks)  # names list of numbers
dat
### with fixed learning rates ###

stanmodel <- rstan::stan_model(model_path)

fit_children <- rstan::sampling(object  = stanmodel,
                       data    = dat,
                       init    = "random",
                       chains  = 2,
                       iter    = 7000,
                       warmup  = 2000,
                       thin    = 1,
                       control = list(adapt_delta   = 0.95,
                                      stepsize      = 1,
                                      max_treedepth = 10),
                       verbose =TRUE)

# save model
#saveRDS(fit, "fit1.rds")
saveRDS(fit_children, "fit_children_2x4.rds")
# fit <- readRDS("fit1.rds") 

parVals <- rstan::extract(fit, permuted = TRUE)

fit_summary <- rstan::summary(fit)

head(fit_summary$summary)
tail(fit_summary$summary)

assoc_active_pair <- rep("NA",dat$T)
for (i in 1:dat$T){
  index <- paste0("assoc_active_pair[",i,"]")
  assoc_active_pair[i] <- fit_summary$summary[index,1]
}
assoc_active_pair <- as.double(assoc_active_pair)

assoc_inactive_pair <- rep("NA",dat$T)
for (i in 1:dat$T){
  index <- paste0("assoc_inactive_pair[",i,"]")
  assoc_inactive_pair[i] <- fit_summary$summary[index,1]
}
assoc_inactive_pair <- as.double(assoc_inactive_pair)

pe_hat <- rep("NA",dat$T)
for (i in 1:dat$T){
  index <- paste0("pe_hat[",i,"]")
  pe_hat[i] <- fit_summary$summary[index,1]
}
pe_hat <- as.double(pe_hat)

deltas <- rep("NA",dat$T)
for (i in 1:dat$T){
  index <- paste0("delta_hat[",i,"]")
  deltas[i] <- fit_summary$summary[index,1]
}
deltas <- as.double(deltas)

v_mod <- rep("NA",dat$N)
for (i in 1:dat$N){
  index <- paste0("v_mod[",i,"]")
  v_mod[i] <- fit_summary$summary[index,1]
}
v_mod <- as.double(v_mod)

raw_data <- raw_data[,-19]
data_out <- cbind(raw_data,assoc_active=as.array(assoc_active_pair),assoc_inactive=as.array(assoc_inactive_pair),deltas=as.array(deltas),drift=rep("NA",dat$T))
for(i in 1:n_subj){
  subj = as.character(subjs[i])
  #data_out[which(data_out$subjID==subj),]$drift = data_out[which(data_out$subjID==subj),]$deltas * v_mod[i]
  data_out[which(data_out$subjID==subj),]$drift = as.double(data_out[which(data_out$subjID==subj),]$drift)
  write.table(data_out[which(data_out$subjID==subj),], paste0(subj,"_params",".csv"), sep=",", row.names=FALSE, col.names=TRUE)
}



write.csv(assoc_active_pair,paste(data_path,"/assoc_active_pair.csv",sep=""),quote=FALSE,row.names=FALSE)
write.csv(assoc_inactive_pair,paste(data_path,"/assoc_inactive_pair.csv",sep=""),quote=FALSE,row.names=FALSE)
write.csv(pe_hat,paste(data_path,"/pe_hat.csv",sep=""),quote=FALSE,row.names=FALSE)
write.csv(deltas,paste(data_path,"/deltas.csv",sep=""),quote=FALSE,row.names=FALSE)

#######################################  fit children data  #####################################################

fit_children <- rstan::sampling(object  = stanmodel,
                                data    = dat,
                                init    = "random",
                                chains  = 2,
                                iter    = 7000,
                                warmup  = 2000,
                                thin    = 1,
                                control = list(adapt_delta   = 0.95,
                                               stepsize      = 1,
                                               max_treedepth = 10),
                                verbose =TRUE)

# save model
#saveRDS(fit, "fit1.rds")
# fit <- readRDS("fit1.rds") 

parValsCh <- rstan::extract(fit_children, permuted = TRUE)

fit_summary_ch <- rstan::summary(fit_children)

fit_summary_ch$summary[70:90,]

alpha <- rep("NA",dat$N)
for (i in 1:dat$N){
  index <- paste0("alpha[",i,"]")
  alpha[i] <- fit_summary_ch$summary[index,1]
}
alpha <- as.double(alpha)

alpha_mod <- rep("NA",dat$N)
for (i in 1:dat$N){
  index <- paste0("a_mod[",i,"]")
  alpha_mod[i] <- fit_summary_ch$summary[index,1]
}
alpha_mod <- as.double(alpha_mod)

drift_mod <- rep("NA",dat$N)
for (i in 1:dat$N){
  index <- paste0("v_mod[",i,"]")
  drift_mod[i] <- fit_summary_ch$summary[index,1]
}
drift_mod <- as.double(drift_mod)

tau <- v_mod <- rep("NA",dat$N)
for (i in 1:dat$N){
  index <- paste0("tau[",i,"]")
  tau[i] <- fit_summary_ch$summary[index,1]
}
tau <- as.double(tau)

deltas <- rep("NA",dat$T)
for (i in 1:dat$T){
  index <- paste0("delta_hat[",i,"]")
  deltas[i] <- fit_summary$summary_ch[index,1]
}
deltas <- as.double(deltas)


assoc_active_pair <- rep("NA",dat$T)
for (i in 1:dat$T){
  index <- paste0("assoc_active_pair[",i,"]")
  assoc_active_pair[i] <- fit_summary$summary_ch[index,1]
}
assoc_active_pair <- as.double(assoc_active_pair)

assoc_inactive_pair <- rep("NA",dat$T)
for (i in 1:dat$T){
  index <- paste0("assoc_inactive_pair[",i,"]")
  assoc_inactive_pair[i] <- fit_summary$summary_ch[index,1]
}
assoc_inactive_pair <- as.double(assoc_inactive_pair)

pe_hat <- rep("NA",dat$T)
for (i in 1:dat$T){
  index <- paste0("pe_hat[",i,"]")
  pe_hat[i] <- fit_summary$summary_ch[index,1]
}
pe_hat <- as.double(pe_hat)

tbtregs <- cbind(raw_data[,1:3],deltas=as.array(deltas),drift=rep("NA",dat$T),assoc_active=as.array(assoc_active_pair),assoc_inactive=as.array(assoc_inactive_pair))
mean_drift <-  rep("NA",dat$N)
for(i in 1:n_subj){
  subj = as.character(subjs[i])
  tbtregs[which(tbtregs$subjID==subj),]$drift = tbtregs[which(tbtregs$subjID==subj),]$deltas * drift_mod[i]
  tbtregs[which(tbtregs$subjID==subj),]$drift = as.double(tbtregs[which(tbtregs$subjID==subj),]$drift)
  mean_drift[i] = mean(tbtregs[which(tbtregs$subjID==subj),]$drift)
  #write.table(data_out[which(data_out$subjID==subj),], paste0(subj,"_params",".csv"), sep=",", row.names=FALSE, col.names=TRUE)
}

tbtregs$drift <- as.double(tbtregs$drift)
mean_drift <- aggregate(list(drift=tbtregs$drift),list(subjID=tbtregs$subjID), mean)
subj_params <- cbind(mean_drift,alpha,alpha_mod,tau,drift_mod)

raw_data <- raw_data[,-19]
data_out <- cbind(raw_data,assoc_active=as.array(assoc_active_pair),assoc_inactive=as.array(assoc_inactive_pair),deltas=as.array(deltas),drift=rep("NA",dat$T))
for(i in 1:n_subj){
  subj = as.character(subjs[i])
  #data_out[which(data_out$subjID==subj),]$drift = data_out[which(data_out$subjID==subj),]$deltas * v_mod[i]
  data_out[which(data_out$subjID==subj),]$drift = as.double(data_out[which(data_out$subjID==subj),]$drift)
  write.table(data_out[which(data_out$subjID==subj),], paste0(subj,"_params",".csv"), sep=",", row.names=FALSE, col.names=TRUE)
}



write.csv(assoc_active_pair,paste(data_path,"/assoc_active_pair.csv",sep=""),quote=FALSE,row.names=FALSE)
write.csv(assoc_inactive_pair,paste(data_path,"/assoc_inactive_pair.csv",sep=""),quote=FALSE,row.names=FALSE)
write.csv(pe_hat,paste(data_path,"/pe_hat.csv",sep=""),quote=FALSE,row.names=FALSE)
write.csv(deltas,paste(data_path,"/deltas.csv",sep=""),quote=FALSE,row.names=FALSE)