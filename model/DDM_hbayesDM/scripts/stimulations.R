library(rstan)
library(RWiener)
library(boot) #needed for inverse logit

path <- dirname(rstudioapi::getActiveDocumentContext()$path)
model_path <- paste0(path,"/rlddm.stan")
data_path <- paste0(path,"/test_input.txt")
setwd(path)
raw_data <- data.table::fread(file = data_path, header = TRUE, sep = "\t", data.table = TRUE,
                              fill = TRUE, stringsAsFactors = TRUE, logical01 = FALSE)

raw_data <- cbind(rep(substr(data_path,1,12),dim(raw_data)[1]),raw_data)
#D<-D[D$resp!=0,] # remove 'too slow ' responses 
### Rename and transform some columns
colnames(raw_data)[1] <- "subjID"

raw_data$rt <- raw_data$rt/1000
names(raw_data)[names(raw_data)=="rt"] <- "RT"
raw_data <- raw_data[which(raw_data$RT > 0.15),]
raw_data$subjID = rep('01',nrow(raw_data))

names(raw_data)[names(raw_data)=="fb"] <- "correct"

## prepare data for jags
#raw_data$row <- seq.int(nrow(raw_data))
DT_trials <- raw_data[, .N, by = "subjID"]
subjs     <- DT_trials$subjID
n_subj    <- length(subjs)
# get minRT
minRT <- with(raw_data, aggregate(RT, by = list(y = subjID), FUN = min)[["x"]])
# assign new trial number for excluded decisions

for (subj in subjs){
  sub <- which(raw_data$subjID==subj)
  raw_data[sub,]$trial <- seq.int(nrow(raw_data[sub,]))
}
# first is Sx1 matrix identifying all first trials of a subject for each choice
first <- which(raw_data$trial==1)
# last is a Sx1 matrix identifying all last trials of a subject for each choice
last <- as.integer(first + DT_trials$N - 1)
# incorrect is the inverse vector of choice and is needed to update the ev for the non-choices
raw_data$incorrect <- as.integer(ifelse(raw_data$correct==1, 0, 1))
# define the values for the rewards
value <- ifelse(raw_data$correct==1, 1, 0)
## all RT with negative choices -> -1
#new_RT <- ifelse(raw_data$correct==1, raw_data$RT*-1, raw_data$RT)
## # obs
n_trials <- nrow(raw_data)
## 
stims <- raw_data$aStim

# encoding for simulation: lower response=1, upper response =2 
raw_data$response = raw_data$correct+1
raw_data$nonresponse = abs(raw_data$correct-2)

dat <- list("N" = n_subj, "T"=n_trials,"RTbound" = 0.15,"minRT" = minRT, "iter" = raw_data$trial, "response" = raw_data$response,"nonresponse" = raw_data$nonresponse,
            "RT" = raw_data$RT, "first" = first, "last" = last, "value"=value, "stims" = stims)  # names list of numbers


############ SIMULATE MODEL ###############

## Hyperparameters (group)
mu_pr <- rnorm(6,0,1)
aa <- rnorm(n = 6, mean = 0, s = 0.2)
while(any(aa<0)) { aa <- rnorm(n = 6, mean = 0, s = 0.2) }
sigma <- aa


## Priors
eta_pr = matrix(data= NA, nrow=dat$N, ncol=2)
eta = matrix(data= NA, nrow=dat$N, ncol=2)
a_mod_pr = matrix(data= NA, nrow=dat$N, ncol=1)
a_mod = matrix(data= NA, nrow=dat$N, ncol=1)
v_mod_pr = matrix(data= NA, nrow=dat$N, ncol=1)
v_mod = matrix(data= NA, nrow=dat$N, ncol=1)
tau_pr = matrix(data= NA, nrow=dat$N, ncol=1)
tau = matrix(data= NA, nrow=dat$N, ncol=1)


for(s in 1:dat$N){
  aa <- rnorm(n = s, mean = 0, s = 1)
  while(any(aa<0)) { aa <- rnorm(n = s, mean = 0, s = 1) }
  alpha_pr <- aa
  alpha = exp(mu_pr[1] + sigma[1] * alpha_pr)
  
  for (i in 1:2){
    eta_pr[s,i] <- rnorm(1,0,1)
    eta[s,i] <- exp(mu_pr[2] + sigma[2] * eta_pr[s,i])
  }
  aa <- rnorm(n = s, mean = 0, s = 1)
  while(any(aa< -0.5 | aa>2)) { aa <- rnorm(n = s, mean = 0, s = 1) }
  a_mod_pr <- aa
  a_mod = exp(mu_pr[4] + sigma[4] * a_mod_pr)
  
  aa <- rnorm(n = s, mean = 0, s = 1)
  while(any(aa< 0 | aa>10)) { aa <- rnorm(n = s, mean = 0, s = 1) }
  v_mod_pr <- aa
  v_mod = exp(mu_pr[5] + sigma[5] * v_mod_pr)
  
  RTbound <- 0.3
  minRT <- 0.35
  tau <- rnorm(s,0,1)
  while(any(tau > RTbound | tau < 0)) {tau <- abs(rnorm(s,0,0.1))}
}

## manually set parameters

v_mod <- rep(3.566, n_subj)
eta[,1] <- rep(0.03, n_subj)
eta[,2] <- rep(0.1, n_subj)


### prediction error (works)
# generate EV with given parameters: V1: updates via weighted prediction error
ev <- matrix(data=0, nrow=dat$T,ncol=2)
pe <- list()
pe <- rep(0,dat$T)
delta <- list()
delta <- rep(0,dat$T)
for(s in 1:dat$N){
  # initialize lower(=1) and upper(=2) bound value
  ev[first[s],1] <- 0.5
  ev[first[s],2] <- 0.5
  for(trial in (first[s]:last[s]-1)){
    delta[trial] <- (ev[trial,2]-ev[trial,1]) * v_mod[s]
    # response=1: lower bound, response=2: upper bound
    pe[trial] <- (dat$value[trial]-ev[trial,dat$response[trial]])
    ev[trial+1,dat$response[trial]] <- ev[trial,dat$response[trial]] + eta[s,dat$response[trial]]*(dat$value[trial]-ev[trial,dat$response[trial]])
    ev[trial+1,dat$nonresponse[trial]] <- ev[trial,dat$nonresponse[trial]]
  }
  delta[last[s]] =  (ev[last[s]-1,2]-ev[last[s]-1,1]) * v_mod[s]
}
learning <- cbind(ev,pe,dat$response,delta)

ggplot() + 
  geom_line(aes(x=1:38, y=learning[,1]),colour="red",linetype="dashed") +
  geom_line(aes(x=1:38, y=learning[,2]), colour="blue",linetype="twodash")

ggplot() + 
  geom_line(aes(x=1:38, y=delta),colour="red",linetype="dashed")

### prediction errors with ilogit(eta) (arrives at upper/lower bound much faster)
ev <- matrix(data=0, nrow=dat$T,ncol=2)
pe <- list()
pe <- rep(0,dat$T)
delta <- list()
delta <- rep(0,dat$T)
for(s in 1:dat$N){
  # initialize lower and upper bound value
  ev[first[s],1] <- 0.5
  ev[first[s],2] <- 0.5
  for(trial in (first[s]:last[s]-1)){
    delta[trial] <- (ev[trial,2]-ev[trial,1]) * v_mod[s]
    # response=1: lower bound, response=2: upper bound
    pe[trial] <- (dat$value[trial]-ev[trial,dat$response[trial]])
    ev[trial+1,dat$response[trial]] <- ev[trial,dat$response[trial]] + inv.logit(eta[s,dat$response[trial]])*(dat$value[trial]-ev[trial,dat$response[trial]])
    ev[trial+1,dat$nonresponse[trial]] <- ev[trial,dat$nonresponse[trial]]
  }
  delta[last[s]] =  (ev[last[s]-1,2]-ev[last[s]-1,1]) * v_mod[s]
}
learning <- cbind(ev,pe,dat$response,delta)

ggplot() + 
  # red is the pe for the lower bound -> goes to 0 because person has learned and doesn't 
  geom_line(aes(x=1:38, y=learning[,1]),colour="red",linetype="dashed") +
  geom_line(aes(x=1:38, y=learning[,2]), colour="blue",linetype="twodash")

ggplot() + 
  geom_line(aes(x=1:38, y=delta),colour="red",linetype="dashed")


### prediction errors where we learn from unsuccessful trials
ev <- matrix(data=0, nrow=dat$T,ncol=2)
pe <- list()
pe <- rep(0,dat$T)
delta <- list()
delta <- rep(0,dat$T)
for(s in 1:dat$N){
  # initialize lower and upper bound value
  ev[first[s],1] <- 0.5
  ev[first[s],2] <- 0.5
  for(trial in (first[s]:last[s]-1)){
    delta[trial] <- (ev[trial,2]-ev[trial,1]) * v_mod[s]
    # response=1: lower bound, response=2: upper bound
    pe[trial] <- (dat$value[trial]-ev[trial,dat$response[trial]])
    ev[trial+1,dat$response[trial]] <- ev[trial,dat$response[trial]] + eta[s,dat$response[trial]]*(dat$value[trial]-ev[trial,dat$response[trial]])
    ev[trial+1,dat$nonresponse[trial]] <- ev[trial,dat$nonresponse[trial]]
  }
  delta[last[s]] =  (ev[last[s]-1,2]-ev[last[s]-1,1]) * v_mod[s]
}
learning <- cbind(ev,pe,dat$response,delta)

ggplot() + 
  geom_line(aes(x=1:38, y=learning[,1]),colour="red",linetype="dashed") +
  geom_line(aes(x=1:38, y=learning[,2]), colour="blue",linetype="twodash")

ggplot() + 
  geom_line(aes(x=1:38, y=delta),colour="red",linetype="dashed")


ev <- matrix(data=0, nrow=dat$T,ncol=2)
pe <- list()
pe <- rep(0,dat$T)
delta <- list()
delta <- rep(0,dat$T)
for(s in 1:dat$N){
  # initialize lower and upper bound value
  ev[first[s],1] <- 0.5
  ev[first[s],2] <- 0.5
  for(trial in (first[s]:last[s]-1)){
    delta[trial] <- (ev[trial,1]-ev[trial,2]) * v_mod[s]
    # response=1: lower bound, response=2: upper bound
    pe[trial] <- (dat$value[trial]-ev[trial,dat$response[trial]])
    ev[trial+1,dat$response[trial]] <- ev[trial,dat$response[trial]] + eta[s,dat$response[trial]]*(dat$value[trial]-ev[trial,dat$response[trial]])
    ev[trial+1,dat$nonresponse[trial]] <- ev[trial,dat$nonresponse[trial]]
  }
  delta[last[s]] =  (ev[last[s]-1,1]-ev[last[s]-1,2]) * v_mod[s]
}
learning <- cbind(ev,pe,dat$response,delta)

ggplot() + 
  geom_line(aes(x=1:38, y=learning[,1]),colour="red",linetype="dashed") +
  geom_line(aes(x=1:38, y=learning[,2]), colour="blue",linetype="twodash") 

## 

ev <- matrix(data = 0, nrow=dat$T, ncol=2)
delta <- rep(0, dat$T)
for(s in 1:dat$N){
  # value for correct trial
  ev[first[s],1] <- 0
  # value for incorrect trial
  ev[first[s],2] <- 0
  for(trial in (first[s]+1:last[s]-1)){
    delta[trial] = (ev[trial,1]-ev[trial,2]) #* v_mod[s]
    ev[trial+1,1] = ev[trial,1] + abs(eta[s,1] * (dat$value[trial]-ev[trial,1]))
    ev[trial+1,2] = ev[trial,2] + abs(eta[s,2] * (dat$value[trial]-ev[trial,2]))
  }
}

# annahme: wenn man sich für eine option entscheidet, entscheidet man sich gleichzeitig gegen die andere. 
# 