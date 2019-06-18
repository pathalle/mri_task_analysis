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


### data loading and preprocessing

path <- dirname(rstudioapi::getActiveDocumentContext()$path)
model_path <- paste0(path,"/rlddm_invlog.stan")
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
# automatically filter missed responses (since RT = 0)
raw_data <- raw_data[which(raw_data$RT > 0.15),]
raw_data$subjID = rep('01',nrow(raw_data))

# fb = 0 incorrect, fb = 1 correct, (fb = 2 missed)
names(raw_data)[names(raw_data)=="fb"] <- "correct"

## prepare data for jags
#raw_data$row <- seq.int(nrow(raw_data))
DT_trials <- raw_data[, .N, by = "subjID"]
subjs     <- DT_trials$subjID
n_subj    <- length(subjs)
# get minRT
minRT <- with(raw_data, aggregate(RT, by = list(y = subjID), FUN = min)[["x"]])
ifelse(is.null(dim(minRT)),minRT<-as.array(minRT))
# assign new trial number for excluded decisions

for (subj in subjs){
  sub <- which(raw_data$subjID==subj)
  raw_data[sub,]$trial <- seq.int(nrow(raw_data[sub,]))
}
# first is Sx1 matrix identifying all first trials of a subject for each choice
first <- which(raw_data$trial==1)
# if N=1 transform int to 1-d array
ifelse(is.null(dim(first)),first<-as.array(first))
# last is a Sx1 matrix identifying all last trials of a subject for each choice
last <- as.integer(first + DT_trials$N - 1)
ifelse(is.null(dim(last)),last<-as.array(last))
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

# encoding for simulation: lower (incorrect) response=1, upper (correct) response =2 
raw_data$response = raw_data$correct+1
raw_data$nonresponse = abs(raw_data$correct-2)

dat <- list("N" = n_subj, "T"=n_trials,"RTbound" = 0.15,"minRT" = minRT, "iter" = raw_data$trial, "response" = raw_data$response,"nonresponse" = raw_data$nonresponse,
            "RT" = raw_data$RT, "first" = first, "last" = last, "value"=value, "stims" = stims)  # names list of numbers


stanmodel_invlog <- rstan::stan_model(model_path)

fit_invlog <- rstan::sampling(object  = stanmodel_invlog,
                       data    = dat,
                       init    = "random",
                       chains  = 2,
                       iter    = 10000,
                       warmup  = 2000,
                       thin    = 1,
                       control = list(adapt_delta   = 0.95,
                                      stepsize      = 1,
                                      max_treedepth = 10),
                       verbose =TRUE)

parValsinvl <- rstan::extract(fit_invlog, permuted = TRUE)

rstan::stan_diag(fit_invlog, info = 'sample') # shows three plots together
rstan::stan_par(fit_invlog, par = "alpha")


## now access the data of the fit 


fit_summary_invlog <- rstan::summary(fit_invlog)

print(fit_summary_invlog$summary)

#########



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
