library(rjags)
library(coda)
load.module("wiener")
## for easy summary statistics of the chains
library(mmcc)

path <- dirname(rstudioapi::getActiveDocumentContext()$path)
#path <- "/home/padraigh/Dokumente/Uni/NSC/Thesis/mri_task_analysis/model/rlddm_jags/"
model_path <- paste0(path,"/rlddm_reduced.jag")
data_path <- paste0(path,"/input_6subj_wtrials_1block.txt")
setwd(path)
raw_data <- data.table::fread(file = data_path, header = TRUE, sep = "\t", data.table = TRUE,
                              fill = TRUE, stringsAsFactors = TRUE, logical01 = FALSE)

raw_data <- raw_data[which(raw_data$RT > 0.3),]
names(raw_data)[names(raw_data)=="choice"] <- "correct"

## prepare data for jags
#raw_data$row <- seq.int(nrow(raw_data))
DT_trials <- raw_data[, .N, by = "subjID"]
subjs     <- DT_trials$subjID
n_subj    <- length(subjs)
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
raw_data$incorrect <- as.integer(ifelse(raw_data$correct==1, 2, 1))
# define the values for the rewards
value <- ifelse(raw_data$correct==1, 0, 1)
## all RT with negative choices -> -1
new_RT <- ifelse(raw_data$correct==1, raw_data$RT*-1, raw_data$RT)


#############  for testing purposes ############
for(s in 1:length(subjs)){
  for (trial in (first[s]):(last[s]-1)) {  
  print(raw_data$RT[trial])
  }
}
# Check if this atually does not print the last trials:
# There should be an error for RT[781] because there are only 780 trials]
raw_data$RT[781] # -> error: correct!
# for S1, the last trial is 0.931. Thus RT[last[1]-1] should be this number
raw_data$RT[last[1]-1] # -> correct!

###############################################

dat <- list("S" = n_subj, "iter" = raw_data$trial, "correct" = raw_data$correct, "incorrect" = raw_data$incorrect,
            "RT" = new_RT, "first" = first, "last" = last, "value"=value)  # names list of numbers

##### Initial values
inits <- list( etag_mu=0.2, ag_mu=1.7,ig_mu=0, mg_mu=3.5, tg_mu=0.3)


jags.m <- jags.model( file = model_path, inits=inits, data=dat, n.chains=4, n.adapt=1000)


### take a look at the posterior distributions
params <- c("eta","i","m","t")
samples <- coda.samples(jags.m, params, n.iter = 4000)
plot(samples_eta)

tidy(samples, chain = TRUE)

# these are the values for the positive and negative learning rates
# eta[s,1] are the negative feedbacks 
plot(samps[,c(1,32)])
plot(samps[,c(2,33)])
plot(samps[,c(3,34)])
# etc. -> compare the positive and negative learning rate
save(samps, file="sampling_rlddm_eta_190423.R")
samps_as_matrix <- as.matrix(samps, iters = FALSE, chains = FALSE)
