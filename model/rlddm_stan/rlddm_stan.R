library(rstan)

path <- dirname(rstudioapi::getActiveDocumentContext()$path)
model_path <- paste0(path,"/rlddm.stan")
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
raw_data$incorrect <- as.integer(ifelse(raw_data$correct==1, 2, 1))
# define the values for the rewards
value <- ifelse(raw_data$correct==1, 0, 1)
## all RT with negative choices -> -1
new_RT <- ifelse(raw_data$correct==1, raw_data$RT*-1, raw_data$RT)
## # obs
n_trials <- nrow(raw_data)


dat <- list("N" = n_subj, "T"=n_trials,"RTbound" = 0.05,"minRT" = minRT, "iter" = raw_data$trial, "correct" = raw_data$correct, "incorrect" = raw_data$incorrect,
            "RT" = new_RT, "first" = first, "last" = last, "value"=value)  # names list of numbers

stanmodel_arg <- rstan::stan_model(model_path)

fit <- rstan::sampling(object  = stanmodel_arg,
                       data    = dat,
                       init    = "random",
                       chains  = 2,
                       iter    = 4000,
                       warmup  = 1000,
                       thin    = 1,
                       control = list(adapt_delta   = 0.95,
                                      stepsize      = 1,
                                      max_treedepth = 10))

parVals <- rstan::extract(fit, permuted = TRUE)