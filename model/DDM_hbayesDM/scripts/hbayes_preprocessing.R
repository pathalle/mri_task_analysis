## 

library(rstan)
library(data.table)

preprocess_func = function(raw_data, general_info, RTbound = 0.1) {
  # Use raw_data as a data.frame
  raw_data <- as.data.frame(raw_data)
  
  # Use general_info of raw_data
  subjs   <- general_info$subjs
  n_subj  <- general_info$n_subj
  
  # Number of upper and lower boundary responses
  Nu <- with(raw_data, aggregate(choice == 2, by = list(y = subjid), FUN = sum)[["x"]])
  Nl <- with(raw_data, aggregate(choice == 1, by = list(y = subjid), FUN = sum)[["x"]])
  
  # Reaction times for upper and lower boundary responses
  RTu <- array(-1, c(n_subj, max(Nu)))
  RTl <- array(-1, c(n_subj, max(Nl)))
  for (i in 1:n_subj) {
    subj <- subjs[i]
    subj_data <- subset(raw_data, raw_data$subjid == subj)
    
    RTu[i, 1:Nu[i]] <- subj_data$rt[subj_data$choice == 2]  # (Nu/Nl[i]+1):Nu/Nl_max will be padded with 0's
    RTl[i, 1:Nl[i]] <- subj_data$rt[subj_data$choice == 1]  # 0 padding is skipped in likelihood calculation
  }
  
  # Minimum reaction time
  minRT <- with(raw_data, aggregate(rt, by = list(y = subjid), FUN = min)[["x"]])
  
  # Wrap into a list for Stan
  data_list <- list(
    N       = n_subj,   # Number of subjects
    Nu_max  = max(Nu),  # Max (across subjects) number of upper boundary responses
    Nl_max  = max(Nl),  # Max (across subjects) number of lower boundary responses
    Nu      = Nu,       # Number of upper boundary responses for each subject
    Nl      = Nl,       # Number of lower boundary responses for each subject
    RTu     = RTu,      # Upper boundary response times
    RTl     = RTl,      # Lower boundary response times
    minRT   = minRT,    # Minimum RT for each subject
    RTbound = RTbound   # Lower bound of RT across all subjects (e.g., 0.1 second)
  )
  
  # Returned data_list will directly be passed to Stan
  return(data_list)
}

path <- "/home/padraigh/Dokumente/Uni/NSC/Thesis/mri_task_analysis/model/rlddm_hbayes/"
setwd(path)
model_path <- paste0(path,"choiceRT_rlddm.stan")
data_path <- paste0(path,"test_input.txt")
stanmodel_arg <- rstan::stan_model(model_path)


data_columns    = c("subjID", "choice", "RT")
raw_data <- data.table::fread(file = data_path, header = TRUE, sep = "\t", data.table = TRUE,
                              fill = TRUE, stringsAsFactors = TRUE, logical01 = FALSE)



colnames_raw_data <- colnames(raw_data)
insensitive_data_columns <- tolower(gsub("_", "", data_columns, fixed = TRUE))
colnames(raw_data) <- tolower(gsub("_", "", colnames(raw_data), fixed = TRUE))
complete_rows       <- complete.cases(raw_data[, insensitive_data_columns, with = FALSE])
sum_incomplete_rows <- sum(!complete_rows)

####################################################
##   Prepare general info about the raw data   #####
####################################################

subjs    <- NULL   # List of unique subjects (1D)
n_subj   <- NULL   # Total number of subjects (0D)

b_subjs  <- NULL   # Number of blocks per each subject (1D)
b_max    <- NULL   # Maximum number of blocks across all subjects (0D)

t_subjs  <- NULL   # Number of trials (per block) per subject (2D or 1D)
t_max    <- NULL   # Maximum number of trials across all blocks & subjects (0D)

# To avoid NOTEs by R CMD check
.N <- NULL
subjid <- NULL

DT_trials <- raw_data[, .N, by = c("subjid", "block")]
DT_blocks <- DT_trials[, .N, by = "subjid"]
subjs     <- DT_blocks$subjid
n_subj    <- length(subjs)
b_subjs   <- DT_blocks$N
b_max     <- max(b_subjs)
t_subjs   <- array(0, c(n_subj, b_max))
for (i in 1:n_subj) {
  subj <- subjs[i]
  b <- b_subjs[i]
  t_subjs[i, 1:b] <- DT_trials[subjid == subj]$N
}
t_max     <- max(t_subjs)


general_info <- list(subjs, n_subj, b_subjs, b_max, t_subjs, t_max)
names(general_info) <- c("subjs", "n_subj", "b_subjs", "b_max", "t_subjs", "t_max")

preprocessed <- preprocess_func(raw_data,general_info,RTbound = 0.1)
