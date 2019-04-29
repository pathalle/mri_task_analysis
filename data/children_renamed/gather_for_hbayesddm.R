library(readr)
library(tibble)
library(tidyr)

# ------------------------------------------------------------------------------------------------
# Gather all data 
# ------------------------------------------------------------------------------------------------
#dataPath = "N:/studies/AllRead/Feedback learning/ABMP performance/Adults/logs/"
dataPath = "N:/studies/AllRead/Feedback learning/ABMP performance/Children/logs_excluded/"
setwd(dataPath)
files <- dir(pattern="*.txt")
datalist <- list()

for (f in 1:length(files)){
  #Read file
  no_col <- max(count.fields(files[f], sep = "\t"))
  D <- read_delim(files[f],"\t", escape_double = FALSE, locale = locale(), trim_ws = TRUE)
  
  #exclude practice trials and unnecessary rows (e.g., with avg_resp).
  D <- D[2:dim(D)[1],]
  D <- D[which(D$block==1|D$block==2|D$block==3|D$block==4),] 
  
  #minor adjustments
  colnames(D)[grep("rt",colnames(D))] <- "RT" 
  D[names(D)[names(D) == "RT"]] <-  D[names(D)[names(D) == "RT"]]/1000 # RTs to seconds
  
  
  # if frame does contain the colname resp -> change
  # don't use grep here, because it will find all instances of 'resp' (even colnames with name respOnset) 
  names(D)[names(D) == ".result"] <- "result"
  names(D)[names(D) == "V-file"] <- "vFile"
  names(D)[names(D) == "V-File"] <- "vFile"
  names(D)[names(D) == "A-file"] <- "aFile"
  names(D)[names(D) == "A-File"] <- "aFile"
  names(D)[names(D) == "Vstim"] <- "vStim"
  names(D)[names(D) == "Astim"] <- "aStim"
  names(D)[names(D) == "fb"] <- "choice"
  
  # choice has to be 1/2 -> add +1 to the choice vector (0/1->1/2)
  D$choice = D$choice + 1
  #add subjects
  D <- as_tibble(cbind(rep(substr(files[f],1,7),dim(D)[1]),D)) # use "as_tibble" to avoid issues with the unicode symbols display
  colnames(D)[1] <- "subjID"
  datalist[[f]] <- D
}

Gather <- as_tibble(data.table::rbindlist(datalist, fill=TRUE)) # combine all data frames in one
gather_to_txt <- Gather[,c(grep("subjID",colnames(Gather)),grep("choice",colnames(Gather)),grep("RT",colnames(Gather)))]
gather_to_txt <- gather_to_txt[which(gather_to_txt$RT>=0.2),]
#write.table(gather_to_txt,file = "model_input_abmp_adults.txt",sep="\t",row.names = FALSE,quote=FALSE)
write.table(gather_to_txt,file = "model_input_abmp_children_190402.txt",sep="\t",row.names = FALSE,quote=FALSE)
