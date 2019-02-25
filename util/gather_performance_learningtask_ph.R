
library(readr)
library(rstudioapi)
library(readr, lib="\\\\idnetapp-homes3.uzh.ch\\phalle$\\Documents\\R\\win-library\\3.5")
library(data.table)

#set inputs
dirinput <- dirname(rstudioapi::getActiveDocumentContext()$path)
# make sure output directory exists already
diroutput <- paste(dirinput,"output",sep="/")

#dirinput <- "O:\studies\allread\models_PH\logs"
#diroutput <- "O:\studies\allread\models_PH\logs\output"
#task = "fbl"

setwd(dirinput)
files <- dir(pattern=".txt")

#read data per file and combine in array
datalist <- list()

for (i in 1:length(files)){
  no_col <- max(count.fields(files[i], sep = "\t"))
  D <- read_delim(
    files[i],"\t", escape_double = FALSE, locale = locale(), trim_ws = TRUE, skip_empty_rows=TRUE)
  D <- D[2:dim(D)[1],] # remove first row (with LR or RL )
  D <- D[which(D$block==1|D$block==2|D$block==3|D$block==4),]  #exclude practice trials and unnecessary rows (e.g., with avg_resp). It should have now 200 x 4 = 800 rows
  D <- cbind(as.integer(rep(substr(files[i],2,4)),dim(D)[1]),D)
  #D<-D[D$resp!=0,] # remove 'too slow ' responses 
  colnames(D[1]) <- "subjID"
  colnames(D)[grep("rt",colnames(D))] <- "RT" 
  D[grep("RT",colnames(D))] <- D[grep("RT",colnames(D))]/1000 # RTs in seconds
  colnames(D)[grep("resp",colnames(D))] <- "choice" 
  #D <- D[complete.cases(D),] # remove NAs
  D <- D[,c(grep("subjID",colnames(D)),grep("RT",colnames(D))),grep("choice",colnames(D))]
  datalist[[i]] <- D
}

Gather <- data.table::rbindlist(datalist) # combine all data frames in one

#Save as CSV
setwd(diroutput)
write.table(Gather,file = paste("performance_all_",task,".txt",sep=""),sep="\t",row.names = FALSE,quote=FALSE)

D <- read_delim(
  files[1],"\t", escape_double = FALSE, locale = locale(), trim_ws = TRUE, skip_empty_rows=TRUE)
D <- D[2:dim(D)[1],] # remove first row (with LR or RL )
D <- cbind(rep(as.integer(substr(files[1],14,15)), dim(D)[1]),D) # put subj id
rep(as.integer(substr(files[1],14,15)), dim(D)[1])
