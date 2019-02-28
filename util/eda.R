
library(readr)
library(rstudioapi)
library(readr)
library(data.table)
library(ggplot2)

#set inputs
dirinput <- dirname(rstudioapi::getActiveDocumentContext()$path)
# make sure output directory exists already
diroutput <- paste(dirinput,"output",sep="/")

#dirinput <- "O:\studies\allread\models_PH\logs"
#diroutput <- "O:\studies\allread\models_PH\logs\output"

setwd(dirinput)
files <- dir(pattern=".txt")

#read data per file and combine in array
datalist <- list()

for (i in 1:length(files)){
  no_col <- max(count.fields(files[i], sep = "\t"))
  D <- read_delim(
    files[i],"\t", escape_double = FALSE, locale = locale(), trim_ws = TRUE)
  D <- D[2:dim(D)[1],] # remove first row (with LR or RL )
  D <- D[which(D$block==1|D$block==2|D$block==3|D$block==4),]  #exclude practice trials and unnecessary rows (e.g., with avg_resp). It should have now 200 x 4 = 800 rows
  D <- cbind(rep(as.integer(substr(files[i],14,15)),dim(D)[1]),D)
  D<-D[D$resp!=0,] # remove 'too slow ' responses 
  colnames(D)[1] <- "subjID"
  colnames(D)[grep("rt",colnames(D))] <- "RT" 
  D[grep("RT",colnames(D))] <- D[grep("RT",colnames(D))]/1000 # RTs in seconds
  # if frame does contain the colname resp -> change
  # don't use grep here, because it will find all instances of 'resp' (even colnames with name respOnset) 
  names(D)[names(D) == "V-file"] <- "vFile"
  names(D)[names(D) == "V-File"] <- "vFile"
  names(D)[names(D) == "A-file"] <- "aFile"
  names(D)[names(D) == "A-File"] <- "aFile"
  names(D)[names(D) == "Vstim"] <- "vStim"
  names(D)[names(D) == "Astim"] <- "aStim"
  D <- D[complete.cases(D), ]
  #D <- D[,c(grep("subjID",colnames(D)),grep("choice",colnames(D)),grep("RT",colnames(D)))]
  datalist[[i]] <- D
}

Gather <- data.table::rbindlist(datalist, fill=TRUE) # combine all data frames in one

#Save as CSV
#setwd(diroutput)
write.table(Gather,file = "performance_all_aduls.txt",sep="\t",row.names = FALSE,quote=FALSE)
rm(new_D1)
rm(new_D2)
rm(new_D3)
new_D <- subset(Gather, subjID == 1 & block==1)
new_D$integer = as.integer(new_D$trial)
symbols <- unique(new_D$vFile)
new_D1 <- subset(Gather, subjID == 1 & block==1 & vFile==symbols[1])
new_col <- cumsum(new_D1$fb)
new_col_trial <- 1:nrow(new_D1)
new_cols <- cbind(new_col, new_col_trial)
new_D1 <- cbind(new_cols,new_D1)
colnames(new_D1)[c(1,2)] <- c("cumsum_fb", "trial_separate")
new_D2 <- subset(Gather, subjID == 1 & block==1 & vFile==symbols[2])
new_col <- cumsum(new_D2$fb)
new_col_trial <- 1:nrow(new_D2)
new_cols <- cbind(new_col, new_col_trial)
new_D2 <- cbind(new_cols,new_D2)
colnames(new_D2)[c(1,2)] <- c("cumsum_fb", "trial_separate")
new_D3 <- subset(Gather, subjID == 1 & block==1 & vFile==symbols[3])
new_col <- cumsum(new_D3$fb)
new_col_trial <- 1:nrow(new_D3)
new_cols <- cbind(new_col, new_col_trial)
new_D3 <- cbind(new_cols,new_D3)
colnames(new_D3)[c(1,2)] <- c("cumsum_fb", "trial_separate")
new_D <- rbind(new_D1,new_D2,new_D3)
ggplot(data=new_D, aes(x=trial_separate, y=cumsum_fb, group=vFile, color=vFile)) +
  geom_line()+
  geom_point()

