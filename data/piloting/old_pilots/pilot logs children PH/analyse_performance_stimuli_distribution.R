
library(readr)
library(rstudioapi)
library(readr, lib="\\\\idnetapp-homes3.uzh.ch\\phalle$\\Documents\\R\\win-library\\3.5")
library(data.table)
library(ggplot2)
library(tibble)

#set inputs
dirinput <- dirname(rstudioapi::getActiveDocumentContext()$path)
# make sure output directory exists already
diroutput <- dirinput

task = "FeedLearn"

setwd(dirinput)
files <- dir(pattern=".txt")

#read data per file and combine in array
datalist <- list()

for (i in 1:length(files)){
  no_col <- max(count.fields(files[i], sep = "\t"))
  D <- read_delim(
    files[i],"\t", escape_double = FALSE, locale = locale(), trim_ws = TRUE, skip_empty_rows=TRUE)
  D <- D[2:dim(D)[1],] # remove first row (with LR or RL )
  miss <- length(D$resp==0)
  #D$trial <- as.numeric(D$trial)
  #D <- D[which(D$block==1|D$block==2|D$block==3|D$block==4),]  #exclude practice trials and unnecessary rows (e.g., with avg_resp). It should have now 200 x 4 = 800 rows
  D <- cbind(rep(as.integer(substr(files[i],1,2)),dim(D)[1]),D)
  #D<-D[D$resp!=0,] # remove 'too slow ' responses 
  colnames(D)[1] <- "subjID"
  colnames(D)[grep("rt",colnames(D))] <- "RT" 
  D[grep("RT",colnames(D))] <- D[grep("RT",colnames(D))]/1000 # RTs in seconds
  D <- as_tibble(cbind(D,paste(D$vFile,D$aFile)))
  colnames(D)[ncol(D)] <- "pair"
  D["hit"] = "NA"
  for (j in 1:max(as.numeric(D$trial))){
    if (D[j,]$fb == 1){
      D[j,]$hit = j
    }
    else D[j,]$hit = "NA"
  }
  datalist[[i]] <- D
}

Gather <- as_tibble(data.table::rbindlist(datalist, fill=TRUE)) # combine all data frames in one
#Save as CSV
setwd(diroutput)
write.table(Gather,file = paste("performance_all_",task,".txt",sep=""),sep="\t",row.names = FALSE,quote=FALSE)

x1 = factor(Gather$hit, levels=1:36)
Gather <- cbind(Gather, x1)

subj <- unique(Gather$subjID)
df_subj <- list()
for (i in 1:length(subj)){
  df_subj[[i]] <- subset(Gather, subjID == subj[i])
  blocks <- unique(df_subj[[i]]$block)
  for (j in 0:(length(blocks)-1)){
    data <- subset(Gather, subjID==i & block==j & hit!="NA")
    p<-ggplot(data, aes(x=x1, color=pair)) + 
      geom_histogram(breaks=seq(1, 36, by=1), fill="white", stat="count")+
      scale_x_discrete(drop = FALSE)+
      geom_density(alpha=.2, fill="#FF6666") 
    ggsave(p, file=paste("Hits Distribution","Subj_",i,"_B",j,".png", sep=""),width = 6, height = 6, scale=1)
  }
}

data <- subset(Gather, subjID=1 & block=1 & hit!="NA")
p<-ggplot(Gather, aes(x=hit, color=pair)) + 
  geom_histogram(color="black", fill="white")+
  geom_density(alpha=.2, fill="#FF6666") 
p

rhat(ddm_model2)
printFit(ddm_model2)
