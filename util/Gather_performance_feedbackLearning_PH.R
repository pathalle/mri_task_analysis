# Gather task performance from txt files (Presentation output)
#--------------------------------------------------------------
library(readr)

#set inputs
dirinput <- "N:/studies/AllRead/Feedback learning/ABMP logs/Children/logs/"
diroutput <- "N:/studies/AllRead/Feedback learning/ABMP logs/Children/"
task = "LSB"

#read files
setwd(dirinput)
files <- dir(pattern=paste(task,".txt",sep=""))

#read data per file and combine in array
datalist <- list()
for (i in 1:length(files)){
  
  D <- read_delim(files[i],"\t", escape_double = FALSE, locale = locale(), trim_ws = TRUE) 
  D <- D[2:dim(D)[1],] # remove first row (with LR or RL )
  D <- D[which(D$block==1|D$block==2|D$block==3|D$block==4),]  #exclude practice trials and unnecessary rows (e.g., with avg_resp). It should have now 200 x 4 = 800 rows
  D <- cbind(as.integer(rep(substr(files[i],2,4)),dim(D)[1]),D) #add subject ID in a first column
  #add group info
  if (D[1,1] < 100) { 
    D <- cbind(rep("Typ",dim(D)[1]),D) #add subject ID in a first column
  }else if (D[1,1] > 100){ D <- cbind(rep("Dys",dim(D)[1]),D) #add subject ID in a first column
  }
  D[grep("RT",colnames(D))] <- D[grep("RT",colnames(D))]/1000 # RTs in seconds
  D<-D[D$resp!=0,] # remove 'too slow ' responses 
  # rename some columns
  colnames(D)[1] <- "group"
  colnames(D)[2] <- "subj_idx"
  colnames(D)[grep("C/I",colnames(D))] <- "Stim" 
  colnames(D)[grep("fb",colnames(D))] <- "response" 
  colnames(D)[grep("RT",colnames(D))] <- "rt" 
  D <- D[complete.cases(D),] # remove NAs
  D <- D[,c(grep("group",colnames(D)),grep("subj_idx",colnames(D)),grep("block",colnames(D)),grep("Stim",colnames(D)),grep("rt",colnames(D)),grep("response",colnames(D)))]
  #Add bin (1,2,3,4 quartiles) 
  D <- cbind(D,rep(0,dim(D)[1]))
  colnames(D)[length(D)]<- "Quartile"
 
  for (ii in 1:max(unique(D$block))){ 
    currblockidx <- which(D$block==ii)
    currblockCs <- currblockidx[which(D$Stim[currblockidx]=='C')]
     for (iii in 1:length(currblockCs)) {
       if (iii <=25){  D$Quartile[currblockCs[iii]]<-1
       } else if (iii > 25 & iii <= 50) { D$Quartile[currblockCs[iii]]<-2
       } else if (iii > 50 & iii <= 75) { D$Quartile[currblockCs[iii]]<-3
       } else if (iii > 75 & iii <= 100) {  D$Quartile[currblockCs[iii]]<-4
       }
     }
  }
   for (ii in 1:max(unique(D$block))){ 
    currblockidx <- which(D$block==ii)
    currblockCs <- currblockidx[which(D$Stim[currblockidx]=='I')] # inconsistent trials
    for (iii in 1:length(currblockCs)) {
      if (iii <=25){  D$Quartile[currblockCs[iii]]<-1
      } else if (iii > 25 & iii <= 50) { D$Quartile[currblockCs[iii]]<-2
      } else if (iii > 50 & iii <= 75) { D$Quartile[currblockCs[iii]]<-3
      } else if (iii > 75 & iii <= 100) {  D$Quartile[currblockCs[iii]]<-4
      }
    }
  }
  
  #Add half (1,2) indices 
  D <- cbind(D,rep(0,dim(D)[1]))
  colnames(D)[length(D)]<- "Half"
  for (ii in 1:max(unique(D$block))){ 
    currblockidx <- which(D$block==ii)
    currblockCs <- currblockidx[which(D$Stim[currblockidx]=='C')]
    for (iii in 1:length(currblockCs)) {
      if (iii <=50){  D$Half[currblockCs[iii]]<-1
      } else if (iii > 50) { D$Half[currblockCs[iii]]<-2
      }
    }
  }
  for (ii in 1:max(unique(D$block))){ 
    currblockidx <- which(D$block==ii)
    currblockCs <- currblockidx[which(D$Stim[currblockidx]=='I')] #inconsistent trials
    for (iii in 1:length(currblockCs)) {
      if (iii <=50){  D$Half[currblockCs[iii]]<-1
      } else if (iii > 50) { D$Half[currblockCs[iii]]<-2
      }
    }
  }
  
  #Concatenate to list
  datalist[[i]] <- D
  
}

Gather <- data.table::rbindlist(datalist) # combine all data frames in one
  
#Save as CSV
setwd(diroutput)
write.csv(Gather,file = paste("Gathered_performance_",task,"_bins.csv",sep=""),row.names = FALSE,quote=FALSE)