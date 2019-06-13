
library(readr)
library(rstudioapi)
library(readr)
library(data.table)
library(ggplot2)
library(dplyr)

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
  D <- cbind(rep((substr(files[i],1,3)),dim(D)[1]),D)
}
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
  #D <- D[complete.cases(D), ]
  #D <- D[,c(grep("subjID",colnames(D)),grep("choice",colnames(D)),grep("RT",colnames(D)))]
  datalist[[i]] <- D
}

Gather <- data.table::rbindlist(datalist, fill=TRUE) # combine all data frames in one

#Save as CSV
#setwd(diroutput)
#write.table(Gather,file = "performance_all_aduls.txt",sep="\t",row.names = FALSE,quote=FALSE)
#rm(new_D1)
#rm(new_D2)
#rm(new_D3)
subj <- unique(Gather$subjID)
df_subj <- list()
plot_list <- list()
Gather$fb = as.integer(Gather$fb)
for (i in 1:length(subj)){
  df_subj[[i]] <- subset(Gather, subjID == i)
  blocks <- unique(df_subj[[i]]$block)
  for (j in 1:length(blocks)){
   df_subj_block <- subset(df_subj[[i]],block==j)
   symbols <- unique(df_subj_block$vFile)
   df_subj_block_symbol_list <- list()
   for (k in 1:length(symbols)){
     if (length(unique(df_subj_block$vFile))>=1){
       df_subj_block_symbol_list[[k]] <- subset(df_subj_block, vFile==symbols[k])
       new_col <- cumsum(df_subj_block_symbol_list[[k]]$fb)
       new_col_trial <- 1:nrow(df_subj_block_symbol_list[[k]])
       new_cols <- cbind(new_col,new_col_trial)
       df_subj_block_symbol_list[[k]] <- cbind(new_cols,df_subj_block_symbol_list[[k]])
       colnames(df_subj_block_symbol_list[[k]])[c(1,2)] <- c("cumsum_fb", "trial_separate")
       df_subj_block_symbol <- bind_rows(df_subj_block_symbol_list)
       #print(df_subj_block_symbol)
     }
     plot <- ggplot(data=df_subj_block_symbol, aes(x=trial_separate, y=cumsum_fb, group=vFile, color=vFile)) +
       geom_line()+ geom_point()
     #print(plot)
     ggsave(plot, file=paste("plot",i,j,".png", sep="_"), scale=2)
   }
  }
}


# # for 1 case  
# new_D <- subset(Gather, subjID == 1 & block==1)
# new_D$integer = as.integer(new_D$trial)
# symbols <- unique(new_D$vFile)
# new_D1 <- subset(Gather, subjID == 1 & block==1 & vFile==symbols[1])
# new_col <- cumsum(new_D1$fb)
# new_col_trial <- 1:nrow(new_D1)
# new_cols <- cbind(new_col, new_col_trial)
# new_D1 <- cbind(new_cols,new_D1)
# colnames(new_D1)[c(1,2)] <- c("cumsum_fb", "trial_separate")
# new_D2 <- subset(Gather, subjID == 1 & block==1 & vFile==symbols[2])
# new_col <- cumsum(new_D2$fb)
# new_col_trial <- 1:nrow(new_D2)
# new_cols <- cbind(new_col, new_col_trial)
# new_D2 <- cbind(new_cols,new_D2)
# colnames(new_D2)[c(1,2)] <- c("cumsum_fb", "trial_separate")
# new_D3 <- subset(Gather, subjID == 1 & block==1 & vFile==symbols[3])
# new_col <- cumsum(new_D3$fb)
# new_col_trial <- 1:nrow(new_D3)
# new_cols <- cbind(new_col, new_col_trial)
# new_D3 <- cbind(new_cols,new_D3)
# colnames(new_D3)[c(1,2)] <- c("cumsum_fb", "trial_separate")
# new_D <- rbind(new_D1,new_D2,new_D3)
# ggplot(data=new_D, aes(x=trial_separate, y=cumsum_fb, group=vFile, color=vFile)) +
#   geom_line()+
#   geom_point()

