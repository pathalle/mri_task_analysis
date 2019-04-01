
library(readr, lib="\\\\idnetapp-homes3.uzh.ch\\phalle$\\Documents\\R\\win-library\\3.5")
library(rstudioapi, lib="\\\\idnetapp-homes3.uzh.ch\\phalle$\\Documents\\R\\win-library\\3.5")
library(readr, lib="\\\\idnetapp-homes3.uzh.ch\\phalle$\\Documents\\R\\win-library\\3.5")
library(data.table, lib="\\\\idnetapp-homes3.uzh.ch\\phalle$\\Documents\\R\\win-library\\3.5")
library(dplyr, lib="\\\\idnetapp-homes3.uzh.ch\\phalle$\\Documents\\R\\win-library\\3.5")
library(ggplot2, lib="\\\\idnetapp-homes3.uzh.ch\\phalle$\\Documents\\R\\win-library\\3.5")
library(gridExtra)

#set inputs
dirinput <- dirname(rstudioapi::getActiveDocumentContext()$path)
# make sure output directory exists already
diroutput <- dirinput

task = "FeedLearn"

setwd(dirinput)
files <- dir(pattern=".txt")

#For each block: accuracy (% correct) every six trials (irrespective of pair). And mean accuracy across blocks every six trials. 

gather_data <- function(files){
  datalist <- list()
  for (i in 1:length(files)){
    no_col <- max(count.fields(files[i], sep = "\t"))
    D <- read_delim(
      files[i],"\t", escape_double = FALSE, locale = locale(), trim_ws = TRUE, skip_empty_rows=TRUE)
    D <- D[2:dim(D)[1],] # remove first row (with LR or RL )
    D <- D[which(D$block==1|D$block==2|D$block==3|D$block==4),]  #exclude practice trials and unnecessary rows (e.g., with avg_resp). It should have now 200 x 4 = 800 rows
    D <- cbind(rep(as.integer(substr(files[i],1,2)),dim(D)[1]),D)
    D<-D[D$resp!=0,] # remove 'too slow ' responses 
    ### Rename and transform some columns
    colnames(D)[1] <- "subjID"
    D[,1] <- as.factor(D[,1])
    D$trial = as.integer(D$trial)
    colnames(D)[grep("rt",colnames(D))] <- "RT" 
    D[grep("RT",colnames(D))] <- D[grep("RT",colnames(D))]/1000 # RTs in seconds
    # if frame does contain the colname resp -> change
    # don't use grep here, because it will find all instances of 'resp' (even colnames with name respOnset) 
    names(D)[names(D) == "resp"] <- "choice"
    D <- D[complete.cases(D), ]
    #D <- D[,c(grep("subjID",colnames(D)),grep("choice",colnames(D)),grep("RT",colnames(D)))]
    D <- as_tibble(cbind(D,paste(D$vFile,D$aFile)))
    colnames(D)[ncol(D)] <- "pair"
    datalist[[i]] <- D
  }
  transformed <- data.table::rbindlist(datalist) # combine all data frames in on
  return(transformed)
}

compute_cumulative_sums <- function(data){
  df_subj <- list()
  data<-data[data$resp!=0,] # remove 'too slow ' responses 
  new_data <- data[FALSE,]
  new_cols <- data.frame(cumsum_fb = integer(0), trial_separate = integer(0))
  new_data <- cbind(new_cols,new_data)
  for(i in unique(data$subjID)){
    df_subj[[i]] <- subset(data, subjID == i)
    for(j in unique(df_subj[[i]]$block)){
      df_subj_block <- subset(df_subj[[i]],block==j)
      # compute cumulative sum for each sound-symbol pair in a given block of a subject
      for(k in unique(df_subj_block$pair)){
        df_subj_block_pair <- list()
        df_subj_block_pair[[k]] <- subset(df_subj_block, pair==k)
        new_col <- cumsum(df_subj_block_pair[[k]]$fb)
        new_col_trial <- 1:nrow(df_subj_block_pair[[k]])
        new_cols <- cbind(new_col,new_col_trial)
        df_subj_block_pair[[k]]<- as_tibble(cbind(new_cols,df_subj_block_pair[[k]]))
        colnames(df_subj_block_pair[[k]])[c(1,2)] <- c("cumsum_fb", "trial_separate")
        new_data <- as_tibble(rbind(new_data,df_subj_block_pair[[k]]))
      }
    }
  }
  return(as_tibble(new_data))
}

get_summary_stats <- function(data){
  ### #How many missing responses per block (fb == 2)
  miss_per_block <- data %>%
    select(subjID,RT, fb,block) %>%
    filter(fb==2)  %>%
    group_by(subjID,block) %>%
    tally()
  
  ### #RTs per block
  RT_per_block <- data %>%
    select(subjID,RT, fb,block) %>%
    group_by(subjID, fb,block) %>%
    summarise(mean_rt = mean(RT))
  
  ### #RTs across blocks
  RT_across_blocks <- data %>%
    select(subjID,RT, fb) %>%
    group_by(subjID, fb) %>%
    summarise(mean_rt = mean(RT))
  
  return(list(
    "miss_per_block"=miss_per_block,"rt_per_block"=RT_per_block,"rt_across_blocks"=RT_across_blocks))
}

data <- gather_data(files)
summary_stats <- get_summary_stats(data)
data_with_cumsum <- compute_cumulative_sums(data)

## misses per block
p <- ggplot(summary_stats$miss_per_block, aes(subjID,n,fill=block)) + 
  geom_bar(stat="identity", position="dodge")

df <- split(data_with_cumsum,f = data_with_cumsum$subjID)

attach(df)
for(i in df){
  cumSumPlot <- ggplot(data=i, aes(x=trial_separate, y=cumsum_fb, group=pair, color=pair)) +
    geom_line()+
    geom_point(aes(fill=pair),colour="black",alpha=.5, shape=21, size=3,position=position_dodge(0.2))+
    facet_grid(cols = vars(block)) +
    scale_x_continuous(breaks = unique(trial_separate),limits=c(1,6.5))  +
    scale_y_continuous(breaks = c(0,1,2,3,4,5,6),limits=c(0,6))   +
    guides(alpha=FALSE)+
    theme(axis.title = element_text(size=12),
          title = element_text(size=14),
          plot.subtitle = element_text(size=14,color="darkblue"),
          legend.text=element_text(size=12),
          panel.grid.major = element_line(colour="white"),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "gray88")) +
    labs(title=paste(subjID," Cumulative sum of hits per pair",sep=""))
  ggsave(cumSumPlot, file=paste("Accuracy","Subj_",i$subjID[1],".png", sep=""),width = 6, height = 6, scale=1)
}



#Save as CSV
#setwd(diroutput)
#write.table(Gather,file = paste("performance_all_",task,".txt",sep=""),sep="\t",row.names = FALSE,quote=FALSE)