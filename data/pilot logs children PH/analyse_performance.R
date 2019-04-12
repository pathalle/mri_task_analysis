
library(readr)
library(rstudioapi)
library(readr)
library(data.table)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(plyr)

#set inputs
dirinput <- dirname(rstudioapi::getActiveDocumentContext()$path)
# make sure output directory exists already
diroutput <- dirinput

task = "FeedLearn"

setwd(dirinput)
files <- dir(pattern=".txt")

gather_data <- function(files){
  datalist <- list()
  for (i in 1:length(files)){
    no_col <- max(count.fields(files[i], sep = "\t"))
    D <- read_delim(
      files[i],"\t", escape_double = FALSE, locale = locale(), trim_ws = TRUE)
    D <- D[2:dim(D)[1],] # remove first row (with LR or RL )
    D <- D[which(D$block==1|D$block==2|D$block==3|D$block==4),]  #exclude practice trials and unnecessary rows (e.g., with avg_resp). It should have now 200 x 4 = 800 rows
    D <- cbind(rep(as.integer(substr(files[i],1,2)),dim(D)[1]),D)
    #D<-D[D$resp!=0,] # remove 'too slow ' responses 
    ### Rename and transform some columns
    colnames(D)[1] <- "subjID"
    D[,1] <- as.factor(D[,1])
    D$trial = as.integer(D$trial)
    D$match = as.factor(D$match)
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
  data<-data[data$choice!=0,] # remove 'too slow ' responses 
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
        #df_subj_block_pair[[k]] <- subset(df_subj_block, pair==k & fb != 2 )
        df_subj_block_pair[[k]] <- subset(df_subj_block, pair==k )
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
    filter(fb!=2)  %>%
    group_by(subjID, fb,block) %>%
    summarise(mean_rt = mean(RT))
  
  ### #RTs across blocks
  RT_across_blocks <- data %>%
    select(subjID,RT, fb) %>%
    filter(fb!=2)  %>%
    group_by(subjID, fb) %>%
    summarise(mean_rt = mean(RT))

  ## hits per sextile
  correct_per_sextile <- data %>%
    select(subjID, fb,block, sextile) %>%
    filter(fb==1) %>%
    group_by(subjID,block,sextile) %>%
    tally()
  
  return(list(
    "miss_per_block"=miss_per_block,"rt_per_block"=RT_per_block,"rt_across_blocks"=RT_across_blocks,"hits_per_sextile"=correct_per_sextile))
}

split_trials <- function(data){
  data$trial <- as.integer(data$trial)
  data$sextile <- 0
  data[which(data$trial <= 6),]$sextile = 1
  data[which(data$trial > 6 & data$trial <= 12),]$sextile = 2
  data[which(data$trial > 12 & data$trial <= 18),]$sextile = 3
  data[which(data$trial > 18 & data$trial <= 24),]$sextile = 4
  data[which(data$trial > 24 & data$trial <= 30),]$sextile = 5
  data[which(data$trial > 25 & data$trial <= 36),]$sextile = 6
  return(data)
  }


########################### prepare data #########################
data <- gather_data(files)
data_sextiles <- split_trials(data)
summary_stats <- get_summary_stats(data_sextiles)
data_with_cumsum <- compute_cumulative_sums(data_sextiles)
# add column with information to which sextile a sequence of trials belongs"


########################### start plotting #########################
## hits per sextile
summary_stats_sextiles <- summary_stats$hits_per_sextile
attach(summary_stats_sextiles)
for(i in unique(summary_stats_sextiles$subjID)){
  for (j in unique(summary_stats_sextiles$block)){
    subset_all = subset(summary_stats_sextiles, subjID==i)
    plots_per_sextiles <- ggplot(data=summary_stats_sextiles, aes(x=sextile, y=n/6, fill=match, colour=match)) +
      geom_line()+
      scale_x_continuous(breaks = unique(sextile),limits=c(0.5,6.5))  +
      scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1),limits=c(0,1))   +
      geom_point(aes(fill=block, colour=block),colour="black",alpha=.5, shape=21, size=3,position=position_dodge(0.2))+
      xlab("Sextiles") +
      ylab("Accuracy") +
      labs(title="Accuracy across pair blocks (6 trials)") +
      facet_grid(rows=vars(block))
     ggsave(plots_per_sextiles, file=paste("accuracy_across_sextiles","_subj_",i,".png", sep=""),width = 6, height = 6, scale=1)
  }
}
    
plots_per_sextiles <- ggplot(data=summary_stats$hits_per_sextile, aes(x=sextile, y=n/6, fill=match, colour=match)) +
  geom_line()+
  scale_x_continuous(breaks = unique(trial_separate),limits=c(0.5,6.5))  +
  scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1),limits=c(0,1))   +
  geom_point(aes(fill=block, colour=block),colour="black",alpha=.5, shape=21, size=3,position=position_dodge(0.2))+
  xlab("Sextiles") +
  ylab("Accuracy") +
  labs(title="Accuracy across pair blocks (6 trials)") +
  facet_grid(rows=vars(block))


hits_per_sextile
#ggsave(hits_per_sextile, file=paste("accuracy_per_6_trials.png",width = 6, height = 6, scale=1))

miss_per_block <- ggplot(summary_stats$miss_per_block, aes(subjID,n,fill=block)) + 
  geom_bar(stat="identity", position="dodge") +
  ylab("N° of misses")
miss_per_block

data_with_cumsum$block <- revalue(data_with_cumsum$block, c("1"="Block 1", "2"="Block 2"))

attach(data_with_cumsum)
for(i in unique(data_with_cumsum$subjID)){
  for (j in unique(data_with_cumsum$block)){
    
    subset_all = subset(data_with_cumsum, subjID==i)
    subset_no_miss = subset(data_with_cumsum, subjID==i & fb!=2)
    cdat <- ddply(subset_no_miss, c("block","fb"), summarise, RT.mean=mean(RT))
    mean_rt_across_blocks <- ggplot(data=subset_no_miss, aes(x=trial, y=RT, fill=fb, colour=fb)) +
      geom_line()+
      geom_hline(data=cdat, aes(yintercept=RT.mean),color=c("black","blue","black","blue"),
                 linetype="dashed", size=0.3) +
      scale_x_continuous(breaks = c(1,7,13,19,25,31),limits=c(-0.5,37))  +
      #scale_y_continuous(breaks = c(0,0.5,1,1.5),limits=c(0,1.75))   +
      #geom_point(aes(fill=block, colour=block),colour="black",alpha=.5, shape=21, size=3,position=position_dodge(0.2))+
      xlab("Trial") +
      ylab("mean RT") +
      labs(title=paste("RT across trials"," Subj ",i,sep="")) +
      #annotate("text", x = c(36,36), y = c(1.8,2), label = cdat$RT.mean) +
      facet_grid(rows = vars(block))
      ggsave(mean_rt_across_blocks, file=paste("rt_across_blocks","_subj_",i,".png", sep=""),width = 6, height = 6, scale=1)
  }
}

cdat <- ddply(subset, c("block","match"), summarise, RT.mean=mean(RT))






#Save as CSV
#setwd(diroutput)
#write.table(Gather,file = paste("performance_all_",task,".txt",sep=""),sep="\t",row.names = FALSE,quote=FALSE)