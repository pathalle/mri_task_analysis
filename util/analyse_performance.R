library(readr)
library(rstudioapi)
library(readr)
library(data.table)
library(plyr)
library(dplyr)
library(ggplot2)
library(gridExtra)


task <- "fbl_kloten"
#set inputs
#dirinput <- "N:/Users/phaller/mri_task_analysis/data/piloting/piloting_kloten/2x3"
dirinput <- "N:/Users/phaller/mri_task_analysis/data/piloting/piloting_kloten/2x4"
# make sure output directory exists already
diroutput <- "N:/Users/phaller/mri_task_analysis/data/piloting/analysis/"

setwd(dirinput)
files <- dir(pattern=".txt", recursive=TRUE)

gather_data <- function(files){
  # summarize all data in 1 data frame
  datalist <- list()
  for (i in 1:length(files)){
    no_col <- max(count.fields(files[i], sep = "\t"))
    D <- read_delim(
      files[i],"\t", escape_double = FALSE, locale = locale(), trim_ws = TRUE)
    D <- cbind(rep(substr(files[i],17,22),dim(D)[1]),D)
    #D<-D[D$resp!=0,] # remove 'too slow ' responses 
    ### Rename and transform some columns
    colnames(D)[1] <- "subj_idx"
    #D[,1] <- strtoi(D[,1])
    D$trial = as.integer(D$trial)
    D$aStim = as.integer(D$aStim)
    #colnames(D)[grep("rt",colnames(D))] <- "RT" 
    D[grep("rt",colnames(D))] <- D[grep("rt",colnames(D))]/1000 # RTs in seconds
    # don't use grep here, because it will find all instances of 'resp' (even colnames with name respOnset) 
    #names(D)[names(D) == "fb"] <- "response"
    #D <- as_tibble(cbind(D,paste(D$vFile,D$aFile)))
    #colnames(D)[ncol(D)] <- "pair"
    datalist[[i]] <- D
  }
  transformed <- data.table::rbindlist(datalist) # combine all data frames in on
  return(transformed)
}
 
compute_cumulative_sums <- function(data){
  df_subj <- list()
  #data<-data[data$choice!=0,] # remove 'too slow ' responses 
  new_data <- data[FALSE,]
  new_cols <- data.frame(cumsum_fb = integer(0), trial_separate = integer(0))
  new_data <- cbind(new_cols,new_data)
  for(i in unique(data$subjID)){
    df_subj[[i]] <- subset(data, subjID == i)
    for(j in unique(df_subj[[i]]$block)){
      df_subj_block <- subset(df_subj[[i]],block==j)
      # compute cumulative sum for each auditory stimulus in a given block of a subject
      for(k in unique(df_subj_block$aStim)){
        df_subj_block_astim <- list()
        df_subj_block_astim[[k]] <- subset(df_subj_block, aStim==k )
        new_col <- cumsum(df_subj_block_astim[[k]]$fb)
        new_col_trial <- 1:nrow(df_subj_block_astim[[k]])
        new_cols <- cbind(new_col,new_col_trial)
        df_subj_block_astim[[k]]<- as_tibble(cbind(new_cols,df_subj_block_astim[[k]]))
        colnames(df_subj_block_astim[[k]])[c(1,2)] <- c("cumulsum_fb", "trial_separate")
        new_data <- as_tibble(rbind(new_data,df_subj_block_astim[[k]]))
        # reorder data
        new_data <- new_data[
          with(new_data, order(subjID, block,trial)),
          ]
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
  correct_per_quartile <- data %>%
    select(subjID, fb,block, quartile) %>%
    filter(fb==1) %>%
    group_by(subjID,block,quartile) %>%
    tally()
  
  return(list(
    "miss_per_block"=miss_per_block,"rt_per_block"=RT_per_block,"rt_across_blocks"=RT_across_blocks,"hits_per_sextile"=correct_per_sextile))
}

split_trials <- function(data){
  data$trial <- as.integer(data$trial)
  data$quartile <- 0
  data[which(data$trial <= 10),]$quartile = 
  data[which(data$trial > 10 & data$trial <= 20),]$quartile = 2
  data[which(data$trial > 20 & data$trial <= 30),]$quartile = 3
  data[which(data$trial > 30 & data$trial <= 40),]$quartile = 4
  data <- data[
    with(data, order(subjID, block,trial)),
    ]
  return(data)
}


## prepare data
data <- gather_data(files)
#data <- data[which(data$rt > 0.15),]
#data <- data[which(data$block ==1)]
#data_with_cumulsum <- compute_cumulative_sums(data)
#data_quartiles <- split_trials(data_with_cumulsum)
#summary_stats <- get_summary_stats(data_sextiles)
#data_with_cumsum <- compute_cumulative_sums(data_sextiles)
# add column with information to which sextile a sequence of trials belongs"
#Save as CSV

## show how many observation (=trials) per subj per block
trials_per_subj_per_block <- data %>%
  select(subj_idx,block) %>%
  group_by(subj_idx,block) %>%
  tally() 

trials_per_subj <- data %>%
  select(subj_idx) %>%
  group_by(subj_idx) %>%
  tally() 

fb_per_subj_per_block <- data %>%
  select(subj_idx, fb,block) %>%
  group_by(subj_idx,block) %>%
  tally() 

plotting <-  data %>%
  select(subj_idx, fb,block) %>%
  filter(fb==1) %>%
  group_by(subj_idx,block) %>%
  tally() 

plotting$n <- plotting$n/32
plotting$block = as.factor(plotting$block)
ggplot(plotting, aes(x=block, y=n, colour=block)) + 
  geom_boxplot(outlier.colour="red", outlier.shape=8,
               outlier.size=4) + 
  scale_y_continuous(limits=0:1) +
  ylab("accuracy") +
  ggtitle("Mean accuracy per block (2x4)") +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5) + 
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))

b1 <- plotting[which(plotting$block==1),]
sum(b1$n)/nrow(b1)

b2 <- plotting[which(plotting$block==2),]
sum(b2$n)/nrow(b2)

b3 <- plotting[which(plotting$block==3),]
sum(b3$n)/nrow(b3)

overall_mean <- sum(plotting$n)/nrow(plotting)

data_for_rt_plot <- data[which(data$fb!=2),]
mean_rts <- data_for_rt_plot %>% 
  summarize(average=mean(rt))


ggplot(data=data_for_rt_plot, aes(x=trial, y=rt, fill=subj_idx, colour=subj_idx)) +
  geom_line()+
  scale_x_continuous(breaks = c(1,4,8,12,16,20,24,28,32),limits=c(-0.5,32))  +
  #scale_y_continuous(breaks = c(0,0.5,1,1.5),limits=c(0,1.75))   +
  #geom_point(aes(fill=block, colour=block),colour="black",alpha=.5, shape=21, size=3,position=position_dodge(0.2))+
  xlab("Trial") +
  ylab("mean RT")


setwd(diroutput)
write_csv(data,path = paste("summarized_performance_",task,".csv",sep=""),col_names = TRUE,quote=FALSE)

#######




attach(data_with_cumulsum)
for(i in unique(data_with_cumulsum$subjID)){
  subset_all = subset(data_with_cumulsum, subjID==i & fb!=2)
  subset_all$fb = factor(subset_all$fb, labels=c("Neg","Pos"))
  p = ggplot(subset_all, aes(x=RT, fill=fb, color=fb)) +
    geom_histogram(fill="white", alpha=0.5, position="dodge") +
    geom_density(alpha=.2) +
    facet_grid(rows = vars(block))
  ggsave(p, file=paste("rt_distribution","_subj_",i,".png", sep=""),path= diroutput,width = 3, height = 6, scale=1)
}


attach(data_with_cumulsum)
for(i in unique(data_with_cumulsum$subjID)){
  subset_all = subset(data_with_cumulsum, subjID==i)
  #subset_no_miss = subset(data_with_cumulsum, subjID==i & fb!=2)
  cdat <- ddply(subset_all, c("block","fb"), summarise, RT.mean=mean(RT))
  mean_rt_across_blocks <- ggplot(data=subset_all, aes(x=trial, y=RT, fill=fb, colour=fb)) +
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
  ggsave(mean_rt_across_blocks, file=paste("rt_across_blocks","_subj_",i,".png", sep=""),path= diroutput,width = 6, height = 6, scale=1)
}

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



cdat <- ddply(subset, c("block","match"), summarise, RT.mean=mean(RT))
#Save as CSV
#setwd(diroutput)
#write.table(Gather,file = paste("performance_all_",task,".txt",sep=""),sep="\t",row.names = FALSE,quote=FALSE)