library(readr)
library(rstudioapi)
library(readr)
library(data.table)
library(plyr)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(tidyr)
library(wesanderson)
library(viridis)

task <- "fbl_kloten"
#set inputs
#dirinput <- "N:/Users/phaller/mri_task_analysis/data/piloting/piloting_kloten/2x3_24"
#dirinput <- "N:/Users/phaller/mri_task_analysis/data/piloting/piloting_kloten/2x3_32"
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
  for(i in unique(data$subj_idx)){
    df_subj[[i]] <- subset(data, subj_idx == i)
    for(j in unique(df_subj[[i]]$block)){
      df_subj_block <- subset(df_subj[[i]],block==j)
      # compute cumulative sum for each auditory stimulus in a given block of a subject
      for(k in unique(df_subj_block$aStim)){
        df_subj_block_astim <- list()
        df_subj_block_astim[[k]] <- subset(df_subj_block, aStim==k )
        new_col <- cumsum(df_subj_block_astim[[k]]$fbprime)
        new_col_trial <- 1:nrow(df_subj_block_astim[[k]])
        new_cols <- cbind(new_col,new_col_trial)
        df_subj_block_astim[[k]]<- as_tibble(cbind(new_cols,df_subj_block_astim[[k]]))
        colnames(df_subj_block_astim[[k]])[c(1,2)] <- c("cumulsum_fb", "trial_separate")
        new_data <- as_tibble(rbind(new_data,df_subj_block_astim[[k]]))
        # reorder data
        new_data <- new_data[
          with(new_data, order(subj_idx, block,trial)),
          ]
      }
    }
  }
  return(as_tibble(new_data))
}

get_summary_stats <- function(data){
  ### #How many missing responses per block (fb == 2)
  miss_per_block <- data %>%
    select(subj_idx,RT, fb,block) %>%
    filter(fb==2)  %>%
    group_by(subjID,block) %>%
    tally()
  
  ### #RTs per block
  RT_per_block <- data %>%
    select(subj_idx,RT, fb,block) %>%
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

split_trials_24 <- function(data){
  data$trial <- as.integer(data$trial)
  data$quartile <- 0
  data[which(data$trial <= 6),]$quartile = 1
  data[which(data$trial > 6 & data$trial <= 12),]$quartile = 2
  data[which(data$trial > 12 & data$trial <= 18),]$quartile = 3
  data[which(data$trial > 18 & data$trial <= 24),]$quartile = 4
  data <- data[
    with(data, order(subj_idx, block,trial)),
    ]
  return(data)
}

split_trials_32 <- function(data){
  data$trial <- as.integer(data$trial)
  data$quartile <- 0
  data[which(data$trial <= 8),]$quartile = 1
    data[which(data$trial > 8 & data$trial <= 16),]$quartile = 2
  data[which(data$trial > 16 & data$trial <= 24),]$quartile = 3
  data[which(data$trial > 24 & data$trial <= 32),]$quartile = 4
  data <- data[
    with(data, order(subj_idx, block,trial)),
    ]
  return(data)
}

split_trials_octile <- function(data){
  data$trial <- as.integer(data$trial)
  data$octile <- 0
  data[which(data$trial <= 8),]$octile = 1
  data[which(data$trial > 8 & data$trial <= 16),]$octile = 2
  data[which(data$trial > 16 & data$trial <= 24),]$octile = 3
  data[which(data$trial > 24 & data$trial <= 32),]$octile = 4
  data <- data[
    with(data, order(subj_idx, block,trial)),
    ]
  return(data)
}

wes_cols= wes_palette("GrandBudapest1", n = 2)

## load data
data <- gather_data(files)

# redefine fb for visualization
data$fbprime <- rep("NA",nrow(data))
data[which(data$fb==0),]$fbprime = -1
data[which(data$fb==1),]$fbprime = 1
data[which(data$fb==2),]$fbprime = 0

# for the visualization, for the subject who completed block 3, b3 will count as b2
data[which(data$block==3),]$block <- rep(2,nrow(data[which(data$block==3),]))
data$block = as.factor(data$block)
levels(data$block) <- c("Block 1", "Block 2")
data$aStim = as.factor(data$aStim)

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

mean_accuracy <-  data %>%
  select(subj_idx, fb,block) %>%
  filter(fb==1) %>%
  group_by(subj_idx,block) %>%
  tally() 

mean_accuracy$n <- mean_accuracy$n/32
ggplot(mean_accuracy, aes(x=block, y=n, colour=block)) + 
  geom_boxplot(outlier.colour="red", outlier.shape=8,
               outlier.size=4) + 
  scale_y_continuous(limits=0:1) +
  ylab("accuracy") +
  ggtitle("Mean accuracy per block (2x4[32])") +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5) + 
  scale_color_manual(values=c("#999999", "#E69F00"))

# compute mean accuracy per block
acc_b1 <- mean_accuracy[which(mean_accuracy$block=="Block 1"),]
mean_acc_b1 <- mean(acc_b1$n)
sd_acc_b1 <- sd(acc_b1$n)

acc_b2 <- mean_accuracy[which(mean_accuracy$block=="Block 2"),]
mean_acc_b2 <- mean(acc_b2$n)
sd_acc_b2 <- sd(acc_b2$n)

overall_mean <- mean(mean_accuracy$n)
overall_sd <- sd(mean_accuracy$n)

# compute mean reaction time for each participant for each block
data_for_rt_plot <- data[which(data$fb!=2),]
mean_rts = aggregate(data_for_rt_plot$rt,
                by = list(subj_idx = data_for_rt_plot$subj_idx, block = data_for_rt_plot$block),
                FUN = mean)

ggplot(mean_rts, aes(x=block, y=x, colour=block)) + 
  geom_boxplot(outlier.colour="red", outlier.shape=8,
               outlier.size=4) + 
  scale_y_continuous(limits=c(0.6,1.8)) +
  ylab("RT") +
  ggtitle("Mean reaction time per block (2x4[32])") +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5) + 
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))

rt_b1 <- mean_rts[which(mean_rts$block=="Block 1"),]
mean_rt_b1 <- mean(rt_b1$x)
sd_rt_b1 <- sd(rt_b1$x)

rt_b2 <- mean_rts[which(mean_rts$block=="Block 2"),]
mean_rt_b2 <- mean(rt_b2$x)
sd_rt_b2 <- sd(rt_b2$x)

overall_mean_rt <- mean(mean_rts$x)
overall_sd_rt <- sd(mean_rts$x)

# compute cumulative sums and quartiles
#data_no_miss <- data[which(data$fb!=2),]
data_with_cumulsum <- compute_cumulative_sums(data)
data_octile <- split_trials_octile(data_with_cumulsum) # make sure you uncomment the last line of the function with 32-version

###### get individual learning trajectories
attach(data_octile)
ggplot(data=data_octile[which(data_octile$block=="Block 1"),], aes(x=trial_separate, y=cumulsum_fb, group=aStim, color=aStim)) +
  geom_line()+
  geom_point(aes(fill=aStim),colour="black",alpha=.5, shape=21, size=3,position=position_dodge(0.2))+
  scale_x_continuous(breaks = unique(trial_separate),limits=c(0.5,13.5))  +
  scale_y_continuous(breaks = -13:13,limits=c(-13,13)) +
  ylab("Cumulative sum of feedback")+
  xlab("Trial per stimulus") +
  guides(alpha=FALSE)+
  theme(axis.title = element_text(size=12),
        title = element_text(size=14),
        plot.subtitle = element_text(size=14,color="darkblue"),
        legend.text=element_text(size=12),
        panel.grid.major = element_line(colour="white"),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "gray88"))+
  facet_wrap( ~ subj_idx, ncol=3) +
  labs(title="Cumulative sum of hits per astim (2x4[32][B1])")



attach(data_octile)
ggplot(data=data_octile[which(data_octile$block==1),], aes(x=trial_separate, y=cumulsum_fb, group=aStim, color=aStim)) +
  geom_line()+
  geom_point(aes(fill=aStim),colour="black",alpha=.5, shape=21, size=3,position=position_dodge(0.2))+
  scale_x_continuous(breaks = unique(trial_separate),limits=c(1,10.5))  +
  scale_y_continuous(breaks = c(0,1,2,3,4,5,6,7,8,9,10),limits=c(0,10)) +
  ylab("Cumulative sum of feedback")+
  xlab("Trial per stimulus") +
guides(alpha=FALSE)+
  theme(axis.title = element_text(size=12),
        title = element_text(size=14),
        plot.subtitle = element_text(size=14,color="darkblue"),
        legend.text=element_text(size=12),
        panel.grid.major = element_line(colour="white"),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "gray88"))+
  facet_wrap( ~ subj_idx, ncol=3) +
  labs(title="Cumulative sum of hits per auditory stimulus (2x3[24][Block 1])")


###### get individual learning trajectories
attach(data_octile)
ggplot(data=data_octile[which(data_octile$block==2),], aes(x=trial_separate, y=cumulsum_fb, group=aStim, color=aStim)) +
  geom_line()+
  geom_point(aes(fill=aStim),colour="black",alpha=.5, shape=21, size=3,position=position_dodge(0.2))+
  scale_x_continuous(breaks = unique(trial_separate),limits=c(1,10.5))  +
  scale_y_continuous(breaks = c(0,1,2,3,4,5,6,7,8,9,10),limits=c(0,10)) +
  ylab("Cumulative sum of feedback")+
  xlab("Trial per stimulus") +
  guides(alpha=FALSE)+
  theme(axis.title = element_text(size=12),
        title = element_text(size=14),
        plot.subtitle = element_text(size=14,color="darkblue"),
        legend.text=element_text(size=12),
        panel.grid.major = element_line(colour="white"),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "gray88"))+
  facet_wrap( ~ subj_idx, ncol=3) +
  labs(title="Cumulative sum of hits per auditory stimulus (2x3[24][Block 2])")
#facet_wrap( ~ subj_idx, ncol=3)

## correct responses per quartile
correct_per_octile <-  data_octile %>%
  select(subj_idx, fb,block,octile) %>%
  filter(fb==1) %>%
  group_by(subj_idx,block,octile) %>%
  tally() 

correct_per_octile$octile = as.factor(correct_per_octile$octile)


ggplot(correct_per_octile, aes(octile,n)) +
  geom_boxplot(fill="lightgray")+
  geom_line(aes(group=subj_idx,col=subj_idx),alpha=0.5,linetype="dashed",position=position_jitter(w=0.05, h=0.05)) +
  scale_color_viridis(discrete = TRUE, option = "C")+
  scale_fill_viridis(discrete = TRUE) +
  geom_point(aes(group=subj_idx,shape=subj_idx,col=subj_idx),size=2) +
  scale_shape_manual(values=1:nlevels(correct_per_octile$subj_idx)) +
  xlab("Quartile") +
  ylab("Number of correct responses") +
  ggtitle("Correct responses per quartile for each block (2x4)[24]") +
  theme_bw()+
  stat_summary(fun.y=median, geom="smooth", aes(group=0),lwd=0.8,col=wes_cols[2])+
  #geom_jitter(aes(col=subj_idx))+
  facet_grid(.~block)

ggplot(correct_per_octile, aes(octile,n,group=block)) +
  geom_line(aes(color=block),linetype="dashed",alpha=1) +
  #stat_summary(fun.y=median, geom="smooth", aes(group=0),lwd=0.8,col=wes_cols[2])+
  #scale_color_viridis(discrete = TRUE, option = "C")+
  #scale_fill_viridis(discrete = TRUE) +
  #geom_point(aes(group=subj_idx,shape=subj_idx,col=subj_idx),size=2) +
  #scale_shape_manual(values=1:nlevels(correct_per_octile$subj_idx)) +
  xlab("Quartile") +
  ylab("Number of correct responses") +
  ggtitle("Correct responses per quartile for each subject per block (2x4)[32]") +
  facet_wrap( ~ subj_idx, ncol=3)

  #geom_jitter(aes(col=subj_idx))+


miss_per_block <- data %>%
  select(subj_idx,rt, fb,block) %>%
  filter(fb==2)  %>%
  group_by(subj_idx,block) %>%
  tally()
miss_per_block$n <- miss_per_block$n/24
sum(miss_per_block$n)/nrow(miss_per_block)

#######

setwd(diroutput)
write_csv(data,path = paste("summarized_performance_",task,".csv",sep=""),col_names = TRUE,quote=FALSE)

######

attach(data_with_cumulsum)
for(i in unique(data_with_cumulsum$subjID)){
  subset_all = subset(data_with_cumulsum, sub_idx==i & fb!=2)
  subset_all$fb = factor(subset_all$fb, labels=c("Neg","Pos"))
  p = ggplot(subset_all, aes(x=rt, fill=fb, color=fb)) +
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