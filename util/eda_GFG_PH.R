library(readr)
library(rstudioapi)
library(data.table)
library(ggplot2)
library(dplyr)
library(wesanderson)



# Clear all variables
rm(list=ls(all=TRUE))
#set inputs
dirinput <- "N:/Users/phaller/mri_task_analysis/data/children_renamed" 
diroutput <- "N:/Users/phaller/mri_task_analysis/data/visualization"  # output dir must exist
# get list of files 
files <- dir(pattern="C*.txt")

# -------------------------------------------------------------------------------------------------------------------------------------------
# Gather all data
# -------------------------------------------------------------------------------------------------------------------------------------------
setwd(dirinput)
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
    colnames(D) <- colnames(D)[c(1:(which(names(D)=="C/I")-1),(which(names(D)=="C/I")+1):length(colnames(D)))]
    D<- D[,1:(length(D)-1)]
    D[names(D)[names(D) == "RT"]] <-  D[names(D)[names(D) == "RT"]]/1000 # RTs to seconds
    
  
 # if frame does contain the colname resp -> change
 # don't use grep here, because it will find all instances of 'resp' (even colnames with name respOnset) 
  names(D)[names(D) == ",result"] <- "result"
  names(D)[names(D) == "V-file"] <- "vFile"
  names(D)[names(D) == "V-File"] <- "vFile"
  names(D)[names(D) == "A-file"] <- "aFile"
  names(D)[names(D) == "A-File"] <- "aFile"
  names(D)[names(D) == "Vstim"] <- "vStim"
  names(D)[names(D) == "Astim"] <- "aStim"
  
  
  #add subjects
  D <- as_tibble(cbind(rep(substr(files[f],1,7),dim(D)[1]),D)) # use "as_tibble" to avoid issues with the unicode symbols display
  colnames(D)[1] <- "subjID"
  datalist[[f]] <- D
}
Gather <- as_tibble(data.table::rbindlist(datalist, fill=TRUE)) # combine all data frames in one


# -------------------------------------------------------------------------------------------------------------------------------------------
# Plots
# -------------------------------------------------------------------------------------------------------------------------------------------

setwd(diroutput)
#write.table(Gather,file = "performance_all_aduls.txt",sep="\t",row.names = FALSE,quote=FALSE)
#rm(new_D1)
#rm(new_D2)
#rm(new_D3)
subj <- unique(Gather$subjID)
df_subj <- list()
Gather$fb = as.integer(Gather$fb)
#loop thru subject,then block, then symbol
for (i in 1:length(subj)){
  df_subj[[i]] <- subset(Gather, subjID == subj[i])
  blocks <- unique(df_subj[[i]]$block)
  for (j in 1:length(blocks)){
    df_subj_block <- subset(df_subj[[i]],block==j)
    symbols <- unique(df_subj_block$vFile)
    phoneme_set <- unique(df_subj_block$aFile)
    df_subj_block_symbol <- list()
    df_subj_block_phonemes <- list()
    for (k in 1:length(symbols)){
       if (length(unique(df_subj_block$vFile))>=2){
         df_subj_block_symbol[[k]] <- subset(df_subj_block, vFile==symbols[k])
         new_col <- cumsum(df_subj_block_symbol[[k]]$fb)
         new_col_trial <- 1:nrow(df_subj_block_symbol[[k]])
         new_cols <- cbind(new_col,new_col_trial)
         df_subj_block_symbol[[k]]<- as_tibble(cbind(new_cols,df_subj_block_symbol[[k]]))
         colnames(df_subj_block_symbol[[k]])[c(1,2)] <- c("cumsum_fb", "trial_separate")
       }
      if (length(unique(df_subj_block$aFile))>=2){
        df_subj_block_phonemes[[k]] <- subset(df_subj_block, aFile==phoneme_set[[k]])
        new_col_phonemes <- cumsum(df_subj_block_phonemes[[k]]$fb)
        df_subj_block_phonemes[[k]]<- as_tibble(cbind(new_col_phonemes,df_subj_block_phonemes[[k]]))
        colnames(df_subj_block_phonemes[[k]])[1] <- "cumsum_fb_phonemes"
      }
    }
    
    data2plot<- as_tibble(data.table::rbindlist(df_subj_block_symbol) )
    data2plot_phonemes <- as_tibble(data.table::rbindlist(df_subj_block_phonemes) )
  #Count errors missing and correct responses and make a summary line (used in plots)
    nErrors <- length(which(df_subj_block$fb==0))
    nMiss <- length(which(df_subj_block$RT==0))
    nCorr <- length(which(df_subj_block$fb==1))
    counts <- rbind(paste("Errors(",nErrors,")",sep=""),paste("Misses(",nMiss,")",sep=""),paste("Correct(",nCorr,")",sep=""))

# #PLOT CUMULATIVE PROBABILITIES PER STIMULI
# #-------------------------------------------------------------------
    attach(data2plot)
    cumSumPlot <- ggplot(data=data2plot, aes(x=trial_separate, y=cumsum_fb, group=vFile, color=vFile)) +
                   geom_line()+
                   geom_point(aes(fill=vFile),colour="black",alpha=.5, shape=21, size=3,position=position_dodge(0.2))+
                   scale_x_continuous(breaks = unique(trial_separate),limits=c(1,10.5))  +
                   scale_y_continuous(breaks = c(0,1,2,3,4,5,6,7,8,9,10),limits=c(0,10))   +
                   guides(alpha=FALSE)+
                   theme(axis.title = element_text(size=12),
                        title = element_text(size=14),
                        plot.subtitle = element_text(size=14,color="darkblue"),
                        legend.text=element_text(size=12),
                        panel.grid.major = element_line(colour="white"),
                        panel.grid.minor = element_blank(),
                        panel.background = element_rect(fill = "gray88"))+
                  labs(title=paste(subjID[i]," Cumulative sum of hits per visual stimulus",sep=""),subtitle=paste("Block ",j,"(",paste(counts,collapse = " | "),")",sep = ""),size=12)
    #save
    ggsave(cumSumPlot, file=paste("SumHits_",subjID[i],"_B",j,".png", sep=""),width = 6, height = 6, scale=1)

#PLOT CUMULATIVE PROBABILITIES PER PHONEME
#-------------------------------------------------------------------
    attach(data2plot_phonemes)
    cumSumPlot <- ggplot(data=data2plot_phonemes, aes(x=trial_separate, y=cumsum_fb_phonemes, group=aFile, color=aFile)) +
      geom_line()+ 
      geom_point(aes(fill=aFile),colour="black",alpha=.5, shape=21, size=3,position=position_dodge(0.2))+
      scale_x_continuous(breaks = unique(trial_separate),limits=c(1,10.5))  +
      scale_y_continuous(breaks = c(0,1,2,3,4,5,6,7,8,9,10),limits=c(0,10))   +
      guides(alpha=FALSE)+
      theme(axis.title = element_text(size=12),
            title = element_text(size=14),
            plot.subtitle = element_text(size=14,color="darkblue"),
            legend.text=element_text(size=12),
            panel.grid.major = element_line(colour="white"),
            panel.grid.minor = element_blank(),
            panel.background = element_rect(fill = "gray88"))+
      labs(title=paste(subjID[i]," Cumulative sum of hits per phoneme",sep=""),subtitle=paste("Block ",j,"(",paste(counts,collapse = " | "),")",sep = ""),size=12)
    #save 
    ggsave(cumSumPlot, file=paste("SumHitsPhoneme_",subjID[i],"_B",j,".png", sep=""),width = 6, height = 6, scale=1)
         
#PLOT HISTOGRAM: different colors for correct/errors + nErrors
#-------------------------------------------------------------------
   #Plot
    df_subj_block$fb <- as.factor(df_subj_block$fb) 
    attach(df_subj_block) 
    hplot <- ggplot(df_subj_block,aes(x=RT,fill=fb))+
                geom_histogram(color="black",position=position_dodge(0.1),alpha=.4,binwidth = .2,show.legend = TRUE)+
                geom_density(aes(y=.2*..count..,color=fb),alpha=.2,show.legend = TRUE) +
                scale_y_continuous(breaks = c(0,1,2,3,4,5,6,7,8,9,10),limits=c(0,6))  +
                scale_x_continuous(breaks = c(0,0.25,0.5,0.75,1,1.25,1.5,1.75,2,2.25,2.5,2.75,3,3.25,3.5,3.75,4),limits=c(0,4.5))+ 
                   guides(alpha=FALSE) +
                  theme(panel.grid.minor =  element_blank(),
                        panel.grid.major.y = element_line(color="white"),
                        panel.grid.major.x = element_blank(),
                        panel.background = element_rect(fill = "gray89"),
                        axis.title = element_text(size=12),
                        axis.text.x = element_text(angle = 60, hjust = 1),
                        title = element_text(size=14),
                        plot.subtitle = element_text(size=12,color="darkblue"))+
                  labs(title=paste(subjID[i]," RTs distribution - Block ",j,sep=""),subtitle=paste(counts,collapse = " | "),size=12)
            #save 
            ggsave(hplot, file=paste("Hist_",subjID[i],"_B",j,".png", sep=""),width = 6, height = 6, scale=1)
      
   }
  }
