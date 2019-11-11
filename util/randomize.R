# Randomizer for feedback learning task
# Patrick Haller, July 2019

library(tidyr)
library(ggplot2)
library(dummies)
library(viridis)
library(plyr)

######### KIDS VERSION #########

# 4 sound-symbol pairs

# all possible combinations for each auditory stimulus
# first digit: sound, second digit: left symbol, third digit: right symbol
a1 <- c(112,113,114,121,131,141)
a2 <- c(221,223,224,212,232,242)
a3 <- c(331,332,334,313,323,343)
a4 <- c(441,442,443,414,424,434)

triples <- list(a1,a2,a3,a4)

shuffle <- function(triples){
  all <- sample(triples)
  seq <- c()
  seq <- c(
    # miniblock 1: 3xa1, 2xa2, 2xa3,3xa4
    # RESTRICTION: miniblock begins with sequence of 3x the same auditory stimulus
    c(sample(all[[1]],3),
      sample(c(sample(all[[2]],2),sample(all[[3]],2),sample(all[[4]],3)))),
    # miniblock 2: 2xa1, 4xa2, 2xa3,2xa4
    c(sample(all[[2]],3),
      sample(c(sample(all[[2]],1),sample(all[[1]],2),sample(all[[3]],2),sample(all[[4]],2)))),
    # miniblock 3: 3xa1, 2xa2, 3xa3,2xa4
    c(sample(all[[3]],3),
      sample(c(sample(all[[2]],2),sample(all[[1]],3),sample(all[[4]],2)))),
    # miniblock 4: 2xa1, 2xa2, 3xa3,3xa4
    # RESTRICTION: miniblock ends with sequence of 3x the same auditory stimulus
    c(sample(c(sample(all[[2]],2),sample(all[[1]],2),sample(all[[3]],3))),
      sample(all[[4]],3))
  )
  
  df <- data.frame(trial=seq.int(1,40),triple=as.character(seq), astim=factor(substr(as.character(seq),1,1)),dummy_y=1)
  #compute identical n-grams for each astim
  identical_seq <- rle(as.character(df$astim))
  ngrams <- as.character(identical_seq$lengths)
  counts <- count(ngrams)
  # RESTRICTION: No stimulus is presented more than 3 times after each other
  if (!(4 %in% identical_seq$lengths|5 %in% identical_seq$lengths|6 %in% identical_seq$lengths)){
    # RESTRICTION: The same stimulus is presented twice in a row at most 4 times
    if (table$freq[2] <= 4){
      return(df)
    }
  }
  else
    shuffle(triples)
}


# sample 20 different pseudo-randomized sequences and write them out
trials <- list()
for (i in 1:20){
  myseq <- shuffle(triples)
  trials[[i]] <- as.character(myseq$triple)
}

lapply(trials, write, "randomized_sequences.txt", append=TRUE, ncolumns=40)

# check if each stimulus is presented the same amount of times
table(myseq$astim)

ggplot(myseq,aes(trial, fill=astim)) +
  geom_bar() + 
  scale_color_viridis(discrete = TRUE, option = "E")+
  scale_fill_viridis(discrete = TRUE) 

ggplot(myseq,aes(trial, color=astim,fill=astim)) +
  stat_density(position="identity",adjust = 1/2,alpha=0.1)+
  scale_color_viridis(discrete = TRUE, option = "E")+
  scale_fill_viridis(discrete = TRUE) 


########## ADULT VERSION ###########

# 5 auditory, 8 visual (3 distractors)

all_1 <- c()

nd_1 <- c(112,113,114,115,121,131,141,151)
nd_2 <- c(221,223,224,225,212,232,242,252)
nd_3 <- c(331,332,334,335,313,323,343,353)
nd_4 <- c(441,442,443,445,414,424,434,454)
nd_5 <- c(551,552,553,554,515,525,535,545)

d_1 <- c(116,117,118,161,171,181)
d_2 <- c(226,227,228,262,272,282)
d_3 <- c(336,337,338,363,373,383)
d_4 <- c(446,447,448,464,474,484)
d_5 <- c(556,557,558,565,575,585)


nd_all <- list(nd_1,nd_2,nd_3,nd_4,nd_5)
d_all <- list(d_1,d_2,d_3,d_4,d_5)

nd_sample <- c()
for(nd in nd_all){
  nd_sample <- c(nd_sample,sample(nd,size=4))
}

d_sample <- c()
for(d in d_all){
  d_sample <- c(d_sample,sample(d,size=4))
}

sample_both <- c(nd_sample,d_sample)
final_sample_a <- sample(sample_both)

all_1 <- c(all_1,final_sample_a)
all_1
