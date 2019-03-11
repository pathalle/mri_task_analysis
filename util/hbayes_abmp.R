# Install all necessary packages into a path that doesn't include the $ sign
# (otherwise you will run into problems with rstan)

 # check current library path (if default is already a valid name, jump the next section) 
.libPaths()
# output may look like this
#[1] "\\\\idnetapp-homes3.uzh.ch/phalle$/Documents/R/win-library/3.5"
#[2] "C:/Program Files/R/R-3.5.2/library" 

assign(".lib.loc", "C:/Program Files/R/R-3.5.2/library", envir = environment(.libPaths))


#install.packages("rstan", lib="C:\\Program Files\\R\\R-3.5.2\\library")
#install.packages("StanHeaders", lib="C:\\Program Files\\R\\R-3.5.2\\library")
#install.packages("rstantools", lib="C:\\Program Files\\R\\R-3.5.2\\library")
#install.packages("hBayesDM", lib="C:\\Program Files\\R\\R-3.5.2\\library")
#install.packages("Rcpp", lib="C:\\Program Files\\R\\R-3.5.2\\library")

# after successful installation, load required packages
library("StanHeaders", lib.loc="C:/Program Files/R/R-3.5.2/library")
library("rstan", lib.loc="C:/Program Files/R/R-3.5.2/library")
library("Rcpp", lib.loc="C:/Program Files/R/R-3.5.2/library")
library("hBayesDM", lib.loc="C:/Program Files/R/R-3.5.2/library")

# get sample files

#dataPath = system.file("extdata/gng_exampleData.txt", package="hBayesDM")
#dataPath = "N:/Users/phaller/R_hBayesDM/gng_exampleData.txt"
#output1 = gng_m1(data=dataPath, niter=2000, nwarmup=500, nchain=4, ncore=4)

# ------------------------------------------------------------------------------------------------
# Gather all data (children)
# ------------------------------------------------------------------------------------------------
dataPath = "N:/Users/phaller/mri_task_analysis/data/Rita"
files <- dir(pattern="*LSB.txt")

setwd(dataPath)
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
  names(D)[names(D) == "resp"] <- "choice"
  
  
  #add subjects
  D <- as_tibble(cbind(rep(substr(files[f],1,7),dim(D)[1]),D)) # use "as_tibble" to avoid issues with the unicode symbols display
  colnames(D)[1] <- "subjID"
  datalist[[f]] <- D
}
Gather <- as_tibble(data.table::rbindlist(datalist, fill=TRUE)) # combine all data frames in one
gather_to_txt <- Gather[,c(grep("subjID",colnames(Gather)),grep("choice",colnames(Gather)),grep("RT",colnames(Gather)))]
gather_to_txt <- gather_to_txt[which(gather_to_txt$RT!=0),]
write.table(gather_to_txt,file = "model_input.txt",sep="\t",row.names = FALSE,quote=FALSE)


# ------------------------------------------------------------------------------------------------
# Fit model
# ------------------------------------------------------------------------------------------------

# exclude RT < 200ms

ddm_model = choiceRT_ddm(
  paste(dataPath, "model_input.txt", sep="/"), saveDir = dataPath, niter=4000, nwarmup=1000, nchain=4, ncore=4, inits="fixed",
  RTbound=0.01)

#save(ddm_model, file = paste(dataPath, "model_children_190305.RData", sep="/"))
load(paste(dataPath2, "model_children_190305.RData", sep="/"))

plot(ddm_model)

ddm_model$allIndPars # all posterior parameters

ddm_model$fit
plot(ddm_model$fit)

plot(ddm_model, type="trace", inc_warmup=T)  
plot(ddm_model, type="trace", fontSize=11)   # traceplot of hyper parameters. Set font size 11.

plotInd(ddm_model, "tau")


# ------------------------------------------------------------------------------------------------
# fit adult data
# ------------------------------------------------------------------------------------------------
dataPath2 = "N:/Users/phaller/mri_task_analysis/model/choiceRT"

ddm_model2 = choiceRT_ddm(
  data=paste(dataPath2, "adults_190305_model_input.txt", sep="/"), niter=4000, nwarmup=1000, nchain=4, ncore=4, inits="random",
  RTbound=0.01)

save(ddm_model2, file = paste(dataPath2, "model_adults_190305.RData", sep="/"))
load(paste(dataPath2, "model_adults_190305.RData", sep="/"))

plot(ddm_model2)

ddm_model2$allIndPars # all posterior parameters

ddm_model2$fit
plot(ddm_model2$fit)

plot(ddm_model2, type="trace", inc_warmup=T)  
plot(ddm_model2, type="trace", fontSize=11)   # traceplot of hyper parameters. Set font size 11.

plotInd(ddm_model2, "tau")

plotInd(ddm_model2, "delta")

# ------------------------------------------------------------------------------------------------
# group comparison
# ------------------------------------------------------------------------------------------------

diffDistAlpha = ddm_model$parVals$mu_alpha - ddm_model2$parVals$mu_alpha
HDIofMCMC(diffDistAlpha)
plotHDI(diffDistAlpha)
plotInd(ddm_model, "alpha")
plotInd(ddm_model2, "alpha")

plotInd(ddm_model, "delta")
ddm_model$rawdata
plotInd(ddm_model2, "delta")

diffDistDelta = ddm_model$parVals$mu_delta - ddm_model2$parVals$mu_delta
HDIofMCMC(diffDistDelta)
plotHDI(diffDistDelta)

# ------------------------------------------------------------------------------------------------
# some sanity checks: is simulated data similar to actual data?
# ------------------------------------------------------------------------------------------------

