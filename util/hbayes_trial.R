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
# Fit model
# ------------------------------------------------------------------------------------------------

# exclude RT < 200ms

dataPath = "N:/Users/phaller/mri_task_analysis/model/choiceRT"

ddm_model = choiceRT_ddm(
  paste(dataPath, "children_190305_model_input.txt", sep="/"), niter=4000, nwarmup=1000, nchain=4, ncore=4, inits="fixed",
  RTbound=0.01)

# save model / if already fitted load.

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

