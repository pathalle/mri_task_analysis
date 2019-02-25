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

dataPath = system.file("extdata/gng_exampleData.txt", package="hBayesDM")
dataPath = "N:/Users/phaller/R_hBayesDM/gng_exampleData.txt"
output1 = gng_m1(data=dataPath, niter=20, nwarmup=5, nchain=4, ncore=4)

output1$allIndPars
output1$fit
plot(output1, type="trace", fontSize=11)   # traceplot of hyper parameters. Set font size 11.

dataPath = "N:/Users/phaller/mri_task_analysis/util/output/performance_all_fbl.txt"
ddm_model = choiceRT_ddm(data=dataPath, niter=4000, nwarmup=1000, nchain=4, ncore=4, inits="random")
