
## --------------------------------------------------------- ##
##
##  Routines to reproduce the results of the paper:
##  "When do parents bury a child? Quantifying uncertainty
##  in the parental age at offspring loss"
## 
##  STEP 6: sensitivity analysis
##
##  NOTE: in order to run the analysis, you should download and 
##  install the MortalitySmooth and svcm packages from the CRAN archive, 
##  available at:
##  https://cran.r-project.org/src/contrib/Archive/MortalitySmooth/
##  https://cran.r-project.org/src/contrib/Archive/svcm/
##  (download the latest .tar.gz files and install them from the
##   Packages/Install/ Install from: Package Archive File)
##
##  Code by Ugofilippo Basellini (2022) unless otherwise stated
##
##  sessionInfo() details:
##
##  R version 4.0.2 (2020-06-22)
##  Platform: x86_64-w64-mingw32/x64 (64-bit)
##  Running under: Windows 10 x64 (build 19044)
##
##  attached base packages:
##  splines   grid    stats   graphics  grDevices utils   datasets 
##  methods   base
## 
##  other attached packages:
##  MortalitySmooth_2.3.4 lattice_0.20-41       svcm_0.1.2           
##  Matrix_1.2-18         fields_11.6           spam_2.5-1           
##  dotCall64_1.0-0       forcats_0.5.1         stringr_1.4.0        
##  dplyr_1.0.7           purrr_0.3.4           readr_1.4.0          
##  tidyr_1.1.4           tibble_3.1.5          ggplot2_3.3.5        
##  tidyverse_1.3.1
##
## --------------------------------------------------------- ##

## cleaning the workspace
rm(list=ls())

## packages
library(tidyverse) 
library(fields)
library(MortalitySmooth)

## saving plots and results?
SAVE.PLOT <- T
SAVE.RES <- T

## start by computing population weights
## for Swedish females in 1900 (for standardization)
cou <- "SWE"
load(paste0("data/clean/",cou,"_arrays_clean.rdata"))

## inputs
ages_mother <- 0:100
cohorts <- as.numeric(colnames(D))[1]:2000

## weights
wei.swe <- POP[,which(cohorts==1900)]/sum(POP[,which(cohorts==1900)])
plot(ages_mother,wei.swe)
sum(wei.swe)

## loading country-specific data and functions
cou <- "SWE"
load(paste0("data/clean/",cou,"_arrays_clean.rdata"))
source("functions/SSEfunctions.R")

## inputs
reprod_ages <- 15:50
ages_child <- 0:(max(ages_mother) - min(reprod_ages))
n_row <- length(ages_mother)
n_col <- length(ages_child)

## age-axis
x <- min(reprod_ages):max(ages_mother)
m <- length(x)

## cohorts
cohorts <- as.numeric(colnames(D))[1]:2000
n <- length(cohorts)

## standardize population counts
## so to have the same distribution of the reference
POP.STD <- matrix(NA,nrow = nrow(POP),ncol = ncol(POP))
i <- 1
for (i in 1:ncol(POP)){
  ## weights current population
  wei.curr <- POP[,i]/sum(POP[,i])
  # plot(ages_mother,wei.swe)
  # lines(ages_mother,wei.curr)
  POP.STD[,i] <- (wei.swe/wei.curr)*POP[,i]
}

## compare original vs standardized distribution
matplot(ages_mother,POP[,1:n],t="l",lty=1,col=rainbow(n))
matplot(ages_mother,POP.STD[,1:n],t="l",lty=1,col=rainbow(n))

## compare population weights in original vs standardized populations
WEI <- WEI.STD <- matrix(NA,nrow = nrow(POP),ncol = ncol(POP))
i <- 1
for (i in 1:ncol(POP)){
  ## weights current population
  WEI[,i] <- POP[,i]/sum(POP[,i])
  WEI.STD[,i] <- POP.STD[,i]/sum(POP.STD[,i])
}

## plotting
matplot(ages_mother,WEI[,1:n],t="l",lty=1,col=rainbow(n))
lines(ages_mother,wei.swe,lty=2,lwd=2)

matplot(ages_mother,WEI.STD[,1:n],t="l",lty=1,col=rainbow(n))
lines(ages_mother,wei.swe,lty=2,lwd=2)


## derive death counts over ages and cohorts 
## (summed up for all child deaths' ages)
## using STANDARDIZED population
Z <- deriveZ(x=x, cohorts = cohorts, H=H, POP=POP.STD, Q=Q,
             reprod_ages = reprod_ages, ages_mother=ages_mother)
image.plot(cohorts,x,t(Z))

## fitting SSE model (decomposition)
fitSSE <- SSEfun(x=x, cohorts = cohorts, Z=Z, cou = cou)

## plotting and saving the SSE outcomes
if (SAVE.PLOT) pdf(file=paste0("figures/sensitivity/",cou,"_fit.pdf"),width = 8, height = 6)
whys <- floor(seq(1,n, length=4))
par(mfrow=c(2,2))
for(i in 1:length(whys)){
  plot(x, Z[,whys[i]], col=8, pch=16,ylim=range(Z[,whys[i]],0),
       main=paste(cohorts[whys[i]]),xlab="age", ylab="counts")
  lines(x, fitSSE$MU.hat[,whys[i]], col=2, lwd=3)
  lines(fitSSE$x1, fitSSE$GAMMA1.hat[,whys[i]], col=3, lwd=2, lty=2)
  lines(fitSSE$x2, fitSSE$GAMMA2.hat[,whys[i]], col=4, lwd=2, lty=2)
}
par(mfrow=c(1,1))
if (SAVE.PLOT) dev.off()

## plotting the two components
matplot(fitSSE$x1, fitSSE$GAMMA1.hat, t="l", col=heat.colors(n), lty=1,
        xlab="age")
matplot(fitSSE$x2, fitSSE$GAMMA2.hat, t="l", col=rainbow(n), lty=1,
        xlab="age")

## computing summary measures
SSEsummary <- SSEsummaryFUN(x=x,cohorts = cohorts,Z=Z,
                            MU.hat = fitSSE$MU.hat,x1=fitSSE$x1,x2=fitSSE$x2,
                            GAMMA1.hat=fitSSE$GAMMA1.hat, GAMMA2.hat=fitSSE$GAMMA2.hat)

## plotting and saving summary measures
if (SAVE.PLOT) pdf(file=paste0("figures/sensitivity/",cou,"_summary.pdf"),width = 8, height = 6)
par(mfrow=c(1,3))
plot(cohorts,SSEsummary$w1,main="Weight Comp 1",t="l",lwd=2)
plot(cohorts,SSEsummary$M1,main="Mean",ylim=range(SSEsummary$M1,SSEsummary$M2),t="l",lwd=2)
lines(cohorts,SSEsummary$M2,col=2,lwd=2)
legend("right",c("Comp 1","Comp 2"),col=1:2,lty=1,lwd=2)
plot(cohorts,SSEsummary$SD1,main="Stand Dev",ylim=range(SSEsummary$SD1,SSEsummary$SD2),t="l",lwd=2)
lines(cohorts,SSEsummary$SD2,col=2,lwd=2)
par(mfrow=c(1,1))
if (SAVE.PLOT) dev.off()

## saving results data
if (SAVE.RES) save.image(file=paste0("results/sensitivity/",cou,"_SSE.Rdata"))

## END