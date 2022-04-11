
## --------------------------------------------------------- ##
##
## Routines to reproduce the results of the paper:
## "When do parents bury a child? Quantifying uncertainty
## in the parental age at offspring loss"
## 
## STEP 3: decomposition analysis
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

## loading data and functions
cou <- "SWE"
load(paste0("data/clean/",cou,"_arrays_clean.rdata"))
source("functions/SSEfunctions.R")

## inputs
ages_mother <- 0:100
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

## derive death counts over ages and cohorts 
## (summed up for all child deaths' ages)
Z <- deriveZ(x=x, cohorts = cohorts, H=H, POP=POP, Q=Q,
             reprod_ages = reprod_ages, ages_mother=ages_mother)
image.plot(cohorts,x,t(Z))

## fitting SSE model (decomposition)
fitSSE <- SSEfun(x=x, cohorts = cohorts, Z=Z, cou = cou)

## plotting and saving the SSE outcomes
if (SAVE.PLOT) pdf(file=paste0("figures/decomp/",cou,"_fit.pdf"),width = 8, height = 6)
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
if (SAVE.PLOT) pdf(file=paste0("figures/decomp/",cou,"_summary.pdf"),width = 8, height = 6)
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
if (SAVE.RES) save.image(file=paste0("results/decomp/",cou,"_SSE.Rdata"))

## END