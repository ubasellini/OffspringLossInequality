

## function to derive death counts over ages and cohorts
deriveZ <- function(x,cohorts,H,POP,Q,
                    reprod_ages=15:50,ages_mother=0:100){

  ## dimensions
  m <- length(x)
  n <- length(cohorts)
  ## child ages and NU dimensions
  ages_child <- 0:(max(ages_mother) - min(reprod_ages))
  n_row <- length(ages_mother)
  n_col <- length(ages_child)
  ## empty matrix to store results
  Z <- matrix(NA,m,n)
  jj <- 1
  for (jj in 1:n){
    my.coh <- cohorts[jj]
    whi.coh <- which(cohorts==my.coh)
    myH <- H[,1:n_col,whi.coh]
    mypop <- POP[,whi.coh]
    myQ <- matrix(NA,n_row,n_col)
    ## adjusting D 
    for(j in 1:n_col){
      for(i in 1:n_row){
        wrong <- 
          (i - 1) < min(reprod_ages) | 
          (i - 1) > max(ages_mother) | 
          (j - 1) < max(c((i - 1) - max(reprod_ages), 0)) | 
          (j - 1) > ((i - 1) - min(reprod_ages))
        if(wrong){
          myQ[i,j] <- NA
        }else{
          myQ[i,j] <- Q[j, whi.coh + i - j]
        }
      } # end rows
    } # end cols
    
    ## pop as matrix
    m <- ncol(myH)
    MYPOP <- mypop %*% t(rep(1,m))
    ## multiplying matrices
    NU <- myH*MYPOP*myQ
    nu <- rowSums(NU,na.rm = T)
    nu2 <- nu[ages_mother%in%x]
    Z[,jj] <- nu2
    
  }
  return(Z)
}

## function for estimating the SSE
SSEfun <- function(x,cohorts,Z,cou,
                   lambdas=NULL, kappa=NULL,max.it=NULL){
  
  ## dimensions
  y <- cohorts
  m <- length(x)
  n <- length(cohorts)
  mn <- m*n
  ## regression weights in case of zero death counts
  WEI <- matrix(1,m,n)
  WEI[Z==0] <- 0
  wei <- c(WEI)
  ## response in column-vector
  z <- c(Z)
  
  ## select ages for each component
  ## infant component
  staINF <- min(x)
  endINF <- max(x) - 20
  ## adult component
  staAGI <- min(x) + 10
  endAGI <- max(x)
  ## select data and create weights for each age-window
  ## infant component
  p1 <- c(staINF, endINF)
  x1 <- x[x>=p1[1] & x<=p1[2]]
  W1 <- matrix(as.numeric(x%in%x1), m, n)
  Z1 <- matrix(Z[which(W1==1)], ncol=n)
  ## adult component
  p2 <- c(staAGI, endAGI)
  x2 <- x[x>=p2[1] & x<=p2[2]]
  W2 <- matrix(as.numeric(x%in%x2), m, n)
  Z2 <- matrix(Z[which(W2==1)], ncol=n)
  
  ## constructing the bases
  
  ## BASIS for the infant component
  ## B-splines over ages
  degx1 <- floor(length(x1)/3)
  xl1 <- min(x)
  xr1 <- max(x)
  xmax1 <- xr1 + 0.01 * (xr1 - xl1)
  xmin1 <- xl1 - 0.01 * (xr1 - xl1)
  Bx1 <- MortSmooth_bbase(x1, xmin1, xmax1, degx1, 3)
  nbx1 <- ncol(Bx1)
  ## B-splines over years
  degy1 <- floor(n/3)
  yl1 <- min(y)
  yr1 <- max(y)
  ymax1 <- yr1 + 0.01 * (yr1 - yl1)
  ymin1 <- yl1 - 0.01 * (yr1 - yl1)
  By1 <- MortSmooth_bbase(y, ymin1, ymax1, degy1, 3)
  nby1 <- ncol(By1)
  ## final basis for the infant component
  B1 <- kronecker(By1, Bx1)
  
  ## BASIS for the adult component
  ## B-splines over ages
  degx2 <- floor(length(x2)/3)
  xl2 <- min(x2)
  xr2 <- max(x2)
  xmax2 <- xr2 + 0.01 * (xr2 - xl2)
  xmin2 <- xl2 - 0.01 * (xr2 - xl2)
  Bx2 <- MortSmooth_bbase(x2, xmin2, xmax2, degx2, 3)
  nbx2 <- ncol(Bx2)
  ## B-splines over years
  degy2 <- floor(n/3)
  yl2 <- min(y)
  yr2 <- max(y)
  ymax2 <- yr2 + 0.01 * (yr2 - yl2)
  ymin2 <- yl2 - 0.01 * (yr2 - yl2)
  By2 <- MortSmooth_bbase(y, ymin2, ymax2, degy2, 3)
  nby2 <- ncol(By2)
  ## final basis for the aging component
  B2 <- kronecker(By2, Bx2)
  
  ## complete model matrix as a list
  ## in their order and with dimensions
  XX <- list(X1=B1, X2=B2)
  nx <- length(XX)
  nc <- unlist(lapply(XX, ncol))
  ## indicators for the coefficients of each basis
  ind2 <- cumsum(nc)
  ind1 <- c(1, ind2[1:(nx-1)]+1)
  ind <- NULL
  for(i in 1:nx){
    ind[[i]] <- ind1[i]:ind2[i]
  }
  ## indicator for the fitted values
  ## rows
  wr1 <- rep(W1[,1], n)
  wr2 <- rep(W2[,1], n)
  indR <- list(wr1=wr1, wr2=wr2)
  ## columns
  wc1 <- rep(1, nc[1])
  wc2 <- rep(1, nc[2])
  indC <- list(wc1=wc1, wc2=wc2)
  
  ## country-specific starting values
  ## for the infant component
  if (cou == "IND"| cou == "USA"){
    ZZ1 <- Z1
    ZZ1[x1>=35] <- 0
    WW1 <- matrix(0,length(x1),n)
    WW1[x1<=35,] <- 1
    fit1 <- Mort2Dsmooth(x=x1, y=y, Z=ZZ1,
                         # offset=log(E1),
                         W=WW1,
                         method=3,
                         lambdas=c(10^1, 10^2),
                         ndx = c(degx1, degy1),
                         deg = c(3,3), pord = c(2,2))
    ## starting coefficients for the childhood component
    coef1.st <- c(fit1$coef)
    ## starting fitted values for the childhood component
    z1.st <- exp(XX[[1]] %*% coef1.st) # *c(E1) ## in vector
    Z1.st <- matrix(z1.st, length(x1), n) ## in matrix
    ## in matrix over all ages
    Z1.stA <- matrix(0, m, n) 
    Z1.stA[W1==1] <- z1.st
  }else if (cou == "AGO" | cou == "GTM" | cou == "BFA" | cou == "TCD" | 
            cou == "SEN" | cou == "MLI" | cou == "GHA"){
    ZZ1 <- Z1
    ZZ1[x1>=45] <- 0
    WW1 <- matrix(0,length(x1),n)
    WW1[x1<=45,] <- 1
    fit1 <- Mort2Dsmooth(x=x1, y=y, Z=ZZ1,
                         # offset=log(E1),
                         W=WW1,
                         method=3,
                         lambdas=c(10^1, 10^1),
                         ndx = c(degx1, degy1),
                         deg = c(3,3), pord = c(2,2))
    ## starting coefficients for the childhood component
    coef1.st <- c(fit1$coef)
    ## starting fitted values for the childhood component
    z1.st <- exp(XX[[1]] %*% coef1.st) # *c(E1) ## in vector
    Z1.st <- matrix(z1.st, length(x1), n) ## in matrix
    ## in matrix over all ages
    Z1.stA <- matrix(0, m, n) 
    Z1.stA[W1==1] <- z1.st
  }else{
    ZZ1 <- Z1
    ZZ1[x1>=40] <- 0
    WW1 <- matrix(0,length(x1),n)
    WW1[x1<=40,] <- 1
    fit1 <- Mort2Dsmooth(x=x1, y=y, Z=ZZ1,
                         # offset=log(E1),
                         W=WW1,
                         method=3,
                         lambdas=c(10^1, 10^2),
                         ndx = c(degx1, degy1),
                         deg = c(3,3), pord = c(2,2))
    ## starting coefficients for the childhood component
    coef1.st <- c(fit1$coef)
    ## starting fitted values for the childhood component
    z1.st <- exp(XX[[1]] %*% coef1.st) # *c(E1) ## in vector
    Z1.st <- matrix(z1.st, length(x1), n) ## in matrix
    ## in matrix over all ages
    Z1.stA <- matrix(0, m, n) 
    Z1.stA[W1==1] <- z1.st
  }
  
  ## for the adult component
  if (cou=="ZWE" | cou == "IND" | cou == "USA" | cou == "CHN" | cou == "MMR"){
    ZZ2 <- Z2
    ZZ2[x2<=40] <- 0
    WW2 <- matrix(0,length(x2),n)
    WW2[x2>=40] <- 1
    fit2 <- Mort2Dsmooth(x=x2, y=y, Z=ZZ2, 
                         W=WW2,
                         method=3, lambdas=c(10^1, 10^2),
                         ndx = c(degx2, degy2),
                         deg = c(3,3), pord = c(2,2))
    ## starting coefficients for the senescence component
    coef2.st <- c(fit2$coef)
    ## starting fitted values for the senescence component
    z2.st <- exp(XX[[2]] %*% coef2.st) #*e ## in vector
    Z2.st <- matrix(z2.st, length(x2), n) ## in matrix
    ## in matrix over all ages
    Z2.stA <- matrix(0, m, n) 
    Z2.stA[W2==1] <- z2.st
  }else if (cou=="AGO"){
    ZZ2 <- Z2
    ZZ2[x2<=55] <- 0
    WW2 <- matrix(0,length(x2),n)
    WW2[x2>=55] <- 1
    fit2 <- Mort2Dsmooth(x=x2, y=y, Z=ZZ2, 
                         W=WW2,
                         method=3, lambdas=c(10^1, 10^1),
                         ndx = c(degx2, degy2),
                         deg = c(3,3), pord = c(2,2))
    ## starting coefficients for the senescence component
    coef2.st <- c(fit2$coef)
    ## starting fitted values for the senescence component
    z2.st <- exp(XX[[2]] %*% coef2.st) #*e ## in vector
    Z2.st <- matrix(z2.st, length(x2), n) ## in matrix
    ## in matrix over all ages
    Z2.stA <- matrix(0, m, n) 
    Z2.stA[W2==1] <- z2.st
  }else{
    ZZ2 <- Z2
    ZZ2[x2<=50] <- 0
    WW2 <- matrix(0,length(x2),n)
    WW2[x2>=50] <- 1
    fit2 <- Mort2Dsmooth(x=x2, y=y, Z=ZZ2, 
                         W=WW2,
                         method=3, lambdas=c(10^1, 10^2),
                         ndx = c(degx2, degy2),
                         deg = c(3,3), pord = c(2,2))
    ## starting coefficients for the senescence component
    coef2.st <- c(fit2$coef)
    ## starting fitted values for the senescence component
    z2.st <- exp(XX[[2]] %*% coef2.st) #*e ## in vector
    Z2.st <- matrix(z2.st, length(x2), n) ## in matrix
    ## in matrix over all ages
    Z2.stA <- matrix(0, m, n) 
    Z2.stA[W2==1] <- z2.st
  }
  
  ## all deaths
  Z.st <- Z1.stA + Z2.stA 
  
  if (any(Z.st<=0)){
    cat("WARNING - starting values < 0", "\n")
    break
  } 
  
  ## concatenating starting coefficients
  coef.st <- as.vector(c(coef1.st,coef2.st))
  
  ## shape penalties
  
  ## penalty terms for the infant component
  ## including log-concaveness
  Dcon1 <- kronecker(diag(nby1), diff(diag(nbx1),
                                      diff=2))
  wcon1 <- rep(0, (nbx1-2)*nby1)
  Wcon1 <- diag(wcon1)
  ## smooth penalty for the age
  Dx1 <- diff(diag(nbx1), diff=3)
  tDDx1 <- t(Dx1) %*% Dx1
  Px1 <- kronecker(diag(nby1), tDDx1)
  ## smooth penalty for the year
  Dy1 <- diff(diag(nby1), diff=2)
  tDDy1 <- t(Dy1) %*% Dy1
  Py1 <- kronecker(tDDy1, diag(nbx1))
  
  ## penalty terms for the adult component
  ## including log-concaveness
  Dcon2 <- kronecker(diag(nby2), diff(diag(nbx2),
                                      diff=2))
  wcon2 <- rep(0, (nbx2-2)*nby2)
  Wcon2 <- diag(wcon2)
  ## smooth penalty for the age
  Dx2 <- diff(diag(nbx2), diff=3)
  tDDx2 <- t(Dx2) %*% Dx2
  Px2 <- kronecker(diag(nby2), tDDx2)
  ## smooth penalty for the year
  Dy2 <- diff(diag(nby2), diff=2)
  tDDy2 <- t(Dy2) %*% Dy2
  Py2 <- kronecker(tDDy2, diag(nbx2))
  
  ## SSE objective function 
  SSEestimationFUN <- function(lambdas, coef.st){
    lambda.x1 <- lambdas[1]
    lambda.x2 <- lambdas[2]
    lambda.y  <- lambdas[3]
    cat(lambdas, "\n")
    ## smoothing each component
    Pxy1 <- lambda.x1 * Px1 + lambda.y * Py1
    Pxy2 <- lambda.x2 * Px2 + lambda.y * Py2
    ## final penalty for smoothness
    P <- bdiag.spam(Pxy1, Pxy2)
    coef <- coef.st
    conv <- TRUE
    ## iteration
    for(it in 1:max.it){
      ## penalty for the concaveness for child
      Pcon10 <- t(Dcon1) %*% Wcon1 %*% Dcon1
      Pcon1 <- kappa * Pcon10
      ## penalty for the concaveness for aging
      Pcon20 <- t(Dcon2) %*% Wcon2 %*% Dcon2
      Pcon2 <- kappa * Pcon20
      ## final shape penalty term
      Psha <- bdiag.spam(Pcon1, Pcon2)
      ## linear predictor
      eta <- numeric(nx*mn)
      for(i in 1:nx){
        eta0 <- rep(0, mn)
        eta0[which(indR[[i]]==1)] <- XX[[i]]%*%coef[ind[[i]]]
        eta[1:mn+(i-1)*mn] <- eta0
      }
      ## components
      gamma <- exp(eta) * c(unlist(indR))
      ## expected values
      mu <- numeric(mn)
      for(i in 1:nx){
        mu <- (1 * gamma[1:mn+(i-1)*mn]) + mu
      }    
      ## weights for the IWLS
      w <- mu
      ## modified model matrix for a CLM
      U <- matrix(NA, mn, sum(nc))
      for(i in 1:nx){
        u <- gamma[1:mn+(i-1)*mn]/mu * 1
        u[wei==0] <- 0
        XXi <- matrix(0, nrow=mn, ncol=nc[i])
        XXi[which(indR[[i]]==1),which(indC[[i]]==1)]<-XX[[i]]
        U0 <- u * XXi
        U[,ind[[i]]] <- U0
      }
      ## regression parts for the P-CLM
      tUWU <- t(U) %*% (w * U)
      tUWUpP <- tUWU + P + Psha + Pr
      r <- z - mu
      tUr <- t(U) %*% r
      ## updating coefficients with a d-step 
      coef.old <- coef
      coefTRY <- try(solve(tUWUpP, tUr + tUWU%*%coef),silent = T)
      if(class(coefTRY)[1]=="try-error" | any(is.na(coefTRY)) | all(coefTRY==0)){
        conv <- FALSE
        break
      }else{
        coef <- coefTRY
      }
      coef <- d*coef.old + (1-d)*coef
      ## update weights for shape constraints
      ## infant, log-concaveness
      Wcon1.old <- Wcon1
      wcon1 <- rep(0, nrow(Dcon1))
      COEF <- matrix(coef[ind[[1]]], nbx1, nby1)
      diff.COEF <- apply(COEF, 2, diff, diff=2)>=0
      wcon1[c(diff.COEF)] <- 1
      Wcon1 <- diag(wcon1)
      ## aging, log-concaveness
      Wcon2.old <- Wcon2
      wcon2 <- rep(0, nrow(Dcon2))
      COEF <- matrix(coef[ind[[2]]], nbx2, nby2)
      diff.COEF <- apply(COEF, 2, diff, diff=2)>=0
      wcon2[c(diff.COEF)] <- 1
      Wcon2 <- diag(wcon2)  
      ## convergence criterion for the coefficients
      dif.coef<-max(abs(coef.old-coef))/max(abs(coef))
      ## stopping loop at convergence
      if(dif.coef < 1e-04 & it > 4) break
    }
    if (conv){
    ## compute devaince
    zz <- z
    zz[z==0] <- 10^-8
    mumu <- mu
    mumu[mu==0] <- 10^-8
    dev <- 2*sum(z * log(zz/mumu))
    ## effective dimensions
    H <- solve(tUWUpP, tUWU)
    diagH <- diag(H)
    ed <- sum(diagH)
    ## BIC
    bic <- dev + log(mn)*ed
    ## monitoring part
    cat(it, signif(dif.coef,8), bic, "\n\n")
    ## returning list with various objects
    out <- list(bic=bic, coef=coef, it=it,
                dif.coef=dif.coef, diagH=diagH,conv=conv)
    }else{
      out <- list(bic=NA, coef=NA, it=NA,
                  dif.coef=NA, diagH=NA,conv=conv)
    }
    return(out)
  }
  
  ## shape parameter
  if (is.null(kappa)) kappa <- 10^6
  ## solving step
  d <- 0.5
  ## small ridge penalty for numerical stability
  Pr <- 10^-4 * diag.spam(length(coef.st))
  ## maximum number of iteration
  if (is.null(max.it)) max.it <- 250
  ## smoothing parameters
  if (is.null(lambdas)) lambdas <- c(10^1.5,10^1.5,10^1.5)
  ## start time
  start_time <- Sys.time()
  ## estimating the SSE
  fitFIN <- SSEestimationFUN(lambdas, coef.st)
  if (!fitFIN$conv){
    ## if not converged, use higher smoothing parameters
    lambdas <- c(10^2,10^2,10^2)
    fitFIN <- SSEestimationFUN(lambdas, coef.st)
  }
  ## end time
  end_time <- Sys.time()
  print(end_time-start_time)
  
  ## estimated coefficients
  if (fitFIN$conv){
    coef.hat <- fitFIN$coef
    
    ## linear predictor for each component
    etas <- NULL
    for(i in 1:nx){
      etas[[i]] <- XX[[i]] %*% coef.hat[ind[[i]]]
    }
    eta1.hat <- etas[[1]]
    eta2.hat <- etas[[2]]
    
    ## linear predictor in a matrix over the whole x
    ETA.hat <- matrix(NA, mn, nx)
    for(i in 1:nx){
      ETA.hat[which(indR[[i]]==1),i] <- etas[[i]]
    }
    ## linear predictor for overall mortality
    eta.hat <- log(apply(exp(ETA.hat), 1, sum, na.rm=TRUE))
    
    ## fitted values for each component
    gamma1.hat <- exp(eta1.hat)
    gamma2.hat <- exp(eta2.hat)
    gamma.hat <- exp(eta.hat) * c(unlist(indR))
    
    GAMMA1.hat <- matrix(gamma1.hat,length(x1),n)
    GAMMA2.hat <- matrix(gamma2.hat,length(x2),n)
    
    ## expected values
    # mu.hat <- exp(eta.hat)*e
    mu.hat <- exp(eta.hat)
    MU.hat <- matrix(mu.hat,m,n)
    
    ## deviance residuals
    res1 <- sign(z - mu.hat)
    res2 <- sqrt(2 * (z * log(z/mu.hat) - z + mu.hat))
    res <- res1 * res2 ## in vector
    RES <- matrix(res, m, n) ## in matrix
  }
  ## output
  out <- list(MU.hat=MU.hat,GAMMA1.hat=GAMMA1.hat,GAMMA2.hat=GAMMA2.hat,
              RES=RES,coef.hat=coef.hat,x1=x1,x2=x2)
  return(out)
  
}

## function to compute summary measures from SSE
SSEsummaryFUN <- function(x,cohorts,Z,MU.hat,x1,x2,GAMMA1.hat,GAMMA2.hat){
  
  ## dimensions
  m <- length(x)
  n <- length(cohorts)
  
  ## weights of distributions
  w1 <- apply(GAMMA1.hat,2,sum)/apply(MU.hat,2,sum)
  w2 <- apply(GAMMA2.hat,2,sum)/apply(MU.hat,2,sum)
  
  ## mean and sd of distibutions
  FX <- FX.hat <- matrix(NA,m,n)
  FX1.hat <- matrix(NA,length(x1),n)
  FX2.hat <- matrix(NA,length(x2),n)
  for (jj in 1:n){
    FX[,jj] <- Z[,jj]/sum(Z[,jj])
    FX.hat[,jj] <- MU.hat[,jj]/sum(MU.hat[,jj])
    FX1.hat[,jj] <- GAMMA1.hat[,jj]/sum(MU.hat[,jj])
    FX2.hat[,jj] <- GAMMA2.hat[,jj]/sum(MU.hat[,jj])
  }
  apply(FX.hat,2,sum)
  apply(FX1.hat,2,sum)+apply(FX2.hat,2,sum)
  
  ## mean and sd
  ax <- 0.5
  M1 <- M2 <- SD1 <- SD2 <- rep(NA,n)
  for (jj in 1:n){
    M1[jj] <- sum((fitSSE$x1+ax)*FX1.hat[,jj])/sum(FX1.hat[,jj])
    M2[jj] <- sum((fitSSE$x2+ax)*FX2.hat[,jj])/sum(FX2.hat[,jj])
    v1 <- sum(FX1.hat[,jj] * ((fitSSE$x1+ax) - M1[jj])^2) / sum(FX1.hat[,jj])
    v2 <- sum(FX2.hat[,jj] * ((fitSSE$x2+ax) - M2[jj])^2) / sum(FX2.hat[,jj])
    SD1[jj] <- sqrt(v1)
    SD2[jj] <- sqrt(v2)
  }
  
  ## output
  out <- list(M1=M1,M2=M2,SD1=SD1,SD2=SD2,
              w1=w1)
  return(out)
  
}

