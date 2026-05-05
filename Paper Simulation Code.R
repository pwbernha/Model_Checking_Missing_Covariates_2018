#####################################################################################
#NOTE: PLEASE DO NOT USE CODE FOR ANY PUBLICATION PURPOSES WITHOUT FIRST CONTACTING #
#      PAUL BERNHARDT AT PAUL.BERNHARDT@VILLANOVA.EDU                               #
#####################################################################################

####################Necessary Packages for Simulation######################
library(mcmc)
library(MCMCpack)
library(msm)
library(mvtnorm)
library(survival)
library(cubature)
library(mice)


###################### Simulation Functions #######################

#Generating Data...somewhat lazily here requires defining parameter values
#within the function rather than upfront
DataGen <- function(N,n){
  Ydata <- Xdata <- n0 <- n1 <- n2 <- n12 <- list()
  BetaX <- c(0.5,-0.01,0.1,-10)
  sigmaX <- matrix(c(225,20,20,225),2,2)
  Beta <- c(-5, 0.03, 0.0005, 0.02, -1)
  ParMiss <- matrix(c(-2,0.05,0.1,0.1,-3,0.05,-0.1,0.1),2,4,byrow=TRUE)
  sigma <- 1
  for(k in 1:N){
    #Generating Covariates
    p <- rbinom(n,1,0.55)
    Age <- p*rtnorm(n,18,40,lower=0,upper=42)+(1-p)*rtnorm(n,53,16,lower=42)
    Gender <- rbinom(n,1,0.488-0.001477778*Age+0.00004259259*Age^2)
    
    #x covariates represent serum LDL and systolic blood pressure (relationships based on Japanese study)
    x <- matrix(0,n,2)
    for(i in 1:n) x[i,] <- rmvnorm(1,c(100+BetaX[1]*Age[i]+BetaX[2]*Gender[i],115+BetaX[3]*Age[i]+BetaX[4]*Gender[i]),sigmaX)			
    Xgen <- cbind(x,Age,Gender)
    
    #Response representing the exponential fraction of an artery that is blocked
    #Ygen <- rnorm(n, Beta[1] + Beta[2]*Xgen[,1]+Beta[3]*Xgen[,2]+Beta[4]*Xgen[,3]+Beta[5]*Xgen[,4], sigma)	
    Ygen <- Beta[1] + Beta[2]*Xgen[,1]+Beta[3]*Xgen[,2]+Beta[4]*Xgen[,3]+Beta[5]*Xgen[,4]+ sqrt(3)*rgamma(n, 3, 3)
    
    #Generating Missing Values for LDL and BP
    r1 <- rbinom(n, 1, pnorm(ParMiss[1,1]+ParMiss[1,2]*Xgen[,3]+ParMiss[1,3]*Xgen[,4]+ParMiss[1,4]*Ygen))
    r2 <- rbinom(n, 1, pnorm(ParMiss[2,1]+ParMiss[2,2]*Xgen[,3]+ParMiss[2,3]*Xgen[,4]+ParMiss[2,4]*Ygen))
    
    Y <- Ygen
    X <- Xgen
    
    X[r1==1,1] <- Inf
    X[r2==1,2] <- Inf
    
    #The following orders the data based on location of missing
    Yobs <- Y[X[,1]!=Inf & X[,2]!=Inf]
    Y1miss <- Y[X[,1]==Inf & X[,2]!=Inf]
    Y2miss <- Y[X[,1]!=Inf & X[,2]==Inf]
    Y12miss <- Y[X[,1]==Inf & X[,2]==Inf]
    Ydata[[k]] <- c(Yobs,Y1miss,Y2miss,Y12miss)
    Xobs <- X[X[,1]!=Inf & X[,2]!=Inf,]
    X1miss <- X[X[,1]==Inf & X[,2]!=Inf,]
    X2miss <- X[X[,1]!=Inf & X[,2]==Inf,]
    X12miss <- X[X[,1]==Inf & X[,2]==Inf,]
    Xdata[[k]] <- rbind(Xobs,X1miss,X2miss,X12miss)
    
    #Number of complete, missing 1, missing 2 or missing both 
    n0[[k]] <- length(Yobs)
    n1[[k]] <- length(Y1miss)
    n2[[k]] <- length(Y2miss)
    n12[[k]] <- length(Y12miss)
  }
  return(list(Ydata,Xdata,n0,n1,n2,n12))
}

#likelihood functions used to obtained initial maximum likelihood estimates
like1 <- function(x,X,Y,par) (dnorm(Y,par[1]+par[2]*x+X[2:4]%*%par[3:5],par[6])*dnorm(x, par[7]+ X[2:4]%*%par[8:10],par[11])*dnorm(X[2], par[12]+X[3:4]%*%par[13:14],par[15]))
like2 <- function(x,X,Y,par) (dnorm(Y,par[1]+par[2]*X[1]+par[3]*x+X[3:4]%*%par[4:5],par[6])*dnorm(X[1], par[7]+ x*par[8] +X[3:4]%*%par[9:10],par[11])*dnorm(x, par[12]+X[3:4]%*%par[13:14],par[15]))
like12 <- function(x,X,Y,par) (dnorm(Y,par[1]+par[2]*x[1]+par[3]*x[2]+X[3:4]%*%par[4:5],par[6])*dnorm(x[1], par[7]+ x[2]*par[8]+X[3:4]%*%par[9:10],par[11])*dnorm(x[2], par[12]+X[3:4]%*%par[13:14],par[15]))

#Maximization function for purpose of initial estimates, not needed unless numerical issues arose
MaxMissing <- function(X,Y,n0,n1,n2,n12){
  like <- function(par) {
    fc<-rep(0,n)
    fc[1:n0] <- -log(dnorm(Y[1:n0],par[1]+par[2]*X[1:n0,1]+X[1:n0,2:4]%*%par[3:5],par[6])*dnorm(X[1:n0,1], par[7]+ X[1:n0,2:4]%*%par[8:10],par[11])*dnorm(X[1:n0,2], par[12]+X[1:n0,3:4]%*%par[13:14],par[15]))
    for(i in (n0+1):(n0+n1)) fc[i] <- -log(integrate(like1,lower=0, upper=200,X=X[i,],Y=Y[i], par=par)$value)
    for(i in (n0+n1+1):(n0+n1+n2)) fc[i] <- -log(integrate(like2,lower=0, upper=200,X=X[i,],Y=Y[i], par=par)$value)
    for(i in (n0+n1+n2+1):n) fc[i] <- -log(adaptIntegrate(like12,lower=c(0,0), upper=c(400,400),maxEval=500, X=X[i,],Y=Y[i], par=par)$integral)
    sum(fc)
  }
  like
}

#Maximization function for purpose of initial estimates based on complete cases only 
MaxMissingSimp <- function(X,Y,n0){
  like <- function(par) {
    fc<-rep(0,n0)
    fc[1:n0] <- -log(dnorm(Y[1:n0],par[1]+par[2]*X[1:n0,1]+X[1:n0,2:4]%*%par[3:5],par[6])*dnorm(X[1:n0,1], par[7]+ X[1:n0,2:4]%*%par[8:10],par[11])*dnorm(X[1:n0,2], par[12]+X[1:n0,3:4]%*%par[13:14],par[15]))
    sum(fc)
  }
  like
}

#log likelihoods for Gibbs sampler with missing 1, missing 2 or missing both
lPost1 <- function(x,X,Y,pars) log(dnorm(Y,pars[1]+pars[2]*x+X[2:4]%*%pars[3:5],pars[6])*dnorm(x, pars[7]+ X[2:4]%*%pars[8:10],pars[11])*dnorm(X[2], pars[12]+X[3:4]%*%pars[13:14],pars[15]))
lPost2 <- function(x,X,Y,pars) log(dnorm(Y,pars[1]+pars[2]*X[1]+pars[3]*x+X[3:4]%*%pars[4:5],pars[6])*dnorm(X[1], pars[7]+ x*pars[8] +X[3:4]%*%pars[9:10],pars[11])*dnorm(x, pars[12]+X[3:4]%*%pars[13:14],pars[15]))
lPost12 <- function(x,X,Y,pars) log(dnorm(Y,pars[1]+pars[2]*x[1]+pars[3]*x[2]+X[3:4]%*%pars[4:5],pars[6])*dnorm(x[1], pars[7]+ x[2]*pars[8]+X[3:4]%*%pars[9:10],pars[11])*dnorm(x[2], pars[12]+X[3:4]%*%pars[13:14],pars[15]))

#Bayesian method for obtaining parameter estimates and imputations
BayesMethods <- function(Y,X,NumImp,Inits,n0,n1,n2,n12){
  #Parameter storage
  keep.Ximp <- list()
  keep.Beta <- matrix(0,NumImp, 5)
  keep.BetaX1 <- matrix(0,NumImp, 4)
  keep.BetaX2 <- matrix(0,NumImp, 3)
  keep.sigma <- rep(0,NumImp)
  keep.sigmaX1 <- rep(0,NumImp)
  keep.sigmaX2 <- rep(0,NumImp)
  keep.Mean <- matrix(0,NumImp,5)
  
  #Initial values/defs:
  BetaI <- Inits[1:5]
  sigmaI <- Inits[6]
  BetaX1I <- Inits[7:10]
  sigmaX1I <- Inits[11]
  BetaX2I <- Inits[12:14]
  sigmaX2I <- Inits[15]
  MeanI <- Inits[1:5]
  
  
  Ximp <- X
  for(j in 1:NumImp){
    if(j%%50==0) print(j)
    
    #####Drawing missing X's#####
    invisible(capture.output(if(j==1){
      if(n1>0){
        for(k in 1:n1){tryCatch(Ximp[(n0+k),1] <- MCMCmetrop1R(lPost1,theta.init=(mean(X[1:n0,1])+rnorm(1,0,sigmaX1I)),burnin=10, mcmc=1, thin=1,X=X[n0+k,],Y=Y[n0+k], pars=c(BetaI,sigmaI,BetaX1I,sigmaX1I,BetaX2I,sigmaX2I)), error=function(...) Ximp[n0+k,1] <<- mean(X[1:n0,1]))}
      }
      
      if(n2>0){
        for(k in 1:n2){
          tryCatch(Ximp[(n0+n1+k),2] <- MCMCmetrop1R(lPost2,theta.init=(mean(X[1:n0,2])+rnorm(1,0,sigmaX2I)),burnin=10, mcmc=1, thin=1, X=X[n0+n1+k,], Y=Y[n0+n1+k], pars=c(BetaI,sigmaI,BetaX1I,sigmaX1I,BetaX2I,sigmaX2I)), error=function(...) Ximp[n0+n1+k,2] <<- mean(X[1:n0,2]))
        }
      }
      
      if(n12>0){
        for(k in 1:n12){
          tryCatch(Ximp[(n0+n1+n2+k),1:2] <- MCMCmetrop1R(lPost12,theta.init=(colMeans(X[1:n0,1:2])+rmvnorm(1,c(0,0),diag(c(sigmaX1I^2,sigmaX2I^2),2))),burnin=300, mcmc=1, thin=1, X=X[n0+n1+n2+k,], Y=Y[n0+n1+n2+k], pars=c(BetaI,sigmaI,BetaX1I,sigmaX1I,BetaX2I,sigmaX2I)), error=function(...) Ximp[n0+n1+n2+k,1:2] <<- colMeans(X[1:n0,1:2]))
        }
      }
    }))
    
    invisible(capture.output(if(j>1){
      if(n1>0){
        for(k in 1:n1){
          Ximp[(n0+k),1] <- MCMCmetrop1R(lPost1,theta.init=(Ximp[n0+k,1]),burnin=10, mcmc=1, thin=1, X=X[n0+k,],Y=Y[n0+k], pars=c(BetaI,sigmaI,BetaX1I,sigmaX1I,BetaX2I,sigmaX2I))
        }
      }
      
      if(n2>0){
        for(k in 1:n2){
          Ximp[(n0+n1+k),2] <- MCMCmetrop1R(lPost2,theta.init=(Ximp[n0+n1+k,2]),burnin=10, mcmc=1, thin=1, X=X[n0+n1+k,], Y=Y[n0+n1+k], pars=c(BetaI,sigmaI,BetaX1I,sigmaX1I,BetaX2I,sigmaX2I))
        }
      }	
      
      if(n12>0){
        for(k in 1:n12){
          Ximp[(n0+n1+n2+k),1:2] <- MCMCmetrop1R(lPost12,theta.init=(Ximp[n0+n1+n2+k,1:2]),burnin=300, mcmc=1, thin=1, X=X[n0+n1+n2+k,], Y=Y[n0+n1+n2+k], pars=c(BetaI,sigmaI,BetaX1I,sigmaX1I,BetaX2I,sigmaX2I))
        }
      }	
    }))
    
    keep.Ximp[[j]] <- Ximp
    
    #Update Beta (Prior, conditional on sigmaI, is N(MeanI, 100*I))
    BetaI <- as.vector(rmvnorm(1,solve(t(cbind(1,Ximp))%*%cbind(1,Ximp)+(1/100)*sigmaI^2*diag(5))%*%(t(cbind(1,Ximp))%*%Y+(1/100)*sigmaI^2*diag(5)%*%MeanI),solve(t(cbind(1,Ximp))%*%cbind(1,Ximp)+(1/100)*sigmaI^2*diag(5))))
    keep.Beta[j,]<-BetaI
    
    #Update MeanI (Prior is N(0,25*I))
    MeanI <- as.vector(rmvnorm(1, (1/5)*BetaI, 20*diag(5)))
    keep.Mean[j,]<-MeanI
    
    #Update sigma (Prior for sigmaI^2 is IG(1,0.1))
    sigmaI <-sqrt(rinvgamma(1,1+n/2,0.1+0.5*(t(Y)%*%Y+t(MeanI)%*%(1/100*diag(5))%*%MeanI-t(solve(t(cbind(1,Ximp))%*%cbind(1,Ximp)+(1/100)*diag(5))%*%(t(cbind(1,Ximp))%*%Y+(1/100)*diag(5)%*%MeanI))%*%(t(cbind(1,Ximp))%*%cbind(1,Ximp)+(1/100)*diag(5))%*%solve(t(cbind(1,Ximp))%*%cbind(1,Ximp)+(1/100)*diag(5))%*%(t(cbind(1,Ximp))%*%Y+(1/100)*diag(5)%*%MeanI))))
    keep.sigma[j] <- sigmaI
    
    #Update BetaX1 (prior is N(0,100*I)
    BetaX1I <- BetaUps <-  rmvnorm(1,solve(t(cbind(1,Ximp[,2:4]))%*%cbind(1,Ximp[,2:4])/sigmaX1I^2+1/100*diag(4))%*%(t(cbind(1,Ximp[,2:4]))%*%Ximp[,1]/sigmaX1I^2),solve(t(cbind(1,Ximp[,2:4]))%*%cbind(1,Ximp[,2:4])/sigmaX1I^2+1/100*diag(4)))
    keep.BetaX1[j,]<-BetaX1I
    
    #Update BetaX2I (prior is N(0,100*I)
    BetaX2I  <- rmvnorm(1,solve(t(cbind(1,Ximp[,3:4]))%*%cbind(1,Ximp[,3:4])/sigmaX2I^2+1/100*diag(3))%*%(t(cbind(1,Ximp[,3:4]))%*%Ximp[,2]/sigmaX2I^2),solve(t(cbind(1,Ximp[,3:4]))%*%cbind(1,Ximp[,3:4])/sigmaX2I^2+1/100*diag(3)))
    keep.BetaX2[j,]<-BetaX2I
    
    #Update sigmaX1 (prior is IG(1,0.1)
    sigmaX1I <- sqrt(rinvgamma(1,1+n/2,0.1+0.5*(t(Ximp[,1])%*%Ximp[,1]-t(solve(t(cbind(1,Ximp[,2:4]))%*%cbind(1,Ximp[,2:4])+(1/100)*diag(4))%*%(t(cbind(1,Ximp[,2:4]))%*%Ximp[,1]))%*%(t(cbind(1,Ximp[,2:4]))%*%cbind(1,Ximp[,2:4])+(1/100)*diag(4))%*%solve(t(cbind(1,Ximp[,2:4]))%*%cbind(1,Ximp[,2:4])+(1/100)*diag(4))%*%(t(cbind(1,Ximp[,2:4]))%*%Ximp[,1]))))
    keep.sigmaX1[j] <- sigmaX1I
    
    #Update sigmaX2 (prior is IG(1,0.1)
    sigmaX2I <-sqrt(rinvgamma(1,1+n/2,0.1+0.5*(t(Ximp[,2])%*%Ximp[,2]-t(solve(t(cbind(1,Ximp[,3:4]))%*%cbind(1,Ximp[,3:4])+(1/100)*diag(3))%*%(t(cbind(1,Ximp[,3:4]))%*%Ximp[,2]))%*%(t(cbind(1,Ximp[,3:4]))%*%cbind(1,Ximp[,3:4])+(1/100)*diag(3))%*%solve(t(cbind(1,Ximp[,3:4]))%*%cbind(1,Ximp[,3:4])+(1/100)*diag(3))%*%(t(cbind(1,Ximp[,3:4]))%*%Ximp[,2]))))
    keep.sigmaX2[j] <- sigmaX2I
  }
  
  return(list(keep.Ximp,keep.Beta,keep.sigma,keep.BetaX1,keep.BetaX2,keep.sigmaX1,keep.sigmaX2, keep.Mean))
}

#Function for imputating X given parameters; not needed in this simulation
Ximputation <- function(Y,X,pars,B){
  
  Imps <- list()
  for(j in 1:B){
    Imps[[j]] <- X
  }
  
  for(k in 1:n){
    SEED <- sample(1:1000000,1) #problem with random number generator for function?
    #vals <- MCMC(lPost12, 50*B, init=X[k,1:2], scale = matrix(c(225,28,28,225),2,2), X=X[k,], Y=Y[k], pars=pars)$samples[seq(50,50*B,50),] #alternative method
    invisible(capture.output(vals <- MCMCmetrop1R(lPost12,theta.init=(X[k,1:2]+rmvnorm(1,c(0,0),matrix(c(225,28,28,225),2,2))),burnin=100, mcmc=100*B, thin=100, V=matrix(c(225,28,28,225),2,2), seed=SEED, X=X[k,], Y=Y[k], pars=pars)))
    for(j in 1:B){
      Imps[[j]][k,1:2] <- vals[j,]
    }
  }
  return(Imps)	
}


#Finding Komogorov-Smirnov p-values using bootstrapping based on Bayesian model imputations
KSstat1 <- function(Y, Ximps, BetaPost, SigmaPost, NumImp,B){
  pvals <- rep(0,NumImp)
  
  for(i in 1:NumImp){
    Fit <- lm(Y ~ Ximps[[i]])
    sigmaobs <- summary(Fit)[[6]]
    KSobs <- ks.test(Fit$residuals,"pnorm",0,sigmaobs)$statistic
    
    KSsim <-  rep(0,B)
    for(j in 1:B){
      Ynew <- rnorm(n, cbind(1,Ximps[[i]])%*%Fit$coefficients, sigmaobs)	
      FIT <- lm(Ynew ~ Ximps[[i]])
      sigmasim <- summary(FIT)[[6]]
      KSsim[j] <- ks.test(FIT$residuals,"pnorm",0,sigmasim)$statistic
    }
    pvals[i] <- sum(KSobs<KSsim)/B
  }
  
  return(pvals)
}

#Finding Komogorov-Smirnov p-values using posterior predictive checking based on Bayesian model imputations
KSstat2 <- function(Y, Ximps, pars, NumImp, B){
  pvals <- rep(0,NumImp)
  for(j in 1:NumImp){
    Fit <- lm(Y ~ Ximps[[j]])
    sigmaobs <- summary(Fit)[[6]]
    KSobs <- ks.test(Fit$residuals,"pnorm",0,sigmaobs)$statistic
    
    KSsim <- rep(0,B)
    for(k in 1:B){
      Ynew <- rnorm(n, cbind(1,Ximps[[j]])%*%pars[j,1:5], pars[j,6])	
      Fit <- lm(Ynew ~ Ximps[[j]])
      sigmasim <- summary(Fit)[[6]]
      KSsim[k] <- ks.test(Fit$residuals,"pnorm",0,sigmasim)$statistic
    }
    
    pvals[j] <- sum(KSobs<KSsim)/B
    
  }
  return(pvals)
}

#Finding Komogorov-Smirnov p-values using bootstrapping based on FCS pmm imputations
KSstat3 <- function(Data,NumImp,B){
  invisible(capture.output(Imputations <- mice(Data, m=NumImp)))
  
  pvals <- rep(0,NumImp)
  for(i in 1:NumImp){
    Ximp <- as.matrix(complete(Imputations,i))
    FIT <-  lm(Ximp[,1] ~ Ximp[,2:5])
    sigmaobs <- summary(FIT)[[6]]
    KSobs <- ks.test(FIT$residuals,"pnorm",0,sigmaobs)$statistic
    
    KSsim <- rep(0,B)
    for(k in 2:(B+1)){
      Ynew <- rnorm(n, cbind(1,Ximp[,2:5])%*%FIT$coefficients, sigmaobs)
      Fit <- lm(Ynew ~ Ximp[,2:5])
      sigmasim <- summary(Fit)[[6]]
      KSsim[k-1] <- ks.test(Fit$residuals,"pnorm",0,sigmasim)$statistic
    }
    pvals[i] <- sum(KSobs < KSsim)/B
  }
  return(pvals)
}



#############################Simulation#####################################

##### Initial parameters and some initialization of vectors #####
Pval1a <- Pval1b <- Pval2a <- Pval2b <- Pval3a <- Pval3b <- rep(0,N)
AllVals <- list()
XYVals <- list()
NumImp <- 100
N <- 200
n <- 100


##### Data Generation #####
DataSets <- DataGen(N,n)

#Simulation over N datasets
for(j in 1:N){
  
  #Defining jth dataset based on generated datasets
  Y <- DataSets[[1]][[j]]
  X <- DataSets[[2]][[j]]
  n0 <- DataSets[[3]][[j]]
  n1 <- DataSets[[4]][[j]]
  n2 <- DataSets[[5]][[j]]
  n12 <- DataSets[[6]][[j]]
  
  #Obtaining initial values
  parest <- optim(c(Beta,sigma,c(120,0,0,0,10,120,0,0,10)), MaxMissingSimp(X,Y,n0))
  pars <- parest$par
  
  #Obtaining Bayesian parameter and imputation draws from posterior/posterior predictive
  Values <- BayesMethods(Y,X,NumImp+50,pars,n0,n1,n2,n12)
  AllVals[[j]] <- Values
  
  ##### Obtaining K-S p-values #####
  
  #p-value for bootstraping based on posterior predictive imputations
  KSvals1 <- KSstat1(Y,Values[[1]][51:(NumImp+50)],Values[[2]][51:(NumImp+50),], Values[[3]][51:(NumImp+50)], NumImp,20)
  Pval1a[j] <- mean(KSvals1)
  KSvals1b <- KSstat1(Y,Values[[1]][51:(NumImp+50)],Values[[2]][51:(NumImp+50),], Values[[3]][51:(NumImp+50)], NumImp,1)
  Pval1b[j] <- mean(KSvals1b)
  
  
  #p-value for posterior predictive checks based on posterior predictive imputations and posterior parameter values
  KSvals2 <- KSstat2(Y,Values[[1]][51:(NumImp+50)], cbind(Values[[2]][50:(NumImp+49),],Values[[3]][50:(NumImp+49)]), NumImp, 20)
  Pval2a[j] <- mean(KSvals2)
  KSvals2b <- KSstat2(Y,Values[[1]][51:(NumImp+50)], cbind(Values[[2]][50:(NumImp+49),],Values[[3]][50:(NumImp+49)]), NumImp, 1)
  Pval2b[j] <- mean(KSvals2b)
  
  #To get p-values with FCS, first obtaining FCS imputations
  Data <- cbind(Y,X)
  Data[Data[,2]==Inf,2] <- NA
  Data[Data[,3]==Inf,3] <- NA
  rownames(Data) <- c()
  
  #p-values based on bootstrapping with FCS imputations
  KSvals3a <- KSstat3(Data,NumImp,20)
  Pval3a[j] <- mean(KSvals3a)
  KSvals3b<- KSstat3(Data,NumImp,1)
  Pval3b[j] <- mean(KSvals3b)
}

#Finding final p-value means
Pvals <- cbind(Pval1a,Pval1b,Pval2a,Pval2b,Pval3a,Pval3b)
M <- colMeans(Pvals)
S <- apply(Pvals,2,sd) 
round(as.vector(matrix(c(M,S),2,6,byrow=TRUE)),3)

###Power plotting code (given all data for p-values across samples sizes###
###Where Po is power matrix for all of the sample sizes in "values" below for all 6 strategies###
par(mfrow=c(1,2))
values <- c(20,30,40,60,80,100,150,200,250,300)
plot(values,(Po[,1]), main="Power Curves", cex.lab=1.2, xlab="Sample Size", ylab="Power",ylim=c(0,1))
points(values,(Po[,3]), col="blue")
points(values,(Po[,5]), col="red")
lo1 <- loess((Po[,1]) ~ values,span=0.6)
lines(values,predict(lo1),lwd=2, col="black",lty=1)
lo2 <- loess((Po[,3]) ~ values,span=0.6)
lines(values,predict(lo2),lwd=2, col="blue",lty=2)
lo3 <- loess((Po[,5]) ~ values,span=0.6)
lines(values,predict(lo3),lwd=2, col="red",lty=3)
legend(145,0.4,c("Bootstrap (BS)","Posterior Predictive (PP)","Reimputed (RI)"), lty=c(1,2,3),col=c("black","blue","red"))

###Mean P-value plotting code (given all data for p-values across samples sizes###
###Where P is mean p-value matrix for all of the sample sizes in "values" below for all 6 strategies###
values <- c(20,30,40,60,80,100,150,200,250,300)
plot(values,(P[,1]), main="Mean P-value Curves", cex.lab=1.2, xlab="Sample Size", ylab="P-value",ylim=c(0,0.5))
points(values,(P[,3]), col="blue")
points(values,(P[,5]), col="red")
lo1 <- loess((P[,1]) ~ values,span=0.7)
lines(values,predict(lo1),lwd=2, col="black",lty=1)
lo2 <- loess((P[,3]) ~ values,span=0.7)
lines(values,predict(lo2),lwd=2, col="blue",lty=2)
lo3 <- loess((P[,5]) ~ values,span=0.7)
lines(values,predict(lo3),lwd=2, col="red",lty=3)
legend(145,0.4,c("Bootstrap (BS)","Posterior Predictive (PP)","Reimputed (RI)"), lty=c(1,2,3),col=c("black","blue","red"))

