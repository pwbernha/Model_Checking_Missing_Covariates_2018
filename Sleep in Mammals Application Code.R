################Necessary Packages for Application###################
library(SDMTools)
library(mice)
library(sfsmisc)
library(survival)
library(ggplot2)


#######################Reading in Data##############################
sleep <- read.csv("sleep.csv", header=TRUE)
#Note: Data is easily obtained at many online locations, for example, http://www.statsci.org/data/general/sleep.html

######################Organizing Data###############################
Sleep <- cbind(log(sleep[,3]), sleep[,6:5], log(sleep[,7:8]), log(sleep[,2]))
colnames(Sleep) <- c("LBrainWgt", "SleepTime", "DreamTime", "LSpan", "LGest","LBodyWgt")
SleepA <- cbind(sleep[,3], sleep[,6:5], log(sleep[,7:8]), log(sleep[,2]))
colnames(SleepA) <- c("BrainWgt", "SleepTime", "DreamTime", "LSpan", "LGest","LBodyWgt")


##########Obtaining Imputations, Estimates, and Final Model#########

#Obtaining Imputations for both potential response cases (log or untransformed)
Imputations <- mice(Sleep,m=300)
ImputationsA <- mice(SleepA,m=300)

#Fitting models based on completed data sets for both potential response cases
lmfit <- with(Imputations, lm(LBrainWgt ~  SleepTime + DreamTime + LSpan + LGest + LBodyWgt))
coefs <- summary(pool(lmfit))[,1]
lmfitA <- with(ImputationsA, lm(BrainWgt ~  SleepTime + DreamTime + LSpan + LGest + LBodyWgt))
coefsA <- summary(pool(lmfitA))[,1]

#Finding an estimate for sigma
sig <- NULL
for(i in 1:300){
  sig <- c(sig, summary(lmfit[[4]][[i]])$sigma)
}
sig <- mean(sig)

#Residual Histograms/QQ plots for both response cases (in example here, for 7th dataset)
par(mfrow=c(1,2))
hist(lmfit[[4]][[7]]$residual, freq=FALSE, main="Log-transformed Brain Weight", xlab="Residuals", cex.lab=1.2, breaks=12,xlim=c(-1.5,1.5))
curve(dnorm(x, mean(lmfit[[4]][[7]]$residual), sd(lmfit[[4]][[7]]$residual)),  col="darkblue", lwd=2, add=TRUE)
hist(lmfitA[[4]][[7]]$residual, freq=FALSE, main="Untransformed Response", xlab="Residuals", cex.lab=1.2, breaks=20,xlim=c(-2000,4000))
curve(dnorm(x, mean(lmfitA[[4]][[7]]$residual), sd(lmfitA[[4]][[7]]$residual)),  col="darkblue", lwd=2, add=TRUE)

par(mfrow=c(1,2))
qqnorm(lmfit[[4]][[7]]$residual,  main="Log-transformed Brain Weight", xlab="Standard Normal Quantiles", ylab="Sample Quantiles", cex.lab=1.2)
qqline(lmfit[[4]][[7]]$residual)
qqnorm(lmfitA[[4]][[7]]$residual,  main="Untransformed Brain Weight", xlab="Standard Normal Quantiles", ylab="Sample Quantiles", cex.lab=1.2)
qqline(lmfitA[[4]][[7]]$residual)



###############Completed/Stacked Data Residual Histograms#############

#Residual histograms for three example data sets (5, 30, 246)
layout(matrix(c(1,1,1,2,2,2,3,3,3,4,4,5,5,5,5,5,6,6), 2, 9, byrow = TRUE))
for(i in  c(5,30,246)){
  hist(lmfit[[4]][[i]]$residual, main=paste("Completed Data Set", i), xlab="Residuals", freq=FALSE, breaks=20, xlim=c(-1.5,1.5))
  curve(dnorm(x, mean(lmfit[[4]][[i]]$residual), sd(lmfit[[4]][[i]]$residual)),  col="darkblue", lwd=2, add=TRUE)
}

#Residual histograms for stacked data set

#Complete Cases
cc <- which(is.na(Sleep[,2])==FALSE &is.na(Sleep[,3])==FALSE & is.na(Sleep[,4])==FALSE &is.na(Sleep[,5])==FALSE)

#Adding blank plot in row 2
frame() 

#Finding residuals for stacked data 
resids <- fits <- NULL
for(i in 1:300){
  resids <- c(resids,(Sleep[,1]-cbind(1,as.matrix(complete(Imputations,i)[,2:6]))%*%coefs))
}
WeightsOne <- rep(1/(300*62),62)
WeightsOne[cc] <- 1/300
Weights <- rep(WeightsOne,300)

#Plotting Stacked Data Histogram
hist(resids, freq=FALSE, main="Stacked Data Set", xlab="Residuals", cex.lab=1.2, breaks=20, xlim=c(-1.5,1.5))
curve(dnorm(x, wt.mean(resids,Weights), wt.sd(resids,Weights)),  col="darkblue", lwd=2, add=TRUE)



#######################Residual vs. Fitted Plots########################
#Finding Stacked Data Set and associated Weights
resids <- fits <- NULL
resids <- (Sleep[cc,1]-cbind(1,as.matrix(Sleep[cc,2:6]))%*%coefs)
fits <- cbind(1,as.matrix(Sleep[cc,2:6]))%*%coefs
for(i in 1:300){
  resids <- c(resids,(Sleep[-cc,1]-cbind(1,as.matrix(complete(Imputations,i)[-cc,2:6]))%*%coefs))
  fits <- c(fits,(cbind(1,as.matrix(complete(Imputations,i)[-cc,2:6]))%*%coefs))
}
Weights <- c(rep(1,length(cc),),rep(1/300,(62-length(cc))*300))

#Plotting Residuals versus predicted (stacked)
par(mfrow=c(1,2))
plot(fits,resids, xlab="Fitted Response Values", ylab="Residuals", main="Residuals vs. Fitted Values",cex.lab=1.2 )
l.fit <- loess(resids ~ fits, weights=Weights, span=0.8)
l.pred <- predict(l.fit,seq(min(fits),max(fits),0.01))
lines(seq(min(fits),max(fits),0.01),l.pred,col="blue",lwd=3)

#Finding average fitted value for each observation (for complete cases, average is of a single estimate)
Mean.fits <- fits[1:42]
for(i in 1:20) Mean.fits <- c(Mean.fits, mean(fits[seq(42+i,6042,20)]))

#Creating blank plot with an appropriate range
plot(fits,resids, xlab="Fitted Response Values", ylab="Residuals", main="Observed and Randomly Generated Loess Curve Fits",cex.lab=1.2, col="white")

#Generating a new response vector, finding imputations and estimates, and fitting a weighted loess curves to residuals
for(i in 1:50){
  y <- rnorm(62,Mean.fits,sig)
  Sleepy <-cbind(y,rbind(Sleep[cc,2:6],Sleep[-cc,2:6]))
  Imps <- mice(Sleepy,m=300)
  lmfit2 <- with(Imps, lm(y ~  SleepTime + DreamTime + LSpan + LGest + LBodyWgt))
  coefs2 <- summary(pool(lmfit2))[,1]
  resids2 <- fits2 <- NULL
  resids2 <- (Sleepy[1:42,1]-cbind(1,as.matrix(Sleepy[1:42,2:6]))%*%coefs2)
  fits2 <- cbind(1,as.matrix(Sleepy[1:42,2:6]))%*%coefs2
  for(i in 1:300){
    resids2 <- c(resids2,(Sleepy[43:62,1]-cbind(1,as.matrix(complete(Imps,i)[43:62,2:6]))%*%coefs2))
    fits2 <- c(fits2,(cbind(1,as.matrix(complete(Imps,i)[43:62,2:6]))%*%coefs2))
  }
  l.fit2 <- loess(resids2 ~ fits2, weights=Weights, span=0.8)
  l.pred2 <- predict(l.fit2,seq(min(fits2),max(fits2),0.01))
  lines(seq(min(fits2),max(fits2),0.01),l.pred2,col="black",lwd=1)
}
lines(seq(min(fits),max(fits),0.01),l.pred,col="blue",lwd=3)

#################Partial Regression Leverage Plots####################
par(mfrow=c(1,2))

#####Partial Regression Plot for Sleep#####
#Fitting model without total sleep
lmfit1 <- with(Imputations, lm(LBrainWgt ~  DreamTime + LSpan + LGest + LBodyWgt))
coefsNoSleep <- summary(pool(lmfit1))[,1]

#Fitting model on total sleep
lmfit2 <- with(Imputations, lm(SleepTime ~  DreamTime + LSpan + LGest + LBodyWgt))
coefsOnSleep <- summary(pool(lmfit2))[,1]

#Finding residuals for both models
residsNoSleep <- residsOnSleep <- NULL
residsNoSleep <- (Sleep[cc,1]-as.matrix(cbind(1,Sleep[cc,3:6]))%*%coefsNoSleep)
residsOnSleep <- (Sleep[cc,2]-as.matrix(cbind(1,Sleep[cc,3:6]))%*%coefsOnSleep)
for(i in 1:300){
  residsNoSleep <- c(residsNoSleep,(Sleep[-cc,1]-cbind(1,as.matrix(complete(Imputations,i)[-cc,3:6]))%*%coefsNoSleep))
  residsOnSleep <- c(residsOnSleep,(as.vector(complete(Imputations,i)[-cc,2])-cbind(1,as.matrix(complete(Imputations,i)[-cc,3:6]))%*%coefsOnSleep))
}

#Plotting both sets of residuals against one another, with loess curve
plot(residsNoSleep,residsOnSleep, xlab="Sleep Residuals", ylab="Log Brain Weight Residuals", main="Partial Regression Leverage Plot for Sleep",cex.lab=1.2 )
l.fit <- loess(residsNoSleep ~ residsOnSleep, weights=Weights, span=0.8)
l.pred <- predict(l.fit,seq(min(residsOnSleep),max(residsOnSleep),0.01))
lines(seq(min(residsOnSleep),max(residsOnSleep),0.01),l.pred,lwd=3)

#####Partial Regression Plot for Span#####
#Fitting model without lifespan
lmfitnew3 <- with(Imputations, lm(LBrainWgt ~  SleepTime + DreamTime  + LGest + LBodyWgt))
coefsNoSpan <- summary(pool(lmfitnew3))[,1]

#Fitting model on lifespan
lmfitnew4 <- with(Imputations, lm(LSpan ~ SleepTime + DreamTime + LGest + LBodyWgt))
coefsOnSpan <- summary(pool(lmfitnew4))[,1]

#Finding residuals for both models
residsNoSpan <- residsOnSpan <- NULL
residsNoSpan <- (Sleep[cc,1]-as.matrix(cbind(1,Sleep[cc,c(2:3,5:6)]))%*%coefsNoSpan)
residsOnSpan <- (Sleep[cc,4]-as.matrix(cbind(1,Sleep[cc,c(2,3,5,6)]))%*%coefsOnSpan)
for(i in 1:300){
  residsNoSpan <- c(residsNoSpan,(Sleep[-cc,1]-cbind(1,as.matrix(complete(Imputations,i)[-cc,c(2:3,5:6)]))%*%coefsNoSpan))
  residsOnSpan <- c(residsOnSpan,(as.vector(complete(Imputations,i)[-cc,4])-cbind(1,as.matrix(complete(Imputations,i)[-cc,c(2:3,5:6)]))%*%coefsOnSpan))
}

#Plotting both sets of residuals against one another, with loess curve
plot(residsNoSpan,residsOnSpan, xlab="Log Lifespan Residuals", ylab="Log Brain Weight Residuals", main="Partial Regression Leverage Plot for Log Lifespan",cex.lab=1.2 )
l.fit <- loess(residsNoSpan ~ residsOnSpan, weights=Weights, span=0.8, degree=1)
l.pred <- predict(l.fit,seq(min(residsOnSpan),max(residsOnSpan),0.01))
lines(seq(min(residsOnSpan),max(residsOnSpan),0.01),l.pred,lwd=3)



##################################Goodness-of-Fit Tests##################################

#########Shapiro-Wilk p-values#########
pvals <- rep(0,300)
for(i in 1:300){
  pvals[i] <- shapiro.test(lmfit[[4]][[i]]$residual)$p
}

#Mean Shapiro-Wilks p-value
mean(pvals)


#########Kolmogorv-Smirnov p-values#######

####Method 1 in paper based on frequentist bootstrapping####
B <- 1000 #Of bootstrap samples
pvalk1 <- rep(0,300)
for(i in 1:300){
  print(i)
  Ximp <- as.matrix(complete(Imputations,i))
  
  #Obtaining K-S Stat for ith imputed data set
  FIT <- lmfit[[4]][[i]]
  sigmaimp <- summary(FIT)[[6]]
  KSobs <- ks.test(FIT$residuals,"pnorm",0,sigmaimp)$statistic
  
  #Obtaining K-S Stat for B bootstrap samples
  KSsim <- rep(0,B)
  for(j in 1:B){
    Ynew <- rnorm(62, cbind(1,Ximp[,2:6])%*%FIT$coefficients,sigmaimp)	
    Fit <- lm(Ynew ~ Ximp[,2:6])
    sigmanew <- summary(Fit)[[6]]
    KSsim[j] <- ks.test(Fit$residuals,"pnorm",0,sigmanew)$statistic
  }
  
  pvalk1[i] <- sum(KSobs < KSsim)/B	#p-value for ith data set
  
}

#Mean Kolmogorv-Smirnov p-value
mean(pvalk1)

#Plotting histogram of p-values for S-W and K-S simultaneously
hist(pvalk1, main="Completed Data Goodness-of-Fit Tests", xlab="P-values", cex.lab=1.2, xlim=c(0,1), ylim=c(0,40), breaks=20, col=rgb(0.1,0.1,0.1,0.5))
hist(pvals, col=rgb(0.8,0.8,0.8,0.5), add=T, breaks=20)

#Adding a fancy legend
box()
legend(0.75,35,c("Kolmogorov-Smirnov", "Shapiro-Wilk"), lty=c(2,1),col="white",cex=1.05)
polygon(x=c(0.78,0.78,0.795),y=c(32,33.5,33.5), col="grey43")
polygon(x=c(0.78,0.795,0.795),y=c(32,32,33.5), col="gray66")
polygon(x=c(0.78,0.78,0.795),y=c(29.5,31,31), col="gray66")
polygon(x=c(0.78,0.795,0.795),y=c(29.5,29.5,31), col="gray88")

#Alternative Method for Plotting
#dat <- data.frame(factor(rep(c("S", "K"), each=300)),c(pvals,pvalk1))
#colnames(dat) <- c("Type","pval")
#ggplot(dat, aes(x=pval, color=Type)) +
#  geom_histogram(fill="grey41", alpha=0.5, position="identity",binwidth=0.03) +scale_color_grey()+scale_fill_grey() +
# theme_classic()+theme_bw(base_size = 12, base_family = "")



##################################Influence Diagnostics##################################

#####Cook's Distances#####
#Calculating Cook's distance
Cooks <- matrix(0,62,300)
for(i in 1:300){
  Cooks[,i] <- cooks.distance(lmfit[[4]][[i]])
}

#Creating vector of names for purpose of plotting boxplots of Cook's distances for each mammal
name <- c("Afr. Elephant","Afr. Giant Rat", "Arctic Fox",  "Arctic Squirrel", "Asian Elephant",          
          "Baboon", "Big Br. Bat"," Braz. Tapir", "Cat","Chimpanzee", "Chinchilla", "Cow",  "Desert Hedgehog", "Donkey",               
          "Eastern Mole", "Echidna", "Euro. Hedgehog", "Galago", "Genet", "Giant Armadillo", "Giraffe", "Goat", "Golden Hamster",
          "Gorilla",  "Gray Seal", "Gray Wolf", "Ground Squirrel", "Guinea Pig", "Horse", "Jaguar", "Kangaroo",
          "Short-tl Shrew", "Little Br. Bat", "Man", "Mole Rat", "Mtn. Beaver", "Mouse", "Muskshrew",             
          "N.A. Opossum", "Banded Armadillo","Okapi", "Owl Monkey", "Patas Monkey", "Phanlanger", "Pig", "Rabbit", "Raccoon", "Rat", "Red Fox",
          "Rhesus Monkey", "Rock Hyrax(H.)", "Rock Hyrax(P.)",  "Roe Deer", "Sheep", "Slow Loris", "Starnosed Mole", "Tenrec","Tree Hyrax", "Tree Shrew",              
          "Vervet", "Water Opossum", "Y.B. Marmot")

par(mar=c(8,4,2.1,1.1))
boxplot.matrix(t(Cooks),las=2,names=name, ylab="Cook's Distance", main="Boxplots of Cook's Distances for 62 Mammal Species") 
abline(h=(4/62))

######dfbetas#####
#Calculating DFbetas based on completed data set averages
DFbetas <- matrix(0,62,6)
for(i in 1:300){
  DFbetas <- DFbetas + dfbetas(lmfit[[4]][[i]])/300
}
DFbetas[c(7,13,33,34),]

#Creating Stacked matrix
Stacked <- NULL
for(i in 1:300){
  Stacked <- rbind(Stacked,as.matrix(complete(Imputations,i)))
}


#Calculating DFbetas based on stacked data set
Fits <- summary(lm(Stacked[,1] ~ Stacked[,2] + Stacked[,3] + Stacked[,4] + Stacked[,5] + Stacked[,6]))
XX <- Fits$cov.unscaled
BETA <- Fits$coefficients[,1]
DFbetas2 <- matrix(0,62,6)
for(i in 1:62){
  Stackedb <- Stacked[-seq(i,18600,62),]
  fit <- summary(lm(Stackedb[,1] ~ Stackedb[,2] + Stackedb[,3] + Stackedb[,4] + Stackedb[,5] + Stackedb[,6]))
  DFbetas2[i,] <- (BETA - fit$coefficients[,1])/(sqrt(fit$sigma^2*diag(XX)*300/(1-24/300)))
}
DFbetas2[c(7,13,33,34),]


#################################Testing Imputation Models######################################
#reimputing data for all covariates with missing values, concatenated to imputed original dataset (for 10000 imputations)
Imputations2 <- mice(rbind(as.matrix(Sleep),cbind(Sleep[,1],matrix(NA,62,4),Sleep[,6])),m=10000)
indMat <- matrix(0,10000,7)
for(i in 1:10000){
  com <- as.matrix(complete(Imputations2,i))
  lmfit2a <-lm(com[1:62,1] ~  com[1:62,2] + com[1:62,3] + com[1:62,4] + com[1:62,5] + com[1:62,6])
  lmfit2b <-lm(com[63:124,1] ~  com[63:124,2] + com[63:124,3] + com[63:124,4] + com[63:124,5] + com[63:124,6])
  indMat[i,1:6] <- as.numeric(summary(lmfit2a)$coefficients[,1] > summary(lmfit2b)$coefficients[,1])
  indMat[i,7] <- summary(lmfit2a)$sigma > summary(lmfit2b)$sigma
}
colMeans(indMat) #proportion of each parameter estimate that is greater for the reimputed data



