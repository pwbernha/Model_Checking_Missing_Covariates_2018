################Necessary Packages for Application###################
library(foreign)
library(ResourceSelection)
library(pROC)
library(mitools)

#######################Reading in Data##############################
SIS <- read.dta(file="SIS.dta")
#Note: Data is available upon request from Columbia Population Reseach Center:
#http://cupop.columbia.edu/research/research-areas/social-indicators-survey-sis

######################Organizing Data###############################

#Only considering those with a job
SIS2 <- SIS[which(SIS[,286]=='1 yes'),]

#Getting/scaling to annual salary (subject to missing)
Salary <- as.vector(SIS2[,299])
Salary[is.na(SIS2[,296])==FALSE] <- as.vector(SIS2[is.na(SIS2[,296])==FALSE,296])*52
Salary[is.na(SIS2[,297])==FALSE] <- as.vector(SIS2[is.na(SIS2[,297])==FALSE,297])*26
Salary[is.na(SIS2[,298])==FALSE] <- as.vector(SIS2[is.na(SIS2[,298])==FALSE,298])*12
Salary <- Salary/1000

#Health Insurance Status (defining missing broadly)
HI <- as.vector(SIS2[,169])
HI[HI=="9. RF"] <- NA
HI[HI=="8. DK"] <- NA
HI[HI=="1 yes"] <- 1
HI[HI=="2. no"] <- 0
HI <- as.numeric(HI)

#Marriage Status (defining missing broadly)
Married <- as.vector(SIS2[,66])
Married[Married=="9. RF"] <- NA
Married[Married!="1. married" & is.na(Married)==FALSE] <- 0
Married[Married=="1. married"] <- 1
Married <- as.numeric(Married)

#Age
Age <- as.vector(SIS2[,68])

#Health Status (defining missing broadly)
HS <- as.vector(SIS2[,154])
HS[HS=="8. DK"] <- NA
HS[HS=="9. RF"] <- NA
HS[HS=="1. excellent" | HS=="2. very good" | HS=="3. good"] <- 1
HS[HS=="4. fair" | HS=="5. poor" ] <- 0
HS <- as.numeric(HS)

#Sex
Sex <- as.vector(SIS2[,69])
Sex[Sex=="1 Male"] <- 0 
Sex[Sex=="2 Female"] <- 1
Sex <- as.numeric(Sex)

#Creating dataframe of all variables
Health <- data.frame(as.factor(HI),Age,as.factor(HS),as.factor(Married),Salary,Sex)
colnames(Health) <- c("HI", "Age", "HS", "Married", "Salary", "Sex")



###################Obtaining imputations, estimates################
Imputations <- mice(Health,300)
fit <- with(Imputations, glm(HI ~ Married + HS + Sex + Salary, family="binomial"))
summary(pool(fit))



###################Hosmer-Lemeshow Test################
HL <- rep(0,300)
for(i in 1:300){
  HL[i] <- hoslem.test(fit[[4]][[i]]$y, fitted(fit[[4]][[i]]), g=10)$p.value
}
mean(HL)
sum(HL<0.05)/300 #Percent of data sets with p-value<0.05



######################ROC Plots#########################
par(mfrow=c(1,2))
AUC <- rep(0,300)

#Calculating & plotting ROC curves
for(i in 1:300){
  preds=predict(fit[[4]][[i]])
  roc1=roc(fit[[4]][[i]]$y ~ preds)
  AUC[i] <- roc1$auc
  if(i==1) plot(roc1, main="ROC Curves for Completed Data Sets", cex.lab=1.2) else if(i<21) plot(roc1,add=TRUE)
}
mean(AUC)			#mean AUC
text(0.45,0.2, expression(paste(bar(" AUC"), " = 0.647")))

#Plotting 10, 50, 90 percentile ROC curves
A1 <- sort(AUC)[30]	#10th percentile curve
A3 <- sort(AUC)[150]	#50th percentile curve
A2 <- sort(AUC)[270]	#90th percentile curve
plot(roc(fit[[4]][[which(AUC==A1)]]$y ~ predict(fit[[4]][[which(AUC==A1)]])), main="ROC Percentile Curves", cex.lab=1.2, lty=2)
plot(roc(fit[[4]][[which(AUC==A2)]]$y ~ predict(fit[[4]][[which(AUC==A2)]])), add=TRUE, lty=3)
plot(roc(fit[[4]][[which(AUC==A3)]]$y ~ predict(fit[[4]][[which(AUC==A3)]])), add=TRUE)
legend(0.67,0.25,c("10th Percentile, AUC=0.612", "50th Percentile, AUC=0.648", "90th Percentile, AUC=0.677"), lty=c(2,1,3))



######################NMAR Toy Example########################
Health2 <- data.frame(Salary,as.factor(HI),Age,as.factor(HS),as.factor(Married),Sex)
colnames(Health2) <- c("Salary", "HI", "Age", "HS", "Married", "Sex")
Imputations2 <- mice(Health2,800)  #larger number imps since some will be rejected

#Missingness probability function
Pr <- function(Salary,Sex) exp(0.5-0.02*Salary-0.2*Sex)/(1+exp(0.5-0.02*Salary-0.2*Sex))
miss <- which(is.na(Health2[,1]==TRUE))

#Obtaining Stacked Data set
FinalImps <- complete(Imputations2,"long")

#Obtaining imputations with adjusted probabilities
s<-0
for(i in 1:1338){
  if(any(miss==i)){
    s <- s+1
    missprobs <- Pr(as.numeric(Imputations2$imp$Salary[s,]), as.numeric(Health2[i,6]))
    pobs <- rbinom(800,1,(1-missprobs))
    FinalImps[seq(i,1338*sum(pobs==0),1338),] <- FinalImps[seq(i,1070400 ,1338),][pobs==0,] 
  }
}

Imps <- list()
for(i in 1:300){
  Imps[[i]] <- as.data.frame(FinalImps[(i*1338-1337):(i*1338),3:8])
}
Imps2 <- imputationList(Imps)

lmfit <- with(Imps2, glm(HI ~  HS + Married + Salary + Sex, family="binomial"))
lmfit <- as.mira(lmfit)
summary(pool(lmfit))

#Finding summary statistics for the imputations under nmar
par(mfrow=c(1,2))
AUC <- rep(0,300)
for(i in 1:300){
  preds=predict(lmfit[[4]][[i]])
  roc1=roc(lmfit[[4]][[i]]$y ~ preds)
  AUC[i] <- roc1$auc
  if(i==1) plot(roc1, main="ROC Curves for Completed Data Sets", cex.lab=1.2) else if(i<21) plot(roc1,add=TRUE)
}
mean(AUC)

HL <- rep(0,300)
for(i in 1:300){
  HL[i] <- hoslem.test(lmfit[[4]][[i]]$y, fitted(lmfit[[4]][[i]]), g=10)$p.value
}
mean(HL)
