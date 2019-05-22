#' Global sensitivity analysis
#' @date April 30, 2019
#' @author Lauren White
#' @email lwhite@sesync.org

#' Adapted from: https://daphnia.ecology.uga.edu/drakelab/wp-content/uploads/2015/07/sensitivity-ebola.pdf
#' load ode functions for each of the models
source('~/JustinianPlague/Plague_model_functions.R')
#' load uniform and non uniform LHS distributions
source('~/JustinianPlague/LHSnonuniform.R')


#install.packages('lhs')
require(lhs) #add the lhs library
library(sensitivity)
require(ggplot2)
library(tidyverse)

set.seed(2718) #set random seed
times <- seq(0, 5000, by= 1) #time sequence to integrate over for all models
h <- 100 #choose number of parameter sets/subdivisions to sample to sample
uniform=FALSE #unfiform or non-uniform distributions

#### Sensitivity Analysis for Pneumonic SIR Model ####----------------------------

#expected parameter values
parameters <- c(beta_p = 0.45, gamma_p = 1/2.5, b_h=1/(25*365), d_h=1/(25*365)) #you can play with transmission and recovery rates here

if(uniform==TRUE){ params.set<-LHS_pSIR.uniform  #uniform LHS distribution from `LHSnonuniform.R`
} else{ params.set<-LHS_pSIR} #non-uniform LHS distribution from `LHSnonuniform.R`

init <- c(S_h = 499999, I_h = 1, D_h=0)
pSIR <- data.frame(matrix(rep(NA,h),nrow=h, ncol= ncol(params.set)+5))
colnames(pSIR)<-c(names(parameters), "MaxInf", "MaxDur", "Thresh1", "Thresh0.1", "Thresh0.01")

for(i in 1:h){
  pSIR[i,1:4] <- params <- as.list(c(params.set[i,]))
  out <- as.data.frame(ode(init, times, pneumonicSIR, params))
  pSIR[i,5] <- max(out$D_h)
  pSIR[i,6] <- which.max(out$D_h)
  
  difference<-as.vector(NA)
  for(j in 2:length(out$D_h)){
    difference[j-1]<-out$D_h[j]-out$D_h[j-1]
  }
  pSIR[i,7]<-ifelse(length(which(difference>1)[length(which(difference>1))])>0, which(difference>1)[length(which(difference>1))], 0) 
  pSIR[i,8]<-ifelse(length(which(difference>0.1)[length(which(difference>0.1))])>0, which(difference>0.1)[length(which(difference>0.1))], 0) 
  pSIR[i,9]<-ifelse(length(which(difference>0.01)[length(which(difference>0.01))])>0, which(difference>0.01)[length(which(difference>0.01))], 0) 
}

#plot boxplots
par(mfrow=c(1,2))
boxplot(pSIR$MaxInf~pSIR$beta_p, main= expression(paste("Effect of ", beta[p], " on Size")))
boxplot(pSIR$Thresh1~pSIR$beta_p, main= expression(paste("Effect of ", beta[p], " on Duration")))
boxplot(pSIR$MaxInf~pSIR$gamma_p,  main= expression(paste("Effect of ", gamma[p], " on Size")))
boxplot(pSIR$Thresh1~pSIR$gamma_p, main= expression(paste("Effect of ", gamma[p], " on Duration")))
boxplot(pSIR$MaxInf~pSIR$b_h, main= expression(paste("Effect of ", b[h], " on Size")))
boxplot(pSIR$Thresh1~pSIR$b_h, main= expression(paste("Effect of ", b[h], " on Duration")))
boxplot(pSIR$MaxInf~pSIR$d_h, main= expression(paste("Effect of ", d[h], " on Size")))
boxplot(pSIR$Thresh1~pSIR$d_h, main= expression(paste("Effect of ", d[h], " on Duration")))

par(mfrow=c(1,2))
boxplot(pSIR$MaxInf, main= "Outbreak Size", ylab= "Number of Dead Humans", ylim=c(0,500000))
boxplot(pSIR$Thresh1, main= "Outbreak Duration", ylab="Time (Days)")

bonferroni.alpha <- 0.05/5
prcc_size <- pcc(pSIR[,1:4], pSIR$MaxInf, nboot = 1000, rank=TRUE, conf=1-bonferroni.alpha)
prcc_duration <- pcc(pSIR[,1:4], pSIR$Thresh1, nboot = 1000, rank=TRUE, conf=1-bonferroni.alpha)

#plot correlation coefficients and confidence intervals for epidemic size and duration
size<-prcc_size$PRCC
size$param<-rownames(size)
colnames(size)[4:5] <- c("maxCI", "minCI")

ggplot(size, aes(x=param, y = original)) +
  geom_point(size = 4)+ 
  geom_errorbar(aes(ymax = maxCI, ymin = minCI))+
  ggtitle("PRCC Outbreak Size: Pneumonic SIR")

duration<-prcc_duration$PRCC
duration$param<-rownames(duration)
colnames(duration)[4:5] <- c("maxCI", "minCI")

ggplot(duration, aes(x=param, y = original)) +
  geom_point(size = 4)+ 
  geom_errorbar(aes(ymax = maxCI, ymin = minCI))+
  ggtitle("PRCC Outbreak Duration: Pneumonic SIR")


#### Sensitivity Analysis for Pneumonic SEIR Model #### ----------------------------

#expected parameter values
init <- c(S_h = 499999, E_h= 0, I_h = 1, D_h=0) #population size, and how many individuals start in each susceptible, infected, or removed category
parameters <- c(beta_p = 0.45, sigma_p= 1/4.3, gamma_p = 1/2.5, b_h=1/(25*365), d_h=1/(25*365)) #you can play with transmission and recovery rates here

if(uniform==TRUE){ params.set<-LHS_pSEIR.uniform  #uniform LHS distribution from `LHSnonuniform.R`
} else{ params.set<-LHS_pSEIR} #non-uniform LHS distribution from `LHSnonuniform.R`

pSEIR <- data.frame(matrix(rep(NA,h),nrow=h, ncol= ncol(params.set)+5))
colnames(pSEIR)<-c(names(parameters), "MaxInf", "MaxDur", "Thresh1", "Thresh0.1", "Thresh0.01")

for(i in 1:h){
  pSEIR[i,1:5] <- params <- as.list(c(params.set[i,]))
  out <- as.data.frame(ode(init, times, pneumonicSEIR, params))
  pSEIR[i,6] <- max(out$D_h)
  pSEIR[i,7] <- which.max(out$D_h)
  
  difference<-as.vector(NA)
  for(j in 2:length(out$D_h)){
    difference[j-1]<-out$D_h[j]-out$D_h[j-1]
  }
  pSEIR[i,8]<-ifelse(length(which(difference>1)[length(which(difference>1))])>0, which(difference>1)[length(which(difference>1))], 0) 
  pSEIR[i,9]<-ifelse(length(which(difference>0.1)[length(which(difference>0.1))])>0, which(difference>0.1)[length(which(difference>0.1))], 0) 
  pSEIR[i,10]<-ifelse(length(which(difference>0.01)[length(which(difference>0.01))])>0, which(difference>0.01)[length(which(difference>0.01))], 0) 
  
}
par(mfrow=c(1,2))
boxplot(pSEIR$MaxInf~pSEIR$beta_p, main= expression(paste("Effect of ", beta[p], " on Size")))
boxplot(pSEIR$Thresh1~pSEIR$beta_p, main= expression(paste("Effect of ", beta[p], " on Duration")))
boxplot(pSEIR$MaxInf~pSEIR$sigma_p, main= expression(paste("Effect of ", sigma[p], " on Size")))
boxplot(pSEIR$Thresh1~pSEIR$sigma_p, main= expression(paste("Effect of ", sigma[p], " on Duration")))
boxplot(pSEIR$MaxInf~pSEIR$gamma_p,  main= expression(paste("Effect of ", gamma[p], " on Size")))
boxplot(pSEIR$Thresh1~pSEIR$gamma_p, main= expression(paste("Effect of ", gamma[p], " on Duration")))
boxplot(pSEIR$MaxInf~pSEIR$b_h, main= expression(paste("Effect of ", b[h], " on Size")))
boxplot(pSEIR$Thresh1~pSEIR$b_h, main= expression(paste("Effect of ", b[h], " on Duration")))
boxplot(pSEIR$MaxInf~pSEIR$d_h, main= expression(paste("Effect of ", d[h], " on Size")))
boxplot(pSEIR$Thresh1~pSEIR$d_h, main= expression(paste("Effect of ", d[h], " on Duration")))

par(mfrow=c(1,2))
boxplot(pSEIR$MaxInf, main= "Outbreak Size", ylab= "Number of Dead Humans", ylim=c(0,500000))
boxplot(pSEIR$Thresh1, main= "Outbreak Duration", ylab="Time (Days)")

bonferroni.alpha <- 0.05/5
prcc_size <- pcc(pSEIR[,1:5], pSEIR[,6], nboot = 1000, rank=TRUE, conf=1-bonferroni.alpha)
prcc_duration <- pcc(pSEIR[,1:5], pSEIR[,7], nboot = 1000, rank=TRUE, conf=1-bonferroni.alpha)

#plot correlation coefficients and confidence intervals for epidemic size and duration
size<-prcc_size$PRCC
size$param<-rownames(size)
colnames(size)[4:5] <- c("maxCI", "minCI")

ggplot(size, aes(x=param, y = original)) +
  geom_point(size = 4)+ 
  geom_errorbar(aes(ymax = maxCI, ymin = minCI))+
  ggtitle("PRCC Outbreak Size: Pneumonic SEIR")

duration<-prcc_duration$PRCC
duration$param<-rownames(duration)
colnames(duration)[4:5] <- c("maxCI", "minCI")

ggplot(duration, aes(x=param, y = original)) +
  geom_point(size = 4)+ 
  geom_errorbar(aes(ymax = maxCI, ymin = minCI))+
  ggtitle("PRCC Outbreak Duration: Pneumonic SEIR")

# Bubonic SIR model -------------------------------------------------------

#Define initial conditions and expected parameter values
init <- c(S_r=499999, I_r=1, R_r=0, D_r=0, H=6, Fl=0, S_h = 500000, I_h = 0, R_h=0, D_h=0) #population size, and how many individuals start in each susceptible, infected, or removed category
parameters <- c(beta_r = 0.09, alpha=3/500000, gamma_r = 1/5.15, g_r=0.1, r_f=0.0084, K_f=6, d_f=1/5, beta_h=0.19, gamma_h=1/10, g_h=0.34, b_h=1/(25*365), d_h=1/(25*365)) #you can play with transmission and recovery rates here

if(uniform==TRUE){ params.set<-LHS_bSIR.uniform  #uniform LHS distribution from `LHSnonuniform.R`
} else{ params.set<-LHS_bSIR} #non-uniform LHS distribution from `LHSnonuniform.R`

bSIR <- data.frame(matrix(rep(NA,h),nrow=h, ncol= ncol(params.set)+5))
colnames(bSIR)<-c(names(parameters), "MaxInf", "MaxDur", "Thresh1", "Thresh0.1", "Thresh0.01")

for(i in 1:h){
  bSIR[i,1:length(parameters)] <- params <- as.list(c(params.set[i,]))
  out <- as.data.frame(ode(init, times, bubonicSIR, params))
  bSIR[i,length(parameters)+1] <- max(out$D_h)
  bSIR[i,length(parameters)+2] <- which.max(out$D_h)
  difference<-as.vector(NA)
  for(j in 2:length(out$D_h)){
    difference[j-1]<-out$D_h[j]-out$D_h[j-1]
  }
  bSIR[i,length(parameters)+3]<-ifelse(length(which(difference>1)[length(which(difference>1))])>0, which(difference>1)[length(which(difference>1))], 0) 
  bSIR[i,length(parameters)+4]<-ifelse(length(which(difference>0.1)[length(which(difference>0.1))])>0, which(difference>0.1)[length(which(difference>0.1))], 0) 
  bSIR[i,length(parameters)+5]<-ifelse(length(which(difference>0.01)[length(which(difference>0.01))])>0, which(difference>0.01)[length(which(difference>0.01))], 0) 
  
}

par(mfrow=c(1,2))
boxplot(bSIR$MaxInf~bSIR$beta_r, main= expression(paste("Effect of ", beta[r], " on Size")))
boxplot(bSIR$Thresh1~bSIR$beta_r, main= expression(paste("Effect of ", beta[r], " on Duration")))
boxplot(bSIR$MaxInf~bSIR$alpha, main= expression(paste("Effect of ", alpha, " on Size")))
boxplot(bSIR$Thresh1~bSIR$alpha, main= expression(paste("Effect of ", alpha, " on Duration")))
boxplot(bSIR$MaxInf~bSIR$gamma_r,  main= expression(paste("Effect of ", gamma[r], " on Size")))
boxplot(bSIR$Thresh1~bSIR$gamma_r, main= expression(paste("Effect of ", gamma[r], " on Duration")))
boxplot(bSIR$MaxInf~bSIR$g_r, main= expression(paste("Effect of ", g[r], " on Size")))
boxplot(bSIR$Thresh1~bSIR$g_r, main= expression(paste("Effect of ", g[r], " on Duration")))
boxplot(bSIR$MaxInf~bSIR$r_f, main= expression(paste("Effect of ", r[f], " on Size")))
boxplot(bSIR$Thresh1~bSIR$r_f, main= expression(paste("Effect of ", r[f], " on Duration")))
boxplot(bSIR$MaxInf~bSIR$K_f, main= expression(paste("Effect of ", K[f], " on Size")))
boxplot(bSIR$Thresh1~bSIR$K_f, main= expression(paste("Effect of ", K[f], " on Duration")))
boxplot(bSIR$MaxInf~bSIR$d_f, main= expression(paste("Effect of ", d[f], " on Size")))
boxplot(bSIR$Thresh1~bSIR$d_f, main= expression(paste("Effect of ", d[f], " on Duration")))
boxplot(bSIR$MaxInf~bSIR$beta_h, main= expression(paste("Effect of ", beta[h], " on Size")))
boxplot(bSIR$Thresh1~bSIR$beta_h, main= expression(paste("Effect of ", beta[h], " on Duration")))
boxplot(bSIR$MaxInf~bSIR$gamma_h, main= expression(paste("Effect of ", gamma[h], " on Size")))
boxplot(bSIR$Thresh1~bSIR$gamma_h, main= expression(paste("Effect of ", gamma[h], " on Duration")))
boxplot(bSIR$MaxInf~bSIR$g_h, main= expression(paste("Effect of ", g[h], " on Size")))
boxplot(bSIR$Thresh1~bSIR$g_h, main= expression(paste("Effect of ", g[h], " on Duration")))
boxplot(bSIR$MaxInf~bSIR$b_h, main= expression(paste("Effect of ", b[h], " on Size")))
boxplot(bSIR$Thresh1~bSIR$b_h, main= expression(paste("Effect of ", b[h], " on Duration")))
boxplot(bSIR$MaxInf~bSIR$d_h, main= expression(paste("Effect of ", d[h], " on Size")))
boxplot(bSIR$Thresh1~bSIR$d_h, main= expression(paste("Effect of ", d[h], " on Duration")))

par(mfrow=c(1,2))
boxplot(bSIR$MaxInf, main= "Outbreak Size", ylab= "Number of Dead Humans", ylim=c(0,500000))
boxplot(bSIR$Thresh1, main= "Outbreak Duration", ylab="Time (Days)")

bonferroni.alpha <- 0.05/5
prcc_size <- pcc(bSIR[,1:length(parameters)], bSIR[,length(parameters)+1], nboot = 1000, rank=TRUE, conf=1-bonferroni.alpha)
prcc_duration <- pcc(bSIR[,1:length(parameters)], bSIR[,length(parameters)+2], nboot = 1000, rank=TRUE, conf=1-bonferroni.alpha)

#plot correlation coefficients and confidence intervals for epidemic size and duration
size<-prcc_size$PRCC
size$param<-rownames(size)
colnames(size)[4:5] <- c("maxCI", "minCI")

ggplot(size, aes(x=param, y = original)) +
  geom_point(size = 4)+ 
  geom_errorbar(aes(ymax = maxCI, ymin = minCI))+
  ggtitle("PRCC Outbreak Size: Bubonic Plague SIR")


duration<-prcc_duration$PRCC
duration$param<-rownames(duration)
colnames(duration)[4:5] <- c("maxCI", "minCI")

ggplot(duration, aes(x=param, y = original)) +
  geom_point(size = 4)+ 
  geom_errorbar(aes(ymax = maxCI, ymin = minCI))+
  ggtitle("PRCC Outbreak Duration: Bubonic Plague SIR")


# Bubonic SEIR model -------------------------------------------------------


#Define initial conditions and expected parameter values
init <- c(S_r=499999, I_r=1, R_r=0, D_r=0, H=6, Fl=0, S_h = 500000, E_h=0, I_h = 0, R_h=0, D_h=0) #population size, and how many individuals start in each susceptible, infected, or removed category
parameters <- c(beta_r = 0.09, alpha=3/500000, gamma_r = 1/5.15, g_r=0.1, r_f=0.0084, K_f=6, d_f=1/5, beta_h=0.19, sigma_h= 1/4, gamma_h=1/10, g_h=0.34, b_h=1/(25*365), d_h=1/(25*365)) #you can play with transmission and recovery rates here

if(uniform==TRUE){ params.set<-LHS_bSEIR.uniform  #uniform LHS distribution from `LHSnonuniform.R`
} else{ params.set<-LHS_bSEIR} #non-uniform LHS distribution from `LHSnonuniform.R`

bSEIR <- data.frame(matrix(rep(NA,h),nrow=h, ncol= ncol(params.set)+5))
colnames(bSEIR)<-c(names(parameters), "MaxInf", "MaxDur", "Thresh1", "Thresh0.1", "Thresh0.01")

for(i in 1:h){
  bSEIR[i,1:length(parameters)] <- params <- as.list(c(params.set[i,]))
  out <- as.data.frame(ode(init, times, bubonicSEIR, params))
  bSEIR[i,length(parameters)+1] <- max(out$D_h)
  bSEIR[i,length(parameters)+2] <- which.max(out$D_h)
  difference<-as.vector(NA)
  for(j in 2:length(out$D_h)){
    difference[j-1]<-out$D_h[j]-out$D_h[j-1]
  }
  bSEIR[i,length(parameters)+3]<-ifelse(length(which(difference>1)[length(which(difference>1))])>0, which(difference>1)[length(which(difference>1))], 0) 
  bSEIR[i,length(parameters)+4]<-ifelse(length(which(difference>0.1)[length(which(difference>0.1))])>0, which(difference>0.1)[length(which(difference>0.1))], 0) 
  bSEIR[i,length(parameters)+5]<-ifelse(length(which(difference>0.01)[length(which(difference>0.01))])>0, which(difference>0.01)[length(which(difference>0.01))], 0) 
  
}

par(mfrow=c(1,2))
boxplot(bSEIR$MaxInf~bSEIR$beta_r, main= expression(paste("Effect of ", beta[r], " on Size")))
boxplot(bSEIR$Thresh1~bSEIR$beta_r, main= expression(paste("Effect of ", beta[r], " on Duration")))
boxplot(bSEIR$MaxInf~bSEIR$alpha, main= expression(paste("Effect of ", alpha, " on Size")))
boxplot(bSEIR$Thresh1~bSEIR$alpha, main= expression(paste("Effect of ", alpha, " on Duration")))
boxplot(bSEIR$MaxInf~bSEIR$gamma_r,  main= expression(paste("Effect of ", gamma[r], " on Size")))
boxplot(bSEIR$Thresh1~bSEIR$gamma_r, main= expression(paste("Effect of ", gamma[r], " on Duration")))
boxplot(bSEIR$MaxInf~bSEIR$g_r, main= expression(paste("Effect of ", g[r], " on Size")))
boxplot(bSEIR$Thresh1~bSEIR$g_r, main= expression(paste("Effect of ", g[r], " on Duration")))
boxplot(bSEIR$MaxInf~bSEIR$r_f, main= expression(paste("Effect of ", r[f], " on Size")))
boxplot(bSEIR$Thresh1~bSEIR$r_f, main= expression(paste("Effect of ", r[f], " on Duration")))
boxplot(bSEIR$MaxInf~bSEIR$K_f, main= expression(paste("Effect of ", K[f], " on Size")))
boxplot(bSEIR$Thresh1~bSEIR$K_f, main= expression(paste("Effect of ", K[f], " on Duration")))
boxplot(bSEIR$MaxInf~bSEIR$d_f, main= expression(paste("Effect of ", d[f], " on Size")))
boxplot(bSEIR$Thresh1~bSEIR$d_f, main= expression(paste("Effect of ", d[f], " on Duration")))
boxplot(bSEIR$MaxInf~bSEIR$beta_h, main= expression(paste("Effect of ", beta[h], " on Size")))
boxplot(bSEIR$Thresh1~bSEIR$beta_h, main= expression(paste("Effect of ", beta[h], " on Duration")))
boxplot(bSEIR$MaxInf~bSEIR$sigma_h, main= expression(paste("Effect of ", sigma[h], " on Size")))
boxplot(bSEIR$Thresh1~bSEIR$sigma_h, main= expression(paste("Effect of ", sigma[h], " on Duration")))
boxplot(bSEIR$MaxInf~bSEIR$gamma_h, main= expression(paste("Effect of ", gamma[h], " on Size")))
boxplot(bSEIR$Thresh1~bSEIR$gamma_h, main= expression(paste("Effect of ", gamma[h], " on Duration")))
boxplot(bSEIR$MaxInf~bSEIR$g_h, main= expression(paste("Effect of ", g[h], " on Size")))
boxplot(bSEIR$Thresh1~bSEIR$g_h, main= expression(paste("Effect of ", g[h], " on Duration")))
boxplot(bSEIR$MaxInf~bSEIR$b_h, main= expression(paste("Effect of ", b[h], " on Size")))
boxplot(bSEIR$Thresh1~bSEIR$b_h, main= expression(paste("Effect of ", b[h], " on Duration")))
boxplot(bSEIR$MaxInf~bSEIR$d_h, main= expression(paste("Effect of ", d[h], " on Size")))
boxplot(bSEIR$Thresh1~bSEIR$d_h, main= expression(paste("Effect of ", d[h], " on Duration")))

par(mfrow=c(1,2))
boxplot(bSEIR$MaxInf, main= "Outbreak Size", ylab= "Number of Dead Humans", ylim=c(0,500000))
boxplot(bSEIR$Thresh1, main= "Outbreak Duration", ylab="Time (Days)")

bonferroni.alpha <- 0.05/5
prcc_size <- pcc(bSEIR[,1:length(parameters)], bSEIR[,length(parameters)+1], nboot = 1000, rank=TRUE, conf=1-bonferroni.alpha)
prcc_duration <- pcc(bSEIR[,1:length(parameters)], bSEIR[,length(parameters)+2], nboot = 1000, rank=TRUE, conf=1-bonferroni.alpha)

#plot correlation coefficients and confidence intervals for epidemic size and duration
size<-prcc_size$PRCC
size$param<-rownames(size)
colnames(size)[4:5] <- c("maxCI", "minCI")

ggplot(size, aes(x=param, y = original)) +
  geom_point(size = 4)+ 
  geom_errorbar(aes(ymax = maxCI, ymin = minCI))+
  ggtitle("PRCC Outbreak Size: Bubonic Plague SEIR")
  

duration<-prcc_duration$PRCC
duration$param<-rownames(duration)
colnames(duration)[4:5] <- c("maxCI", "minCI")

ggplot(duration, aes(x=param, y = original)) +
  geom_point(size = 4)+ 
  geom_errorbar(aes(ymax = maxCI, ymin = minCI))+
  ggtitle("PRCC Outbreak Duration: Bubonic Plague SEIR")


# Bubonic/Pneumonic SEIR --------------------------------------------------

#Define initial conditions and expected parameter values
init <- c(S_r=499999, I_r=1, R_r=0, D_r=0, H=6, Fl=0, S_h = 500000, E_b=0, E_p=0, I_b = 0, I_p =0, R_h=0, D_h=0) #population size, and how many individuals start in each susceptible, infected, or removed category
parameters <- c(beta_r = 0.09, alpha=3/500000, gamma_r = 1/5.15, g_r=0.1, r_f=0.0084, K_f=6, d_f=1/5, beta_b=0.19, beta_p = 0.45, sigma_b= 1/6, sigma_p=1/4.3, gamma_b=1/10, gamma_p=1/2.5, p=0.2, g_h=0.34, b_h=1/(25*365), d_h=1/(25*365)) #you can play with transmission and recovery rates here

if(uniform==TRUE){ params.set<-LHS_bpSEIR.uniform  #uniform LHS distribution from `LHSnonuniform.R`
} else{ params.set<-LHS_bpSEIR} #non-uniform LHS distribution from `LHSnonuniform.R`

bpSEIR <- data.frame(matrix(rep(NA,h),nrow=h, ncol= ncol(params.set)+5))
colnames(bpSEIR)<-c(names(parameters), "MaxInf", "MaxDur", "Thresh1", "Thresh0.1", "Thresh0.01")

for(i in 1:h){
  bpSEIR[i,1:length(parameters)] <- params <- as.list(c(params.set[i,]))
  out <- as.data.frame(ode(init, times, bubonic_pneumonicSEIR, params))
  bpSEIR[i,length(parameters)+1] <- max(out$D_h)
  bpSEIR[i,length(parameters)+2] <- which.max(out$D_h)
  difference<-as.vector(NA)
  for(j in 2:length(out$D_h)){
    difference[j-1]<-out$D_h[j]-out$D_h[j-1]
  }
  bpSEIR[i,length(parameters)+3]<-ifelse(length(which(difference>1)[length(which(difference>1))])>0, which(difference>1)[length(which(difference>1))], 0) 
  bpSEIR[i,length(parameters)+4]<-ifelse(length(which(difference>0.1)[length(which(difference>0.1))])>0, which(difference>0.1)[length(which(difference>0.1))], 0) 
  bpSEIR[i,length(parameters)+5]<-ifelse(length(which(difference>0.01)[length(which(difference>0.01))])>0, which(difference>0.01)[length(which(difference>0.01))], 0) 
  
}

par(mfrow=c(1,2))
boxplot(bpSEIR$MaxInf~bpSEIR$beta_r, main= expression(paste("Effect of ", beta[r], " on Size")))
boxplot(bpSEIR$Thresh1~bpSEIR$beta_r, main= expression(paste("Effect of ", beta[r], " on Duration")))
boxplot(bpSEIR$MaxInf~bpSEIR$alpha, main= expression(paste("Effect of ", alpha, " on Size")))
boxplot(bpSEIR$Thresh1~bpSEIR$alpha, main= expression(paste("Effect of ", alpha, " on Duration")))
boxplot(bpSEIR$MaxInf~bpSEIR$gamma_r,  main= expression(paste("Effect of ", gamma[r], " on Size")))
boxplot(bpSEIR$Thresh1~bpSEIR$gamma_r, main= expression(paste("Effect of ", gamma[r], " on Duration")))
boxplot(bpSEIR$MaxInf~bpSEIR$g_r, main= expression(paste("Effect of ", g[r], " on Size")))
boxplot(bpSEIR$Thresh1~bpSEIR$g_r, main= expression(paste("Effect of ", g[r], " on Duration")))
boxplot(bpSEIR$MaxInf~bpSEIR$r_f, main= expression(paste("Effect of ", r[f], " on Size")))
boxplot(bpSEIR$Thresh1~bpSEIR$r_f, main= expression(paste("Effect of ", r[f], " on Duration")))
boxplot(bpSEIR$MaxInf~bpSEIR$K_f, main= expression(paste("Effect of ", K[f], " on Size")))
boxplot(bpSEIR$Thresh1~bpSEIR$K_f, main= expression(paste("Effect of ", K[f], " on Duration")))
boxplot(bpSEIR$MaxInf~bpSEIR$d_f, main= expression(paste("Effect of ", d[f], " on Size")))
boxplot(bpSEIR$Thresh1~bpSEIR$d_f, main= expression(paste("Effect of ", d[f], " on Duration")))
boxplot(bpSEIR$MaxInf~bpSEIR$beta_b, main= expression(paste("Effect of ", beta[b], " on Size")))
boxplot(bpSEIR$Thresh1~bpSEIR$beta_b, main= expression(paste("Effect of ", beta[b], " on Duration")))
boxplot(bpSEIR$MaxInf~bpSEIR$sigma_b, main= expression(paste("Effect of ", sigma[b], " on Size")))
boxplot(bpSEIR$Thresh1~bpSEIR$sigma_b, main= expression(paste("Effect of ", sigma[b], " on Duration")))
boxplot(bpSEIR$MaxInf~bpSEIR$gamma_b, main= expression(paste("Effect of ", gamma[b], " on Size")))
boxplot(bpSEIR$Thresh1~bpSEIR$gamma_b, main= expression(paste("Effect of ", gamma[b], " on Duration")))
boxplot(bpSEIR$MaxInf~bpSEIR$beta_p, main= expression(paste("Effect of ", beta[p], " on Size")))
boxplot(bpSEIR$Thresh1~bpSEIR$beta_p, main= expression(paste("Effect of ", beta[p], " on Duration")))
boxplot(bpSEIR$MaxInf~bpSEIR$sigma_p, main= expression(paste("Effect of ", sigma[p], " on Size")))
boxplot(bpSEIR$Thresh1~bpSEIR$sigma_p, main= expression(paste("Effect of ", sigma[p], " on Duration")))
boxplot(bpSEIR$MaxInf~bpSEIR$gamma_p,  main= expression(paste("Effect of ", gamma[p], " on Size")))
boxplot(bpSEIR$Thresh1~bpSEIR$gamma_p, main= expression(paste("Effect of ", gamma[p], " on Duration")))
boxplot(bpSEIR$MaxInf~bpSEIR$g_h, main= expression(paste("Effect of ", g[h], " on Size")))
boxplot(bpSEIR$Thresh1~bpSEIR$g_h, main= expression(paste("Effect of ", g[h], " on Duration")))
boxplot(bpSEIR$MaxInf~bpSEIR$p, main= expression(paste("Effect of ", p, " on Size")))
boxplot(bpSEIR$Thresh1~bpSEIR$p, main= expression(paste("Effect of ", p, " on Duration")))
boxplot(bpSEIR$MaxInf~bpSEIR$b_h, main= expression(paste("Effect of ", b[h], " on Size")))
boxplot(bpSEIR$Thresh1~bpSEIR$b_h, main= expression(paste("Effect of ", b[h], " on Duration")))
boxplot(bpSEIR$MaxInf~bpSEIR$d_h, main= expression(paste("Effect of ", d[h], " on Size")))
boxplot(bpSEIR$Thresh1~bpSEIR$d_h, main= expression(paste("Effect of ", d[h], " on Duration")))

par(mfrow=c(1,2))
boxplot(bpSEIR$MaxInf, main= "Outbreak Size", ylab= "Number of Dead Humans", ylim=c(0,500000))
boxplot(bpSEIR$Thresh1, main= "Outbreak Duration", ylab="Time (Days)")

bonferroni.alpha <- 0.05/5
prcc_size <- pcc(bpSEIR[,1:length(parameters)], bpSEIR[,length(parameters)+1], nboot = 1000, rank=TRUE, conf=1-bonferroni.alpha)
prcc_duration <- pcc(bpSEIR[,1:length(parameters)], bpSEIR[,length(parameters)+2], nboot = 1000, rank=TRUE, conf=1-bonferroni.alpha)

#plot correlation coefficients and confidence intervals for epidemic size and duration
size<-prcc_size$PRCC
size$param<-rownames(size)
colnames(size)[4:5] <- c("maxCI", "minCI")

ggplot(size, aes(x=param, y = original)) +
  geom_point(size = 4)+ 
  geom_errorbar(aes(ymax = maxCI, ymin = minCI))+
  ggtitle("PRCC Outbreak Size: Bubonic/Pneumonic Plague SEIR")

duration<-prcc_duration$PRCC
duration$param<-rownames(duration)
colnames(duration)[4:5] <- c("maxCI", "minCI")

ggplot(duration, aes(x=param, y = original)) +
  geom_point(size = 4)+ 
  geom_errorbar(aes(ymax = maxCI, ymin = minCI))+
  ggtitle("PRCC Outbreak Duration: Bubonic/Pneumonic Plague SEIR")

# Human Ectoparasite Model --------------------------------------------------
#Define initial conditions and parameter values
init <- c(S_l=15*499999, I_l=0, S_h = 499999, I_low=0, I_high=1, R_h=0, D_h=0) #population size, and how many individuals start in each susceptible, infected, or removed category
parameters <- c(r_l = 0.11, K_l=15, beta_low=0.04, beta_high=0.3, beta_l=0.05, gamma_lice = 1/3, gamma_low=1/2, gamma_high=1/8, g_h=0.34, b_h=1/(25*365), d_h=1/(25*365)) #you can play with transmission and recovery rates here

if(uniform==TRUE){ params.set<-LHS_eSIR.uniform  #uniform LHS distribution from `LHSnonuniform.R`
} else{ params.set<-LHS_eSIR} #non-uniform LHS distribution from `LHSnonuniform.R`


eSIR <- data.frame(matrix(rep(NA,h),nrow=h, ncol= ncol(params.set)+5))
colnames(eSIR)<-c(names(parameters), "MaxInf", "MaxDur", "Thresh1", "Thresh0.1", "Thresh0.01")

for(i in 1:h){
  eSIR[i,1:length(parameters)] <- params <- as.list(c(params.set[i,]))
  out <- as.data.frame(ode(init, times, liceSIR, params))
  eSIR[i,length(parameters)+1] <- max(out$D_h)
  eSIR[i,length(parameters)+2] <- which.max(out$D_h)
  difference<-as.vector(NA)
  for(j in 2:length(out$D_h)){
    difference[j-1]<-out$D_h[j]-out$D_h[j-1]
  }
  eSIR[i,length(parameters)+3]<-ifelse(length(which(difference>1)[length(which(difference>1))])>0, which(difference>1)[length(which(difference>1))], 0) 
  eSIR[i,length(parameters)+4]<-ifelse(length(which(difference>0.1)[length(which(difference>0.1))])>0, which(difference>0.1)[length(which(difference>0.1))], 0) 
  eSIR[i,length(parameters)+5]<-ifelse(length(which(difference>0.01)[length(which(difference>0.01))])>0, which(difference>0.01)[length(which(difference>0.01))], 0) 
  
}

par(mfrow=c(1,2))
boxplot(eSIR$MaxInf~eSIR$r_l, main= expression(paste("Effect of ", r[l], " on Size")))
boxplot(eSIR$Thresh1~eSIR$r_l, main= expression(paste("Effect of ", r[l], " on Duration")))
boxplot(eSIR$MaxInf~eSIR$K_l, main= expression(paste("Effect of ", K[l], " on Size")))
boxplot(eSIR$Thresh1~eSIR$K_l, main= expression(paste("Effect of ", K[l], " on Duration")))
boxplot(eSIR$MaxInf~eSIR$beta_low, main= expression(paste("Effect of ", beta[low], " on Size")))
boxplot(eSIR$Thresh1~eSIR$beta_low, main= expression(paste("Effect of ", beta[low], " on Duration")))
boxplot(eSIR$MaxInf~eSIR$beta_high, main= expression(paste("Effect of ", beta[high], " on Size")))
boxplot(eSIR$Thresh1~eSIR$beta_high, main= expression(paste("Effect of ", beta[high], " on Duration")))
boxplot(eSIR$MaxInf~eSIR$beta_l, main= expression(paste("Effect of ", beta[l], " on Size")))
boxplot(eSIR$Thresh1~eSIR$beta_l, main= expression(paste("Effect of ", beta[l], " on Duration")))
boxplot(eSIR$MaxInf~eSIR$gamma_lice,  main= expression(paste("Effect of ", gamma[lice], " on Size")))
boxplot(eSIR$Thresh1~eSIR$gamma_lice, main= expression(paste("Effect of ", gamma[lice], " on Duration")))
boxplot(eSIR$MaxInf~eSIR$gamma_low,  main= expression(paste("Effect of ", gamma[low], " on Size")))
boxplot(eSIR$Thresh1~eSIR$gamma_low, main= expression(paste("Effect of ", gamma[low], " on Duration")))
boxplot(eSIR$MaxInf~eSIR$gamma_high,  main= expression(paste("Effect of ", gamma[high], " on Size")))
boxplot(eSIR$Thresh1~eSIR$gamma_high, main= expression(paste("Effect of ", gamma[high], " on Duration")))
boxplot(eSIR$MaxInf~eSIR$g_h, main= expression(paste("Effect of ", g[h], " on Size")))
boxplot(eSIR$Thresh1~eSIR$g_h, main= expression(paste("Effect of ", g[h], " on Duration")))
boxplot(eSIR$MaxInf~eSIR$b_h, main= expression(paste("Effect of ", b[h], " on Size")))
boxplot(eSIR$Thresh1~eSIR$b_h, main= expression(paste("Effect of ", b[h], " on Duration")))
boxplot(eSIR$MaxInf~eSIR$d_h, main= expression(paste("Effect of ", d[h], " on Size")))
boxplot(eSIR$Thresh1~eSIR$d_h, main= expression(paste("Effect of ", d[h], " on Duration")))

par(mfrow=c(1,2))
boxplot(eSIR$MaxInf, main= "Outbreak Size", ylab= "Number of Dead Humans", ylim=c(0,500000))
boxplot(eSIR$Thresh1, main= "Outbreak Duration", ylab="Time (Days)")

bonferroni.alpha <- 0.05/5
prcc_size <- pcc(eSIR[,1:length(parameters)], eSIR[,length(parameters)+1], nboot = 1000, rank=TRUE, conf=1-bonferroni.alpha)
prcc_duration <- pcc(eSIR[,1:length(parameters)], eSIR[,length(parameters)+2], nboot = 1000, rank=TRUE, conf=1-bonferroni.alpha)

#plot correlation coefficients and confidence intervals for epidemic size and duration
size<-prcc_size$PRCC
size$param<-rownames(size)
colnames(size)[4:5] <- c("maxCI", "minCI")

ggplot(size, aes(x=param, y = original)) +
  geom_point(size = 4)+ 
  geom_errorbar(aes(ymax = maxCI, ymin = minCI))+
  ggtitle("PRCC Outbreak Size: Ectoparasite SIR")

duration<-prcc_duration$PRCC
duration$param<-rownames(duration)
colnames(duration)[4:5] <- c("maxCI", "minCI")

ggplot(duration, aes(x=param, y = original)) +
  geom_point(size = 4)+ 
  geom_errorbar(aes(ymax = maxCI, ymin = minCI))+
  ggtitle("PRCC Outbreak Duration: Ectoparasite SIR")

# Comparative Figure for All Models ---------------------------------------

comp_size<-data.frame(pSIR=pSIR$MaxInf, pSEIR= pSEIR$MaxInf, bSIR=bSIR$MaxInf, bSEIR= bSEIR$MaxInf, bpSEIR= bpSEIR$MaxInf, eSIR=eSIR$MaxInf)
comp_dur<- data.frame(pSIR=pSIR$Thresh1, pSEIR= pSEIR$Thresh1, bSIR=bSIR$Thresh1, bSEIR= bSEIR$Thresh1, bpSEIR= bpSEIR$Thresh1, eSIR=eSIR$Thresh1)

long_DFsize <- comp_size %>% gather(Model, NumberDead, c(pSIR, pSEIR, bSIR, bSEIR, bpSEIR, eSIR))

long_DFdur <- comp_dur %>% gather(Model, Duration, c(pSIR, pSEIR, bSIR, bSEIR, bpSEIR, eSIR))

ggplot(long_DFsize, aes(Model, NumberDead)) +  geom_boxplot()+ geom_jitter(alpha=0.5) +
  ylab("Number Dead")+
  scale_x_discrete(labels = c(bpSEIR="Bubonic/Pneumonic (SEIR)", bSIR="Bubonic Plague (SIR)", bSEIR="Bubonic Plague (SEIR)", pSEIR="Pneumonic Plague (SEIR)", pSIR="Pneumonic Plague (SIR)", eSIR="Ectoparasites"))

ggplot(long_DFdur, aes(Model, Duration)) +  geom_boxplot()+ geom_jitter(alpha=0.5) +
  ylab("Time (Days)")+
  scale_x_discrete(labels = c("bpSEIR"= "Bubonic/Pneumonic (SEIR)", "bSEIR"= "Bubonic Plague (SEIR)", "pSEIR"="Pneumonic Plague (SEIR)", "pSIR"="Pneumonic Plague (SIR)"))
