#' Global sensitivity analysis
#' @date April 30, 2019
#' @author Lauren White
#' @email lwhite@sesync.org

#' Adapted from: https://daphnia.ecology.uga.edu/drakelab/wp-content/uploads/2015/07/sensitivity-ebola.pdf
#' load ode functions for each of the models
source('~/JustinianPlague/Plague_model_functions.R')

#install.packages('lhs')
require(lhs) #add the lhs library
library(sensitivity)
require(ggplot2)
library(tidyverse)

h <- 1000 #choose number of parameter sets to sample to sample
set.seed(2718) #set random seed



# Sensitivity Analysis for Pneumonic SIR Model ----------------------------

#expected parameter values
parameters <- c(beta_p = 0.45, gamma_p = 1/2.5, b_h=1/(25*365), d_h=1/(25*365)) #you can play with transmission and recovery rates here
lhs<-maximinLHS(h, length(parameters)) #here 4 is the number of parameters required for the ODE system

beta_p.min <- 0.42
beta_p.max <- 0.48
gamma_p.min <- 1/3.7
gamma_p.max <- 1/1.3
b_h.min <- 1/(30*365)
b_h.max <- 1/(20*365)
d_h.min <- 1/(30*365)
d_h.max <- 1/(20*365)

params.set <- cbind( beta_p = lhs[,1]*(beta_p.max-beta_p.min)+beta_p.min,
                     gamma_p = lhs[,2]*(gamma_p.max-gamma_p.min)+gamma_p.min,
                     b_h = lhs[,3]*(b_h.max-b_h.min)+b_h.min,
                     d_h = lhs[,4]*(d_h.max-d_h.min)+d_h.min)


h2 <-250

init <- c(S_h = 499999, I_h = 1, D_h=0)
times <- seq(0, 2000, by= 1)

pSIR <- data.frame(matrix(rep(NA,h2),nrow=h2, ncol= ncol(params.set)+5))
colnames(pSIR)<-c(names(parameters), "MaxInf", "MaxDur", "Thresh1", "Thresh0.1", "Thresh0.01")
for(i in 1:h2){
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

#names(pSIR) <- c(names(params),'MaxInf', 'Thresh1')
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

#save(pSIR, file='pSIR.Rdata')
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

# Sensitivity Analysis for Pneumonic SEIR Model ----------------------------


#expected parameter values
init <- c(S_h = 499999, E_h= 0, I_h = 1, D_h=0) #population size, and how many individuals start in each susceptible, infected, or removed category
parameters <- c(beta_p = 0.45, sigma_p= 1/4.3, gamma_p = 1/2.5, b_h=1/(25*365), d_h=1/(25*365)) #you can play with transmission and recovery rates here
times <- seq(0, 5000, by= 1)
lhs<-maximinLHS(h,length(parameters)) #here 5 is the number of parameters required for the ODE system
h2 <-250

beta_p.min <- 0.42
beta_p.max <- 0.48
sigma_p.min <- 1/6.1
sigma_p.max <- 1/2.5
gamma_p.min <- 1/3.7
gamma_p.max <- 1/1.3
b_h.min <- 1/(30*365)
b_h.max <- 1/(20*365)
d_h.min <- 1/(30*365)
d_h.max <- 1/(20*365)

params.set <- cbind( beta_p = lhs[,1]*(beta_p.max-beta_p.min)+beta_p.min,
                     sigma_p = lhs[,2]*(sigma_p.max-sigma_p.min)+sigma_p.min,
                     gamma_p = lhs[,3]*(gamma_p.max-gamma_p.min)+gamma_p.min,
                     b_h = lhs[,4]*(b_h.max-b_h.min)+b_h.min,
                     d_h = lhs[,5]*(d_h.max-d_h.min)+d_h.min)


pSEIR <- data.frame(matrix(rep(NA,h2),nrow=h2, ncol= ncol(params.set)+5))
colnames(pSEIR)<-c(names(parameters), "MaxInf", "MaxDur", "Thresh1", "Thresh0.1", "Thresh0.01")

for(i in 1:h2){
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

lhs<-maximinLHS(h, length(parameters)) 

times <- seq(0, 1000, by= 1)
h2 <-250

beta_r.min <- 0.04
beta_r.max <- 0.14
alpha.min <- 0.39/500000
alpha.max <- 20/500000
gamma_r.min <- 1/4.71
gamma_r.max <- 1/5.59
g_r.min <- 0
g_r.max <- 0.37
r_f.min <- 0.0084
r_f.max <- 0.055
K_f.min <-3.29
K_f.max <-11.17
d_f.min <- 1/11.66
d_f.max <- 1/1
beta_h.min <- 0.18
beta_h.max <- 0.20
gamma_h.min <- 1/26
gamma_h.max <- 1/10
g_h.min <- 0.3
g_h.max <- 0.4
b_h.min <- 1/(30*365)
b_h.max <- 1/(20*365)
d_h.min <- 1/(30*365)
d_h.max <- 1/(20*365)


# create_parameters<-function(lhs_matrix, parameters){
#   min_max<-data.frame(min_val=paste0(names(parameters),".min"), max_val=paste0(names(parameters), ".max"))
#   for(i in 1:length(parameters)){
#     param_matrix[i]<-lhs_matrix[i]*(max_val[i]-min_val[i])+min_val[i]
#   }
# }

params.set <- cbind( beta_r = lhs[,1]*(beta_r.max-beta_r.min)+beta_r.min,
                     alpha = lhs[,2]*(alpha.max-alpha.min)+alpha.min,
                     gamma_r = lhs[,3]*(gamma_r.max-gamma_r.min)+gamma_r.min,
                     g_r = lhs[,4]*(g_r.max-g_r.min)+g_r.min,
                     r_f = lhs[,5]*(r_f.max-r_f.min)+r_f.min,
                     K_f = lhs[,6]*(K_f.max-K_f.min)+K_f.min,
                     d_f = lhs[,7]*(d_f.max-d_f.min)+d_f.min,
                     beta_h= lhs[,8]*(beta_h.max-beta_h.min)+beta_h.min,
                     gamma_h = lhs[,9]*(gamma_h.max-gamma_h.min)+gamma_h.min,
                     g_h = lhs[,10]*(g_h.max-g_h.min)+g_h.min,
                     b_h = lhs[,11]*(b_h.max-b_h.min)+b_h.min,
                     d_h = lhs[,12]*(d_h.max-d_h.min)+d_h.min)


bSIR <- data.frame(matrix(rep(NA,h2),nrow=h2, ncol= ncol(params.set)+5))
colnames(bSIR)<-c(names(parameters), "MaxInf", "MaxDur", "Thresh1", "Thresh0.1", "Thresh0.01")

for(i in 1:h2){
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

lhs<-maximinLHS(h, length(parameters)) #here 5 is the number of parameters required for the ODE system

times <- seq(0, 1000, by= 1)
h2 <-250

beta_r.min <- 0.04
beta_r.max <- 0.14
alpha.min <- 0.39/500000
alpha.max <- 20/500000
gamma_r.min <- 1/4.71
gamma_r.max <- 1/5.59
g_r.min <- 0
g_r.max <- 0.37
r_f.min <- 0.0084
r_f.max <- 0.055
K_f.min <-3.29
K_f.max <-11.17
d_f.min <- 1/11.66
d_f.max <- 1/1
beta_h.min <- 0.18
beta_h.max <- 0.20
sigma_h.min <- 1/6
sigma_h.max <- 1/2
gamma_h.min <- 1/26
gamma_h.max <- 1/10
g_h.min <- 0.3
g_h.max <- 0.4
b_h.min <- 1/(30*365)
b_h.max <- 1/(20*365)
d_h.min <- 1/(30*365)
d_h.max <- 1/(20*365)


# create_parameters<-function(lhs_matrix, parameters){
#   min_max<-data.frame(min_val=paste0(names(parameters),".min"), max_val=paste0(names(parameters), ".max"))
#   for(i in 1:length(parameters)){
#     param_matrix[i]<-lhs_matrix[i]*(max_val[i]-min_val[i])+min_val[i]
#   }
# }

params.set <- cbind( beta_r = lhs[,1]*(beta_r.max-beta_r.min)+beta_r.min,
                     alpha = lhs[,2]*(alpha.max-alpha.min)+alpha.min,
                     gamma_r = lhs[,3]*(gamma_r.max-gamma_r.min)+gamma_r.min,
                     g_r = lhs[,4]*(g_r.max-g_r.min)+g_r.min,
                     r_f = lhs[,5]*(r_f.max-r_f.min)+r_f.min,
                     K_f = lhs[,6]*(K_f.max-K_f.min)+K_f.min,
                     d_f = lhs[,7]*(d_f.max-d_f.min)+d_f.min,
                     beta_h= lhs[,8]*(beta_h.max-beta_h.min)+beta_h.min,
                     sigma_h = lhs[,9]*(sigma_h.max-sigma_h.min)+sigma_h.min,
                     gamma_h = lhs[,10]*(gamma_h.max-gamma_h.min)+gamma_h.min,
                     g_h = lhs[,11]*(g_h.max-g_h.min)+g_h.min,
                     b_h = lhs[,12]*(b_h.max-b_h.min)+b_h.min,
                     d_h = lhs[,13]*(d_h.max-d_h.min)+d_h.min)


bSEIR <- data.frame(matrix(rep(NA,h2),nrow=h2, ncol= ncol(params.set)+5))
colnames(bSEIR)<-c(names(parameters), "MaxInf", "MaxDur", "Thresh1", "Thresh0.1", "Thresh0.01")

for(i in 1:h2){
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

lhs<-maximinLHS(h, length(parameters)) 

times <- seq(0, 1000, by= 1)
h2 <-1000

beta_r.min <- 0.04
beta_r.max <- 0.14
alpha.min <- 0.39/500000
alpha.max <- 20/500000
gamma_r.min <- 1/4.71
gamma_r.max <- 1/5.59
g_r.min <- 0
g_r.max <- 0.37
r_f.min <- 0.0084
r_f.max <- 0.055
K_f.min <-3.29
K_f.max <-11.17
d_f.min <- 1/11.66
d_f.max <- 1/1
beta_b.min <- 0.18
beta_b.max <- 0.20
sigma_b.min <- 1/6
sigma_b.max <- 1/2
gamma_b.min <- 1/26
gamma_b.max <- 1/10
beta_p.min <- 0.42
beta_p.max <- 0.48
sigma_p.min <- 1/6.1
sigma_p.max <- 1/2.5
gamma_p.min <- 1/3.7
gamma_p.max <- 1/1.3
p.min<- 0
p.max<- 0.4
g_h.min <- 0.3
g_h.max <- 0.4
b_h.min <- 1/(30*365)
b_h.max <- 1/(20*365)
d_h.min <- 1/(30*365)
d_h.max <- 1/(20*365)

# create_parameters<-function(lhs_matrix, parameters){
#   min_max<-data.frame(min_val=paste0(names(parameters),".min"), max_val=paste0(names(parameters), ".max"))
#   for(i in 1:length(parameters)){
#     param_matrix[i]<-lhs_matrix[i]*(max_val[i]-min_val[i])+min_val[i]
#   }
# }

params.set <- cbind( beta_r = lhs[,1]*(beta_r.max-beta_r.min)+beta_r.min,
                     alpha = lhs[,2]*(alpha.max-alpha.min)+alpha.min,
                     gamma_r = lhs[,3]*(gamma_r.max-gamma_r.min)+gamma_r.min,
                     g_r = lhs[,4]*(g_r.max-g_r.min)+g_r.min,
                     r_f = lhs[,5]*(r_f.max-r_f.min)+r_f.min,
                     K_f = lhs[,6]*(K_f.max-K_f.min)+K_f.min,
                     d_f = lhs[,7]*(d_f.max-d_f.min)+d_f.min,
                     beta_b= lhs[,8]*(beta_b.max-beta_b.min)+beta_b.min,
                     sigma_b = lhs[,9]*(sigma_b.max-sigma_b.min)+sigma_b.min,
                     gamma_b = lhs[,10]*(gamma_b.max-gamma_b.min)+gamma_b.min,
                     beta_p = lhs[,11]*(beta_p.max-beta_p.min)+beta_p.min,
                     sigma_p = lhs[,12]*(sigma_p.max-sigma_p.min)+sigma_p.min,
                     gamma_p = lhs[,13]*(gamma_p.max-gamma_p.min)+gamma_p.min,
                     g_h = lhs[,14]*(g_h.max-g_h.min)+g_h.min,
                     p = lhs[,15]*(p.max-p.min)+p.min,
                     b_h = lhs[,16]*(b_h.max-b_h.min)+b_h.min,
                     d_h = lhs[,17]*(d_h.max-d_h.min)+d_h.min)


bpSEIR <- data.frame(matrix(rep(NA,h2),nrow=h2, ncol= ncol(params.set)+5))
colnames(bpSEIR)<-c(names(parameters), "MaxInf", "MaxDur", "Thresh1", "Thresh0.1", "Thresh0.01")

for(i in 1:h2){
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
lhs<-maximinLHS(h, length(parameters)) 

times <- seq(0, 1000, by= 1)
h2 <-1000

r_l.min<-0.10
r_l.max<-0.12
K_l.min<-10.5
K_l.max<-67.7
beta_low.min<- 0
beta_low.max<- 0.05
beta_high.min<-0
beta_high.max<-1
beta_l.min<-0
beta_l.max<-0.1
gamma_lice.min<-2
gamma_lice.max<-4
gamma_low.min<-0
gamma_low.max<-4
gamma_high.min<-0
gamma_high.max<-4
g_h.min <- 0.3
g_h.max <- 0.4
b_h.min <- 1/(30*365)
b_h.max <- 1/(20*365)
d_h.min <- 1/(30*365)
d_h.max <- 1/(20*365)

# create_parameters<-function(lhs_matrix, parameters){
#   min_max<-data.frame(min_val=paste0(names(parameters),".min"), max_val=paste0(names(parameters), ".max"))
#   for(i in 1:length(parameters)){
#     param_matrix[i]<-lhs_matrix[i]*(max_val[i]-min_val[i])+min_val[i]
#   }
# }

params.set <- cbind( r_l = lhs[,1]*(r_l.max-r_l.min)+r_l.min,
                     K_l = lhs[,2]*(K_l.max-K_l.min)+K_l.min,
                     beta_low = lhs[,3]*(beta_low.max-beta_low.min)+beta_low.min,
                     beta_high = lhs[,4]*(beta_high.max-beta_high.min)+beta_high.min,
                     beta_l = lhs[,5]*(beta_l.max-beta_l.min)+beta_l.min,
                     gamma_lice = lhs[,6]*(gamma_lice.max-gamma_lice.min)+gamma_lice.min,
                     gamma_low = lhs[,7]*(gamma_low.max-gamma_low.min)+gamma_low.min,
                     gamma_high= lhs[,8]*(gamma_high.max-gamma_high.min)+gamma_high.min,
                     g_h = lhs[,9]*(g_h.max-g_h.min)+g_h.min,
                     b_h = lhs[,10]*(b_h.max-b_h.min)+b_h.min,
                     d_h = lhs[,11]*(d_h.max-d_h.min)+d_h.min)


eSIR <- data.frame(matrix(rep(NA,h2),nrow=h2, ncol= ncol(params.set)+5))
colnames(eSIR)<-c(names(parameters), "MaxInf", "MaxDur", "Thresh1", "Thresh0.1", "Thresh0.01")

for(i in 1:h2){
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
  ggtitle("PRCC Outbreak Size: Bubonic/Pneumonic Plague SEIR")

duration<-prcc_duration$PRCC
duration$param<-rownames(duration)
colnames(duration)[4:5] <- c("maxCI", "minCI")

ggplot(duration, aes(x=param, y = original)) +
  geom_point(size = 4)+ 
  geom_errorbar(aes(ymax = maxCI, ymin = minCI))+
  ggtitle("PRCC Outbreak Duration: Bubonic/Pneumonic Plague SEIR")

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
