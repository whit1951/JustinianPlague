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


h <- 1000 #choose number of parameter sets to sample to sample
set.seed(6242015) #set random seed



# Sensitivity Analysis for Pneumonic SIR Model ----------------------------

#expected parameter values
parameters <- c(beta_p = 0.45, gamma_p = 1/2.5, b_h=1/(25*365), d_h=1/(25*365)) #you can play with transmission and recovery rates here
lhs<-maximinLHS(h,4) #here 4 is the number of parameters required for the ODE system

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
times <- seq(0, 1000, by= 1)

pSIR <- data.frame(matrix(rep(NA,h2),nrow=h2, ncol= ncol(params.set)+2))
for(i in 1:h2){
  pSIR[i,1:4] <- params <- as.list(c(params.set[i,]))
  out <- as.data.frame(ode(init, times, pneumonicSIR, params))
  pSIR[i,5] <- max(out$D_h)
  pSIR[i,6] <- which.max(out$D_h)
}
difference<-as.vector(NA)
for(i in 2:length(out$D_h)){
  difference[i-1]<-out$D_h[i]-out$D_h[i-1]
}
which(difference<0.0001)[1]
names(pSIR) <- c(names(params),'outbreak.size', 'outbreak.duration')
par(mfrow=c(1,2))
boxplot(pSIR$outbreak.size~pSIR$beta_p, main= expression(paste("Effect of ", beta[p], " on Size")))
boxplot(pSIR$outbreak.duration~pSIR$beta_p, main= expression(paste("Effect of ", beta[p], " on Duration")))
boxplot(pSIR$outbreak.size~pSIR$gamma_p,  main= expression(paste("Effect of ", gamma[p], " on Size")))
boxplot(pSIR$outbreak.duration~pSIR$gamma_p, main= expression(paste("Effect of ", gamma[p], " on Duration")))
boxplot(pSIR$outbreak.size~pSIR$b_h, main= expression(paste("Effect of ", b[h], " on Size")))
boxplot(pSIR$outbreak.duration~pSIR$b_h, main= expression(paste("Effect of ", b[h], " on Duration")))
boxplot(pSIR$outbreak.size~pSIR$d_h, main= expression(paste("Effect of ", d[h], " on Size")))
boxplot(pSIR$outbreak.duration~pSIR$d_h, main= expression(paste("Effect of ", d[h], " on Duration")))

par(mfrow=c(1,2))
boxplot(pSIR$outbreak.size, main= "Outbreak Size", ylab= "Number of Dead Humans", ylim=c(0,500000))
boxplot(pSIR$outbreak.duration, main= "Outbreak Duration", ylab="Time (Days)")

#save(pSIR, file='pSIR.Rdata')
bonferroni.alpha <- 0.05/5
prcc_size <- pcc(pSIR[,1:4], pSIR[,5], nboot = 1000, rank=TRUE, conf=1-bonferroni.alpha)
prcc_duration <- pcc(pSIR[,1:4], pSIR[,6], nboot = 1000, rank=TRUE, conf=1-bonferroni.alpha)

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

lhs<-maximinLHS(h,5) #here 5 is the number of parameters required for the ODE system

#expected parameter values
init <- c(S_h = 499999, E_h= 0, I_h = 1, D_h=0) #population size, and how many individuals start in each susceptible, infected, or removed category
parameters <- c(beta_p = 0.45, sigma_p= 1/4.3, gamma_p = 1/2.5, b_h=1/(25*365), d_h=1/(25*365)) #you can play with transmission and recovery rates here
times <- seq(0, 1000, by= 1)
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


pSEIR <- data.frame(matrix(rep(NA,h2),nrow=h2, ncol= ncol(params.set)+2))
for(i in 1:h2){
  pSEIR[i,1:5] <- params <- as.list(c(params.set[i,]))
  out <- as.data.frame(ode(init, times, pneumonicSEIR, params))
  pSEIR[i,6] <- max(out$D_h)
  pSEIR[i,7] <- which.max(out$D_h)
}
names(pSEIR) <- c(names(params),'outbreak.size', 'outbreak.duration')
par(mfrow=c(1,2))
boxplot(pSEIR$outbreak.size~pSEIR$beta_p, main= expression(paste("Effect of ", beta[p], " on Size")))
boxplot(pSEIR$outbreak.duration~pSEIR$beta_p, main= expression(paste("Effect of ", beta[p], " on Duration")))
boxplot(pSEIR$outbreak.size~pSEIR$sigma_p, main= expression(paste("Effect of ", sigma[p], " on Size")))
boxplot(pSEIR$outbreak.duration~pSEIR$sigma_p, main= expression(paste("Effect of ", sigma[p], " on Duration")))
boxplot(pSEIR$outbreak.size~pSEIR$gamma_p,  main= expression(paste("Effect of ", gamma[p], " on Size")))
boxplot(pSEIR$outbreak.duration~pSEIR$gamma_p, main= expression(paste("Effect of ", gamma[p], " on Duration")))
boxplot(pSEIR$outbreak.size~pSEIR$b_h, main= expression(paste("Effect of ", b[h], " on Size")))
boxplot(pSEIR$outbreak.duration~pSEIR$b_h, main= expression(paste("Effect of ", b[h], " on Duration")))
boxplot(pSEIR$outbreak.size~pSEIR$d_h, main= expression(paste("Effect of ", d[h], " on Size")))
boxplot(pSEIR$outbreak.duration~pSEIR$d_h, main= expression(paste("Effect of ", d[h], " on Duration")))

par(mfrow=c(1,2))
boxplot(pSEIR$outbreak.size, main= "Outbreak Size", ylab= "Number of Dead Humans", ylim=c(0,500000))
boxplot(pSEIR$outbreak.duration, main= "Outbreak Duration", ylab="Time (Days)")

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


# Bubonic SEIR model -------------------------------------------------------


#Define initial conditions and expected parameter values
init <- c(S_r=499999, I_r=1, R_r=0, D_r=0, H=6, Fl=0, S_h = 500000, E_h=0, I_h = 0, R_h=0, D_h=0) #population size, and how many individuals start in each susceptible, infected, or removed category
parameters <- c(beta_r = 0.09, alpha=3/500000, gamma_r = 1/5.15, g_r=0.1, r_f=0.0084, K_f=6, d_f=1/5, beta_h=0.19, sigma_h= 1/4, gamma_h=1/10, g_h=0.34) #you can play with transmission and recovery rates here

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
                     g_h = lhs[,11]*(g_h.max-g_h.min)+g_h.min)


bSEIR <- data.frame(matrix(rep(NA,h2),nrow=h2, ncol= ncol(params.set)+2))
for(i in 1:h2){
  bSEIR[i,1:length(parameters)] <- params <- as.list(c(params.set[i,]))
  out <- as.data.frame(ode(init, times, ratfleaSIR, params))
  bSEIR[i,length(parameters)+1] <- max(out$D_h)
  bSEIR[i,length(parameters)+2] <- which.max(out$D_h)
}
names(bSEIR) <- c(names(params),'outbreak.size', 'outbreak.duration')
par(mfrow=c(1,2))
boxplot(bSEIR$outbreak.size~bSEIR$beta_r, main= expression(paste("Effect of ", beta[r], " on Size")))
boxplot(bSEIR$outbreak.duration~bSEIR$beta_r, main= expression(paste("Effect of ", beta[r], " on Duration")))
boxplot(bSEIR$outbreak.size~bSEIR$alpha, main= expression(paste("Effect of ", alpha, " on Size")))
boxplot(bSEIR$outbreak.duration~bSEIR$alpha, main= expression(paste("Effect of ", alpha, " on Duration")))
boxplot(bSEIR$outbreak.size~bSEIR$gamma_r,  main= expression(paste("Effect of ", gamma[r], " on Size")))
boxplot(bSEIR$outbreak.duration~bSEIR$gamma_r, main= expression(paste("Effect of ", gamma[r], " on Duration")))
boxplot(bSEIR$outbreak.size~bSEIR$g_r, main= expression(paste("Effect of ", g[r], " on Size")))
boxplot(bSEIR$outbreak.duration~bSEIR$g_r, main= expression(paste("Effect of ", g[r], " on Duration")))
boxplot(bSEIR$outbreak.size~bSEIR$r_f, main= expression(paste("Effect of ", r[f], " on Size")))
boxplot(bSEIR$outbreak.duration~bSEIR$r_f, main= expression(paste("Effect of ", r[f], " on Duration")))
boxplot(bSEIR$outbreak.size~bSEIR$K_f, main= expression(paste("Effect of ", K[f], " on Size")))
boxplot(bSEIR$outbreak.duration~bSEIR$K_f, main= expression(paste("Effect of ", K[f], " on Duration")))
boxplot(bSEIR$outbreak.size~bSEIR$d_f, main= expression(paste("Effect of ", d[f], " on Size")))
boxplot(bSEIR$outbreak.duration~bSEIR$d_f, main= expression(paste("Effect of ", d[f], " on Duration")))
boxplot(bSEIR$outbreak.size~bSEIR$beta_h, main= expression(paste("Effect of ", beta[h], " on Size")))
boxplot(bSEIR$outbreak.duration~bSEIR$beta_h, main= expression(paste("Effect of ", beta[h], " on Duration")))
boxplot(bSEIR$outbreak.size~bSEIR$sigma_h, main= expression(paste("Effect of ", sigma[h], " on Size")))
boxplot(bSEIR$outbreak.duration~bSEIR$sigma_h, main= expression(paste("Effect of ", sigma[h], " on Duration")))
boxplot(bSEIR$outbreak.size~bSEIR$gamma_h, main= expression(paste("Effect of ", gamma[h], " on Size")))
boxplot(bSEIR$outbreak.duration~bSEIR$gamma_h, main= expression(paste("Effect of ", gamma[h], " on Duration")))
boxplot(bSEIR$outbreak.size~bSEIR$g_h, main= expression(paste("Effect of ", g[h], " on Size")))
boxplot(bSEIR$outbreak.duration~bSEIR$g_h, main= expression(paste("Effect of ", g[h], " on Duration")))

par(mfrow=c(1,2))
boxplot(bSEIR$outbreak.size, main= "Outbreak Size", ylab= "Number of Dead Humans", ylim=c(0,500000))
boxplot(bSEIR$outbreak.duration, main= "Outbreak Duration", ylab="Time (Days)")

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
parameters <- c(beta_r = 0.09, alpha=3/500000, gamma_r = 1/5.15, g_r=0.1, r_f=0.0084, K_f=6, d_f=1/5, beta_b=0.19, beta_p = 0.45, sigma_b= 1/6, sigma_p=1/4.3, gamma_b=1/10, gamma_p=1/2.5, p=0.2, g_h=0.34) #you can play with transmission and recovery rates here

lhs<-maximinLHS(h, length(parameters)) #here 5 is the number of parameters required for the ODE system

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
                     p = lhs[,15]*(p.max-p.min)+p.min)


bpSEIR <- data.frame(matrix(rep(NA,h2),nrow=h2, ncol= ncol(params.set)+2))
for(i in 1:h2){
  bpSEIR[i,1:length(parameters)] <- params <- as.list(c(params.set[i,]))
  out <- as.data.frame(ode(init, times, bubonic_pneumonicSEIR, params))
  bpSEIR[i,length(parameters)+1] <- max(out$D_h)
  bpSEIR[i,length(parameters)+2] <- which.max(out$D_h)
}
names(bpSEIR) <- c(names(params),'outbreak.size', 'outbreak.duration')
par(mfrow=c(1,2))
boxplot(bpSEIR$outbreak.size~bpSEIR$beta_r, main= expression(paste("Effect of ", beta[r], " on Size")))
boxplot(bpSEIR$outbreak.duration~bpSEIR$beta_r, main= expression(paste("Effect of ", beta[r], " on Duration")))
boxplot(bpSEIR$outbreak.size~bpSEIR$alpha, main= expression(paste("Effect of ", alpha, " on Size")))
boxplot(bpSEIR$outbreak.duration~bpSEIR$alpha, main= expression(paste("Effect of ", alpha, " on Duration")))
boxplot(bpSEIR$outbreak.size~bpSEIR$gamma_r,  main= expression(paste("Effect of ", gamma[r], " on Size")))
boxplot(bpSEIR$outbreak.duration~bpSEIR$gamma_r, main= expression(paste("Effect of ", gamma[r], " on Duration")))
boxplot(bpSEIR$outbreak.size~bpSEIR$g_r, main= expression(paste("Effect of ", g[r], " on Size")))
boxplot(bpSEIR$outbreak.duration~bpSEIR$g_r, main= expression(paste("Effect of ", g[r], " on Duration")))
boxplot(bpSEIR$outbreak.size~bpSEIR$r_f, main= expression(paste("Effect of ", r[f], " on Size")))
boxplot(bpSEIR$outbreak.duration~bpSEIR$r_f, main= expression(paste("Effect of ", r[f], " on Duration")))
boxplot(bpSEIR$outbreak.size~bpSEIR$K_f, main= expression(paste("Effect of ", K[f], " on Size")))
boxplot(bpSEIR$outbreak.duration~bpSEIR$K_f, main= expression(paste("Effect of ", K[f], " on Duration")))
boxplot(bpSEIR$outbreak.size~bpSEIR$d_f, main= expression(paste("Effect of ", d[f], " on Size")))
boxplot(bpSEIR$outbreak.duration~bpSEIR$d_f, main= expression(paste("Effect of ", d[f], " on Duration")))
boxplot(bpSEIR$outbreak.size~bpSEIR$beta_b, main= expression(paste("Effect of ", beta[b], " on Size")))
boxplot(bpSEIR$outbreak.duration~bpSEIR$beta_b, main= expression(paste("Effect of ", beta[b], " on Duration")))
boxplot(bpSEIR$outbreak.size~bpSEIR$sigma_b, main= expression(paste("Effect of ", sigma[b], " on Size")))
boxplot(bpSEIR$outbreak.duration~bpSEIR$sigma_b, main= expression(paste("Effect of ", sigma[b], " on Duration")))
boxplot(bpSEIR$outbreak.size~bpSEIR$gamma_b, main= expression(paste("Effect of ", gamma[b], " on Size")))
boxplot(bpSEIR$outbreak.duration~bpSEIR$gamma_b, main= expression(paste("Effect of ", gamma[b], " on Duration")))
boxplot(bpSEIR$outbreak.size~bpSEIR$beta_p, main= expression(paste("Effect of ", beta[p], " on Size")))
boxplot(bpSEIR$outbreak.duration~bpSEIR$beta_p, main= expression(paste("Effect of ", beta[p], " on Duration")))
boxplot(bpSEIR$outbreak.size~bpSEIR$sigma_p, main= expression(paste("Effect of ", sigma[p], " on Size")))
boxplot(bpSEIR$outbreak.duration~bpSEIR$sigma_p, main= expression(paste("Effect of ", sigma[p], " on Duration")))
boxplot(bpSEIR$outbreak.size~bpSEIR$gamma_p,  main= expression(paste("Effect of ", gamma[p], " on Size")))
boxplot(bpSEIR$outbreak.duration~bpSEIR$gamma_p, main= expression(paste("Effect of ", gamma[p], " on Duration")))
boxplot(bpSEIR$outbreak.size~bpSEIR$g_h, main= expression(paste("Effect of ", g[h], " on Size")))
boxplot(bpSEIR$outbreak.duration~bpSEIR$g_h, main= expression(paste("Effect of ", g[h], " on Duration")))
boxplot(bpSEIR$outbreak.size~bpSEIR$p, main= expression(paste("Effect of ", p, " on Size")))
boxplot(bpSEIR$outbreak.duration~bpSEIR$p, main= expression(paste("Effect of ", p, " on Duration")))

par(mfrow=c(1,2))
boxplot(bpSEIR$outbreak.size, main= "Outbreak Size", ylab= "Number of Dead Humans", ylim=c(0,500000))
boxplot(bpSEIR$outbreak.duration, main= "Outbreak Duration", ylab="Time (Days)")

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


# Comparative Figure for All Models ---------------------------------------

comp_size<-data.frame(pSIR=pSIR$outbreak.size, pSEIR= pSEIR$outbreak.size , bSEIR= bSEIR$outbreak.size, bpSEIR= bpSEIR$outbreak.size)
comp_dur<- data.frame(pSIR=pSIR$outbreak.duration, pSEIR= pSEIR$outbreak.duration , bSEIR= bSEIR$outbreak.duration, bpSEIR= bpSEIR$outbreak.duration)

long_DFsize <- comp_size %>% gather(Model, NumberDead, c(pSIR, pSEIR, bSEIR, bpSEIR))

long_DFdur <- comp_dur %>% gather(Model, Duration, c(pSIR, pSEIR, bSEIR, bpSEIR))

ggplot(long_DFsize, aes(Model, NumberDead)) + geom_jitter() +
  ylab("Number Dead")
  #scale_x_discrete(labels = c("Bubonic/Pneumonic (SEIR)", "Bubonic Plague (SIR)", "Pneumonic Plague (SEIR)", "Pneumonic Plague (SIR)")

ggplot(long_DFdur, aes(Model, Duration)) + geom_jitter() +
  ylab("Time (Days)")+
  scale_x_discrete(labels = c("bpSEIR"= "Bubonic/Pneumonic (SEIR)", "bSEIR"= "Bubonic Plague (SEIR)", "pSEIR"="Pneumonic Plague (SEIR)", "pSIR"="Pneumonic Plague (SIR)"))
