#' Global sensitivity analysis
#' @date April 30, 2019
#' @author Lauren White
#' @email lwhite@sesync.org


#### Sensitivity Analysis for Pneumonic SIR Model ####----------------------------

#expected parameter values
parameters <- c(beta_p = 0.45, gamma_p = 1/2.5, b_h=1/(25*365), d_h=1/(25*365)) #you can play with transmission and recovery rates here

if(uniform==TRUE){ params.set<-LHS_pSIR.uniform  #uniform LHS distribution from `LHSnonuniform.R`
} else{ params.set<-LHS_pSIR} #non-uniform LHS distribution from `LHSnonuniform.R`

init <- c(S_h = 499999, I_h = 1, D_h=0)
pSIR <- data.frame(matrix(rep(NA,h),nrow=h, ncol= ncol(params.set)+6))
colnames(pSIR)<-c(names(parameters), "MaxInf", "MaxDur", "Thresh100", "Thresh10", "Thresh1", "Thresh250")

for(i in 1:h){
  pSIR[i,1:4] <- params <- as.list(c(params.set[i,]))
  out <- as.data.frame(ode(init, times, pneumonicSIR, params))
  pSIR[i,5] <- max(out$D_h)
  pSIR[i,6] <- which.max(out$D_h)
  
  difference<-as.vector(NA)
  for(j in 2:length(out$D_h)){
    difference[j-1]<-out$D_h[j]-out$D_h[j-1]
  }
  pSIR[i,7]<-ifelse(length(which(difference>100)[length(which(difference>100))])>0, which(difference>100)[length(which(difference>100))], 0) 
  pSIR[i,8]<-ifelse(length(which(difference>10)[length(which(difference>10))])>0, which(difference>10)[length(which(difference>10))], 0) 
  pSIR[i,9]<-ifelse(length(which(difference>1)[length(which(difference>1))])>0, which(difference>1)[length(which(difference>1))], 0) 
  pSIR[i,10]<-ifelse(length(which(difference>250)[length(which(difference>250))])>0, which(difference>250)[length(which(difference>250))], 0) 
  
  }




#### Sensitivity Analysis for Pneumonic SEIR Model #### ----------------------------

#expected parameter values
init <- c(S_h = 499999, E_h= 0, I_h = 1, D_h=0) #population size, and how many individuals start in each susceptible, infected, or removed category
parameters <- c(beta_p = 0.45, sigma_p= 1/4.3, gamma_p = 1/2.5, b_h=1/(25*365), d_h=1/(25*365)) #you can play with transmission and recovery rates here

if(uniform==TRUE){ params.set<-LHS_pSEIR.uniform  #uniform LHS distribution from `LHSnonuniform.R`
} else{ params.set<-LHS_pSEIR} #non-uniform LHS distribution from `LHSnonuniform.R`

pSEIR <- data.frame(matrix(rep(NA,h),nrow=h, ncol= ncol(params.set)+6))
colnames(pSEIR)<-c(names(parameters), "MaxInf", "MaxDur", "Thresh100", "Thresh10", "Thresh1", "Thresh250")

for(i in 1:h){
  pSEIR[i,1:5] <- params <- as.list(c(params.set[i,]))
  out <- as.data.frame(ode(init, times, pneumonicSEIR, params))
  pSEIR[i,6] <- max(out$D_h)
  pSEIR[i,7] <- which.max(out$D_h)
  
  difference<-as.vector(NA)
  for(j in 2:length(out$D_h)){
    difference[j-1]<-out$D_h[j]-out$D_h[j-1]
  }
  pSEIR[i,8]<-ifelse(length(which(difference>100)[length(which(difference>100))])>0, which(difference>100)[length(which(difference>100))], 0) 
  pSEIR[i,9]<-ifelse(length(which(difference>10)[length(which(difference>10))])>0, which(difference>10)[length(which(difference>10))], 0) 
  pSEIR[i,10]<-ifelse(length(which(difference>1)[length(which(difference>1))])>0, which(difference>1)[length(which(difference>1))], 0) 
  pSEIR[i,11]<-ifelse(length(which(difference>250)[length(which(difference>250))])>0, which(difference>250)[length(which(difference>250))], 0) 
  
}


# Bubonic SIR model -------------------------------------------------------

#Define initial conditions and expected parameter values
init <- c(S_r=499999, I_r=1, R_r=0, D_r=0, H=6, Fl=0, S_h = 500000, I_h = 0, R_h=0, D_h=0) #population size, and how many individuals start in each susceptible, infected, or removed category
parameters <- c(beta_r = 0.09, alpha=3/500000, gamma_r = 1/5.15, g_r=0.1, r_f=0.0084, K_f=6, d_f=1/5, beta_h=0.19, gamma_h=1/10, g_h=0.34, b_h=1/(25*365), d_h=1/(25*365)) #you can play with transmission and recovery rates here

if(uniform==TRUE){ params.set<-bSIR.uniform  #uniform LHS distribution from `LHSnonuniform.R`
} else{ params.set<-LHS_bSIR} #non-uniform LHS distribution from `LHSnonuniform.R`

bSIR <- data.frame(matrix(rep(NA,h),nrow=h, ncol= ncol(params.set)+6))
colnames(bSIR)<-c(names(parameters), "MaxInf", "MaxDur", "Thresh100", "Thresh10", "Thresh1", "Thresh250")

for(i in 1:h){
  bSIR[i,1:length(parameters)] <- params <- as.list(c(params.set[i,]))
  out <- as.data.frame(ode(init, times, bubonicSIR, params))
  bSIR[i,length(parameters)+1] <- max(out$D_h)
  bSIR[i,length(parameters)+2] <- which.max(out$D_h)
  difference<-as.vector(NA)
  for(j in 2:length(out$D_h)){
    difference[j-1]<-out$D_h[j]-out$D_h[j-1]
  }
  bSIR[i,length(parameters)+3]<-ifelse(length(which(difference>100)[length(which(difference>100))])>0, which(difference>100)[length(which(difference>100))], 0) 
  bSIR[i,length(parameters)+4]<-ifelse(length(which(difference>10)[length(which(difference>10))])>0, which(difference>10)[length(which(difference>10))], 0) 
  bSIR[i,length(parameters)+5]<-ifelse(length(which(difference>1)[length(which(difference>1))])>0, which(difference>1)[length(which(difference>1))], 0) 
  bSIR[i,length(parameters)+6]<-ifelse(length(which(difference>250)[length(which(difference>250))])>0, which(difference>250)[length(which(difference>250))], 0) 
}



# Bubonic SEIR model -------------------------------------------------------

#Define initial conditions and expected parameter values
init <- c(S_r=499999, I_r=1, R_r=0, D_r=0, H=6, Fl=0, S_h = 500000, E_h=0, I_h = 0, R_h=0, D_h=0) #population size, and how many individuals start in each susceptible, infected, or removed category
parameters <- c(beta_r = 0.09, alpha=3/500000, gamma_r = 1/5.15, g_r=0.1, r_f=0.0084, K_f=6, d_f=1/5, beta_h=0.19, sigma_h= 1/4, gamma_h=1/10, g_h=0.34, b_h=1/(25*365), d_h=1/(25*365)) #you can play with transmission and recovery rates here

if(uniform==TRUE){ params.set<-LHS_bSEIR.uniform  #uniform LHS distribution from `LHSnonuniform.R`
} else{ params.set<-LHS_bSEIR} #non-uniform LHS distribution from `LHSnonuniform.R`

bSEIR <- data.frame(matrix(rep(NA,h),nrow=h, ncol= ncol(params.set)+6))
colnames(bSEIR)<-c(names(parameters), "MaxInf", "MaxDur", "Thresh100", "Thresh10", "Thresh1", "Thresh250")

for(i in 1:h){
  bSEIR[i,1:length(parameters)] <- params <- as.list(c(params.set[i,]))
  out <- as.data.frame(ode(init, times, bubonicSEIR, params))
  bSEIR[i,length(parameters)+1] <- max(out$D_h)
  bSEIR[i,length(parameters)+2] <- which.max(out$D_h)
  difference<-as.vector(NA)
  for(j in 2:length(out$D_h)){
    difference[j-1]<-out$D_h[j]-out$D_h[j-1]
  }
  bSEIR[i,length(parameters)+3]<-ifelse(length(which(difference>100)[length(which(difference>100))])>0, which(difference>100)[length(which(difference>100))], 0) 
  bSEIR[i,length(parameters)+4]<-ifelse(length(which(difference>10)[length(which(difference>10))])>0, which(difference>10)[length(which(difference>10))], 0) 
  bSEIR[i,length(parameters)+5]<-ifelse(length(which(difference>1)[length(which(difference>1))])>0, which(difference>1)[length(which(difference>1))], 0) 
  bSEIR[i,length(parameters)+6]<-ifelse(length(which(difference>250)[length(which(difference>250))])>0, which(difference>250)[length(which(difference>250))], 0) 
  
}


# Bubonic SIR: Flea/Rat/Human Model with Rat Carrying Capacity & Resistance ----------------------------------------------------

#Define initial conditions and expected parameter values
init <- c(S_r=499999, I_r=1, R_r=0, D_r=0, H=6, Fl=0, S_h = 500000, I_h = 0, R_h=0, D_h=0) #population size, and how many individuals start in each susceptible, infected, or removed category
parameters <- c(r_r=0.014, K_r=499999, p_r=0.975, d_r=0.00055, beta_r = 0.09, alpha=3/500000, gamma_r = 1/5.15, g_r=0.1, r_f=0.0084, K_f=6, d_f=1/5, beta_h=0.19, gamma_h=1/10, g_h=0.34, b_h=1/(25*365), d_h=1/(25*365)) #you can play with transmission and recovery rates here

if(uniform==TRUE){ params.set<-bSIRrK.uniform  #uniform LHS distribution from `LHSnonuniform.R`
} else{ params.set<-LHS_bSIRrK} #non-uniform LHS distribution from `LHSnonuniform.R`

bSIRrK <- data.frame(matrix(rep(NA,h),nrow=h, ncol= ncol(params.set)+6))
colnames(bSIRrK)<-c(names(parameters), "MaxInf", "MaxDur", "Thresh100", "Thresh10", "Thresh1", "Thresh250")

for(i in 1:h){
  bSIRrK[i,1:length(parameters)] <- params <- as.list(c(params.set[i,]))
  out <- as.data.frame(ode(init, times, bubonicSIRratK, params))
  bSIRrK[i,length(parameters)+1] <- max(out$D_h)
  bSIRrK[i,length(parameters)+2] <- which.max(out$D_h)
  difference<-as.vector(NA)
  for(j in 2:length(out$D_h)){
    difference[j-1]<-out$D_h[j]-out$D_h[j-1]
  }
  bSIRrK[i,length(parameters)+3]<-ifelse(length(which(difference>100)[length(which(difference>100))])>0, which(difference>100)[length(which(difference>100))], 0) 
  bSIRrK[i,length(parameters)+4]<-ifelse(length(which(difference>10)[length(which(difference>10))])>0, which(difference>10)[length(which(difference>10))], 0) 
  bSIRrK[i,length(parameters)+5]<-ifelse(length(which(difference>1)[length(which(difference>1))])>0, which(difference>1)[length(which(difference>1))], 0) 
  bSIRrK[i,length(parameters)+6]<-ifelse(length(which(difference>250)[length(which(difference>250))])>0, which(difference>250)[length(which(difference>250))], 0) 
  
}



# Bubonic SEIR: Flea/Rat/Human Model with Rat Carrying Capacity & Resistance ----------------------------------------------------

#Define initial conditions and expected parameter values
init <- c(S_r=249950, I_r=1, R_r=0, D_r=0, H=6, Fl=0, S_h = 500000, E_h=0, I_h = 0, R_h=0, D_h=0) #population size, and how many individuals start in each susceptible, infected, or removed category
parameters <- c(r_r=0.014, K_r=499999, p_r=0.975, d_r=0.00055, beta_r = 0.09, alpha=3/500000, gamma_r = 1/5.15, g_r=0.1, r_f=0.0084, K_f=6, d_f=1/5, beta_h=0.19, sigma_h= 1/4, gamma_h=1/10, g_h=0.34, b_h=1/(25*365), d_h=1/(25*365)) #you can play with transmission and recovery rates here

if(uniform==TRUE){ params.set<-LHS_bSEIRrK.uniform  #uniform LHS distribution from `LHSnonuniform.R`
} else{ params.set<-LHS_bSEIRrK} #non-uniform LHS distribution from `LHSnonuniform.R`

bSEIRrK <- data.frame(matrix(rep(NA,h),nrow=h, ncol= ncol(params.set)+6))
colnames(bSEIRrK)<-c(names(parameters), "MaxInf", "MaxDur", "Thresh100", "Thresh10", "Thresh1", "Thresh250")

for(i in 1:h){
  bSEIRrK[i,1:length(parameters)] <- params <- as.list(c(params.set[i,]))
  out <- as.data.frame(ode(init, times, bubonicSEIRratsK, params))
  bSEIRrK[i,length(parameters)+1] <- max(out$D_h)
  bSEIRrK[i,length(parameters)+2] <- which.max(out$D_h)
  difference<-as.vector(NA)
  for(j in 2:length(out$D_h)){
    difference[j-1]<-out$D_h[j]-out$D_h[j-1]
  }
  bSEIRrK[i,length(parameters)+3]<-ifelse(length(which(difference>100)[length(which(difference>100))])>0, which(difference>100)[length(which(difference>100))], 0) 
  bSEIRrK[i,length(parameters)+4]<-ifelse(length(which(difference>10)[length(which(difference>10))])>0, which(difference>10)[length(which(difference>10))], 0) 
  bSEIRrK[i,length(parameters)+5]<-ifelse(length(which(difference>1)[length(which(difference>1))])>0, which(difference>1)[length(which(difference>1))], 0) 
  bSEIRrK[i,length(parameters)+6]<-ifelse(length(which(difference>250)[length(which(difference>250))])>0, which(difference>250)[length(which(difference>250))], 0) 
  
}




# Bubonic/Pneumonic SEIR --------------------------------------------------

#Define initial conditions and expected parameter values
init <- c(S_r=499999, I_r=1, R_r=0, D_r=0, H=6, Fl=0, S_h = 500000, E_b=0, E_p=0, I_b = 0, I_p =0, R_h=0, D_h=0) #population size, and how many individuals start in each susceptible, infected, or removed category
parameters <- c(beta_r = 0.09, alpha=3/500000, gamma_r = 1/5.15, g_r=0.1, r_f=0.0084, K_f=6, d_f=1/5, beta_b=0.19, beta_p = 0.45, sigma_b= 1/6, sigma_p=1/4.3, gamma_b=1/10, gamma_p=1/2.5, p=0.2, g_h=0.34, b_h=1/(25*365), d_h=1/(25*365)) #you can play with transmission and recovery rates here

if(uniform==TRUE){ params.set<-bpSEIR.uniform  #uniform LHS distribution from `LHSnonuniform.R`
} else{ params.set<-LHS_bpSEIR} #non-uniform LHS distribution from `LHSnonuniform.R`

bpSEIR <- data.frame(matrix(rep(NA,h),nrow=h, ncol= ncol(params.set)+6))
colnames(bpSEIR)<-c(names(parameters), "MaxInf", "MaxDur", "Thresh100", "Thresh10", "Thresh1", "Thresh250")

for(i in 1:h){
  bpSEIR[i,1:length(parameters)] <- params <- as.list(c(params.set[i,]))
  out <- as.data.frame(ode(init, times, bubonic_pneumonicSEIR, params))
  bpSEIR[i,length(parameters)+1] <- max(out$D_h)
  bpSEIR[i,length(parameters)+2] <- which.max(out$D_h)
  difference<-as.vector(NA)
  for(j in 2:length(out$D_h)){
    difference[j-1]<-out$D_h[j]-out$D_h[j-1]
  }
  bpSEIR[i,length(parameters)+3]<-ifelse(length(which(difference>100)[length(which(difference>100))])>0, which(difference>100)[length(which(difference>100))], 0) 
  bpSEIR[i,length(parameters)+4]<-ifelse(length(which(difference>10)[length(which(difference>10))])>0, which(difference>10)[length(which(difference>10))], 0) 
  bpSEIR[i,length(parameters)+5]<-ifelse(length(which(difference>1)[length(which(difference>1))])>0, which(difference>1)[length(which(difference>1))], 0) 
  bpSEIR[i,length(parameters)+6]<-ifelse(length(which(difference>250)[length(which(difference>250))])>0, which(difference>250)[length(which(difference>250))], 0) 
  
}



# Human Ectoparasite Model --------------------------------------------------
#Define initial conditions and parameter values
# init <- c(S_l=15*499999, I_l=0, S_h = 499999, I_low=0, I_high=1, R_h=0, D_h=0) #population size, and how many individuals start in each susceptible, infected, or removed category
# parameters <- c(r_l = 0.11, K_l=15, beta_low=0.04, beta_high=0.3, beta_l=0.05, gamma_lice = 1/3, gamma_low=1/2, gamma_high=1/8, g_h=0.34, b_h=1/(25*365), d_h=1/(25*365)) #you can play with transmission and recovery rates here
# 
# if(uniform==TRUE){ params.set<-LHS_eSIR.uniform  #uniform LHS distribution from `LHSnonuniform.R`
# } else{ params.set<-LHS_eSIR} #non-uniform LHS distribution from `LHSnonuniform.R`
# 
# 
# eSIR <- data.frame(matrix(rep(NA,h),nrow=h, ncol= ncol(params.set)+6))
# colnames(eSIR)<-c(names(parameters), "MaxInf", "MaxDur", "Thresh1", "Thresh0.1", "Thresh0.01")
# 
# for(i in 1:h){
#   eSIR[i,1:length(parameters)] <- params <- as.list(c(params.set[i,]))
#   out <- as.data.frame(ode(init, times, liceSIR, params))
#   eSIR[i,length(parameters)+1] <- max(out$D_h)
#   eSIR[i,length(parameters)+2] <- which.max(out$D_h)
#   difference<-as.vector(NA)
#   for(j in 2:length(out$D_h)){
#     difference[j-1]<-out$D_h[j]-out$D_h[j-1]
#   }
#   eSIR[i,length(parameters)+3]<-ifelse(length(which(difference>1)[length(which(difference>1))])>0, which(difference>1)[length(which(difference>1))], 0) 
#   eSIR[i,length(parameters)+4]<-ifelse(length(which(difference>0.1)[length(which(difference>0.1))])>0, which(difference>0.1)[length(which(difference>0.1))], 0) 
#   eSIR[i,length(parameters)+5]<-ifelse(length(which(difference>0.01)[length(which(difference>0.01))])>0, which(difference>0.01)[length(which(difference>0.01))], 0) 
#   
# }
# 
# par(mfrow=c(1,2))
# plot(eSIR$MaxInf~eSIR$r_l, main= expression(paste("Effect of ", r[l], " on Size")))
# plot(eSIR$Thresh1~eSIR$r_l, main= expression(paste("Effect of ", r[l], " on Duration")))
# plot(eSIR$MaxInf~eSIR$K_l, main= expression(paste("Effect of ", K[l], " on Size")))
# plot(eSIR$Thresh1~eSIR$K_l, main= expression(paste("Effect of ", K[l], " on Duration")))
# plot(eSIR$MaxInf~eSIR$beta_low, main= expression(paste("Effect of ", beta[low], " on Size")))
# plot(eSIR$Thresh1~eSIR$beta_low, main= expression(paste("Effect of ", beta[low], " on Duration")))
# plot(eSIR$MaxInf~eSIR$beta_high, main= expression(paste("Effect of ", beta[high], " on Size")))
# plot(eSIR$Thresh1~eSIR$beta_high, main= expression(paste("Effect of ", beta[high], " on Duration")))
# plot(eSIR$MaxInf~eSIR$beta_l, main= expression(paste("Effect of ", beta[l], " on Size")))
# plot(eSIR$Thresh1~eSIR$beta_l, main= expression(paste("Effect of ", beta[l], " on Duration")))
# plot(eSIR$MaxInf~eSIR$gamma_lice,  main= expression(paste("Effect of ", gamma[lice], " on Size")))
# plot(eSIR$Thresh1~eSIR$gamma_lice, main= expression(paste("Effect of ", gamma[lice], " on Duration")))
# plot(eSIR$MaxInf~eSIR$gamma_low,  main= expression(paste("Effect of ", gamma[low], " on Size")))
# plot(eSIR$Thresh1~eSIR$gamma_low, main= expression(paste("Effect of ", gamma[low], " on Duration")))
# plot(eSIR$MaxInf~eSIR$gamma_high,  main= expression(paste("Effect of ", gamma[high], " on Size")))
# plot(eSIR$Thresh1~eSIR$gamma_high, main= expression(paste("Effect of ", gamma[high], " on Duration")))
# plot(eSIR$MaxInf~eSIR$g_h, main= expression(paste("Effect of ", g[h], " on Size")))
# plot(eSIR$Thresh1~eSIR$g_h, main= expression(paste("Effect of ", g[h], " on Duration")))
# plot(eSIR$MaxInf~eSIR$b_h, main= expression(paste("Effect of ", b[h], " on Size")))
# plot(eSIR$Thresh1~eSIR$b_h, main= expression(paste("Effect of ", b[h], " on Duration")))
# plot(eSIR$MaxInf~eSIR$d_h, main= expression(paste("Effect of ", d[h], " on Size")))
# plot(eSIR$Thresh1~eSIR$d_h, main= expression(paste("Effect of ", d[h], " on Duration")))
# 
# par(mfrow=c(1,2))
# boxplot(eSIR$MaxInf, main= "Outbreak Size", ylab= "Number of Dead Humans", ylim=c(0,500000))
# boxplot(eSIR$Thresh1, main= "Outbreak Duration", ylab="Time (Days)")
# 
# bonferroni.alpha <- 0.05/5
# prcc_size <- pcc(eSIR[,1:length(parameters)], eSIR[,length(parameters)+1], nboot = 1000, rank=TRUE, conf=1-bonferroni.alpha)
# prcc_duration <- pcc(eSIR[,1:length(parameters)], eSIR[,length(parameters)+2], nboot = 1000, rank=TRUE, conf=1-bonferroni.alpha)
# 
# #plot correlation coefficients and confidence intervals for epidemic size and duration
# size<-prcc_size$PRCC
# size$param<-rownames(size)
# colnames(size)[4:5] <- c("maxCI", "minCI")
# 
# ggplot(size, aes(x=param, y = original)) +
#   geom_point(size = 4)+ 
#   geom_errorbar(aes(ymax = maxCI, ymin = minCI))+
#   ggtitle("PRCC Outbreak Size: Ectoparasite SIR")
# 
# duration<-prcc_duration$PRCC
# duration$param<-rownames(duration)
# colnames(duration)[4:5] <- c("maxCI", "minCI")
# 
# ggplot(duration, aes(x=param, y = original)) +
#   geom_point(size = 4)+ 
#   geom_errorbar(aes(ymax = maxCI, ymin = minCI))+
#   ggtitle("PRCC Outbreak Duration: Ectoparasite SIR")

