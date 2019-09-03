#' Plotting plague model functions
#' April 30, 2019
#' Lauren White @email: lwhite@sesync.org

#install.packages("deSolve")
library(deSolve)
library(tidyr)
library(ggplot2)

source('~/JustinianPlague/Plague_model_functions.R')


# Plotting Bubonic SIR: Flea/Rat/Human Model----------------------------------------------------------------

#Define initial conditions and parameter values
init <- c(S_r=100000, I_r=1, R_r=0, D_r=0, H=6, Fl=0, S_h = 500000, I_h = 0, R_h=0, D_h=0) #population size, and how many individuals start in each susceptible, infected, or removed category
parameters <- c(beta_r = 0.09, alpha=3/500000, gamma_r = 1/5.15, g_r=0.1, r_f=0.0084, K_f=6, d_f=1/5, beta_h=0.19, gamma_h=1/10, g_h=0.34, b_h=1/(25*365), d_h=1/(25*365)) #you can play with transmission and recovery rates here
time <- seq(0, 1000, by = 1) #how long to integrate over [time interval]?

#Run ordinary differential equation solver
bubonicSIR_out <- as.data.frame(ode(y = init, times = time, func = bubonicSIR, parms = parameters))
bubonicSIR_out$time<-NULL

#Plot the output
matplot(time, bubonicSIR_out[,1:4], type = "l", main= "Rats", xlab = "Time (days)", ylab = "Number of Individuals", lwd = 1, lty = 1, bty = "l", 
        col = 1:4)
legend("right", c("Susceptible", "Infected", "Recovered", "Dead"), pch = 1, col = 1:4)

matplot(time, bubonicSIR_out[,5:6], type = "l", main= "Fleas", xlab = "Time (days)", ylab = "Number of Individuals", lwd = 1, lty = 1, bty = "l", 
        col = 1:2)
legend("right", c("Fleas/rat", "Free fleas"), pch = 1, col = 1:2)

matplot(time, bubonicSIR_out[,7:10], type = "l", main= "Humans", xlab = "Time (days)", ylab = "Number of Individuals", lwd = 1, lty = 1, bty = "l", 
        col = 1:5)
legend("right", c("Susceptible", "Infected", "Recovered", "Dead"), pch = 1, col = 1:5)


# Plotting Bubonic SEIR: Flea/Rat/Human Model----------------------------------------------------------------

#Define initial conditions and parameter values
init <- c(S_r=249950, I_r=1, R_r=0, D_r=0, H=6, Fl=0, S_h = 500000, E_h=0, I_h = 0, R_h=0, D_h=0) #population size, and how many individuals start in each susceptible, infected, or removed category
parameters <- c(beta_r = 0.09, alpha=3/500000, gamma_r = 1/5.15, g_r=0.1, r_f=0.0084, K_f=6, d_f=1/5, beta_h=0.19, sigma_h= 1/4, gamma_h=1/10, g_h=0.34, b_h=1/(25*365), d_h=1/(25*365)) #you can play with transmission and recovery rates here
time <- seq(0, 1000, by = 1) #how long to integrate over [time interval]?

#Run ordinary differential equation solver
bubonicSEIR_out <- as.data.frame(ode(y = init, times = time, func = bubonicSEIR, parms = parameters))
bubonicSEIR_out$time<-NULL

#Plot the output
matplot(time, bubonicSEIR_out[,1:4], type = "l", main= "Rats", xlab = "Time (days)", ylab = "Number of Individuals", lwd = 1, lty = 1, bty = "l", 
        col = 1:4)
legend("right", c("Susceptible", "Infected", "Recovered", "Dead"), pch = 1, col = 1:4)

matplot(time, bubonicSEIR_out[,5:6], type = "l", main= "Fleas", xlab = "Time (days)", ylab = "Number of Individuals", lwd = 1, lty = 1, bty = "l", 
        col = 1:2)
legend("right", c("Fleas/rat", "Free fleas"), pch = 1, col = 1:2)

matplot(time, bubonicSEIR_out[,7:11], type = "l", main= "Humans", xlab = "Time (days)", ylab = "Number of Individuals", lwd = 1, lty = 1, bty = "l", 
        col = 1:5)
legend("right", c("Susceptible", "Exposed", "Infected", "Recovered", "Dead"), pch = 1, col = 1:5)


# Plotting Bubonic SIR: Flea/Rat/Human Model with Rat Carrying Capacity & Resistance----------------------------------------------------------------

#Define initial conditions and parameter values
init <- c(S_r=499999, I_r=1, R_r=0, D_r=0, H=6, Fl=0, S_h = 500000, I_h = 0, R_h=0, D_h=0) #population size, and how many individuals start in each susceptible, infected, or removed category
parameters <- c(r_r=0.014, K_r=499999, p_r=0.975, d_r=0.00055, beta_r = 0.09, alpha=3/500000, gamma_r = 1/5.15, g_r=0.1, r_f=0.0084, K_f=6, d_f=1/5, beta_h=0.19, gamma_h=1/10, g_h=0.34, b_h=1/(25*365), d_h=1/(25*365)) #you can play with transmission and recovery rates here
time <- seq(0, 1000, by = 1) #how long to integrate over [time interval]?

#Run ordinary differential equation solver
bSIRrK_out <- as.data.frame(ode(y = init, times = time, func = bubonicSIRratK, parms = parameters))
bSIRrK_out$time<-NULL

#Plot the output
matplot(time, bSIRrK_out[,1:4], type = "l", main= "Rats", xlab = "Time (days)", ylab = "Number of Individuals", lwd = 1, lty = 1, bty = "l", 
        col = 1:4)
legend("right", c("Susceptible", "Infected", "Recovered", "Dead"), pch = 1, col = 1:4)

matplot(time, bSIRrK_out[,5:6], type = "l", main= "Fleas", xlab = "Time (days)", ylab = "Number of Individuals", lwd = 1, lty = 1, bty = "l", 
        col = 1:2)
legend("right", c("Fleas/rat", "Free fleas"), pch = 1, col = 1:2)

matplot(time, bSIRrK_out[,7:10], type = "l", main= "Humans", xlab = "Time (days)", ylab = "Number of Individuals", lwd = 1, lty = 1, bty = "l", 
        col = 1:5)
legend("right", c("Susceptible", "Infected", "Recovered", "Dead"), pch = 1, col = 1:5)


# Plotting Bubonic SEIR: Flea/Rat/Human Model with Rat Carrying Capacity & Resistance----------------------------------------------------------------

#Define initial conditions and parameter values
init <- c(S_r=499999, I_r=1, R_r=0, D_r=0, H=6, Fl=0, S_h = 500000, E_h=0, I_h = 0, R_h=0, D_h=0) #population size, and how many individuals start in each susceptible, infected, or removed category
parameters <- c(r_r=0.014, K_r=499999, p_r=0.975, d_r=0.00055, beta_r = 0.09, alpha=3/500000, gamma_r = 1/5.15, g_r=0.1, r_f=0.0084, K_f=6, d_f=1/5, beta_h=0.19, sigma_h= 1/4, gamma_h=1/10, g_h=0.34, b_h=1/(25*365), d_h=1/(25*365)) #you can play with transmission and recovery rates here
time <- seq(0, 1000, by = 1) #how long to integrate over [time interval]?

#Run ordinary differential equation solver
bSEIRrK_out <- as.data.frame(ode(y = init, times = time, func = bubonicSEIRratsK, parms = parameters))
bSEIRrK_out$time<-NULL

#Plot the output
matplot(time, bSEIRrK_out[,1:4], type = "l", main= "Rats", xlab = "Time (days)", ylab = "Number of Individuals", lwd = 1, lty = 1, bty = "l", 
        col = 1:4)
legend("right", c("Susceptible", "Infected", "Recovered", "Dead"), pch = 1, col = 1:4)

matplot(time, bSEIRrK_out[,5:6], type = "l", main= "Fleas", xlab = "Time (days)", ylab = "Number of Individuals", lwd = 1, lty = 1, bty = "l", 
        col = 1:2)
legend("right", c("Fleas/rat", "Free fleas"), pch = 1, col = 1:2)

matplot(time, bSEIRrK_out[,7:11], type = "l", main= "Humans", xlab = "Time (days)", ylab = "Number of Individuals", lwd = 1, lty = 1, bty = "l", 
        col = 1:5)
legend("right", c("Susceptible", "Exposed", "Infected", "Recovered", "Dead"), pch = 1, col = 1:5)



# Plotting Pneumonic plague SIR -------------------------------------------

#Define initial conditions and parameter values
init <- c(S_h = 499999, I_h = 1, D_h=0) #population size, and how many individuals start in each susceptible, infected, or removed category
parameters <- c(beta_p = 0.45, gamma_p = 1/2.5, b_h=1/(25*365), d_h=1/(25*365)) #you can play with transmission and recovery rates here
time <- seq(0, 1000, by= 1) #how long to integrate over [time interval]?

#Run ordinary differential equation solver
pneumonic_out <- as.data.frame(ode(y = init, times = time, func = pneumonicSIR, parms = parameters))
pneumonic_out$time<-NULL

#Plot the output
matplot(time, pneumonic_out, type = "l", xlab = "Time (days)", ylab = "Number of Individuals", main = "Pneumonic Plague in People", lwd = 1, lty = 1, bty = "l", 
        col = 1:3)
# legend("topright", c("Susceptible", "Infected", "Dead from Plague"), pch = 1, col = 1:4)


# Plotting Pneumonic Plague SEIR ------------------------------------------

#Define initial conditions and parameter values
init <- c(S_h = 499999, E_h= 0, I_h = 1, D_h=0) #population size, and how many individuals start in each susceptible, infected, or removed category
parameters <- c(beta_p = 0.45, sigma_p= 1/4.3, gamma_p = 1/2.5, b_h=1/(25*365), d_h=1/(25*365)) #you can play with transmission and recovery rates here
time <- seq(0, 1000, by= 1) #how long to integrate over [time interval]?

#Run ordinary differential equation solver
pneumonicSEIR_out <- as.data.frame(ode(y = init, times = time, func = pneumonicSEIR, parms = parameters))
pneumonicSEIR_out$time<-NULL

#Plot the output
matplot(time, pneumonicSEIR_out, type = "l", xlab = "Time (days)", ylab = "Number of Individuals", main = "Pneumonic Plague in People", lwd = 1, lty = 1, bty = "l", 
        col = 1:4)
#legend("right", c("Susceptible", "Exposed", "Infected", "Dead from Plague"), pch = 1, col = 1:4)


# Plotting Bubonic/Pneumonic SEIR -----------------------------------------

#Define initial conditions and parameter values
init <- c(S_r=499999, I_r=1, R_r=0, D_r=0, H=6, Fl=0, S_h = 500000, E_b=0, E_p=0, I_b = 0, I_p =0, R_h=0, D_h=0) #population size, and how many individuals start in each susceptible, infected, or removed category
parameters <- c(beta_r = 0.09, alpha=3/500000, gamma_r = 1/5.15, g_r=0.1, r_f=0.0084, K_f=6, d_f=1/5, beta_b=0.19, beta_p = 0.45, sigma_b= 1/6, sigma_p=1/4.3, gamma_b=1/10, gamma_p=1/2.5, p=0.2, g_h=0.34, b_h=1/(25*365), d_h=1/(25*365)) #you can play with transmission and recovery rates here
time <- seq(0, 1000, by = 1) #how long to integrate over [time interval]?

#Run ordinary differential equation solver
bp_out <- as.data.frame(ode(y = init, times = time, func = bubonic_pneumonicSEIR, parms = parameters))
bp_out$time<-NULL

#Plot the output
matplot(time, bp_out[,1:4], type = "l", main= "Rats", xlab = "Time (days)", ylab = "Number of Individuals", lwd = 1, lty = 1, bty = "l", 
        col = 1:4)
legend("right", c("Susceptible", "Infected", "Recovered", "Dead"), pch = 1, col = 1:4)

matplot(time, bp_out[,5:6], type = "l", main= "Fleas", xlab = "Time (days)", ylab = "Number of Individuals", lwd = 1, lty = 1, bty = "l", 
        col = 1:2)
legend("right", c("Fleas/rat", "Free fleas"), pch = 1, col = 1:2)

matplot(time, bp_out[,7:13], type = "l", main= "Humans", xlab = "Time (days)", ylab = "Number of Individuals", lwd = 1, lty = 1, bty = "l", 
        col = 1:7)
legend("right", c("Susceptible", "Exposed Bubonic", "Exposed Pneumonic",  "Infected Bubonic", "Infected Pneumonic", "Recovered", "Dead"), pch = 1, col = 1:7)

#summarize
humans<-data.frame(S=bp_out$S_h, E=bp_out$E_b + bp_out$E_b, I=bp_out$I_b + bp_out$I_p, R=bp_out$R_h, D=bp_out$D_h)
matplot(time, humans[,1:5], type = "l", main= "Humans", xlab = "Time (days)", ylab = "Number of Individuals", lwd = 1, lty = 1, bty = "l", 
        col = 1:5)
legend("right", c("Susceptible", "Exposed",  "Infected", "Recovered", "Dead"), pch = 1, col = 1:5)


# Human ectoparasites -----------------------------------------------------
#Define initial conditions and parameter values
init <- c(S_l=15*499999, I_l=0, S_h = 499999, I_low=0, I_high=1, R_h=0, D_h=0) #population size, and how many individuals start in each susceptible, infected, or removed category
parameters <- c(r_l = 0.11, K_l=15, beta_low=0.04, beta_high=0.3, beta_l=0.05, gamma_lice = 1/3, gamma_low=1/2, gamma_high=1/8, g_h=0.34, b_h=1/(25*365), d_h=1/(25*365)) #you can play with transmission and recovery rates here
time <- seq(0, 1000, by = 1) #how long to integrate over [time interval]?

#Run ordinary differential equation solver
lice_out <- as.data.frame(ode(y = init, times = time, func = liceSIR, parms = parameters))
lice_out$time<-NULL

#Plot the output
matplot(time, lice_out[,1:2], type = "l", main= "Lice", xlab = "Time (days)", ylab = "Number of Individuals", lwd = 1, lty = 1, bty = "l", 
        col = 1:4)
legend("right", c("Susceptible", "Infected"), pch = 1, col = 1:4)

matplot(time, lice_out[,3:7], type = "l", main= "Humans", xlab = "Time (days)", ylab = "Number of Individuals", lwd = 1, lty = 1, bty = "l", 
        col = 1:5)
#legend("right", c("S_h"="Susceptible", "I_low"="Low Infectious", "I_high"="High Infectious",  "R_h"="Recovered", "D_h"="Dead"), pch = 1, col = 1:5)
legend("right", colnames(lice_out)[3:7], pch=1, col=1:5)


# Compare time series of model types --------------------------------------
comp<- data.frame(time= time, pSIR=pneumonic_out$D_h, pSEIR=pneumonicSEIR_out$D_h, bSIR=bubonicSIR_out$D_h, bSEIR=bubonicSEIR_out$D_h, bSIRrK=bSIRrK_out$D_h, bSEIRrK=bSEIRrK_out$D_h, bpSEIR=bp_out$D_h, lice=lice_out$D_h)
# write.csv(comp, "NumberDead.csv")
long_DF <- comp %>% gather(Model, NumberDead, c(pSIR, pSEIR, bSIR, bSEIR, bpSEIR, lice, bSIRrK, bSEIRrK))
ggplot(long_DF, aes(time, NumberDead, col=Model)) + geom_line() +
  xlab("Time (Days)") + ylab("Number Dead")+
  scale_color_discrete(labels = c(lice="Human Ectoparasites", bpSEIR="Bubonic/Pneumonic (SEIR)", bSIR="Bubonic Plague (SIR)", bSEIR="Bubonic Plague (SEIR)", pSEIR="Pneumonic Plague (SEIR)", pSIR="Pneumonic Plague (SIR)", bSIRrK="Bubonic SIR (K & Rest.)", bSEIRrK="Bubonic SEIR (K & Rest.)"))

