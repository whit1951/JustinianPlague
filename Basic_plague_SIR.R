#install.packages("deSolve")
library(deSolve)

#Set up system of ordinary differential equations
sir <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dS <- -beta * S * I
    dI <- beta * S * I - gamma * I
    dR <- gamma * I
    
    return(list(c(dS, dI, dR)))
  })
}

#Define initial conditions and parameter values
init <- c(S = 500000, I = 1, R=0) #population size, and how many individuals start in each susceptible, infected, or removed category
parameters <- c(beta = 0.0734, gamma = 0.5) #you can play with transmission and recovery rates here
time <- seq(0, 15, by = .01) #how long to integrate over [time interval]?

#Run ordinary differential equation solver
out <- as.data.frame(ode(y = init, times = time, func = sir, parms = parameters))
out$time<-NULL

#Plot the output
matplot(time, out, type = "l", xlab = "Time (days)", ylab = "Number of Individuals", main = "SIR Model", lwd = 1, lty = 1, bty = "l", 
        col = 2:4)
legend("right", c("Susceptible", "Infected", "Recovered"), pch = 1, col = 2:4)



# Flea/Rat/Human Model ----------------------------------------------------

#Set up system of ordinary differential equations
ratfleaSIR <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    #Rats
    N_r  <- S_r + I_r + R_r #Total number of rats
    dS_r <- -beta_r * S_r * Fl *(1-exp(-1*alpha*N_r))/N_r #Susceptible rats
    dI_r <- beta_r * S_r * Fl *(1-exp(-1*alpha*N_r))/N_r- gamma_r * I_r #Infected rats
    dR_r <- g_r*gamma_r * I_r #Recovered rats
    dD_r <- (1-g_r)*gamma_r * I_r #Dead rats

    #Fleas
    dH  <- r_f*H*(1-H/K_f) #Fleas/rat
    dFl <- (1-g_r)*gamma_r*I_r*H - d_f*Fl #Fleas in environment
    
    #Humans
    N_h<- S_h + I_h + R_h #Total number of humans
    dS_h<- -beta_h * S_h * Fl *(exp(-1*alpha*N_r))/N_h #susceptible humans
    dE_h<- beta_h * S_h * Fl *(exp(-1*alpha*N_r))/N_h - sigma_h*E_h #exposed humans
    dI_h<- sigma_h*E_h - gamma_h*I_h #infected humans
    dR_h<- g_h*gamma_h * I_h #recovered humans
    dD_h<- (1-g_h)*gamma_h * I_h #dead humans
    
    return(list(c(dS_r, dI_r, dR_r, dD_r, dH, dFl, dS_h, dE_h, dI_h, dR_h, dD_h)))
  })
}

#Define initial conditions and parameter values
init <- c(S_r=499999, I_r=1, R_r=0, D_r=0, H=6, Fl=0, S_h = 500000, E_h=0, I_h = 0, R_h=0, D_h=0) #population size, and how many individuals start in each susceptible, infected, or removed category
parameters <- c(beta_r = 0.09, alpha=3/500000, gamma_r = 1/5.15, g_r=0.1, r_f=0.0084, K_f=6, d_f=1/5, beta_h=0.19, sigma_h= 1/4, gamma_h=1/10, g_h=0.34) #you can play with transmission and recovery rates here
time <- seq(0, 1000, by = .1) #how long to integrate over [time interval]?

#Run ordinary differential equation solver
ratflea_out <- as.data.frame(ode(y = init, times = time, func = ratfleaSIR, parms = parameters))
ratflea_out$time<-NULL

#Plot the output
matplot(time, ratflea_out[,1:4], type = "l", main= "Rats", xlab = "Time (days)", ylab = "Number of Individuals", lwd = 1, lty = 1, bty = "l", 
        col = 1:4)
legend("right", c("Susceptible", "Infected", "Recovered", "Dead"), pch = 1, col = 1:4)

matplot(time, ratflea_out[,5:6], type = "l", main= "Fleas", xlab = "Time (days)", ylab = "Number of Individuals", lwd = 1, lty = 1, bty = "l", 
        col = 1:2)
legend("right", c("Fleas/rat", "Free fleas"), pch = 1, col = 1:2)

matplot(time, ratflea_out[,7:11], type = "l", main= "Humans", xlab = "Time (days)", ylab = "Number of Individuals", lwd = 1, lty = 1, bty = "l", 
        col = 1:5)
legend("right", c("Susceptible", "Exposed", "Infected", "Recovered", "Dead"), pch = 1, col = 1:5)


# Pneumonic plague SIR ----------------------------------------------------

#Set up system of ordinary differential equations
pneumonicSIR <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    N_h <- S_h + I_h + D_h
    dS_h <- b_h*S_h -beta_p * S_h * I_h/N_h -d_h*S_h
    dI_h <- beta_p * S_h * I_h/N_h - gamma_p*I_h -d_h*I_h
    dD_h <- gamma_p*I_h #Deaths due to disease related mortality
    
    return(list(c(dS_h, dI_h, dD_h)))
  })
}

#Define initial conditions and parameter values
init <- c(S_h = 499999, I_h = 1, D_h=0) #population size, and how many individuals start in each susceptible, infected, or removed category
parameters <- c(beta_p = 0.45, gamma_p = 1/2.5, b_h=1/(25*365), d_h=1/(25*365)) #you can play with transmission and recovery rates here
time <- seq(0, 1000, by= .1) #how long to integrate over [time interval]?

#Run ordinary differential equation solver
pneumonic_out <- as.data.frame(ode(y = init, times = time, func = pneumonicSIR, parms = parameters))
pneumonic_out$time<-NULL

#Plot the output
matplot(time, pneumonic_out, type = "l", xlab = "Time (days)", ylab = "Number of Individuals", main = "Pneumonic Plague in People", lwd = 1, lty = 1, bty = "l", 
        col = 1:4)
# legend("topright", c("Susceptible", "Infected", "Dead from Plague"), pch = 1, col = 1:4)



# Pneumonic plague SEIR --------------------------------------------------------

#Set up system of ordinary differential equations
pneumonicSEIR <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    N_h<- S_h + E_h+ I_h + D_h
    dS_h <- b_h*S_h -beta_p * S_h * I_h/N_h -d_h*S_h
    dE_h<-  beta_p * S_h * I_h/N_h - sigma_p*E_h - d_h*E_h
    dI_h <- sigma_p*E_h - gamma_p*I_h -d_h*I_h
    dD_h <- gamma_p*I_h #Deaths due to disease related mortality
    
    return(list(c(dS_h, dE_h, dI_h, dD_h)))
  })
}

#Define initial conditions and parameter values
init <- c(S_h = 499999, E_h= 0, I_h = 1, D_h=0) #population size, and how many individuals start in each susceptible, infected, or removed category
parameters <- c(beta_p = 0.45, sigma_p= 1/4.3, gamma_p = 1/2.5, b_h=1/(25*365), d_h=1/(25*365)) #you can play with transmission and recovery rates here
time <- seq(0, 1000, by= .1) #how long to integrate over [time interval]?

#Run ordinary differential equation solver
pneumonicSEIR_out <- as.data.frame(ode(y = init, times = time, func = pneumonicSEIR, parms = parameters))
pneumonicSEIR_out$time<-NULL

#Plot the output
matplot(time, pneumonicSEIR_out, type = "l", xlab = "Time (days)", ylab = "Number of Individuals", main = "Pneumonic Plague in People", lwd = 1, lty = 1, bty = "l", 
        col = 1:4)
#legend("right", c("Susceptible", "Exposed", "Infected", "Dead from Plague"), pch = 1, col = 1:4)



# Flea/Rat/Human Model + Pneumonic Plague Transmission Route ----------------------------------------------------

#Set up system of ordinary differential equations
bubonic_pneumonicSEIR <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    #Rats
    N_r  <- S_r + I_r + R_r #Total number of rats
    dS_r <- -beta_r * S_r * Fl *(1-exp(-1*alpha*N_r))/N_r #Susceptible rats
    dI_r <- beta_r * S_r * Fl *(1-exp(-1*alpha*N_r))/N_r- gamma_r * I_r #Infected rats
    dR_r <- g_r*gamma_r * I_r #Recovered rats
    dD_r <- (1-g_r)*gamma_r * I_r #Dead rats
    
    #Fleas
    dH  <- r_f*H*(1-H/K_f) #Fleas/rat
    dFl <- (1-g_r)*gamma_r*I_r*H - d_f*Fl #Fleas in environment
    
    #Humans
    N_h<- S_h + E_b + E_p + I_b + I_p + R_h #Total number of humans
    dS_h<- -beta_b * S_h * Fl *(exp(-1*alpha*N_r))/N_h - beta_p * S_h * I_p/N_h #susceptible humans
    dE_b<- beta_b * S_h * Fl *(exp(-1*alpha*N_r))/N_h - sigma_b*E_b #exposed humans (bubonic)
    dE_p<- beta_p * S_h * I_p/N_h - sigma_p*E_p #exposed humans (pneumonic)
    dI_b<- sigma_b*E_b - gamma_b*I_b #infected humans (bubonic)
    dI_p<- sigma_p*E_p + p*gamma_b*I_b - gamma_p*I_p #infected humans (pneumonic)
    dR_h<- g_h*gamma_b * I_b #recovered humans (just bubonic)
    dD_h<- (1-p-g_h)*gamma_b * I_b + gamma_p*I_p  #dead humans (bubonic + pneumonic)
    
    return(list(c(dS_r, dI_r, dR_r, dD_r, dH, dFl, dS_h, dE_b, dE_p, dI_b, dI_p, dR_h, dD_h)))
  })
}

#Define initial conditions and parameter values
init <- c(S_r=499999, I_r=1, R_r=0, D_r=0, H=6, Fl=0, S_h = 500000, E_b=0, E_p=0, I_b = 0, I_p =0, R_h=0, D_h=0) #population size, and how many individuals start in each susceptible, infected, or removed category
parameters <- c(beta_r = 0.09, alpha=3/500000, gamma_r = 1/5.15, g_r=0.1, r_f=0.0084, K_f=6, d_f=1/5, beta_b=0.19, beta_p = 0.45, sigma_b= 1/6, sigma_p=1/4.3, gamma_b=1/10, gamma_p=1/2.5, p=0.2, g_h=0.34) #you can play with transmission and recovery rates here
time <- seq(0, 1000, by = .1) #how long to integrate over [time interval]?

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



# Compare time series of model types --------------------------------------
library(tidyr)
library(ggplot2)
comp<- data.frame(time= time, pSIR=pneumonic_out$D_h, pSEIR=pneumonicSEIR_out$D_h, bSIR=ratflea_out$D_h, bpSEIR=bp_out$D_h)
long_DF <- comp %>% gather(Model, NumberDead, c(pSIR, pSEIR, bSIR, bpSEIR))
ggplot(long_DF, aes(time, NumberDead, col=Model)) + geom_line() +
  xlab("Time (Days)") + ylab("Number Dead")+
  scale_color_discrete(labels = c("Pneumonic Plague (SIR)", "Pneumonic Plague (SEIR)", "Bubonic Plague (SIR)", "Bubonic/Pneumonic SEIR"))
