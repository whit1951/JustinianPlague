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
    N_r<-S_r+I_r+R_r
    dS_r <- -beta_r * S_r * Fl *(1-exp(-1*alpha*N_r))/N_r #Susceptible rats
    dI_r <- -beta_r * S_r * Fl *(1-exp(-1*alpha*N_r))/N_r- gamma_r * I_r #Infected rats
    dR_r <- g_r*gamma_r * I_r #Recovered rats
    dD_r <- (1-g_r)*gamma_r * I_r #Dead rats

    
    dH<-r_f*H*(1-H/K_f) #Fleas/rat
    dFl<-(1-g_r)*gamma_r*I_r*H-d_f*Fl #Fleas in environment
    
    N_h<-S_h+I_h+R_h
    dS_h<- -beta_h * S_h * Fl *(exp(-1*alpha*N_r))/N_h #susceptible humans
    dI_h<- beta_h * S_h * Fl *(exp(-1*alpha*N_r))/N_h- gamma_h * I_h #infected humans
    dR_h<- g_h*gamma_h * I_h #recovered humans
    dD_h<- (1-g_h)*gamma_h * I_h #dead humans
    
    return(list(c(dS_r, dI_r, dR_r, dD_r, dH, dFl, dS_h, dI_h, dR_h, dD_h)))
  })
}

#Define initial conditions and parameter values
init <- c(S_r<-500000, I_r=1, R_r=0, D_r=0, H=5, Fl=0, S_h = 500000, I_h = 0, R_h=0, D_h=0) #population size, and how many individuals start in each susceptible, infected, or removed category
parameters <- c(beta_r = 0.0734, alpha=3/S_r, gamma_r = 1/5.2, g_r=0.1, r_f=0.0084, K_f=6, d_f=1/5, beta_h=0.0734, gamma_h=0.1, g_h=0.2) #you can play with transmission and recovery rates here
time <- seq(0, 15, by = .01) #how long to integrate over [time interval]?

#Run ordinary differential equation solver
out <- as.data.frame(ode(y = init, times = time, func = ratfleaSIR, parms = parameters))
out$time<-NULL

#Plot the output
matplot(time, out, type = "l", xlab = "Time (days)", ylab = "Number of Individuals", main = "SIR Model", lwd = 1, lty = 1, bty = "l", 
        col = 2:4)
legend("right", c("Susceptible", "Infected", "Recovered"), pch = 1, col = 2:4)
