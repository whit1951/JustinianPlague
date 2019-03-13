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
