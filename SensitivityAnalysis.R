# Sensitivity analysis to evaluate Plague SIR ODE models

# Sensitivity Analysis ----------------------------------------------------
install.packages("ODEsensitivity")
library(ODEsensitivity)
library(parallel)
#lower bound for parameters
binf= c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
bsup= c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
parameters <- c("beta_r", "alpha", "gamma_r", "g_r", "r_f", "K_f", "d_f", "beta_b", "beta_p", "sigma_b", "sigma_p", "gamma_b", "gamma_p", "p", "g_h")
init <- c(S_r=500000, I_r=1, R_r=0, D_r=0, H=5, Fl=0, S_h = 500000, E_b=0, E_p=0, I_b = 0, I_p =0, R_h=0, D_h=0) #population size, and how many individuals start in each susceptible, infected, or removed category

ODEmorris(mod=bubonic_pneumonicSEIR, pars=parameters, state_init= init, times= SIRtimes, binf=binf, bsup=bsup, r=500) 
bubonic_pneumonic <- ODEmorris(mod = bubonic_pneumonicSEIR,
                               pars = parameters,
                               state_init = init,
                               times = SIRtimes,
                               binf = binf,
                               bsup = bsup,
                               r = 500,
                               design = list(type = "oat", 
                                             levels = 10, grid.jump = 1),
                               scale = TRUE,
                               ode_method = "ode45",
                               parallel_eval = TRUE,
                               parallel_eval_ncores = parallel::detectCores())


bubonic_pneumonic <- ODEsobol(mod = bubonic_pneumonicSEIR,
                               pars = parameters,
                               state_init = init,
                               times = SIRtimes, n=10,
                              ode_method = "ode45",
                               parallel_eval = TRUE,
                               parallel_eval_ncores = parallel::detectCores())
# ODEsensitivity example --------------------------------------------------
# https://d-nb.info/1160443556/34
                
         
##### Lotka-Volterra equations #####
# The model function:
LVmod <- function(Time, State, Pars) {
 with(as.list(c(State, Pars)), {
   Ingestion    <- rIng  * Prey * Predator
   GrowthPrey   <- rGrow * Prey * (1 - Prey/K)
   MortPredator <- rMort * Predator
   
   dPrey        <- GrowthPrey - Ingestion
   dPredator    <- Ingestion * assEff - MortPredator
   
   return(list(c(dPrey, dPredator)))
 })
}
# The parameters to be included in the sensitivity analysis and their lower 
# and upper boundaries:
LVpars  <- c("rIng", "rGrow", "rMort", "assEff", "K")
LVbinf <- c(0.05, 0.05, 0.05, 0.05, 1)
LVbsup <- c(1.00, 3.00, 0.95, 0.95, 20)
# The initial values of the state variables:
LVinit  <- c(Prey = 1, Predator = 2)
# The timepoints of interest:
LVtimes <- c(0.01, seq(1, 50, by = 1))
# Morris screening:
set.seed(7292)
# Warning: The following code might take very long!

LVres_morris <- ODEmorris(mod = LVmod,
                         pars = LVpars,
                         state_init = LVinit,
                         times = LVtimes,
                         binf = LVbinf,
                         bsup = LVbsup,
                         r = 500,
                         design = list(type = "oat", 
                                       levels = 10, grid.jump = 1),
                         scale = TRUE,
                         ode_method = "lsoda",
                         parallel_eval = TRUE,
                         parallel_eval_ncores = 2)

plot(LVres_morris, pars_plot = c("rIng", "rMort", "assEff"), state_plot = "Predator")

#The model function
sir <- function(time, state, parameters) {
 with(as.list(c(state, parameters)), {
   dS <- -beta * S * I
   dI <- beta * S * I - gamma * I
   dR <- gamma * I
   
   return(list(c(dS, dI, dR)))
 })
}

# The parameters to be included in the sensitivity analysis and their lower 
# and upper boundaries:
SIRpars  <- c("beta", "gamma")
SIRbinf <- c(0, 0)
SIRbsup <- c(1, 1)
         
# The initial values of the state variables:
SIRinit  <- c(S = 1000, I=1, R=0)
# The timepoints of interest:
SIRtimes <- c(0.01, seq(1, 50, by = 1))
# Morris screening:
set.seed(7292)
# Warning: The following code might take very long!

SIRres_morris <- ODEmorris(mod = sir,
                          pars = SIRpars,
                          state_init = SIRinit,
                          times = SIRtimes,
                          binf = SIRbinf,
                          bsup = SIRbsup,
                          r = 500,
                          design = list(type = "oat", 
                                        levels = 10, grid.jump = 1),
                          scale = TRUE,
                          ode_method = "lsoda",
                          parallel_eval = TRUE,
                          parallel_eval_ncores = parallel::detectCores())

plot(SIRres_morris, pars_plot = c("beta", "gamma"),
     state_plot = "I")
plot(SIRres_morris, pars_plot = c("beta", "gamma"),
     state_plot = "S")


#Set up system of ordinary differential equations
pneumonicSIR <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    N_h<- S_h + E_h+ I_h 
    dS_h <- b_h*S_h -beta_p * S_h * I_h/N_h -d_h*S_h
    dE_h<-  beta_p * S_h * I_h/N_h - sigma_p*E_h - d_h*E_h
    dI_h <- sigma_p*E_h - gamma_p*I_h -d_h*I_h
    dD_h <- gamma_p*I_h #Deaths due to disease related mortality
    
    return(list(c(dS_h, dE_h, dI_h, dD_h)))
  })
}
         


# Try manually with `lhs` package -----------------------------------------
# https://daphnia.ecology.uga.edu/drakelab/wp-content/uploads/2015/07/sensitivity-ebola.pdf
#install.packages('lhs')
require(lhs) #add the lhs library
h <- 1000 #choose number of parameter sets to sample to sample
set.seed(6242015)
lhs<-maximinLHS(h,5)

#expected parameter values
parameters <- c(beta_p = 0.0734, sigma_p= 1/4.3, gamma_p = 1/2.5, b_h=1/(25*365), d_h=1/(25*365)) #you can play with transmission and recovery rates here


beta_p.min <- 0.001
beta_p.max <- 1
sigma_p.min <- 1/6
sigma_p.max <- 1/2
gamma_p.min <- 1/3
gamma_p.max <- 1/4
b_h.min <- 1/(30*365)
b_h.max <- 1/(20*365)
d_h.min <- 1/(30*365)
d_h.max <- 1/(20*365)

params.set <- cbind( beta_p = lhs[,1]*(beta_p.max-beta_p.min)+beta_p.min,
                     sigma_p = lhs[,2]*(sigma_p.max-sigma_p.min)+sigma_p.min,
                     gamma_p = lhs[,3]*(gamma_p.max-gamma_p.min)+gamma_p.min,
                     b_h = lhs[,4]*(b_h.max-b_h.min)+b_h.min,
                     d_h = lhs[,5]*(d_h.max-d_h.min)+d_h.min)


h2 <-250

init <- c(S_h = 500000, E_h= 0, I_h = 1, D_h=0)
times <- seq(0, 500, by= .1)

data <- data.frame(matrix(rep(NA,h2),nrow=h2, ncol= ncol(params.set)+2))
for(i in 1:h2){
       data[i,1:5] <- params <- as.list(c(params.set[i,]))
       out <- as.data.frame(ode(init, times, pneumonicSIR, params))
       data[i,6] <- max(out$D_h)
       data[i,7] <- which.max(out$D_h)
}
names(data) <- c(names(params),'outbreak.size', 'outbreak.duration')
par(mfrow=c(2,5))
boxplot(data$outbreak.size~data$beta_p)
boxplot(data$outbreak.duration~data$beta_p)
boxplot(data$outbreak.size~data$sigma_p)
boxplot(data$outbreak.duration~data$sigma_p)
boxplot(data$outbreak.size~data$gamma_p)
boxplot(data$outbreak.duration~data$gamma_p)
boxplot(data$outbreak.size~data$b_h)
boxplot(data$outbreak.duration~data$b_h)
boxplot(data$outbreak.size~data$d_h)
boxplot(data$outbreak.duration~data$d_h)

par(mfrow=c(1,2))
boxplot(data$outbreak.size)
boxplot(data$outbreak.duration)

#save(data, file='data.Rdata')

library(sensitivity)
bonferroni.alpha <- 0.05/5
prcc_size <- pcc(data[,1:5], data[,6], nboot = 1000, rank=TRUE, conf=1-bonferroni.alpha)
prcc_duration <- pcc(data[,1:5], data[,7], nboot = 1000, rank=TRUE, conf=1-bonferroni.alpha)

#plot correlation coefficients and confidence intervals
duration<-prcc_duration$PRCC
duration$param<-rownames(duration)
colnames(duration)[4:5] <- c("maxCI", "minCI")

require(ggplot2)
ggplot(duration, aes(x=param, y = original)) +
  geom_point(size = 4)+ 
  geom_errorbar(aes(ymax = maxCI, ymin = minCI))

