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
         