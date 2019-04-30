#install.packages("deSolve")
library(deSolve)


# Basic SIR ---------------------------------------------------------------

sir <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dS <- -beta * S * I
    dI <- beta * S * I - gamma * I
    dR <- gamma * I
    
    return(list(c(dS, dI, dR)))
  })
}





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


# Bubonic SIR: Flea/Rat/Human Model ----------------------------------------------------

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

