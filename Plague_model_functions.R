#install.packages("deSolve")
library(deSolve)


# Pneumonic plague SIR ----------------------------------------------------

#Set up system of ordinary differential equations
pneumonicSIR <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    N_h <- S_h + I_h
    dS_h <- b_h*S_h -beta_p * S_h * I_h/N_h -d_h*S_h
    dI_h <- beta_p * S_h * I_h/N_h - gamma_p*I_h
    dD_h <- gamma_p*I_h #Deaths due to disease related mortality
    
    return(list(c(dS_h, dI_h, dD_h)))
  })
}




# Pneumonic plague SEIR --------------------------------------------------------

#Set up system of ordinary differential equations
pneumonicSEIR <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    N_h<- S_h + E_h+ I_h
    dS_h <- b_h*S_h -beta_p * S_h * I_h/N_h -d_h*S_h
    dE_h<-  beta_p * S_h * I_h/N_h - sigma_p*E_h
    dI_h <- sigma_p*E_h - gamma_p*I_h
    dD_h <- gamma_p*I_h #Deaths due to disease related mortality
    
    return(list(c(dS_h, dE_h, dI_h, dD_h)))
  })
}


# Bubonic SIR: Flea/Rat/Human Model ----------------------------------------------------
bubonicSIR <- function(time, state, parameters) {
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
    dS_h<- -beta_h * S_h * Fl *(exp(-1*alpha*N_r))/N_h + b_h*(S_h+R_h) - d_h*S_h #susceptible humans
    dI_h<- beta_h * S_h * Fl *(exp(-1*alpha*N_r))/N_h - gamma_h*I_h #infected humans
    dR_h<- g_h*gamma_h * I_h- d_h*I_h #recovered humans
    dD_h<- (1-g_h)*gamma_h * I_h #dead humans from disease
    
    return(list(c(dS_r, dI_r, dR_r, dD_r, dH, dFl, dS_h, dI_h, dR_h, dD_h)))
  })
}

# Bubonic SEIR: Flea/Rat/Human Model ----------------------------------------------------

bubonicSEIR <- function(time, state, parameters) {
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
    N_h<- S_h + E_h+ I_h + R_h #Total number of humans
    dS_h<- -beta_h * S_h * Fl *(exp(-1*alpha*N_r))/N_h + b_h*(S_h+R_h) - d_h*S_h #susceptible humans
    dE_h<- beta_h * S_h * Fl *(exp(-1*alpha*N_r))/N_h - sigma_h*E_h #exposed humans
    dI_h<- sigma_h*E_h - gamma_h*I_h #infected humans
    dR_h<- g_h*gamma_h * I_h - d_h*I_h #recovered humans
    dD_h<- (1-g_h)*gamma_h * I_h #dead humans from disease
    
    return(list(c(dS_r, dI_r, dR_r, dD_r, dH, dFl, dS_h, dE_h, dI_h, dR_h, dD_h)))
  })
}

# Bubonic SIR: Flea/Rat/Human Model with Rat Carrying Capacity & Resistance ----------------------------------------------------
bubonicSIRratK <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    #Rats
    N_r  <- S_r + I_r + R_r #Total number of rats
    dS_r <- r_r*S_r*(1-N_r/K_r)+r_r*R_r*(1-p_r)-d_r*S_r -beta_r * S_r * Fl *(1-exp(-1*alpha*N_r))/N_r #Susceptible rats
    dI_r <- beta_r * S_r * Fl *(1-exp(-1*alpha*N_r))/N_r - gamma_r * I_r #Infected rats
    dR_r <- r_r*R_r*(p_r-N_r/K_r)+g_r*gamma_r * I_r-d_r*R_r #Recovered rats
    dD_r <- (1-g_r)*gamma_r * I_r #Dead rats from plague only
    
    #Fleas
    dH  <- r_f*H*(1-H/K_f) #Fleas/rat
    dFl <- (1-g_r)*gamma_r*I_r*H - d_f*Fl #Fleas in environment
    
    #Humans
    N_h<- S_h + I_h + R_h #Total number of humans
    dS_h<- -beta_h * S_h * Fl *(exp(-1*alpha*N_r))/N_h + b_h*(S_h+R_h) - d_h*S_h #susceptible humans
    dI_h<- beta_h * S_h * Fl *(exp(-1*alpha*N_r))/N_h - gamma_h*I_h #infected humans
    dR_h<- g_h*gamma_h * I_h- d_h*I_h #recovered humans
    dD_h<- (1-g_h)*gamma_h * I_h #dead humans from disease
    
    return(list(c(dS_r, dI_r, dR_r, dD_r, dH, dFl, dS_h, dI_h, dR_h, dD_h)))
  })
}

# Bubonic SEIR: Flea/Rat/Human Model with Rat Carrying Capacity & Resistance ----------------------------------------------------

bubonicSEIRratsK <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    #Rats
    N_r  <- S_r + I_r + R_r #Total number of rats
    dS_r <- r_r*S_r*(1-N_r/K_r)+r_r*R_r*(1-p_r)-d_r*S_r -beta_r * S_r * Fl *(1-exp(-1*alpha*N_r))/N_r #Susceptible rats
    dI_r <- beta_r * S_r * Fl *(1-exp(-1*alpha*N_r))/N_r - gamma_r * I_r #Infected rats
    dR_r <- r_r*R_r*(p_r-N_r/K_r)+g_r*gamma_r * I_r-d_r*R_r #Recovered rats
    dD_r <- (1-g_r)*gamma_r * I_r #Dead rats from plague only
    
    #Fleas
    dH  <- r_f*H*(1-H/K_f) #Fleas/rat
    dFl <- (1-g_r)*gamma_r*I_r*H - d_f*Fl #Fleas in environment
    
    #Humans
    N_h<- S_h + E_h+ I_h + R_h #Total number of humans
    dS_h<- -beta_h * S_h * Fl *(exp(-1*alpha*N_r))/N_h + b_h*(S_h+R_h) - d_h*S_h #susceptible humans
    dE_h<- beta_h * S_h * Fl *(exp(-1*alpha*N_r))/N_h - sigma_h*E_h #exposed humans
    dI_h<- sigma_h*E_h - gamma_h*I_h #infected humans
    dR_h<- g_h*gamma_h * I_h - d_h*I_h #recovered humans
    dD_h<- (1-g_h)*gamma_h * I_h #dead humans from disease
    
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
    dS_h<- -beta_b * S_h * Fl *(exp(-1*alpha*N_r))/N_h - beta_p * S_h * I_p/N_h + b_h*(S_h+R_h)-d_h*S_h #susceptible humans
    dE_b<- beta_b * S_h * Fl *(exp(-1*alpha*N_r))/N_h - sigma_b*E_b #exposed humans (bubonic)
    dE_p<- beta_p * S_h * I_p/N_h - sigma_p*E_p #exposed humans (pneumonic)
    dI_b<- sigma_b*E_b - gamma_b*I_b #infected humans (bubonic)
    dI_p<- sigma_p*E_p + p*gamma_b*I_b - gamma_p*I_p #infected humans (pneumonic)
    dR_h<- g_h*gamma_b * I_b- d_h*R_h #recovered humans (just bubonic)
    dD_h<- (1-p-g_h)*gamma_b * I_b + gamma_p*I_p  #dead humans (bubonic + pneumonic)
    
    return(list(c(dS_r, dI_r, dR_r, dD_r, dH, dFl, dS_h, dE_b, dE_p, dI_b, dI_p, dR_h, dD_h)))
  })
}

# Lice SIR: Ectoparasite/Human Model ----------------------------------------------------
liceSIR <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {

    #Humans
    N_h<- S_h + I_low+ I_high + R_h #Total number of humans
    dS_h<- -beta_l * S_h * I_l/N_h + b_h*(S_h+R_h) - d_h*S_h #susceptible humans
    dI_low<- beta_l * S_h * I_l/N_h- gamma_low*I_low
    dI_high<- (1-g_h)*gamma_low*I_low - gamma_high*I_high #infected humans
    dR_h<- g_h*gamma_low*I_low- d_h*I_high #recovered humans
    dD_h<- gamma_high * I_high #dead humans from disease
    
    #Lice
    N_l=S_l+I_l
    dS_l  <- r_l*S_l*(1-S_l/(K_l*N_h))- (beta_low*I_low+beta_high*I_high)*S_l/N_h #Susceptile lice
    dI_l <- (beta_low*I_low+beta_high*I_high)*S_l/N_h- gamma_lice*I_l #Infected lice
    
    return(list(c(dS_l, dI_l, dS_h, dI_low, dI_high, dR_h, dD_h)))
  })
}

