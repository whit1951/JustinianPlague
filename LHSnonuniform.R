#' Create uniform or non-uniform LHS distributions for PRCC analysis 
#' @author Lauren White
#' @date May 22, 2019

#install.packages('lhs')
#install.packages('triangle')
library(lhs)
library(triangle)

set.seed(2718) #set random seed
h <- 100 #choose number of parameter sets/subdivisions to sample to sample

# Pneumonic Plague SIR ----------------------------------------------------

parameters <- c(beta_p = 0.45, gamma_p = 1/2.5, b_h=1/(25*365), d_h=1/(25*365)) 

#uniform distribution
lhs<-maximinLHS(h, length(parameters)) 

beta_p.min <- 0.42
beta_p.max <- 0.48
gamma_p.min <- 1/3.7
gamma_p.max <- 1/1.3
b_h.min <- 1/(30*365)
b_h.max <- 1/(20*365)
d_h.min <- 1/(30*365)
d_h.max <- 1/(20*365)

LHS_pSIR.uniform <- cbind( beta_p = lhs[,1]*(beta_p.max-beta_p.min)+beta_p.min,
                     gamma_p = lhs[,2]*(gamma_p.max-gamma_p.min)+gamma_p.min,
                     b_h = lhs[,3]*(b_h.max-b_h.min)+b_h.min,
                     d_h = lhs[,4]*(d_h.max-d_h.min)+d_h.min)

#non uniform distributions
LHS<-optimumLHS(n=h, k=length(parameters), maxSweeps=2, eps=.1, verbose=FALSE)
beta_p <-qtriangle(LHS[,1], a=0.42, b=0.48, c=0.45)
gamma_p <-1/qnorm(LHS[,2], mean=2.5, sd=1.2)
b_h <-qunif(LHS[,3], min=1/(30*365), max=1/(20*365))
d_h <-qunif(LHS[,4], min= 1/(30*365), max=1/(20*365))
LHS_pSIR<-data.frame(beta_p, gamma_p, b_h, d_h)

neg<-which(LHS_pSIR<0, arr.ind=T) #check for negative parameter values
for (i in 1:dim(neg)[1])
{
  row<-neg[i,1]
  col<-neg[i,2]
  LHS_pSIR[row,col]<-0
}

# Pneumonic Plague SEIR ----------------------------------------------------
parameters <- c(beta_p = 0.45, sigma_p= 1/4.3, gamma_p = 1/2.5, b_h=1/(25*365), d_h=1/(25*365)) #you can play with transmission and recovery rates here
lhs<-maximinLHS(h,length(parameters)) #here 5 is the number of parameters required for the ODE system

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

LHS_pSEIR.uniform <- cbind( beta_p = lhs[,1]*(beta_p.max-beta_p.min)+beta_p.min,
                     sigma_p = lhs[,2]*(sigma_p.max-sigma_p.min)+sigma_p.min,
                     gamma_p = lhs[,3]*(gamma_p.max-gamma_p.min)+gamma_p.min,
                     b_h = lhs[,4]*(b_h.max-b_h.min)+b_h.min,
                     d_h = lhs[,5]*(d_h.max-d_h.min)+d_h.min)


#non uniform distributions
LHS<-optimumLHS(n=h, k=length(parameters), maxSweeps=2, eps=.1, verbose=FALSE)
beta_p <-qtriangle(LHS[,1], a=0.42, b=0.48, c=0.45)
sigma_p<-1/qnorm(LHS[,2], mean=4.3, sd=1.8)
gamma_p <-1/qnorm(LHS[,3], mean=2.5, sd=1.2)
b_h <-qunif(LHS[,4], min=1/(30*365), max=1/(20*365))
d_h <-qunif(LHS[,5], min= 1/(30*365), max=1/(20*365))
LHS_pSEIR<-data.frame(beta_p, sigma_p, gamma_p, b_h, d_h)

neg<-which(LHS_pSEIR<0, arr.ind=T) #check for negative parameter values
for (i in 1:dim(neg)[1])
{
  row<-neg[i,1]
  col<-neg[i,2]
  LHS_pSEIR[row,col]<-0
}

# Bubonic SIR model -------------------------------------------------------
parameters <- c(beta_r = 0.09, alpha=3/500000, gamma_r = 1/5.15, g_r=0.1, r_f=0.0084, K_f=6, d_f=1/5, beta_h=0.19, gamma_h=1/10, g_h=0.34, b_h=1/(25*365), d_h=1/(25*365))
lhs<-maximinLHS(h, length(parameters)) 

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
gamma_h.min <- 1/26
gamma_h.max <- 1/10
g_h.min <- 0.3
g_h.max <- 0.4
b_h.min <- 1/(30*365)
b_h.max <- 1/(20*365)
d_h.min <- 1/(30*365)
d_h.max <- 1/(20*365)

bSIR.uniform <- cbind( beta_r = lhs[,1]*(beta_r.max-beta_r.min)+beta_r.min,
                     alpha = lhs[,2]*(alpha.max-alpha.min)+alpha.min,
                     gamma_r = lhs[,3]*(gamma_r.max-gamma_r.min)+gamma_r.min,
                     g_r = lhs[,4]*(g_r.max-g_r.min)+g_r.min,
                     r_f = lhs[,5]*(r_f.max-r_f.min)+r_f.min,
                     K_f = lhs[,6]*(K_f.max-K_f.min)+K_f.min,
                     d_f = lhs[,7]*(d_f.max-d_f.min)+d_f.min,
                     beta_h= lhs[,8]*(beta_h.max-beta_h.min)+beta_h.min,
                     gamma_h = lhs[,9]*(gamma_h.max-gamma_h.min)+gamma_h.min,
                     g_h = lhs[,10]*(g_h.max-g_h.min)+g_h.min,
                     b_h = lhs[,11]*(b_h.max-b_h.min)+b_h.min,
                     d_h = lhs[,12]*(d_h.max-d_h.min)+d_h.min)


#non uniform distributions
LHS<-optimumLHS(n=h, k=length(parameters), maxSweeps=2, eps=.1, verbose=FALSE)
beta_r <-qtriangle(LHS[,1], a=0.0, b=0.14, c=0.09)
alpha<- qunif(LHS[,2], min=0.39/500000, max=20/500000) #500000= K_r or S_r(t=0)
gamma_r <-1/qnorm(LHS[,3], mean=5.15, sd=0.44)
g_r<-qtriangle(LHS[,4], a=0.0, b=0.37, c=0.1)
r_f <- qunif(LHS[,5], min=0.0084, max=0.055)
K_f<-qtriangle(LHS[,6], a= 3.29, b=11.17, c=6.57)
d_f<-1/qtriangle(LHS[,7], a=1, b=11.66, c=5)
beta_h<-qnorm(LHS[,8], mean=0.19, sd=0.01)
gamma_h<-1/qtriangle(LHS[,9], a=10, b=26, c=10)
g_h<-qtriangle(LHS[,10], a=0.30, b=0.40, c=0.34)
b_h <-qunif(LHS[,11], min=1/(30*365), max=1/(20*365))
d_h <-qunif(LHS[,12], min= 1/(30*365), max=1/(20*365))
LHS_bSIR<-data.frame(beta_r, alpha, gamma_r, g_r, r_f, K_f, d_f, beta_h, gamma_h, g_h, b_h, d_h)

neg<-which(LHS_bSIR<0, arr.ind=T) #check for negative parameter values
if(length(neg>0)){
  for (i in 1:dim(neg)[1])
  {
    row<-neg[i,1]
    col<-neg[i,2]
    LHS_bSIR[row,col]<-0
  }
}


# Bubonic SEIR model -------------------------------------------------------

parameters <- c(beta_r = 0.09, alpha=3/500000, gamma_r = 1/5.15, g_r=0.1, r_f=0.0084, K_f=6, d_f=1/5, beta_h=0.19, sigma_h= 1/4, gamma_h=1/10, g_h=0.34, b_h=1/(25*365), d_h=1/(25*365)) 
lhs<-maximinLHS(h, length(parameters)) 

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
b_h.min <- 1/(30*365)
b_h.max <- 1/(20*365)
d_h.min <- 1/(30*365)
d_h.max <- 1/(20*365)

LHS_bSEIR.uniform <- cbind( beta_r = lhs[,1]*(beta_r.max-beta_r.min)+beta_r.min,
                     alpha = lhs[,2]*(alpha.max-alpha.min)+alpha.min,
                     gamma_r = lhs[,3]*(gamma_r.max-gamma_r.min)+gamma_r.min,
                     g_r = lhs[,4]*(g_r.max-g_r.min)+g_r.min,
                     r_f = lhs[,5]*(r_f.max-r_f.min)+r_f.min,
                     K_f = lhs[,6]*(K_f.max-K_f.min)+K_f.min,
                     d_f = lhs[,7]*(d_f.max-d_f.min)+d_f.min,
                     beta_h= lhs[,8]*(beta_h.max-beta_h.min)+beta_h.min,
                     sigma_h = lhs[,9]*(sigma_h.max-sigma_h.min)+sigma_h.min,
                     gamma_h = lhs[,10]*(gamma_h.max-gamma_h.min)+gamma_h.min,
                     g_h = lhs[,11]*(g_h.max-g_h.min)+g_h.min,
                     b_h = lhs[,12]*(b_h.max-b_h.min)+b_h.min,
                     d_h = lhs[,13]*(d_h.max-d_h.min)+d_h.min)


#non uniform distributions
LHS<-optimumLHS(n=h, k=length(parameters), maxSweeps=2, eps=.1, verbose=FALSE)
beta_r <-qtriangle(LHS[,1], a=0.0, b=0.14, c=0.09)
alpha<- qunif(LHS[,2], min=0.39/500000, max=20/500000) #500000= K_r or S_r(t=0)
gamma_r <-1/qnorm(LHS[,3], mean=5.15, sd=0.44)
g_r<-qtriangle(LHS[,4], a=0.0, b=0.37, c=0.1)
r_f <- qunif(LHS[,5], min=0.0084, max=0.055)
K_f<-qtriangle(LHS[,6], a= 3.29, b=11.17, c=6.57)
d_f<-1/qtriangle(LHS[,7], a=1, b=11.66, c=5)
beta_h<-qnorm(LHS[,8], mean=0.19, sd=0.01)
sigma_h<-1/qtriangle(LHS[,9], a=2, b=6, c=4)
gamma_h<-1/qtriangle(LHS[,9], a=10, b=26, c=10)
g_h<-qtriangle(LHS[,10], a=0.30, b=0.40, c=0.34)
b_h <-qunif(LHS[,11], min=1/(30*365), max=1/(20*365))
d_h <-qunif(LHS[,12], min= 1/(30*365), max=1/(20*365))
LHS_bSEIR<-data.frame(beta_r, alpha, gamma_r, g_r, r_f, K_f, d_f, beta_h, sigma_h, gamma_h, g_h, b_h, d_h)

neg<-which(LHS_bSEIR<0, arr.ind=T) #check for negative parameter values
if(length(neg>0)){
  for (i in 1:dim(neg)[1])
  {
    row<-neg[i,1]
    col<-neg[i,2]
    LHS_bSEIR[row,col]<-0
  }
}


# Bubonic/Pneumonic SEIR --------------------------------------------------

#Define initial conditions and expected parameter values
parameters <- c(beta_r = 0.09, alpha=3/500000, gamma_r = 1/5.15, g_r=0.1, r_f=0.0084, K_f=6, d_f=1/5, beta_b=0.19, beta_p = 0.45, sigma_b= 1/6, sigma_p=1/4.3, gamma_b=1/10, gamma_p=1/2.5, p=0.2, g_h=0.34, b_h=1/(25*365), d_h=1/(25*365))
lhs<-maximinLHS(h, length(parameters)) 

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
b_h.min <- 1/(30*365)
b_h.max <- 1/(20*365)
d_h.min <- 1/(30*365)
d_h.max <- 1/(20*365)

bpSEIR.uniform <- cbind( beta_r = lhs[,1]*(beta_r.max-beta_r.min)+beta_r.min,
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
                     p = lhs[,15]*(p.max-p.min)+p.min,
                     b_h = lhs[,16]*(b_h.max-b_h.min)+b_h.min,
                     d_h = lhs[,17]*(d_h.max-d_h.min)+d_h.min)


#non uniform distributions
LHS<-optimumLHS(n=h, k=length(parameters), maxSweeps=2, eps=.1, verbose=FALSE)

beta_r <-qtriangle(LHS[,1], a=0.0, b=0.14, c=0.09)
alpha<- qunif(LHS[,2], min=0.39/500000, max=20/500000) #500000= K_r or S_r(t=0)
gamma_r <-1/qnorm(LHS[,3], mean=5.15, sd=0.44)
g_r<-qtriangle(LHS[,4], a=0.0, b=0.37, c=0.1)
r_f <- qunif(LHS[,5], min=0.0084, max=0.055)
K_f<-qtriangle(LHS[,6], a= 3.29, b=11.17, c=6.57)
d_f<-1/qtriangle(LHS[,7], a=1, b=11.66, c=5)
beta_b<-qnorm(LHS[,8], mean=0.19, sd=0.01)
sigma_b<-1/qtriangle(LHS[,9], a=2, b=6, c=4)
gamma_b<-1/qtriangle(LHS[,10], a=10, b=26, c=10)
beta_p <-qtriangle(LHS[,11], a=0.42, b=0.48, c=0.45)
sigma_p<-1/qnorm(LHS[,12], mean=4.3, sd=1.8)
gamma_p <-1/qnorm(LHS[,13], mean=2.5, sd=1.2)
g_h<-qtriangle(LHS[,14], a=0.30, b=0.40, c=0.34)
p<-qunif(LHS[,15], min=0, max=0.40)
b_h <-qunif(LHS[,16], min=1/(30*365), max=1/(20*365))
d_h <-qunif(LHS[,17], min= 1/(30*365), max=1/(20*365))
LHS_bpSEIR<-data.frame(beta_r, alpha, gamma_r, g_r, r_f, K_f, d_f, beta_b, sigma_b, gamma_b,
                       beta_p, sigma_p, gamma_p, g_h, p, b_h, d_h)

neg<-which(LHS_bpSEIR<0, arr.ind=T) #check for negative parameter values
if(length(neg>0)){
  for (i in 1:dim(neg)[1])
  {
    row<-neg[i,1]
    col<-neg[i,2]
    LHS_bpSEIR[row,col]<-0
  }
}

# Human Ectoparasite Model --------------------------------------------------
#Define initial conditions and parameter values
parameters <- c(r_l = 0.11, K_l=15, beta_low=0.04, beta_high=0.3, beta_l=0.05, gamma_lice = 1/3, gamma_low=1/2, gamma_high=1/8, g_h=0.34, b_h=1/(25*365), d_h=1/(25*365)) #you can play with transmission and recovery rates here
lhs<-maximinLHS(h, length(parameters)) 

r_l.min<-0.10
r_l.max<-0.12
K_l.min<-10.5
K_l.max<-67.7
beta_low.min<- 0
beta_low.max<- 0.05
beta_high.min<-0
beta_high.max<-1
beta_l.min<-0
beta_l.max<-0.1
gamma_lice.min<-2
gamma_lice.max<-4
gamma_low.min<-0
gamma_low.max<-4
gamma_high.min<-0
gamma_high.max<-4
g_h.min <- 0.3
g_h.max <- 0.4
b_h.min <- 1/(30*365)
b_h.max <- 1/(20*365)
d_h.min <- 1/(30*365)
d_h.max <- 1/(20*365)


eSIR.uniform <- cbind( r_l = lhs[,1]*(r_l.max-r_l.min)+r_l.min,
                     K_l = lhs[,2]*(K_l.max-K_l.min)+K_l.min,
                     beta_low = lhs[,3]*(beta_low.max-beta_low.min)+beta_low.min,
                     beta_high = lhs[,4]*(beta_high.max-beta_high.min)+beta_high.min,
                     beta_l = lhs[,5]*(beta_l.max-beta_l.min)+beta_l.min,
                     gamma_lice = lhs[,6]*(gamma_lice.max-gamma_lice.min)+gamma_lice.min,
                     gamma_low = lhs[,7]*(gamma_low.max-gamma_low.min)+gamma_low.min,
                     gamma_high= lhs[,8]*(gamma_high.max-gamma_high.min)+gamma_high.min,
                     g_h = lhs[,9]*(g_h.max-g_h.min)+g_h.min,
                     b_h = lhs[,10]*(b_h.max-b_h.min)+b_h.min,
                     d_h = lhs[,11]*(d_h.max-d_h.min)+d_h.min)

#non uniform distributions
LHS<-optimumLHS(n=h, k=length(parameters), maxSweeps=2, eps=.1, verbose=FALSE)

r_l <- qunif(LHS[,1], min=0.01, max=0.20) 
K_l <- qtriangle(LHS[,2], a= 10.5, b=67.7, c=15)
beta_low <- qtriangle(LHS[,3], a=0, b=0.05, c=0.04)
beta_high <- qtriangle(LHS[,4], a=0, b=1, c=0.39)
beta_l <- qunif(LHS[,5],min=0, max=0.1)
gamma_lice <- 1/qtriangle(LHS[,6], a=2, b=4, c=3) 
gamma_low <- 1/qtriangle(LHS[,7], a=0, b=4, c=2) 
gamma_high<- 1/qtriangle(LHS[,8], a=6, b=10, c=8) 
g_h<-qtriangle(LHS[,9], a=0.30, b=0.40, c=0.34)
b_h <-qunif(LHS[,10], min=1/(30*365), max=1/(20*365))
d_h <-qunif(LHS[,11], min= 1/(30*365), max=1/(20*365))

LHS_eSIR<-data.frame(r_l, K_l, beta_low, beta_high, beta_l, gamma_lice, gamma_low, gamma_high, 
                        g_h, b_h, d_h)

neg<-which(LHS_eSIR<0, arr.ind=T) #check for negative parameter values
if(length(neg>0)){
  for (i in 1:dim(neg)[1])
  {
    row<-neg[i,1]
    col<-neg[i,2]
    LHS_eSIR[row,col]<-0
  }
}
