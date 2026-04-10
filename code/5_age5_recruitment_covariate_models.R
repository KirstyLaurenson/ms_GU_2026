#### Code to run multistate CMR models estimating covariate impacts on immature 
#### guillemot recruitment at age 5

# Parameters:
# phi.1: first year survival probability
# phi.2: second year survival probability (= phi.1)
# phi.3: third year survival probability
# phi.4: fourth year survival probability
# phi.5: fifth year non survival probability 
# phi.6: 6+ non-breeder survival probability
# phi.ad: adult survival probability
#
# alpha.1: probability to start breeding when 1 year old (0)
# alpha.2: probability to start breeding when 2 years old (0)
# alpha.3: probability to start breeding when 3 years old (0)
# alpha.4: probability to start breeding when 4 years old
# alpha.5: probability to start breeding when 5 years old
# alpha.6: probability to start breeding when 6 years old
#
# p.NB: recapture probability of individuals that have never bred before
# p.B: recapture probability of breeders

library(dplyr)
library(nimble)

####  data  ####
amos <- read.csv("reduced_guille_eh.csv")
y <- amos %>%
  as.matrix()

head(y)
tail(y)

df <- read.csv("reduced_model_cov2.csv")

# there are no immature guillemots aged 5 in the first four years of the model as
# all guillemots enter the model as chicks - so remove first four years of covariates
# change covariate as appropriate for model
df$lag_wsst[1] <- 0
df$lag_wsst[2] <- 0
df$lag_wsst[3] <- 0
df$lag_wsst[4] <- 0

####  estimate impact of single covariate on age 5 recruitment  ####

### define constants and data
# Definition of constants to be used in model
K <- ncol(y) # number of occasions
N <- nrow(y) # Number of individuals

# first capture occasion 
first <- apply(y, 1, function(x) min(which(x !=0)))

# constants 
my.constants <- list(N = nrow(y),   # individuals
                     K = ncol(y),   # occasions (t)
                     first = first, 
                     cov1 = df$lag_wsst)
# data 
my.data <- list(y = y + 1) # non-detection becomes 1

# initial values
# Zs must go 1, 2, 3, 4 - model puts all birds through these states, and the only states you can remain in
# are NB>3 (5) or B (6). After 4th re-encounter, can remain 5, or become 6
# there is probably a much more streamlined way to do this - but this works for me!!

zinits <- matrix(1, ncol = ncol(y), nrow = nrow(y)) # everyone starts as fledgling

get.last.juv <- function(x) max(which(x==1)) # define last time an individual was seen as a fledgling

last <- apply(y, 1, get.last.juv)

for(i in 1:nrow(zinits)) {
  if((last[i] + 1)<=ncol(zinits)) {
    for(j in (last[i] + 1):ncol(zinits)){
      zinits[i, j] <- 2 # all occasions past fledgling are assigned NB1 
    }
  }
}

# make all occasions from fifth recapture (6th occasion) onwards NB>4 (6)
for(i in 1:nrow(zinits)) {
  if((first[i]+2)==ncol(zinits)) next
  if((first[i]+3)==ncol(zinits)) next  
  if((first[i]+4)==ncol(zinits)) next
  zinits[i, (first[i]+5):ncol(zinits)] <- 6
}

# replace second resighting with 3 (NB2)
for(i in 1:nrow(zinits)) {
  zinits[i, (first[i]+2)] <- 3
}

# replace third resighting with 4 (NB3)
for(i in 1:nrow(zinits)) {
  if((first[i]+2)==ncol(zinits)) next    
  zinits[i, (first[i]+3)] <- 4
}

# replace fourth resighting with 5 (NB4)
for(i in 1:nrow(zinits)) {
  if((first[i]+2)==ncol(zinits)) next   
  if((first[i]+3) ==ncol(zinits)) next
  zinits[i, (first[i]+4)] <- 5
}

# add back breeding sightings
zinits <- ifelse(y == 3, 7, zinits) # seen as breeder

# replace 6s (NB>4) after 7s (B) - once a breeder you remain a breeder - transitions unidirectional
replace_6_after_7 <- function(row) {
  found_7 <- FALSE
  for (i in 1:length(row)) {
    if (row[i] == 7) {
      found_7 <- TRUE
    }
    if (found_7 && row[i] == 6) {
      row[i] <- 7
    }
  }
  return(row)
}

# apply function
for (i in 1:nrow(zinits)) {
  zinits[i, ] <- replace_6_after_7(zinits[i, ])
}

### nimble model:
myCode <- nimbleCode ({
  # ------------------------------------------------------
  # Parameters:
  # phi.1: first year survival probability
  # phi.2: second year survival probability
  # phi.3: third year survival probability
  # phi.4: fourth year survival probability
  # phi.5: fifth year survival probability (NB)
  # phi.6: 6+ non-breeder survival probability
  # phi.ad: adult survival probability
  #
  # alpha.1: probability to start breeding when 1 year old (0)
  # alpha.2: probability to start breeding when 2 years old (0)
  # alpha.3: probability to start breeding when 3 years old (0)
  # alpha.4: probability to start breeding when 4 years old
  # alpha.5: probability to start breeding when 5 years old
  # alpha.6: probability to start breeding when 6 years old
  #
  # p.NB: recapture probability of individuals that have never bred before
  # p.B: recapture probability of breeders
  #
  # -------------------------------------------------------
  # 8 States: 
  # Fledgling
  # NB 1 year old
  # NB 2 years old
  # NB 3 years old
  # NB 4 years old
  # NB >4 years old (older than 4, 5+)
  # Breeder (4+)
  # Dead
  #
  # Observations:
  # 1: not seen
  # 2: juvenile/ringing as chick
  # 3: seen as non-breeder (NB)
  # 4: seen as breeder (B)
  # ------------------------------------------------------
  
  # Priors and constraints:
  
  phi.ad ~ dunif(0, 1)                            # adult survival is constant to avoid trap-dependence problems
  phi.3 ~ dunif(0, 1)
  beta[1] ~ dnorm(0, sd = 1.5)                    # prior intercept
  beta[2] ~ dnorm(0, sd = 1.5)                    # prior slope 
  sdeps ~ dunif(0, 10)                            # sd of variance 
  
  for (t in 1:(K-1)) { # loop over time
    
    phi.1[t] ~ dunif(0, 1)
    phi.2[t] <- phi.1[t]                          # make assumption that first and second year survival are identical
    phi.4[t] ~ dunif(0, 1)    
    phi.5[t] ~ dunif(0, 1)
    phi.6[t] <- phi.5[t]                          # assume that survival in adult non-breeders doesn't vary with age
    
    alpha.1[t] <- 0                               # birds cannot begin breeding before age 4
    alpha.2[t] <- 0
    alpha.3[t] <- 0
    alpha.4[t] ~ dunif(0, 1)
    alpha.6[t] ~ dunif(0, 1)
    
    p.NB[t] ~ dunif(0, 1)
    p.B[t] ~ dunif(0, 1)
    
    logit(alpha.5[t]) <- beta[1] + beta[2] * cov1[t] + eps[t]    # derived parameter 
    
    eps[t] ~ dnorm(mean = 0, sd = sdeps)
    
    # Probabilities of state z(t+1) given state z(t)
    gamma[1,1,t] <- 0
    gamma[1,2,t] <- phi.1[t] * (1 - alpha.1[t])
    gamma[1,3,t] <- 0
    gamma[1,4,t] <- 0
    gamma[1,5,t] <- 0
    gamma[1,6,t] <- 0
    gamma[1,7,t] <- phi.1[t] * alpha.1[t]
    gamma[1,8,t] <- 1 - phi.1[t]
    
    gamma[2,1,t] <- 0
    gamma[2,2,t] <- 0
    gamma[2,3,t] <- phi.2[t] * (1 - alpha.2[t])
    gamma[2,4,t] <- 0
    gamma[2,5,t] <- 0
    gamma[2,6,t] <- 0
    gamma[2,7,t] <- phi.2[t] * alpha.2[t]
    gamma[2,8,t] <- 1 - phi.2[t]
    
    gamma[3,1,t] <- 0
    gamma[3,2,t] <- 0
    gamma[3,3,t] <- 0
    gamma[3,4,t] <- phi.3 * (1 - alpha.3[t])
    gamma[3,5,t] <- 0
    gamma[3,6,t] <- 0
    gamma[3,7,t] <- phi.3 * alpha.3[t]
    gamma[3,8,t] <- 1 - phi.3
    
    gamma[4,1,t] <- 0
    gamma[4,2,t] <- 0
    gamma[4,3,t] <- 0
    gamma[4,4,t] <- 0
    gamma[4,5,t] <- phi.4[t] * (1 - alpha.4[t])
    gamma[4,6,t] <- 0
    gamma[4,7,t] <- phi.4[t] * alpha.4[t]
    gamma[4,8,t] <- 1 - phi.4[t]
    
    gamma[5,1,t] <- 0
    gamma[5,2,t] <- 0
    gamma[5,3,t] <- 0
    gamma[5,4,t] <- 0
    gamma[5,5,t] <- 0
    gamma[5,6,t] <- phi.5[t] * (1 - alpha.5[t])
    gamma[5,7,t] <- phi.5[t] * alpha.5[t]
    gamma[5,8,t] <- 1 - phi.5[t] 
    
    gamma[6,1,t] <- 0
    gamma[6,2,t] <- 0
    gamma[6,3,t] <- 0
    gamma[6,4,t] <- 0
    gamma[6,5,t] <- 0
    gamma[6,6,t] <- phi.6[t] * (1 - alpha.6[t])
    gamma[6,7,t] <- phi.6[t] * alpha.6[t]
    gamma[6,8,t] <- 1 - phi.6[t]
    
    gamma[7,1,t] <- 0
    gamma[7,2,t] <- 0
    gamma[7,3,t] <- 0
    gamma[7,4,t] <- 0
    gamma[7,5,t] <- 0
    gamma[7,6,t] <- 0
    gamma[7,7,t] <- phi.ad
    gamma[7,8,t] <- 1 - phi.ad
    
    gamma[8,1,t] <- 0
    gamma[8,2,t] <- 0
    gamma[8,3,t] <- 0
    gamma[8,4,t] <- 0
    gamma[8,5,t] <- 0
    gamma[8,6,t] <- 0
    gamma[8,7,t] <- 0
    gamma[8,8,t] <- 1
    
    # probabilities of y(t) given z(t)
    omega[1,1,t] <- 1
    omega[1,2,t] <- 0
    omega[1,3,t] <- 0
    omega[1,4,t] <- 0
    
    omega[2,1,t] <- 1
    omega[2,2,t] <- 0
    omega[2,3,t] <- 0
    omega[2,4,t] <- 0
    
    omega[3,1,t] <- 1 - p.NB[t]
    omega[3,2,t] <- 0
    omega[3,3,t] <- p.NB[t]
    omega[3,4,t] <- 0
    
    omega[4,1,t] <- 1 - p.NB[t]
    omega[4,2,t] <- 0
    omega[4,3,t] <- p.NB[t]
    omega[4,4,t] <- 0
    
    omega[5,1,t] <- 1 - p.NB[t]
    omega[5,2,t] <- 0
    omega[5,3,t] <- p.NB[t]
    omega[5,4,t] <- 0
    
    omega[6,1,t] <- 1 - p.NB[t]
    omega[6,2,t] <- 0
    omega[6,3,t] <- p.NB[t]
    omega[6,4,t] <- 0
    
    omega[7,1,t] <- 1 - p.B[t] 
    omega[7,2,t] <- 0
    omega[7,3,t] <- 0
    omega[7,4,t] <- p.B[t]
    
    omega[8,1,t] <- 1 
    omega[8,2,t] <- 0
    omega[8,3,t] <- 0
    omega[8,4,t] <- 0
  }          # t
  
  # likelihood 
  for (i in 1:N) {
    # latent state at first capture - is known
    z[i, first[i]] <- y[i, first[i]] - 1
    for (t in (first[i]+1):K){
      # draw z(t) given z(t-1)
      z[i,t] ~ dcat(gamma[z[i,t-1], 1:8, t-1])
      # draw y(t) given z(t)
      y[i,t] ~ dcat(omega[z[i,t], 1:4, t-1])
    } # t 
  } #i
  
})     # Nimble model

# initial values
initial.values <- function() list(beta = rnorm(2,0,1),
                                  sdeps = runif(1,0,3),
                                  phi.1 = runif(my.constants$K-1,0,1),
                                  phi.3 = runif(1,0,1),
                                  phi.4 = runif(my.constants$K-1,0,1),
                                  phi.5 = runif(my.constants$K-1,0,1),
                                  phi.ad = runif(1,0,1),
                                  alpha.4 = runif(my.constants$K-1,0,1),
                                  alpha.6 = runif(my.constants$K-1,0,1),
                                  p.NB = runif(my.constants$K-1,0,1),
                                  p.B = runif(my.constants$K-1,0,1), 
                                  z = zinits)


# parameters to save
parameters.to.save <- c("beta", "sdeps")
parameters.to.save

# MCMC details 
n.iter <- 150000
n.burnin <- 50000
n.chains <- 3
n.thin <- 20

# run nimble 
system.time({
  mcmc.alpha5.lwsst <- nimbleMCMC(code = myCode,
                                  constants = my.constants,
                                  data = my.data,
                                  inits = initial.values,
                                  monitors = parameters.to.save,
                                  niter = n.iter,
                                  nburnin = n.burnin,
                                  nchains = n.chains, 
                                  thin = n.thin,
                                  samplesAsCodaMCMC = T, 
                                  WAIC = T)
})

save(mcmc.alpha5.lwsst, file = "alpha5_lwsst.RData")