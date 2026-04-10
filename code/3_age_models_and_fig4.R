#### Code to run multistate CMR models and identify best model structure to use for
#### models investigating covariate impacts on juvenile and immature survival and recruitment

# model numbers indicate the order models were constructed and match numbers in
# supplementary table S6

# see individual models for details of structure

# code to plot figure 4 can be found at the end of this script

##############
library(dplyr)
library(nimble)

####  data  ####
amos <- read.csv("reduced_guille_eh.csv")

y <- amos %>%
  as.matrix()

### define constants and data
# Definition of constants to be used in model
K <- ncol(y) # number of occasions
N <- nrow(y) # Number of individuals

# first capture occasion 
first <- apply(y, 1, function(x) min(which(x !=0)))

# constants 
my.constants <- list(N = nrow(y),   # individuals
                     K = ncol(y),   # occasions (t)
                     first = first)

# data 
my.data <- list(y = y + 1) # non-detection becomes 1

#### model no 1 ####
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
  
  for (t in 1:(K-1)) { # loop over time
    
    phi.1[t] ~ dunif(0, 1)
    phi.2[t] <- phi.1[t]                          # make assumption that first and second year survival are identical (juvenile)
    phi.3[t] ~ dunif(0, 1)  
    phi.4[t] <- phi.3[t]                          # third and fourth year survival are identical (immature)
    phi.5[t] <- phi.ad                            # survival doesn't change with age from age 5 (adult)
    phi.6[t] <- phi.ad
    
    alpha.1[t] <- 0                               # birds cannot begin breeding before age 4
    alpha.2[t] <- 0
    alpha.3[t] <- 0
    alpha.4[t] ~ dunif(0, 1)
    alpha.5[t] ~ dunif(0, 1)
    alpha.6[t] ~ dunif(0, 1)
    
    p.NB[t] ~ dunif(0, 1)
    p.B[t] ~ dunif(0, 1)
    
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
    gamma[3,4,t] <- phi.3[t] * (1 - alpha.3[t])
    gamma[3,5,t] <- 0
    gamma[3,6,t] <- 0
    gamma[3,7,t] <- phi.3[t] * alpha.3[t]
    gamma[3,8,t] <- 1 - phi.3[t] 
    
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
initial.values <- function() list(phi.1 = runif(my.constants$K-1,0,1),
                                  phi.3 = runif(my.constants$K-1,0,1),
                                  phi.ad = runif(1,0,1),
                                  alpha.4 = runif(my.constants$K-1,0,1),
                                  alpha.5 = runif(my.constants$K-1,0,1),
                                  alpha.6 = runif(my.constants$K-1,0,1),
                                  p.NB = runif(my.constants$K-1,0,1),
                                  p.B = runif(my.constants$K-1,0,1), 
                                  z = zinits)

initial.values()

# parameters to save
parameters.to.save <- c("phi.1", "phi.3", "phi.ad", "alpha.4", "alpha.5", "alpha.6", "p.NB", "p.B")
parameters.to.save

# MCMC details 
n.iter <- 150000
n.burnin <- 50000
n.chains <- 3
n.thin <- 20

# run nimble 
system.time({
  mcmc.guille.agem1 <- nimbleMCMC(code = myCode,
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

save(mcmc.guille.agem1, file = "age-model-v1.RData")
#### model no 2 ####

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

# make all occasions from fourth recapture (5th occasion) onwards NB>3 (5)
for(i in 1:nrow(zinits)) {
  if((first[i]+2)==ncol(zinits)) next
  if((first[i]+3)==ncol(zinits)) next    
  zinits[i, (first[i]+4):ncol(zinits)] <- 5
}

zinits <- ifelse(y == 2, 5, zinits) 

# replace second resighting with 3 (NB2)
for(i in 1:nrow(zinits)) {
  zinits[i, (first[i]+2)] <- 3
}

# replace third resighting with 4 (NB3)
for(i in 1:nrow(zinits)) {
  if((first[i]+2)==ncol(zinits)) next    
  zinits[i, (first[i]+3)] <- 4
}

# add back breeding sightings
zinits <- ifelse(y == 3, 6, zinits) # seen as breeder

# replace 5s (NB>2) after 6s (B) - once a breeder you remain a breeder - transitions unidirectional
replace_5_after_6 <- function(row) {
  found_6 <- FALSE
  for (i in 1:length(row)) {
    if (row[i] == 6) {
      found_6 <- TRUE
    }
    if (found_6 && row[i] == 5) {
      row[i] <- 6
    }
  }
  return(row)
}

# apply function
for (i in 1:nrow(zinits)) {
  zinits[i, ] <- replace_5_after_6(zinits[i, ])
}

myCode <- nimbleCode ({
  # ------------------------------------------------------
  # Parameters:
  # phi.1: first year survival probability
  # phi.2: second year survival probability
  # phi.3: third year survival probability
  # phi.4: fourth year survival probability
  # phi.5: 5+ non-breeder survival probability
  # phi.ad: adult survival probability
  #
  # alpha.1: probability to start breeding when 1 year old
  # alpha.2: probability to start breeding when 2 years old
  # alpha.3: probability to start breeding when 3 years old
  # alpha.4: probability to start breeding when 4 years old
  # alpha.5: probability to start breeding when 5 years old
  #
  # p.NB: recapture probability of individuals that have never bred before
  # p.B: recapture probability of breeders
  #
  # -------------------------------------------------------
  # 7 States: 
  # Fledgling
  # NB 1 year old
  # NB 2 years old
  # NB 3 years old
  # NB 4+ years old
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
  
  for (t in 1:(K-1)) { # loop over time
    
    phi.1[t] ~ dunif(0, 1)
    phi.2[t] <- phi.1[t]                          # make assumption that first and second year survival are identical
    phi.3[t] ~ dunif(0, 1)
    phi.4[t] <- phi.3[t]                          # survival is the same in age3 and 4 immatures
    phi.5[t] ~ dunif(0, 1)
    
    alpha.1[t] <- 0                               # birds cannot begin breeding before age 4
    alpha.2[t] <- 0
    alpha.3[t] <- 0
    alpha.4[t] ~ dunif(0, 1)
    alpha.5[t] ~ dunif(0, 1)
    
    p.NB[t] ~ dunif(0, 1)
    p.B[t] ~ dunif(0, 1)
    
    # Probabilities of state z(t+1) given state z(t)
    gamma[1,1,t] <- 0
    gamma[1,2,t] <- phi.1[t] * (1 - alpha.1[t])
    gamma[1,3,t] <- 0
    gamma[1,4,t] <- 0
    gamma[1,5,t] <- 0
    gamma[1,6,t] <- phi.1[t] * alpha.1[t]
    gamma[1,7,t] <- 1 - phi.1[t]
    
    gamma[2,1,t] <- 0
    gamma[2,2,t] <- 0
    gamma[2,3,t] <- phi.2[t] * (1 - alpha.2[t])
    gamma[2,4,t] <- 0
    gamma[2,5,t] <- 0
    gamma[2,6,t] <- phi.2[t] * alpha.2[t]
    gamma[2,7,t] <- 1 - phi.2[t]
    
    gamma[3,1,t] <- 0
    gamma[3,2,t] <- 0
    gamma[3,3,t] <- 0
    gamma[3,4,t] <- phi.3[t] * (1 - alpha.3[t])
    gamma[3,5,t] <- 0
    gamma[3,6,t] <- phi.3[t] * alpha.3[t]
    gamma[3,7,t] <- 1 - phi.3[t]
    
    gamma[4,1,t] <- 0
    gamma[4,2,t] <- 0
    gamma[4,3,t] <- 0
    gamma[4,4,t] <- 0
    gamma[4,5,t] <- phi.4[t] * (1 - alpha.4[t])
    gamma[4,6,t] <- phi.4[t] * alpha.4[t]
    gamma[4,7,t] <- 1 - phi.4[t]
    
    gamma[5,1,t] <- 0
    gamma[5,2,t] <- 0
    gamma[5,3,t] <- 0
    gamma[5,4,t] <- 0
    gamma[5,5,t] <- phi.5[t] * (1 - alpha.5[t])
    gamma[5,6,t] <- phi.5[t] * alpha.5[t]
    gamma[5,7,t] <- 1 - phi.5[t]
    
    gamma[6,1,t] <- 0
    gamma[6,2,t] <- 0
    gamma[6,3,t] <- 0
    gamma[6,4,t] <- 0
    gamma[6,5,t] <- 0
    gamma[6,6,t] <- phi.ad
    gamma[6,7,t] <- 1 - phi.ad
    
    gamma[7,1,t] <- 0
    gamma[7,2,t] <- 0
    gamma[7,3,t] <- 0
    gamma[7,4,t] <- 0
    gamma[7,5,t] <- 0
    gamma[7,6,t] <- 0
    gamma[7,7,t] <- 1
    
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
    
    omega[6,1,t] <- 1 - p.B[t]
    omega[6,2,t] <- 0
    omega[6,3,t] <- 0
    omega[6,4,t] <- p.B[t]
    
    omega[7,1,t] <- 1 
    omega[7,2,t] <- 0
    omega[7,3,t] <- 0
    omega[7,4,t] <- 0
    
  }          # t
  
  # likelihood 
  for (i in 1:N) {
    # latent state at first capture - is known
    z[i, first[i]] <- y[i, first[i]] - 1
    for (t in (first[i]+1):K){
      # draw z(t) given z(t-1)
      z[i,t] ~ dcat(gamma[z[i,t-1], 1:7, t-1])
      # draw y(t) given z(t)
      y[i,t] ~ dcat(omega[z[i,t], 1:4, t-1])
    } # t 
  } #i
  
})     # Nimble model

# initial values
initial.values <- function() list(phi.1 = runif(my.constants$K-1,0,1),
                                  phi.3 = runif(my.constants$K-1,0,1),
                                  phi.5 = runif(my.constants$K-1,0,1),
                                  phi.ad = runif(1,0,1),
                                  alpha.4 = runif(my.constants$K-1,0,1),
                                  alpha.5 = runif(my.constants$K-1,0,1),
                                  p.NB = runif(my.constants$K-1,0,1),
                                  p.B = runif(my.constants$K-1,0,1), 
                                  z = zinits)

# parameters to save
parameters.to.save <- c("phi.1", "phi.3", "phi.5", "phi.ad", "alpha.4", "alpha.5", "p.NB", "p.B")

# MCMC details 
n.iter <- 150000
n.burnin <- 50000
n.chains <- 3
n.thin <- 20

# run nimble 
system.time({
  mcmc.guille.agem2 <- nimbleMCMC(code = myCode,
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

save(mcmc.guille.agem2, file = "age-model-v2.RData")

#### model no 3 ####

# initial values

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

# make all occasions from fourth recapture (5th occasion) onwards NB>3 (5)
for(i in 1:nrow(zinits)) {
  if((first[i]+2)==ncol(zinits)) next
  if((first[i]+3)==ncol(zinits)) next    
  zinits[i, (first[i]+4):ncol(zinits)] <- 5
}

zinits <- ifelse(y == 2, 5, zinits) 

# replace second resighting with 3 (NB2)
for(i in 1:nrow(zinits)) {
  zinits[i, (first[i]+2)] <- 3
}

# replace third resighting with 4 (NB3)
for(i in 1:nrow(zinits)) {
  if((first[i]+2)==ncol(zinits)) next    
  zinits[i, (first[i]+3)] <- 4
}

# add back breeding sightings
zinits <- ifelse(y == 3, 6, zinits) # seen as breeder

# replace 5s (NB>2) after 6s (B) - once a breeder you remain a breeder - transitions unidirectional
replace_5_after_6 <- function(row) {
  found_6 <- FALSE
  for (i in 1:length(row)) {
    if (row[i] == 6) {
      found_6 <- TRUE
    }
    if (found_6 && row[i] == 5) {
      row[i] <- 6
    }
  }
  return(row)
}

# apply function
for (i in 1:nrow(zinits)) {
  zinits[i, ] <- replace_5_after_6(zinits[i, ])
}

myCode <- nimbleCode ({
  # ------------------------------------------------------
  # Parameters:
  # phi.1: first year survival probability
  # phi.2: second year survival probability
  # phi.3: third year survival probability
  # phi.4: fourth year survival probability
  # phi.5: 5+ non-breeder survival probability
  # phi.ad: adult survival probability
  #
  # alpha.1: probability to start breeding when 1 year old
  # alpha.2: probability to start breeding when 2 years old
  # alpha.3: probability to start breeding when 3 years old
  # alpha.4: probability to start breeding when 4 years old
  # alpha.5: probability to start breeding when 5 years old
  #
  # p.NB: recapture probability of individuals that have never bred before
  # p.B: recapture probability of breeders
  #
  # -------------------------------------------------------
  # 7 States: 
  # Fledgling
  # NB 1 year old
  # NB 2 years old
  # NB 3 years old
  # NB 4+ years old
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
  
  for (t in 1:(K-1)) { # loop over time
    
    phi.1[t] ~ dunif(0, 1)
    phi.2[t] <- phi.1[t]                          # make assumption that first and second year survival are identical
    phi.3[t] ~ dunif(0, 1)
    phi.4[t] ~ dunif(0, 1)    
    phi.5[t] ~ dunif(0, 1)
    
    alpha.1[t] <- 0                               # birds cannot begin breeding before age 4
    alpha.2[t] <- 0
    alpha.3[t] <- 0
    alpha.4[t] ~ dunif(0, 1)
    alpha.5[t] ~ dunif(0, 1)
    
    p.NB[t] ~ dunif(0, 1)
    p.B[t] ~ dunif(0, 1)
    
    # Probabilities of state z(t+1) given state z(t)
    gamma[1,1,t] <- 0
    gamma[1,2,t] <- phi.1[t] * (1 - alpha.1[t])
    gamma[1,3,t] <- 0
    gamma[1,4,t] <- 0
    gamma[1,5,t] <- 0
    gamma[1,6,t] <- phi.1[t] * alpha.1[t]
    gamma[1,7,t] <- 1 - phi.1[t]
    
    gamma[2,1,t] <- 0
    gamma[2,2,t] <- 0
    gamma[2,3,t] <- phi.2[t] * (1 - alpha.2[t])
    gamma[2,4,t] <- 0
    gamma[2,5,t] <- 0
    gamma[2,6,t] <- phi.2[t] * alpha.2[t]
    gamma[2,7,t] <- 1 - phi.2[t]
    
    gamma[3,1,t] <- 0
    gamma[3,2,t] <- 0
    gamma[3,3,t] <- 0
    gamma[3,4,t] <- phi.3[t] * (1 - alpha.3[t])
    gamma[3,5,t] <- 0
    gamma[3,6,t] <- phi.3[t] * alpha.3[t]
    gamma[3,7,t] <- 1 - phi.3[t]
    
    gamma[4,1,t] <- 0
    gamma[4,2,t] <- 0
    gamma[4,3,t] <- 0
    gamma[4,4,t] <- 0
    gamma[4,5,t] <- phi.4[t] * (1 - alpha.4[t])
    gamma[4,6,t] <- phi.4[t] * alpha.4[t]
    gamma[4,7,t] <- 1 - phi.4[t]
    
    gamma[5,1,t] <- 0
    gamma[5,2,t] <- 0
    gamma[5,3,t] <- 0
    gamma[5,4,t] <- 0
    gamma[5,5,t] <- phi.5[t] * (1 - alpha.5[t])
    gamma[5,6,t] <- phi.5[t] * alpha.5[t]
    gamma[5,7,t] <- 1 - phi.5[t]
    
    gamma[6,1,t] <- 0
    gamma[6,2,t] <- 0
    gamma[6,3,t] <- 0
    gamma[6,4,t] <- 0
    gamma[6,5,t] <- 0
    gamma[6,6,t] <- phi.ad
    gamma[6,7,t] <- 1 - phi.ad
    
    gamma[7,1,t] <- 0
    gamma[7,2,t] <- 0
    gamma[7,3,t] <- 0
    gamma[7,4,t] <- 0
    gamma[7,5,t] <- 0
    gamma[7,6,t] <- 0
    gamma[7,7,t] <- 1
    
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
    
    omega[6,1,t] <- 1 - p.B[t]
    omega[6,2,t] <- 0
    omega[6,3,t] <- 0
    omega[6,4,t] <- p.B[t]
    
    omega[7,1,t] <- 1 
    omega[7,2,t] <- 0
    omega[7,3,t] <- 0
    omega[7,4,t] <- 0
    
  }          # t
  
  # likelihood 
  for (i in 1:N) {
    # latent state at first capture - is known
    z[i, first[i]] <- y[i, first[i]] - 1
    for (t in (first[i]+1):K){
      # draw z(t) given z(t-1)
      z[i,t] ~ dcat(gamma[z[i,t-1], 1:7, t-1])
      # draw y(t) given z(t)
      y[i,t] ~ dcat(omega[z[i,t], 1:4, t-1])
    } # t 
  } #i
  
})     # Nimble model

# initial values
initial.values <- function() list(phi.1 = runif(my.constants$K-1,0,1),
                                  phi.3 = runif(my.constants$K-1,0,1),
                                  phi.4 = runif(my.constants$K-1,0,1),
                                  phi.5 = runif(my.constants$K-1,0,1),
                                  phi.ad = runif(1,0,1),
                                  alpha.4 = runif(my.constants$K-1,0,1),
                                  alpha.5 = runif(my.constants$K-1,0,1),
                                  p.NB = runif(my.constants$K-1,0,1),
                                  p.B = runif(my.constants$K-1,0,1), 
                                  z = zinits)


# parameters to save
parameters.to.save <- c("phi.1", "phi.3", "phi.4", "phi.5", "phi.ad", "alpha.4", "alpha.5", "p.NB", "p.B")

# MCMC details 
n.iter <- 150000
n.burnin <- 50000
n.chains <- 3
n.thin <- 20

# run nimble 
system.time({
  mcmc.guille.agem3 <- nimbleMCMC(code = myCode,
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

save(mcmc.guille.agem3, file = "age-model-v3.RData")

#### model no 4 ####
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
  
  for (t in 1:(K-1)) { # loop over time
    
    phi.1[t] ~ dunif(0, 1)
    phi.2[t] <- phi.1[t]                          # make assumption that first and second year survival are identical
    phi.3[t] ~ dunif(0, 1)
    phi.4[t] ~ dunif(0, 1)    
    phi.5[t] ~ dunif(0, 1)
    phi.6[t] ~ dunif(0, 1)
    
    alpha.1[t] <- 0                               # birds cannot begin breeding before age 4
    alpha.2[t] <- 0
    alpha.3[t] <- 0
    alpha.4[t] ~ dunif(0, 1)
    alpha.5[t] ~ dunif(0, 1)
    alpha.6[t] ~ dunif(0, 1)
    
    p.NB[t] ~ dunif(0, 1)
    p.B[t] ~ dunif(0, 1)
    
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
    gamma[3,4,t] <- phi.3[t] * (1 - alpha.3[t])
    gamma[3,5,t] <- 0
    gamma[3,6,t] <- 0
    gamma[3,7,t] <- phi.3[t] * alpha.3[t]
    gamma[3,8,t] <- 1 - phi.3[t] 
    
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
initial.values <- function() list(phi.1 = runif(my.constants$K-1,0,1),
                                  phi.3 = runif(my.constants$K-1,0,1),
                                  phi.4 = runif(my.constants$K-1,0,1),
                                  phi.5 = runif(my.constants$K-1,0,1),
                                  phi.6 = runif(my.constants$K-1,0,1),
                                  phi.ad = runif(1,0,1),
                                  alpha.4 = runif(my.constants$K-1,0,1),
                                  alpha.5 = runif(my.constants$K-1,0,1),
                                  alpha.6 = runif(my.constants$K-1,0,1),
                                  p.NB = runif(my.constants$K-1,0,1),
                                  p.B = runif(my.constants$K-1,0,1), 
                                  z = zinits)

# parameters to save
parameters.to.save <- c("phi.1", "phi.3", "phi.4", "phi.5", "phi.6", "phi.ad", "alpha.4", "alpha.5", "alpha.6", "p.NB", "p.B")

# MCMC details 
n.iter <- 150000
n.burnin <- 50000
n.chains <- 3
n.thin <- 20

# run nimble 
system.time({
  mcmc.guille.agem4 <- nimbleMCMC(code = myCode,
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

save(mcmc.guille.agem4, file = "age-model-v4.RData")

#### model no 5 ####

# initial values
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
  
  for (t in 1:(K-1)) { # loop over time
    
    phi.1[t] ~ dunif(0, 1)
    phi.2[t] <- phi.1[t]                          # make assumption that first and second year survival are identical
    phi.3[t] ~ dunif(0, 1)
    phi.4[t] ~ dunif(0, 1)    
    phi.5[t] ~ dunif(0, 1)
    phi.6[t] ~ dunif(0, 1)
    
    alpha.1[t] <- 0                               # birds cannot begin breeding before age 4
    alpha.2[t] <- 0
    alpha.3[t] <- 0
    alpha.4[t] ~ dunif(0, 1)
    alpha.5[t] <- alpha.4[t]
    alpha.6[t] ~ dunif(0, 1)
    
    p.NB[t] ~ dunif(0, 1)
    p.B[t] ~ dunif(0, 1)
    
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
    gamma[3,4,t] <- phi.3[t] * (1 - alpha.3[t])
    gamma[3,5,t] <- 0
    gamma[3,6,t] <- 0
    gamma[3,7,t] <- phi.3[t] * alpha.3[t]
    gamma[3,8,t] <- 1 - phi.3[t] 
    
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
initial.values <- function() list(phi.1 = runif(my.constants$K-1,0,1),
                                  phi.3 = runif(my.constants$K-1,0,1),
                                  phi.4 = runif(my.constants$K-1,0,1),
                                  phi.5 = runif(my.constants$K-1,0,1),
                                  phi.6 = runif(my.constants$K-1,0,1),
                                  phi.ad = runif(1,0,1),
                                  alpha.4 = runif(my.constants$K-1,0,1),
                                  alpha.6 = runif(my.constants$K-1,0,1),
                                  p.NB = runif(my.constants$K-1,0,1),
                                  p.B = runif(my.constants$K-1,0,1), 
                                  z = zinits)

# parameters to save
parameters.to.save <- c("phi.1", "phi.3", "phi.4", "phi.5", "phi.6", "phi.ad", "alpha.4", "alpha.6", "p.NB", "p.B")

# MCMC details 
n.iter <- 150000
n.burnin <- 50000
n.chains <- 3
n.thin <- 20

# run nimble 
system.time({
  mcmc.guille.agem5 <- nimbleMCMC(code = myCode,
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

save(mcmc.guille.agem5, file = "age-model-v5.RData")

#### model no 6 ####

# initial values
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
  
  for (t in 1:(K-1)) { # loop over time
    
    phi.1[t] ~ dunif(0, 1)
    phi.2[t] <- phi.1[t]                          # make assumption that first and second year survival are identical
    phi.3[t] ~ dunif(0, 1)
    phi.4[t] ~ dunif(0, 1)    
    phi.5[t] ~ dunif(0, 1)
    phi.6[t] <- phi.5[t]                          # assume that survival in adult non-breeders doesn't vary with age
    
    alpha.1[t] <- 0                               # birds cannot begin breeding before age 4
    alpha.2[t] <- 0
    alpha.3[t] <- 0
    alpha.4[t] ~ dunif(0, 1)                      # but recruitment varies between ages 4, 5, and 6+
    alpha.5[t] ~ dunif(0, 1)
    alpha.6[t] ~ dunif(0, 1)
    
    p.NB[t] ~ dunif(0, 1)
    p.B[t] ~ dunif(0, 1)
    
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
    gamma[3,4,t] <- phi.3[t] * (1 - alpha.3[t])
    gamma[3,5,t] <- 0
    gamma[3,6,t] <- 0
    gamma[3,7,t] <- phi.3[t] * alpha.3[t]
    gamma[3,8,t] <- 1 - phi.3[t] 
    
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
initial.values <- function() list(phi.1 = runif(my.constants$K-1,0,1),
                                  phi.3 = runif(my.constants$K-1,0,1),
                                  phi.4 = runif(my.constants$K-1,0,1),
                                  phi.5 = runif(my.constants$K-1,0,1),
                                  phi.ad = runif(1,0,1),
                                  alpha.4 = runif(my.constants$K-1,0,1),
                                  alpha.5 = runif(my.constants$K-1,0,1),
                                  alpha.6 = runif(my.constants$K-1,0,1),
                                  p.NB = runif(my.constants$K-1,0,1),
                                  p.B = runif(my.constants$K-1,0,1), 
                                  z = zinits)

# parameters to save
parameters.to.save <- c("phi.1", "phi.3", "phi.4", "phi.5", "phi.ad", "alpha.4", "alpha.5", "alpha.6", "p.NB", "p.B")

# MCMC details 
n.iter <- 150000
n.burnin <- 50000
n.chains <- 3
n.thin <- 20

# run nimble 
system.time({
  mcmc.guille.agem6 <- nimbleMCMC(code = myCode,
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

save(mcmc.guille.agem6, file = "age-model-v6.RData")

#### model no 7 ####

# initial values
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
  
  for (t in 1:(K-1)) { # loop over time
    
    phi.1[t] ~ dunif(0, 1)
    phi.2[t] <- phi.1[t]                          # make assumption that first and second year survival are identical
    phi.3[t] ~ dunif(0, 1)
    phi.4[t] ~ dunif(0, 1)    
    phi.5[t] <- phi.4[t]
    phi.6[t] <- phi.4[t]                          # assume that survival in adult non-breeders doesn't vary with age
    
    alpha.1[t] <- 0                               # birds cannot begin breeding before age 4
    alpha.2[t] <- 0
    alpha.3[t] <- 0
    alpha.4[t] ~ dunif(0, 1)                      # but recruitment varies between ages 4, 5, and 6+
    alpha.5[t] ~ dunif(0, 1)
    alpha.6[t] ~ dunif(0, 1)
    
    p.NB[t] ~ dunif(0, 1)
    p.B[t] ~ dunif(0, 1)
    
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
    gamma[3,4,t] <- phi.3[t] * (1 - alpha.3[t])
    gamma[3,5,t] <- 0
    gamma[3,6,t] <- 0
    gamma[3,7,t] <- phi.3[t] * alpha.3[t]
    gamma[3,8,t] <- 1 - phi.3[t] 
    
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
initial.values <- function() list(phi.1 = runif(my.constants$K-1,0,1),
                                  phi.3 = runif(my.constants$K-1,0,1),
                                  phi.4 = runif(my.constants$K-1,0,1),
                                  phi.ad = runif(1,0,1),
                                  alpha.4 = runif(my.constants$K-1,0,1),
                                  alpha.5 = runif(my.constants$K-1,0,1),
                                  alpha.6 = runif(my.constants$K-1,0,1),
                                  p.NB = runif(my.constants$K-1,0,1),
                                  p.B = runif(my.constants$K-1,0,1), 
                                  z = zinits)

# parameters to save
parameters.to.save <- c("phi.1", "phi.3", "phi.4", "phi.ad", "alpha.4", "alpha.5", "alpha.6", "p.NB", "p.B")

# MCMC details 
n.iter <- 150000
n.burnin <- 50000
n.chains <- 3
n.thin <- 20

# run nimble 
system.time({
  mcmc.guille.agem7 <- nimbleMCMC(code = myCode,
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

save(mcmc.guille.agem7, file = "age-model-v7.RData")

#### model no 8 ####

# initial values
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
  
  for (t in 1:(K-1)) { # loop over time
    
    phi.1[t] ~ dunif(0, 1)
    phi.2[t] <- phi.1[t]                          # make assumption that first and second year survival are identical
    phi.4[t] ~ dunif(0, 1)    
    phi.5[t] ~ dunif(0, 1)
    phi.6[t] <- phi.5[t]                          # assume that survival in adult non-breeders doesn't vary with age
    
    alpha.1[t] <- 0                               # birds cannot begin breeding before age 4
    alpha.2[t] <- 0
    alpha.3[t] <- 0
    alpha.4[t] ~ dunif(0, 1)                      # but recruitment varies between ages 4, 5, and 6+
    alpha.5[t] ~ dunif(0, 1)
    alpha.6[t] ~ dunif(0, 1)
    
    p.NB[t] ~ dunif(0, 1)
    p.B[t] ~ dunif(0, 1)
    
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
initial.values <- function() list(phi.1 = runif(my.constants$K-1,0,1),
                                  phi.3 = runif(1,0,1),
                                  phi.4 = runif(my.constants$K-1,0,1),
                                  phi.5 = runif(my.constants$K-1,0,1),
                                  phi.ad = runif(1,0,1),
                                  alpha.4 = runif(my.constants$K-1,0,1),
                                  alpha.5 = runif(my.constants$K-1,0,1),
                                  alpha.6 = runif(my.constants$K-1,0,1),
                                  p.NB = runif(my.constants$K-1,0,1),
                                  p.B = runif(my.constants$K-1,0,1), 
                                  z = zinits)

# parameters to save
parameters.to.save <- c("phi.1", "phi.3", "phi.4", "phi.5", "phi.ad", "alpha.4", "alpha.5", "alpha.6", "p.NB", "p.B")

# MCMC details 
n.iter <- 150000
n.burnin <- 50000
n.chains <- 3
n.thin <- 20

# run nimble 
system.time({
  mcmc.guille.agem8 <- nimbleMCMC(code = myCode,
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

save(mcmc.guille.agem8, file = "age-model-v8.RData")

#### model no 9 ####

# initial values
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
  alpha.4 ~ dunif(0, 1)
  
  for (t in 1:(K-1)) { # loop over time
    
    phi.1[t] ~ dunif(0, 1)
    phi.2[t] <- phi.1[t]                          # make assumption that first and second year survival are identical
    phi.4[t] ~ dunif(0, 1)    
    phi.5[t] ~ dunif(0, 1)
    phi.6[t] <- phi.5[t]                          # assume that survival in adult non-breeders doesn't vary with age
    
    alpha.1[t] <- 0                               # birds cannot begin breeding before age 4
    alpha.2[t] <- 0
    alpha.3[t] <- 0                               # but recruitment varies between ages 4, 5, and 6+
    alpha.5[t] ~ dunif(0, 1)
    alpha.6[t] ~ dunif(0, 1)
    
    p.NB[t] ~ dunif(0, 1)
    p.B[t] ~ dunif(0, 1)
    
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
    gamma[4,5,t] <- phi.4[t] * (1 - alpha.4)
    gamma[4,6,t] <- 0
    gamma[4,7,t] <- phi.4[t] * alpha.4
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
initial.values <- function() list(phi.1 = runif(my.constants$K-1,0,1),
                                  phi.3 = runif(1,0,1),
                                  phi.4 = runif(my.constants$K-1,0,1),
                                  phi.5 = runif(my.constants$K-1,0,1),
                                  phi.ad = runif(1,0,1),
                                  alpha.4 = runif(1,0,1),
                                  alpha.5 = runif(my.constants$K-1,0,1),
                                  alpha.6 = runif(my.constants$K-1,0,1),
                                  p.NB = runif(my.constants$K-1,0,1),
                                  p.B = runif(my.constants$K-1,0,1), 
                                  z = zinits)

# parameters to save
parameters.to.save <- c("phi.1", "phi.3", "phi.4", "phi.5", "phi.ad", "alpha.4", "alpha.5", "alpha.6", "p.NB", "p.B")

# MCMC details 
n.iter <- 150000
n.burnin <- 50000
n.chains <- 3
n.thin <- 20

# run nimble 
system.time({
  mcmc.guille.agem9 <- nimbleMCMC(code = myCode,
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

save(mcmc.guille.agem9, file = "age-model-v9.RData")

#### compare with WAIC ####
mcmc.guille.agem1$WAIC$WAIC
mcmc.guille.agem2$WAIC$WAIC
mcmc.guille.agem3$WAIC$WAIC
mcmc.guille.agem4$WAIC$WAIC
mcmc.guille.agem5$WAIC$WAIC
mcmc.guille.agem6$WAIC$WAIC
mcmc.guille.agem7$WAIC$WAIC
mcmc.guille.agem8$WAIC$WAIC
mcmc.guille.agem9$WAIC$WAIC

#### code to plot figure 4 ####
library(tidyverse)
library(cowplot)
library(MCMCvis)
library(magick)

res8 <- as.data.frame(MCMCsummary(mcmc.guille.agem8$samples, round = 3))
res8 <- res8 %>% rename(q2.5 = `2.5%`, q50 = `50%`, q97.5 = `97.5%`)

# save parameters separately
juvphi <- as.data.frame(MCMCsummary(mcmc.ms.agev8.constant3$samples, round = 3, params = "phi.1"))  # save summary
age4 <- as.data.frame(MCMCsummary(mcmc.ms.agev8.constant3$samples, round = 3, params = "phi.4"))  
age5plus <- as.data.frame(MCMCsummary(mcmc.ms.agev8.constant3$samples, round = 3, params = "phi.5")) 
breed4 <- as.data.frame(MCMCsummary(mcmc.ms.agev8.constant3$samples, round = 3, params = "alpha.4"))  
breed5 <- as.data.frame(MCMCsummary(mcmc.ms.agev8.constant3$samples, round = 3, params = "alpha.5"))  
breed6 <- as.data.frame(MCMCsummary(mcmc.ms.agev8.constant3$samples, round = 3, params = "alpha.6")) 
pnb <- as.data.frame(MCMCsummary(mcmc.ms.agev8.constant3$samples, round = 3, params = "p.NB"))  
pb <- as.data.frame(MCMCsummary(mcmc.ms.agev8.constant3$samples, round = 3, params = "p.B")) 

# guillemot population
pop <- read.csv("covariates_for_models_1995-2020.csv")
pop <- pop %>% select(1,9)

### juveniles
juvphi <- juvphi %>% rename(q2.5 = `2.5%`, q50 = `50%`, q97.5 = `97.5%`)
juvphi$year <- seq(1995, 2020, 1)

range(juvphi$mean)
mean(juvphi$mean)

# 95% credible intervals
juvphi$lower <- juvphi$mean - juvphi$q2.5
juvphi$upper <- juvphi$q97.5 - juvphi$mean

jphi <- function() {
  par(mar = c(1.5, 2.4, 0.5, 0.5))
  plot(juvphi$year, juvphi$mean, type = "l", col = "#4B0055", ylim=range(c(0, 1)), pch=19, xlab="", ylab="", 
       xlim = c(1995, 2020), cex.axis = 0.65, las = 1, bty = "n", axes = F)
  axis(1, ylim = c(1995, 2020), tck = 1, mgp = c(0, 0.005, 0), cex.axis = 0.6, col = "grey90")  
  axis(2, ylim = c(0,1), las = 1, tck = 1, mgp = c(0, 0.2, 0), cex.axis = 0.6, col = "grey90") 
  axis(1, ylim = c(1995, 2020), tck = 1, mgp = c(0, 0.005, 0), cex.axis = 0.6, col = "grey90", 
       at = c(1997.5, 2002.5, 2007.5, 2012.5, 2017.5), labels = F, lwd = 0.5)  
  axis(2, ylim = c(0,1), las = 1, tck = 1, mgp = c(0, 0.2, 0), cex.axis = 0.6, col = "grey90", 
       at = c(0.1, 0.3, 0.5, 0.7, 0.9), labels = F, lwd = 0.5) 
  par(new = T)
  plot(juvphi$year, juvphi$mean, type = "l", col = "#4B0055", ylim=range(c(0, 1)), pch=19, xlab="", ylab="", 
       xlim = c(1995, 2020), cex.axis = 0.65, las = 1, bty = "n", axes = F)
  mtext("Juvenile survival", side = 2, line = 1.2, cex = 0.65)        # y axis title
  axis(2, ylim = c(0,1), las = 1, tck = -0.01, mgp = c(0, 0.2, 0), cex.axis = 0.6)   # y axis
  axis(1, ylim = c(1995, 2020), tck = -0.01, mgp = c(0, 0.005, 0), cex.axis = 0.6)     # x-axis
  # errorbars
  arrows(juvphi$year, juvphi$mean-juvphi$lower, juvphi$year, juvphi$mean+juvphi$upper, length=0.0, angle=90, code=3, col = "#4B0055")
  points(juvphi$year, juvphi$mean, pch = 19, col = "#4B0055", cex = 0.5)
  box(bty = "o")
}

ggdraw(jphi)

### age 4
age4 <- age4 %>% rename(q2.5 = `2.5%`, q50 = `50%`, q97.5 = `97.5%`)
age4$year <- seq(1995, 2020, 1)
age4 <- age4 %>% slice(4:26)

# 95% credible intervals
age4$lower <- age4$mean - age4$q2.5
age4$upper <- age4$q97.5 - age4$mean

mean(age4$mean)

phi4 <- function() {
  par(mar = c(1.5, 2.4, 0.5, 0.5))
  plot(age4$year, age4$mean, type = "l", col = "#4B0055", ylim=range(c(0, 1)), pch=19, xlab="", ylab="", 
       xlim = c(1995, 2020), cex.axis = 0.65, las = 1, bty = "n", axes = F)
  axis(1, ylim = c(1995, 2020), tck = 1, mgp = c(0, 0.005, 0), cex.axis = 0.6, col = "grey90")  
  axis(2, ylim = c(0,1), las = 1, tck = 1, mgp = c(0, 0.2, 0), cex.axis = 0.6, col = "grey90") 
  axis(1, ylim = c(1995, 2020), tck = 1, mgp = c(0, 0.005, 0), cex.axis = 0.6, col = "grey90", 
       at = c(1997.5, 2002.5, 2007.5, 2012.5, 2017.5), labels = F, lwd = 0.5)  
  axis(2, ylim = c(0,1), las = 1, tck = 1, mgp = c(0, 0.2, 0), cex.axis = 0.6, col = "grey90", 
       at = c(0.1, 0.3, 0.5, 0.7, 0.9), labels = F, lwd = 0.5) 
  par(new = T)
  plot(age4$year, age4$mean, type = "l", col = "#4B0055", ylim=range(c(0, 1)), pch=19, xlab="", ylab="", 
       xlim = c(1995, 2020), cex.axis = 0.65, las = 1, bty = "n", axes = F)
  axis(2, ylim = c(0,1), las = 1, tck = -0.01, mgp = c(0, 0.2, 0), cex.axis = 0.6)   # y axis
  mtext("Age 4 immature survival", side = 2, col = "black", line = 1.4, cex = 0.65)        # y axis title
  axis(1, ylim = c(1995, 2020), tck = -0.01, mgp = c(0, 0, 0), cex.axis = 0.6)
  # errorbars
  arrows(age4$year, age4$mean-age4$lower, age4$year, age4$mean+age4$upper, length=0.0, angle=90, code=3, col = "#4B0055")
  points(age4$year, age4$mean, pch = 19, col = "#4B0055", cex = 0.5)
  box(bty = "o")
  
}

ggdraw(phi4)

### age 5
age5plus <- age5plus %>% rename(q2.5 = `2.5%`, q50 = `50%`, q97.5 = `97.5%`)
age5plus$year <- seq(1995, 2020, 1)
age5plus <- age5plus %>% slice(5:26)

# 95% credible intervals
age5plus$lower <- age5plus$mean - age5plus$q2.5
age5plus$upper <- age5plus$q97.5 - age5plus$mean

mean(age5plus$mean)

phi5 <- function() { 
  par(mar = c(1.5, 2.4, 0.5, 0.5))
  plot(age5plus$year, age5plus$mean, type = "l", col = "#4B0055", ylim=range(c(0, 1)), pch=19, xlab="", ylab="", 
       xlim = c(1995, 2020), cex.axis = 0.65, las = 1, bty = "n", axes = F)
  axis(1, ylim = c(1995, 2020), tck = 1, mgp = c(0, 0.005, 0), cex.axis = 0.6, col = "grey90")  
  axis(2, ylim = c(0,1), las = 1, tck = 1, mgp = c(0, 0.2, 0), cex.axis = 0.6, col = "grey90") 
  axis(1, ylim = c(1995, 2020), tck = 1, mgp = c(0, 0.005, 0), cex.axis = 0.6, col = "grey90", 
       at = c(1997.5, 2002.5, 2007.5, 2012.5, 2017.5), labels = F, lwd = 0.5)  
  axis(2, ylim = c(0,1), las = 1, tck = 1, mgp = c(0, 0.2, 0), cex.axis = 0.6, col = "grey90", 
       at = c(0.1, 0.3, 0.5, 0.7, 0.9), labels = F, lwd = 0.5) 
  par(new = T)
  plot(age5plus$year, age5plus$mean, type = "l", col = "#4B0055", ylim=range(c(0, 1)), pch=19, xlab="", ylab="", 
       xlim = c(1995, 2020), cex.axis = 0.65, las = 1, bty = "n", axes = F)
  axis(2, ylim = c(0,1), las = 1, tck = -0.01, mgp = c(0, 0.2, 0), cex.axis = 0.6)   # y axis
  mtext("Age 5+ NB survival", side = 2, col = "black", line = 1.2, cex = 0.65)        # y axis title
  axis(1, ylim = c(1995, 2020), tck = -0.01, mgp = c(0, 0.008, 0), cex.axis = 0.6)
  # errorbars
  arrows(age5plus$year, age5plus$mean-age5plus$lower, age5plus$year, age5plus$mean+age5plus$upper, length=0.0, angle=90, code=3, col = "#4B0055")
  points(age5plus$year, age5plus$mean, pch = 19, col = "#4B0055", cex = 0.5)
  box(bty = "o")
}

ggdraw(phi5)

### recruitment age 4
breed4 <- breed4 %>% rename(q2.5 = `2.5%`, q50 = `50%`, q97.5 = `97.5%`)
breed4$year <- seq(1995, 2020, 1)
breed4 <- breed4 %>% slice(4:26)

# 95% credible intervals
breed4$lower <- breed4$mean - breed4$q2.5
breed4$upper <- breed4$q97.5 - breed4$mean

### recruitment age 5
breed5 <- breed5 %>% rename(q2.5 = `2.5%`, q50 = `50%`, q97.5 = `97.5%`)
breed5$year <- seq(1995, 2020, 1)
breed5 <- breed5 %>% slice(5:26)

# 95% credible intervals
breed5$lower <- breed5$mean - breed5$q2.5
breed5$upper <- breed5$q97.5 - breed5$mean

### recruitment age 6+
breed6 <- breed6 %>% rename(q2.5 = `2.5%`, q50 = `50%`, q97.5 = `97.5%`)
breed6$year <- seq(1995, 2020, 1)
breed6 <- breed6 %>% slice(6:26)

# 95% credible intervals
breed6$lower <- breed6$mean - breed6$q2.5
breed6$upper <- breed6$q97.5 - breed6$mean

### all recruitment
allpsi <- function() {
  par(mar = c(1.5, 2.4, 0.5, 0.5))
  # recruitment age 5
  plot(breed5$year, breed5$mean, type = "l", col = "#4B0055", ylim=range(c(0, 1)), pch=19, xlab="", ylab="", 
       xlim = c(1995, 2020), cex.axis = 0.65, las = 1, bty = "n", axes = F)
  axis(1, ylim = c(1995, 2020), tck = 1, mgp = c(0, 0.005, 0), cex.axis = 0.6, col = "grey90")  
  axis(2, ylim = c(0,1), las = 1, tck = 1, mgp = c(0, 0.2, 0), cex.axis = 0.6, col = "grey90") 
  axis(1, ylim = c(1995, 2020), tck = 1, mgp = c(0, 0.005, 0), cex.axis = 0.6, col = "grey90", 
       at = c(1997.5, 2002.5, 2007.5, 2012.5, 2017.5), labels = F, lwd = 0.5)  
  axis(2, ylim = c(0,1), las = 1, tck = 1, mgp = c(0, 0.2, 0), cex.axis = 0.6, col = "grey90", 
       at = c(0.1, 0.3, 0.5, 0.7, 0.9), labels = F, lwd = 0.5) 
  par(new = T)
  plot(breed5$year, breed5$mean, type = "l", col = "#4B0055", ylim=range(c(0, 1)), pch=19, xlab="", ylab="", 
       xlim = c(1995, 2020), cex.axis = 0.65, las = 1, bty = "n", axes = F)
  # errorbars
  arrows(breed5$year, breed5$mean-breed5$lower, breed5$year, breed5$mean+breed5$upper, length=0.0, angle=90, code=3, col = "#4B0055")
  points(breed5$year, breed5$mean, pch = 15, col = "#4B0055", cex = 0.6)
  # recruitment age 6
  points(breed6$year, breed6$mean, type = "l", col = "#009B95", ylim=range(c(0, 1)), pch=19, xlab="", ylab="", 
         xlim = c(1995, 2020), cex.axis = 0.65, las = 1, bty = "n", axes = F)
  # errorbars
  arrows(breed6$year, breed6$mean-breed5$lower, breed6$year, breed6$mean+breed5$upper, length=0.0, angle=90, code=3, col = "#009B95")
  points(breed6$year, breed6$mean, pch = 17, col = "#009B95", cex = 0.6)
  # recruitment age 4
  points(breed4$year, breed4$mean, type = "l", col = "#FDE333", ylim=range(c(0, 1)), pch=19, xlab="", ylab="", 
         xlim = c(1995, 2020), cex.axis = 0.65, las = 1, bty = "n", axes = F)
  # errorbars
  arrows(breed4$year, breed4$mean-breed4$lower, breed4$year, breed4$mean+breed4$upper, length=0.0, angle=90, code=3, col = "#FDE333")
  points(breed4$year, breed4$mean, pch = 19, col = "#FDE333", cex = 0.5)
  # axes
  axis(2, ylim = c(0,1), las = 1, tck = -0.01, mgp = c(0, 0.2, 0), cex.axis = 0.6)   # y axis
  mtext("Recruitment probability", side = 2, col = "black", line = 1.4, cex = 0.65)        # y axis title
  axis(1, ylim = c(1995, 2020), tck = -0.01, mgp = c(0, 0.008, 0), cex.axis = 0.6)
  box(bty = "o")
  
  legend(2015, 0.97, legend = c("age 4", "age 5", "age 6+"), pch = c(19, 15, 17), lty = 1, col = c("#FDE333", "#4B0055", "#009B95"), 
         cex = 0.6, bty = "n")
  
}

ggdraw(allpsi)

### NB resighting
pnb <- pnb %>% rename(q2.5 = `2.5%`, q50 = `50%`, q97.5 = `97.5%`)
pnb$year <- seq(1995, 2020, 1)
pnb <- pnb %>% slice(2:26)

range(pnb$mean)
mean(pnb$mean)

# 95% credible intervals
pnb$lower <- pnb$mean - pnb$q2.5
pnb$upper <- pnb$q97.5 - pnb$mean

### B resighting
pb <- pb %>% rename(q2.5 = `2.5%`, q50 = `50%`, q97.5 = `97.5%`)
pb$year <- seq(1995, 2020, 1)
pb <- pb %>% slice(6:26)

range(pb$mean)
mean(pb$mean)

# 95% credible intervals
pb$lower <- pb$mean - pb$q2.5
pb$upper <- pb$q97.5 - pb$mean

### resighting together

sighting <- function() {
  par(mar = c(1.5, 2.4, 0.5, 0.5))
  plot(pnb$year, pnb$mean, type = "l", col = "#4B0055", ylim=range(c(0, 1)), pch=19, xlab="", ylab="", 
       xlim = c(1995, 2020), cex.axis = 0.65, las = 1, bty = "n", axes = F)
  axis(1, ylim = c(1995, 2020), tck = 1, mgp = c(0, 0.005, 0), cex.axis = 0.6, col = "grey90")  
  axis(2, ylim = c(0,1), las = 1, tck = 1, mgp = c(0, 0.2, 0), cex.axis = 0.6, col = "grey90") 
  axis(1, ylim = c(1995, 2020), tck = 1, mgp = c(0, 0.005, 0), cex.axis = 0.6, col = "grey90", 
       at = c(1997.5, 2002.5, 2007.5, 2012.5, 2017.5), labels = F, lwd = 0.5)  
  axis(2, ylim = c(0,1), las = 1, tck = 1, mgp = c(0, 0.2, 0), cex.axis = 0.6, col = "grey90", 
       at = c(0.1, 0.3, 0.5, 0.7, 0.9), labels = F, lwd = 0.5) 
  par(new = T)
  plot(pnb$year, pnb$mean, type = "l", col = "#4B0055", ylim=range(c(0, 1)), pch=19, xlab="", ylab="", 
       xlim = c(1995, 2020), cex.axis = 0.65, las = 1, bty = "n", axes = F)
  axis(side = 1, xlim = c(1996, 2020), cex.axis = 0.6, tck = -0.01, mgp = c(0, 0.008, 0))
  mtext("Year", side = 1, col = "black", line = 0.5, cex = 0.65)
  axis(side = 2, ylim = c(0, 1), cex.axis = 0.6, tck = -0.01, mgp = c(0, 0.2, 0), las = 1)
  mtext("Resighting", side = 2, col = "black", line = 1.2, cex = 0.65,)
  # errorbars
  arrows(pnb$year, pnb$mean-pnb$lower, pnb$year, pnb$mean+pnb$upper, length=0.0, angle=90, code=3, col = "#4B0055")
  points(pnb$year, pnb$mean, pch = 19, col = "#4B0055", cex = 0.5)
  box(bty = "o")
  
  lines(pb$year, pb$mean, type = "l", col = "#009B95")
  arrows(pb$year, pb$mean-pb$lower, pb$year, pb$mean+pb$upper, length=0.0, angle=90, code=3, col = "#009B95")
  points(pb$year, pb$mean, pch = 17, col = "#009B95", cex = 0.6)
  
  legend(2016, 0.25, legend = c("B", "NB"), pch = c(17, 19), lty = 1, col = c("#009B95", "#4B0055"), cex = 0.6, 
         bty = "n")
  
}

ggdraw(sighting)

### plot population trend

population <- function() {
  par(mar = c(1.5, 2.4, 0.5, 0.5))
  plot(pop$year, pop$count, type = "l", col = "#4B0055", pch=19, xlab="", ylab="", 
       xlim = c(1995, 2020), ylim = c(6000,30000), cex.axis = 0.65, las = 1, bty = "n", axes = F)
  axis(1, ylim = c(1995, 2020), tck = 1, mgp = c(0, 0.005, 0), cex.axis = 0.6, col = "grey90")  
  axis(2, ylim = c(0,1), las = 1, tck = 1, mgp = c(0, 0.2, 0), cex.axis = 0.6, col = "grey90") 
  axis(1, ylim = c(1995, 2020), tck = 1, mgp = c(0, 0.005, 0), cex.axis = 0.6, col = "grey90", 
       at = c(1997.5, 2002.5, 2007.5, 2012.5, 2017.5), labels = F, lwd = 0.5) 
  axis(2, ylim = c(6000,30000), las = 1, tck = 1, mgp = c(0, 0.2, 0), cex.axis = 0.6, col = "grey90", 
       at = c(7500, 12500, 17500, 22500,27500), labels = F, lwd = 0.5) 
  par(new = T)
  plot(pop$year, pop$count, type = "l", col = "#4B0055", pch=19, xlab="", ylab="", 
       xlim = c(1995, 2020), ylim = c(6000,30000), cex.axis = 0.65, las = 1, bty = "n", axes = F)
  axis(side = 1, xlim = c(1995, 2020), cex.axis = 0.6, tck = -0.01, mgp = c(0, 0.008, 0))
  mtext("Year", side = 1, col = "black", line = 0.5, cex = 0.65)
  axis(side = 2, cex.axis = 0.6, tck = -0.01, mgp = c(0, 0.2, 0), las = 1, ylim=c(6000,30000), labels = F)
  mtext("Population size", side = 2, col = "black", line = 1.7, cex = 0.65,)
  box(bty = 'o')
  
}

ggdraw(population)

(all <- plot_grid(jphi, phi4, phi5, allpsi, sighting,population, nrow = 3, labels = c("A", "B", "C", "D", "E", "F"), label_size = 10))

(all_final <- ggdraw() +
    draw_plot(all) +
    theme(panel.background = element_rect(fill = "white"))+
    panel_border(color = "white"))

save_plot("figure4_2026.pdf", all_final, base_height = 6, base_width = 7, dpi = 600)
save_plot("figure4_2026.png", all_final, base_height = 6, base_width = 7, dpi = 600)