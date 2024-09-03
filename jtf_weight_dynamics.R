#---------------------------------------
# Ultrasonography as a non-invasive technique for assessing diet effects
#   on fish gonad dynamics
#
# Joaquim Tomas-Ferrer, Irene Moro-Martinez, Enrique Massuti-Pascual,
#   Amalia Grau, Miquel Palmer
#---------------------------------------

# Last update: 21 August 2024

# Data an code for estimating parameters of equations 1 to 3

#-----
# 1: Loading data and libraries
#-----

remove(list=ls())
library(cmdstanr)
load("input.RData")

# n_fish: Number of fish (32)
# n_obs: Number of ultrasonographies (335)
# gonad_weight [n_obs] gonad weight (gr) as calculated from ultrasonographies
# ID_fish [n_obs] ID for defining which ultrasonographies correspond to which fish
# time [n_obs] days after Jan 1th each ultrasonography has been done
# n_tank: number of cages (6)
# ID_tank [n_fish] ID for defining which cage correspond to which fish
# food_intake [n_fish]: Average food intake of each fish

#-----       
# 2: STAN model
#-----
sink("model.stan")
cat(" // first line
functions {
}

data{
  int n_obs;                          // number of observations
  int n_fish;                         // number of fish
  int n_tank;                         // number of cages
  array [n_fish] real food_intake;    // avergage food intake of each fish
  array [n_fish] int ID_tank;         // cage ID for each fish
  array [n_obs] real gonad_weight;    // calculated gonad_weight from ultrasonography 
  array [n_obs] real time;            // days from Jan 1th
  array [n_obs] int ID_fish;          // fish ID
  
  // priors
  array [2] real prior_gamma0;
  array [2] real prior_gamma1;
  array [2] real prior_mu_sd;
  array [2] real prior_mu_cage_sd;
  array [2] real prior_beta0;
  array [2] real prior_beta1;
  array [2] real prior_sigma_sd;
  array [2] real prior_sigma_cage_sd;
  array [2] real prior_alpha0;
  array [2] real prior_alpha1;
  array [2] real prior_h_sd;
  array [2] real prior_h_cage_sd;
  array [2] real prior_sd_within;
}

transformed data {
}

parameters {

  // parameters for mu (date at wich the maximum gonad weigth is attained)
  real gamma0;
  real gamma1;                      
  real <lower=0> mu_sd;   
  real <lower=0> mu_cage_sd;
  array [n_tank] real mu_cage;
  array [n_fish] real mu;
  
  // parameters for sigma (spread of the reproductive period)
  real beta0;
  real beta1;
  real <lower=0> sigma_sd;
  real <lower=0> sigma_cage_sd;
  array [n_tank] real sigma_cage;
  array [n_fish] real  <lower=0> sigma;

  // parameters for h (maximum gonad weight)
  real alpha0;
  real alpha1;
  real <lower=0> h_sd;
  real <lower=0> h_cage_sd;
  array [n_tank] real h_cage;
  array [n_fish] real h;

  real <lower=0> sd_within;  // obervation-level error
}

transformed parameters {
  array [n_fish] real mu_hat;
  array [n_fish] real sigma_hat;
  array [n_fish] real h_hat;
  for (i in 1:n_fish){
    mu_hat[i] = gamma0+gamma1*food_intake[i];   // expected date at wich the maximum gonad weigth is attained
    sigma_hat[i] = beta0+beta1*food_intake[i];  // expected spread of the reproductive period
    h_hat[i] = alpha0+alpha1*food_intake[i];    // expected maximum gonad weight
  }

}   

model {

  // mu (date at wich the maximum gonad weigth is attained)
  gamma0 ~ normal(prior_gamma0[1] , prior_gamma0[2]);
  gamma1 ~ normal(prior_gamma1[1] , prior_gamma1[2]);
  mu_sd ~ normal(prior_mu_sd[1] , prior_mu_sd[2]);
  mu_cage_sd ~ normal(prior_mu_cage_sd[1] , prior_mu_cage_sd[2]);
  
  // sigma (spread of the reproductive period)
  beta0 ~ normal(prior_beta0[1] , prior_beta0[2]);
  beta1 ~ normal(prior_beta1[1] , prior_beta1[2]);
  sigma_sd ~ normal(prior_sigma_sd[1] , prior_sigma_sd[2]);
  sigma_cage_sd ~ normal(prior_sigma_cage_sd[1] , prior_sigma_cage_sd[2]);
  
  // h (maximum gonad weight)
  alpha0 ~ normal(prior_alpha0[1] , prior_alpha0[2]);
  alpha1 ~ normal(prior_alpha1[1] , prior_alpha1[2]);
  h_sd ~ normal(prior_h_sd[1] , prior_h_sd[2]);
  h_cage_sd ~ normal(prior_h_cage_sd[1] , prior_h_cage_sd[2]);

  // error at the obervation level
  sd_within ~ normal(prior_sd_within[1] , prior_sd_within[2]);
  
  // random effects at the cage level
  for (i in 1:n_tank){ // cage level
    mu_cage[i] ~ normal(0 , mu_cage_sd);
    sigma_cage[i] ~ normal(0 , sigma_cage_sd);
    h_cage[i] ~ normal(0 , h_cage_sd); 
  } 
 
  // combining fixed and random effects
  for (i in 1:n_fish){ // fish level 
    mu[i] ~ normal(mu_hat[i]+mu_cage[ID_tank[i]] , mu_sd);  
    sigma[i] ~ normal(sigma_hat[i]+sigma_cage[ID_tank[i]] , sigma_sd);
    h[i] ~ normal(h_hat[i]+h_cage[ID_tank[i]] , h_sd); 
  } 
  
  // likelihood
  for (i in 1:n_obs){ //observation level
     gonad_weight[i] ~ normal(h[ID_fish[i]]*exp(-0.5*((time[i]-mu[ID_fish[i]])/sigma[ID_fish[i]])^2), sd_within);

  }
}

generated quantities {
}

" # end of model code
,fill = TRUE)
sink()

#-----
# 3: Compiling model
#-----
mod = cmdstan_model("model.stan")

#-----
# 4: Initializing the model
#-----
# Setting priors
prior_gamma0=c(0,100)
prior_gamma1=c(0,100)
prior_mu_sd=c(0,100)
prior_mu_cage_sd=c(0,100)

prior_beta0=c(0,100)
prior_beta1=c(0,100)
prior_sigma_sd=c(0,100)
prior_sigma_cage_sd=c(0,100)

prior_alpha0=c(100,100)
prior_alpha1=c(0,100)
prior_h_sd=c(100,100)
prior_h_cage_sd=c(100,100)

prior_sd_within=c(100,100)

# Initial values
n.chains = 3 # number of chains
initializer = function() list(
  "gamma0"=30,
  "gamma1"=1,
  "mu"=rep(50,n_fish),
  "mu_cage"=rep(0,n_tank),
  "mu_sd"=10,
  "mu_cage_sd"=10,
  
  "beta0"=30,
  "beta1"=0.1,
  "sigma"=rep(50,n_fish),
  "sigma_cage"=rep(0,n_tank),
  "sigma_sd"=10,
  "sigma_cage_sd"=10,
  
  "alpha0"=100,
  "alpha1"=10,
  "h"=rep(1000,n_fish),
  "h_cage"=rep(0,n_tank),
  "h_sd"=100,
  "h_cage_sd"=100,
  
  "sd_within"=50
)
inits = list()
for (chain in 1:n.chains) inits[[chain]] = initializer()

#-----
# 5: Running
#-----
# predictions
fit = mod$sample(
  data =list (
    n_obs=n_obs,
    n_fish=n_fish,
    n_tank=n_tank,
    ID_tank=ID_tank,
    food_intake=food_intake,
    ID_fish=ID_fish,
    time=time,
    gonad_weight=gonad_weight,
    
    #priors
    prior_gamma0=prior_gamma0,
    prior_gamma1=prior_gamma1,
    prior_mu_sd=prior_mu_sd,
    prior_mu_cage_sd=prior_mu_cage_sd,
    
    prior_beta0=prior_beta0,
    prior_beta1=prior_beta1,
    prior_sigma_sd=prior_sigma_sd,
    prior_sigma_cage_sd=prior_sigma_cage_sd,
    
    prior_alpha0=prior_alpha0,
    prior_alpha1=prior_alpha1,
    prior_h_sd=prior_h_sd,
    prior_h_cage_sd=prior_h_cage_sd,
    
    prior_sd_within=prior_sd_within
  ),
  chains = n.chains,
  parallel_chains = n.chains,
  iter_warmup = 2000,
  iter_sampling = 10000,
  init = inits,
  max_treedepth = 12,
  adapt_delta = 0.9999
)

#-----
# 6: Saving
#-----
# fit$save_object(file = "out.RDS")
# sink(file = "diagnose.txt")
# fit$cmdstan_diagnose()
# sink()
# file.rename("model.stan","model.R")
#-----

#-----
# 7: Results summary (Table 1)
#-----
# fit = readRDS("out.RDS")# Reloading existing definitive results
# table=fit$summary(c("gamma0","gamma1","mu_sd","mu_cage_sd",
#                     "beta0","beta1","sigma_sd","sigma_cage_sd",
#                     "alpha0","alpha1","h_sd","h_cage_sd",
#                     "sd_within"))
# table=data.frame(table)[,c(1,6,3,7,8,9)]
# table[,c(2,3,4,5)]=round(table[,c(2,3,4,5)],5)
# table[,6]=round(table[,6],0)
# table[,2:5]=round(table[,2:5],4)
# table
#--------------

