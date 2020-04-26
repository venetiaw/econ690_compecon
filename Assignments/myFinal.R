setwd("Dropbox/comp_econ/Assignments/")
library(sjmisc)
library(AER)
library(stats)
library(dplyr)
library(nloptr)
data("GSOEP9402")
dat = GSOEP9402
inc = dat$income
#=================================
# [1]: KERNEL DENSITY ESTIMATION
#=================================
#1
band = function(X) {
  n = length(X)
  sigma = sd(X)
  return(((4*sigma^5)/(3*n))^0.2)
}

#2
gaussian = function(X) {
  return(exp(-(X^2)/2)/sqrt(2*pi))
}

#3
kdens = function(x0, X) {
  h = band(X)
  tmp = (x0-X)/h
  return(sum(gaussian(tmp))/(length(tmp)*h))
}

#=====================
# [2]: EMPIRICAL CDF 
#=====================

myEcdf = function(x0, X) {
  return(mean(X < x0))
}

#=============================
# [3]: NUMERICAL INTEGRATION
#=============================

trapzf = function(c, d, n, f, myX) {
  h  = (d-c)/n
  inner = c(c, seq(1, n-1) * h + c, d)
  inner = sapply(X=inner, FUN=f, myX)
  inner[1] = inner[1]/2
  inner[length(inner)] = inner[length(inner)]/2
  return(sum(inner)*h)
}

#========================
# [4]: RESERVATION WAGES
#========================

lambda_u = 2
lambda_e = 0.5
delta = 0.02
r = 0.01

#1
fcdf = function(x0, income_dist) {
  G = myEcdf(x0, income_dist)
  Fbar = ((delta/G)-delta)/lambda_e
  myF = 1-Fbar
  myF[myF<0] = 0
  return(myF)
}

fpdf = function(x0, income_dist) {
  g = kdens(x0, income_dist)
  G = myEcdf(x0, income_dist)
  return(delta*lambda_e*g/((lambda_e*G)^2))
}


#2
get_inner = function(x0, X) {
  F_bar = 1-fcdf(x0, X)
  return(F_bar/(r+delta+lambda_e*F_bar))
}

get_nonlinear = function(psi, b, X=inc, n=100) {
  myInt = trapzf(psi, max(inc), n, get_inner, X)
  return(b + (lambda_u-lambda_e)*myInt - psi)
}

get_reservationWage = function(b, lower=0, upper=max(inc)) {
  return(uniroot(get_nonlinear, b, interval=c(lower, upper))$root)
}

#==========================
# [5]: INDIRECT INFERENCE
#==========================
#1
marriage = dat %>% to_dummy(marital, suffix = 'label') %>% select(-marital_married)
covariates = dat %>% select(income, kids, meducation) 
covariates = cbind(covariates, marriage)
myModel = lm(income~., data=covariates)
out = summary(myModel)
theta_data = out$coefficients[,1:2]

#2.1
observed_attributes = covariates %>% select(-income) 
observed_attributes = cbind(1, observed_attributes)
mat = as.matrix(observed_attributes)

get_eps = function(sigma, num_sim) {
  #each row in eps_mat corresponds to one eps vector
  eps_mat = matrix(rnorm(nrow(mat)*num_sim, mean=0, sd=sigma), nrow = num_sim)
  return(eps_mat)
}

do_sim = function(S, beta_guess, sigma_guess) {
  # S: number of simulations
  # betas: vector of length m, where:
  # - m: number of covariates
  # sigma: scalar
  # Returns: reservation wage for each individual
  m = ncol(mat)
  n = nrow(mat)
  eps_mat = get_eps(sigma_guess, S) 
  myStore = matrix(nrow=n, ncol = S)
  for (i in 1:S) {
    eps = eps_mat[i,]
    #recovery
    b0 = mat %*% beta_guess + eps
    phi = sapply(X=b0, FUN=get_reservationWage)
    myStore[,i] = phi
  }
  return(apply(myStore, FUN=mean, MARGIN=1))
}


#2.2

indirect_inference = function(S, beta_guess, sigma_guess) {
  myPhi = do_sim(S, beta_guess, sigma_guess)
  tmp = cbind(observed_attributes, myPhi)
  myModel = lm(myPhi~., data=tmp[,2:dim(tmp)[2]])
  out = summary(myModel)
  theta_model = out$coefficients[,1:2]
  return(theta_model)
}
# indirect_inference(2, 1:7, 1)

#3
to_min = function(curr_theta_model) {
  diff = (curr_theta_model[,1] - theta_data[,1])^2
  return(sum(diff))
}

main = function(S, main_guess) {
  beta_guess = main_guess[1:length(main_guess)-1]
  sigma_guess = main_guess[length(main_guess)]
  theta_model = indirect_inference(S, beta_guess, sigma_guess)
  val = to_min(theta_model)
  return(val)
}

get_start = function() {
  x = runif(8)
  x[1:3] = x[1:3] * 1000
  return(x)
}

#sometimes you get errors because of the starting points
#just rerun it, will be fine
main_guess = get_start()
optimized_params = bobyqa(main_guess, fn=main, S=2)$par

#4
eachBoostrap = function(s) {
  #1
  base_sample = covariates[sample(1:nrow(covariates), nrow(covariates), replace=T),]
  myModel = lm(income~., data=base_sample)
  out = summary(myModel)
  theta_data = out$coefficients[,1:2]
  #2-3
  observed_attributes = base_sample %>% select(-income)
  observed_attributes = cbind(1, observed_attributes)
  mat = as.matrix(observed_attributes)
  theta_model = bobyqa(get_start(), fn=main, S=s)$par
  return(theta_model[1:7]) #note: last item in array is sigma
}

myBoostrap = function(N, S) {
  #each row corresponds to a boostrap iteration
  myStore = matrix(nrow=N, ncol = 7)
  i=1
  while (i <= N) {
    supposed_output = try(eachBoostrap(S))
    if (class(supposed_output) == 'numeric') {
      myStore[i,] = supposed_output
      i = i + 1
    }
    else {next}
  }
  return(myStore)
}
testThis = myBoostrap(N=100,S=1000)
sd = apply(X=testThis, FUN=mean, MARGIN = 2)

