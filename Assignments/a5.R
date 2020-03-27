# setwd("Dropbox/comp_econ/Assignments/")

dat = read.csv("dat_choices.csv")[,2:26]
l11 = cbind(seq(48,16,-8),seq(48, 112, 16))
l12 = cbind(seq(48,24,-6), seq(48,72,6))
l13 = cbind(seq(48,8,-10), seq(48,104,14))
l14 = cbind(seq(42,18,-6), seq(42,114,18))
l15 = cbind(seq(54, 14, -10), seq(54,110,14))
l1 = rbind(l11, l12, l13, l14, l15)
l21 = c(40,32,24,16,8,42,36,30,24,16,38,28,18,8,0,36,30,24,18,10,44,34,24,14,6)
l22 = c(64,80,96,112,120,66,84,102,120,128,62,76,90,104,112,60,78,96,114,122,68,82,96,110,118)
l2 = cbind(l21,l22)
choice = cbind(l1,l2)

##################
# EX 1: PART 1.1 #
##################
sigma = 0.5
grid = seq(0.1, 2, 0.1)
TOTAL = 1000
# dim(theta_draws): TOTAL*length(grid)
# Initialize the 1000 guesses for each of the guessed thetas
theta_draws = sapply(X=grid, FUN=rnorm, n=TOTAL, sd=0.5)

crra = function(c, theta, eps = 0.01) {
  ifelse(abs(theta-1)<eps,log(c),(c^(1-theta))/(1-theta))
}
EU = function(bundle, theta) {
  return(0.5*(crra(bundle[1], theta)+crra(bundle[2], theta)))
}
freq_sim = function(choicerow, total, guesses) {
  bundle1 = choicerow[1:2]
  bundle2 = choicerow[3:4]
  simulated_probabilities = vector(length = total)
  for (i in 1:total) {
    curr_theta = guesses[i]
    delta = EU(bundle1, theta=curr_theta) -EU(bundle2, theta=curr_theta)
    simulated_probabilities[i] = ifelse(delta<0, 0, 1)
  }
  prob = mean(simulated_probabilities)
  return(prob)
  
}

curr_choice = choice[20,] #do this with all rows of curr_choice
prob = apply(X=theta_draws, FUN=freq_sim, choicerow=curr_choice, total=TOTAL, MARGIN=2)
plot(grid, prob)

##################
# EX 1: PART 1.2 #
##################
library(nloptr)
i120 = dat[120,]
i280 = dat[280,]
i1200 = dat[1200,]
#draw from ~N(0,1)
theta_draws = rnorm(TOTAL)
mod_draw = function(mods, initial_draw=theta_draws) {
  #mods = c(theta_mod, sigma_mod)
  return(mods[1] + initial_draw*mods[2])
}
freq_sim = function(choicerow, total, guesses) {
  bundle1 = choicerow[1:2]
  bundle2 = choicerow[3:4]
  simulated_probabilities = vector(length = total)
  for (i in 1:total) {
    curr_theta = guesses[i]
    delta = EU(bundle1, theta=curr_theta) -EU(bundle2, theta=curr_theta)
    simulated_probabilities[i] = ifelse(delta<0, 0, 1)
  }
  prob = mean(simulated_probabilities)
  return(prob)
}
getLikelihood = function(guesses, total, individual_decisions) {
  prob = apply(X=choice, FUN=freq_sim, total=TOTAL, guesses, MARGIN = 1)
  prob[prob<0.000001] = 0.000001
  prob[prob>0.999999] = 0.999999
  loglik = individual_decisions * log(prob) + (1-individual_decisions)*log(1-prob)
  return(-sum(loglik))
}
likelihood_wrapper = function(guess_modifier, total, individual_decisions) {
  #guess_modifier = c(theta_mod, sigma_mod)
  return(getLikelihood(guesses = mod_draw(guess_modifier), total, individual_decisions))
}
#Individual 120
isres(runif(2), fn=likelihood_wrapper, total=TOTAL, individual_decisions=i120, lower=c(-5,0.0000001), upper=c(5,5))
bobyqa(runif(2), fn=likelihood_wrapper, total=TOTAL, individual_decisions=i120, lower=c(-5,0.0000001))
# This returned:
# -> (0.2961638, 0.8221510), 48.98868 <-
# -> (-1.367849, 4.698385), 19.12612
# -> (0.6808150, 0.8813539), 47.25325

#Individual 280
isres(runif(2), fn=likelihood_wrapper, total=TOTAL, individual_decisions=i280, lower=c(-5,0.0000001), upper=c(5,5))
bobyqa(runif(2), fn=likelihood_wrapper, total=TOTAL, individual_decisions=i280, lower=c(-5,0.0000001))
# This returned:
# -> (0.2277925, 0.6152068), 52.04506 <-
# -> (-1.519861, 3.986422), 17.20104
# -> (-1.016562, 2.616728), 19.80454

#Individual 1200
isres(runif(2), fn=likelihood_wrapper, total=TOTAL, individual_decisions=i1200, lower=c(-5,0.0000001), upper=c(5,5))
bobyqa(runif(2), fn=likelihood_wrapper, total=TOTAL, individual_decisions=i1200, lower=c(-5,0.0000001))
# This returned:
# -> (0.199710, 0.408939), 53.80315 <-
# -> (-2.457150, 6.158071), 14.57715
# -> (0.1103986, 0.5446302), 48.47889
