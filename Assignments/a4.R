library(pracma)
setwd("Dropbox/comp_econ/assignments")
dat = read.csv("dat_choices.csv")

##########
# PART 1 #
##########
#1.1
crra = function(c, theta, eps = 0.01) {
  ifelse(abs(theta-1)<eps,log(c),(c^(1-theta))/(1-theta))
}

#1.2
GRIDSPACE = 0.0001
myTheta = seq(-2,2,GRIDSPACE)
myU = sapply(myTheta, crra, c=10)
plot(myTheta, myU)

##########
# PART 2 #
##########
#2.1
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

#2.2
EU = function(bundle, theta) {
  return(0.5*(crra(bundle[1], theta)+crra(bundle[2], theta)))
}
find_theta = function(choiceSet, a=-6, b=6, tolerance = 10^(-10), max_iter = 100) {
  bundle1 = choiceSet[1:2]
  bundle2 = choiceSet[3:4]
  fa = EU(bundle1, theta=a) -EU(bundle2, theta=a) 
  fb = EU(bundle1, theta=b) - EU(bundle2, theta=b)
  stopifnot(fa*fb <=0)
  count = 0
  while (count <= max_iter) {
    c = 0.5*(a+b)
    fc = EU(bundle1, theta=c) - EU(bundle2, theta=c) 
    if (abs(fc) < tolerance) return(c)
    if (abs(c-a) < tolerance) return(c)
    if (fa*fc < 0) {
      b = c
      fb = fc
    } else {
      a = c
      fa = fc
    }
    count = count + 1
  }
}
thetas = apply(FUN = find_theta, X = choice, MARGIN = 1)

#2.3: identified sets for each individual
identify_individual = function(individualDecision, choiceSet=choice) {
  #to choose 0, theta is below the indifference
  #to choose 1, theta is above the indifference
  i0 = which(myDes == 0)
  i1 = which(myDes == 1)
  if (length(i0) == 0) {
    myMin = -Inf
  } else {myMin = min(thetas[i0])}
  if (length(i1) == 0) {
    myMax = Inf
  } else {myMax = max(thetas[i1])}
  return(c(myMin, myMax))
}
#get distribution for each individual
dist = c()
for (i in 1:dim(dat)[[1]]) {
  myDes = dat[i,2:26]
  dist = cbind(dist, identify_individual(myDes, choice))
}


##########
# PART 3 #
##########
library(evd)
wealth=20
#3.1
getDist = function(choicerow, theta, w=wealth) {
  eu = crra(choicerow + w, rep(theta, 4))
  #expected utility of choice 0
  e0 = EU(choicerow[1:2]+w, theta)
  #expected utility of choice 1
  e1 = EU(choicerow[3:4]+w, theta)
  #p(choose choice 1) = p(e1>e0)
  p = pgumbel(e1-e0)
  return(p)
}

likelihood = function(ydat, choiceSet, theta, w=wealth) {
  p = apply(X = choiceSet, MARGIN = 1, FUN = getDist, theta=theta)
  p[p<0.000001] = 0.000001
  p[p>0.999999] = 0.999999
  loglik = ydat * log(p) + (1-ydat)*log(1-p)
  return(sum(loglik))
}
# likelihood(ydat=dat[900,], theta=0.5, choiceSet=choice)
#3.2
grid = linspace(-10,10,2000)
#for individual 900
i900 = dat[900,]
vals900 = sapply(X=grid, FUN = likelihood, choiceSet=choice, ydat=i900)
grid[which(vals900==min(vals900))]
#for individual 115
i115 = dat[115,]
vals115 = sapply(X=grid, FUN = likelihood, choiceSet=choice, ydat=i115)
grid[which(vals115==min(vals115))]


##########
# PART 4 #
##########

beale = function(x,y) {
  return((1.5-x+x*y)^2 + (2.25-x+x*y^2)^2 + (2.625-x+x*y^3)^2)
}
#4.1
xgrid = linspace(-5,5,10/0.01)
ygrid = linspace(-5,5,10/0.01)
#structure such that beale(xgrid[2], ygrid[1]) = bruteforce[2,1] and so on
bruteforce = sapply(X = ygrid, FUN = beale, x = xgrid)
x0y0 = c(which(bruteforce==min(bruteforce), arr.ind = T)[1,1], which(bruteforce==min(bruteforce), arr.ind = T)[1,2])
ix = x0y0[1]
iy = x0y0[2]
x = xgrid[ix]
y = ygrid[iy]
#4.2
beale = function(vals) {
  x = vals[1]
  y=vals[2]
  return((1.5-x+x*y)^2 + (2.25-x+x*y^2)^2 + (2.625-x+x*y^3)^2)
}
fprime_beale = function(vals, eps=10^(-10)) {
  out = NULL
  for (i in 1:length(vals)) {
    z1=z2=vals
    z1[i] = z1[i] + eps
    z2[i] = z2[i] - eps
    out[i] = beale(z1) - beale(z2)
  }
  return(out/2/eps)
}
fprime_beale(c(x,y))

##########
# PART 5 #
##########
f = function(vals) {
  x = vals[1]
  y = vals[2]
  return((1-x)^2+5*(y-x^2)^2)
}

fprime = function(vals, eps=10^(-6)) {
  out = NULL
  for (i in 1:length(vals)) {
    z1=z2=vals
    z1[i] = z1[i] + eps
    z2[i] = z2[i] - eps
    out[i] = f(z1) - f(z2)
  }
  return(out/2/eps)
}

line_search = function(vals, alpha) {
  x = vals[1]
  y = vals[2]
  primes = fprime(vals)
  xp = primes[1]
  yp = primes[2]
  newf = f(c(x-alpha*xp, y-yp))
  return(newf)
}

steepest_descent1 = function(start=c(0,0), tol=10^(-10), max_iter = 100000) {
  i=1
  alphaVec = c(1, 0.1, 0.01, 0.0001)
  difference = tol + 1
  while (i<max_iter) {
    deriv = fprime(start)
    #find alpha that minimizes value
    val = sapply(X = alphaVec, FUN = line_search, vals=start)
    alpha = alphaVec[which(val==min(val))]
    oldStart=start
    start = start - alpha*deriv
    i = i+1
    difference = f(start) - f(oldStart)
    if (abs(difference) < tol) {return(start)}
  }
  print("Reached max-iter: not optimized")
  return(start)
}




##########
# PART 6 #
##########
library(dplyr)
library(AER)
library(stats)
library(mclogit)
library(sjmisc)
#6.1
data("SmokeBan")
smoke = SmokeBan
smoker = as.integer(smoke$smoker == "yes")
smoke$smoker = smoker
ban = as.integer(smoke$ban == "yes")
smoke$ban = ban
education = to_dummy(smoke$education)[,1:4]
smoke$education=NULL
smoke = cbind(smoke, education)
afam = as.integer(smoke$afam == "yes")
smoke$afam = afam
hispanic = as.integer(smoke$hispanic == "yes")
smoke$hispanic = hispanic
gender = as.integer(smoke$gender == "male")
smoke$gender = gender

probit = glm(smoker ~ ., family = binomial(link='probit'), data = smoke)
betaProbit = as.matrix(coefficients(probit))
xProbit = smoke %>% select(-smoker) %>% mutate(coef = 1)
xProbit = as.matrix(subset(xProbit, select=c(length(names(xProbit)), 1:length(names(xProbit))-1)))
yProbit = as.matrix(smoke$smoker)

smoking_likelihood = function(beta,y=yProbit, x=xProbit) {
  p = pnorm(x%*% beta)
  p[p<0.000001] = 0.000001
  p[p>0.999999] = 0.999999
  L = sum(y*log(p)+(1-y)*log(1-p))
  return(L)
}

out = optim(rnorm(10), fn=smoking_likelihood, method="BFGS",
            control = list(maxit = 30000,trace = TRUE,REPORT = 1),
            y=yProbit,x=xProbit)
myBeta = out$par

#6.2

fisher = -hessian(smoking_likelihood, myBeta)
se = 1/sqrt(fisher)


##########
# PART 7 #
##########
library(nloptr)
individual = dat[5,]
#7.1
getDist = function(s, t, choicerow, w=wealth) {
  #t: theta
  #s: sigma
  #expected utility of choice 0
  e0 = EU(choicerow[1:2]+w, t)
  #expected utility of choice 1
  e1 = EU(choicerow[3:4]+w, t)
  #p(choose choice 1) = p(e1>e0)
  p = pnorm(e1-e0, sd = s)
  return(p)
}

likelihood_norm = function(params, ydat, choiceSet=choice, w=wealth) {
  sigma = params[1]
  theta = params[2]
  p = apply(X = choiceSet, MARGIN = 1, FUN = getDist, s=sigma, t=theta)
  p[p<0.000001] = 0.000001
  p[p>0.999999] = 0.999999
  loglik = ydat * log(p) + (1-ydat)*log(1-p)
  return(sum(loglik))
}

#7.2
mySigma=0.8 #recall: sigma cannot be negative
myTheta=2.1
myParams = c(mySigma,myTheta)
myLower = c(0,-2)
myUpper = c(3,3)
isres(myParams, likelihood_norm, lower=myLower, upper=myUpper, ydat=individual, maxeval = 500000) #doesn't converge: [1] 0.01025624 0.10739206
lbfgs(myParams, likelihood_norm, lower=myLower, ydat=individual) #stuck at local
bobyqa(myParams, likelihood_norm, ydat=individual, lower=myLower) #stuck at local but gives best of locals
neldermead(myParams, likelihood_norm, ydat=individual, lower=myLower, ydat=individual) #stuck at local

#bobyqa for full sample
est = function(person, start=myParams, choiceSet = choice, w=wealth) {
  params = isres(start, likelihood_norm, ydat=person, lower=myLower, upper=myUpper, maxeval=500000)
  return(params$par)
}
dist = apply(X = dat, MARGIN = 1, FUN = est)
hist(dist[2,])

##########
# PART 8 #
##########
#8.1
riskGrid = seq(0.01,2,0.01)
getDist = function(t, s, choicerow, w=wealth) {
  #t: theta
  #s: sigma
  eu = crra(choicerow + w, rep(t, 4))
  #expected utility of choice 0
  e0 = EU(choicerow[1:2]+w, t)
  #expected utility of choice 1
  e1 = EU(choicerow[3:4]+w, t)
  #p(choose choice 1) = p(e1>e0)
  p = pnorm(e1-e0, sd = s)
  return(p)
}
prob = c() #each column corresponds to the distribution for a risk grid
for (i in 1:nrow(choice)) {
  prob = cbind(prob,sapply(FUN=getDist, X=riskGrid, s=0.5, choicerow=choice[i,]))
}
plot(riskGrid,prob[,5])

#8.2
set.seed(123)

simulateData = function(t, s, choicerow, w=wealth) {
  eps = rnorm(2,0,s)
  u1 = EU(choicerow[1:2]+w, t) +eps[1]
  u2 = EU(choicerow[3:4]+w, t) +eps[2]
  ysim = ifelse(u1<=u2,0,1)
  return(ysim)
}

myLikelihood = function(params, ydat, choiceSet=choice, w=wealth) {
  theta = params[1]
  sigma = params[2]
  p = apply(X = choiceSet, MARGIN = 1, FUN = getDist, s=sigma, t=theta)
  p[p<0.000001] = 0.000001
  p[p>0.999999] = 0.999999
  loglik = sum(ydat * log(p) + (1-ydat)*log(1-p))
  return(loglik)
}

est = function(person, start=myParams, choiceSet = choice, w=wealth) {
  params = bobyqa(start, myLikelihood, ydat=person, lower=c(-10,0))
  # print(params)
  return(params$par)
}


monte_carlo = function(test) {
  #@param test=c(theta, sigma)
  #ysim: simulated data for an individual
  #5% CI
  start = c(0,0)
  allThetas = c()
  allSigmas = c()
  for (count in 1:100) {
    ysim = c()
    for (i in 1:25) {
      ysim = c(ysim, simulateData(test[1], test[2], choice[i,]))
    }
    myParams = est(person = ysim, start)
    # myParams = est(ysim, test)
    allThetas = c(allThetas, myParams[1])
    allSigmas = c(allSigmas, myParams[2])
  }
  mean_theta = mean(allThetas)
  sd_theta = sd(allThetas)
  ci_theta = c(mean_theta-1.96*sd_theta, mean_theta+1.96*sd_theta)
  mean_sigma = mean(allSigmas)
  sd_sigma = sd(allSigmas)
  ci_sigma = c(mean_sigma-1.96*sd_sigma, mean_sigma+1.96*sd_sigma)
  return(data.frame(ci_theta, ci_sigma))
}


#8.3
grid=1:8
r = c(0.5,0.5,0.8,0.8,1.5,1.5,2.8,2.8)
s = rep(c(0.1,0.9),4)
param_grid = data.frame(grid, r, s)
monte_grid = data.frame(matrix(ncol = 5, nrow = 0))
colnames(monte_grid) =c('grid', 'theta_lower', 'theta_upper', 'sigma_lower', 'sigma_upper')

for (i in 1:nrow(param_grid)) {
  pair = c(param_grid[i,]$r, param_grid[i,]$s)
  myMonte = monte_carlo(pair)
  monte_grid[i,] = c(i, myMonte[1,1], myMonte[2,1], myMonte[1,2], myMonte[2,2])
}

##############
# EXERCISE 9 #
##############

p1 = rep(75,6*8)
t1 = rep(c(rep(1,6), rep(8,6), rep(31,6), rep(91,6)),2)
p2 = rep(c(75.31, 75.63, 76.25, 78.13, 81.25, 87.50), 8)
t2 = c(rep(31,6*2), rep(61, 6), rep(121, 6), rep(361,6), rep(368,6), rep(391,6), rep(451,6))
table2 = data.frame(p1, t1, p2, t2)


#1
getDist = function(beta, delta, sigma, theta, choicerow, w=wealth) {
  choice1 = choicerow[[1]]
  time1 = choicerow[[2]]
  choice2 = choicerow[[3]]
  time2 = choicerow[[4]]
  #expected utility of choice 0
  if (time1 == 0) {
    util = crra(w, theta) - crra(w+choice1, theta) + beta*theta^(time2)*(crra(w+choice2, theta)-crra(w, theta))
  }
  else {
    util = beta*theta^(time1)*(crra(w, theta)-crra(w+choice1, theta))+ beta*theta^(time2)*(crra(w+choice2, theta)-crra(w, theta))
  }
  p = pnorm(rnorm(1, 0, sigma)-util)
  return(p)
}

myLikelihood = function(params, th, ydat, choiceSet=table2, w=wealth) {
  b = params[1] #beta
  d = params[2] #delta
  s = params[3] #sigma
  #th: theta
  p = apply(X = choiceSet, MARGIN = 1, FUN = getDist, beta=b, delta=d, sigma=s, theta=th)
  p[p<0.000001] = 0.000001
  p[p>0.999999] = 0.999999
  loglik = ydat * log(p) + (1-ydat)*log(1-p)
  return(sum(loglik))
}

#2
thetas = dist[2,]
dat_time = read.csv("dat_time.csv")
estimates = data.frame(matrix(ncol = 4, nrow = 0))
colnames(estimates) =c('individual', 'beta', 'delta', 'sigma')
for (i in 1:dim(dat_time)[[1]]) {
  print(i)
  individual = dat_time[i, 2:49]
  estimates[i,] = c(i, bobyqa(c(0.1,0.1,0.1), myLikelihood, lower=c(0.0001,0.0001,0.0001), upper=c(1,1,5), th=thetas[i], ydat=individual)$par)
}
standardDeviations = apply(X=estimates[,2:4], FUN=sd, MARGIN = 2)

#3
cor(thetas, estimates$beta) #-0.5247573
cor(thetas, estimates$delta) #-0.02230964
cor(thetas, estimates$sigma) #-0.4540421
