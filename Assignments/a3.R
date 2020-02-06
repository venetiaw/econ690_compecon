library(AER)
library(ggplot2)
library(sjmisc)
library(mclogit)
library(dplyr)
data("Affairs")
dat = Affairs

#########################
# [1] BASIC DESCRIPTION #
#########################
#1
isFemale = dat %>% to_dummy(gender, suffix ='label') %>% select(female = gender_female)
dat = dat %>% mutate(female = isFemale[[1]])
hasChildren = dat %>% to_dummy(children, suffix ='label') %>% select(children = children_yes)
dat = dat %>% mutate(child = hasChildren[[1]])
dat %>% select(-c(gender, children)) %>% cor

#2
ggplot(dat) + geom_histogram(aes(x=affairs, color=gender, fill=gender))
ggplot(dat, aes(x=affairs, color=as.factor(age), fill=as.factor(age))) + geom_histogram(binwidth = 1)
ggplot(dat, aes(x=affairs, color=as.factor(religiousness), fill=as.factor(religiousness))) + geom_histogram(binwidth = 1)
ggplot(dat, aes(x=affairs, color=as.factor(education), fill=as.factor(education))) + geom_histogram(binwidth = 1)

########################
# [2] LINEAR REGRESSION#
########################
#1
lm(affairs ~ education+female + occupation+age+child+religiousness+yearsmarried, data = dat)
#2
lm(affairs ~ education+female + occupation + age*female+age+child+religiousness+yearsmarried, data = dat)
lm(affairs ~ education+female + occupation+child*age+age+child+religiousness+yearsmarried, data = dat)

##############
# [3] PROBIT #
##############
#1
dat = dat %>% mutate(affairs_dummy = as.numeric(affairs > 0)) 

#2
affairs_probit = function(y, beta, x) {
  L = prod(pnorm(x%*%beta, 0,1)^y * (1-(pnorm(x%*%beta, 0,1)))^(1-y))
  return(L)
}

#3
educDummies = dat %>% to_dummy(education) %>% select(V1, V2, V3, V4, V5, V6)
names(educDummies) = c('e1', 'e2', 'e3', 'e4', 'e5', 'e6')
religionDummy = dat %>% to_dummy(religiousness) %>% select(-V5)
names(religionDummy) = c('r1', 'r2', 'r3', 'r4')
satisfaction = dat %>% to_dummy(rating) %>% select(-V5)
names(satisfaction) = c('s1', 's2', 's3', 's4')
d = cbind(affair = dat$affairs_dummy, educDummies, female = dat$female, child = dat$child, religionDummy, years = dat$yearsmarried, satisfaction)

#4
probit = glm(affair ~ ., family = binomial(link='probit'), data = d)
betaProbit = as.matrix(coefficients(probit))
xProbit = d %>% select(-affair) %>% mutate(coef = 1)
xProbit = as.matrix(subset(xProbit, select=c(length(names(xProbit)), 1:length(names(xProbit))-1)))
yProbit = as.matrix(d$affair)

#5
affairs_probit(yProbit, betaProbit, xProbit)
exp(logLik(probit))
#the two values are equal

#############
# [4] LOGIT #
#############

logitCDF = function(x) {
  return(exp(x)/(1+exp(x)))
}

affairs_logit = function(y, beta, x) {
  L = prod(logitCDF(x%*%beta)^y * (1-(logitCDF(x%*%beta)))^(1-y))
  return(L)
}

logit = glm(affair ~ ., family = binomial(link='logit'), data = d)
betaLogit = as.matrix(coefficients(logit))
xLogit = xProbit
yLogit = yProbit
affairs_logit(yLogit, betaLogit, xLogit)
exp(logLik(logit))
#they are equal

#########################
# [5] MULTINOMIAL LOGIT #
#########################
dat = dat %>% select(-affairs_dummy,-gender, -children) 
nAffairs = to_dummy(dat$affairs)
names(nAffairs) = c('A0', 'A1', 'A2', 'A3', 'A7', 'A12')
dat = cbind(nAffairs, dat) 
dat$affairs=NULL
dat = cbind(dat, educDummies, religionDummy, satisfaction)
dat = dat %>% select(-religiousness,-education, -rating) 
occ = to_dummy(dat$occupation)[,1:6]
names(occ) = c('o1', 'o2', 'o3', 'o4', 'o5', 'o6')
dat = cbind(dat, occ) %>% select(-occupation)

f1 = as.formula(paste("cbind(A0, A1, A2, A3, A7, A12) ~", paste(names(dat)[c(7:30)], collapse="+")))
x = as.matrix(dat[,7:ncol(dat)])
x = cbind("coef"=1, x)
y = as.matrix(nAffairs)

#1: Conditional logit
#no solution

#2: Multinomial logit
mnlogit = mblogit(f1, data = dat)
betalogit = matrix(coefficients(mnlogit), ncol=5, byrow = T)
betalogit = cbind(0, betalogit)

multinomial_prob = function(y, beta, x) {
  p = exp(x%*%beta)/sum(exp(x%*%beta))
  L = sum(y*log(p))
  return(L)
}

multinomial_prob(y, betalogit, x)


#3: mixed logit
#no solution




