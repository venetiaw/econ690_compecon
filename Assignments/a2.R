library(readstata13)
library(dplyr)
setwd("Dropbox/comp_econ/Assignments")
dat = read.dta13("guiide_master_file_anonymized_full.dta")

#############################
# [1] BASIC CHARACTERISTICS #
#############################
#1
dim(dat)
sink("summary.txt")
summary(dat)
sink()

#2 
varnames = names(dat)
varnames[grep("id", varnames)]
subs = dat$anon_studentid
length(unique(subs)) == length(subs)
#unique id: anon_studentid 

###################
# [2] SMALL BASES #
###################
#1
cond1 = dat %>% select(anon_studentid, varnames[grep("stud_base", varnames)])

#2
cond2 = dat %>% select(anon_studentid, varnames[grep("alladmin", varnames)])

#3
cond3 = dat %>% select(-starts_with("stud_base"), -starts_with("alladmin"))

###################
# [3] CONSISTENCY #
###################
#1
gender1 = names(cond1)[grep("gender", names(cond1), ignore.case = T)]
gender2 = names(cond2)[grep("gender", names(cond2), ignore.case = T)]
gender3 = names(cond3)[grep("gender", names(cond3), ignore.case = T)]

sub1 = cond1 %>% select(stud_base_GENDER)
sub2 = cond2 %>% select(gender2)
sub3 = cond3 %>% select(gender3)
genders = data.frame(sub1, sub2, sub3) #no they are not consistent

#2
names(cond1)[grep("age", names(cond1), ignore.case = T)] #stud_base_age
age1 = cond1 %>% select(stud_base_age)
names(cond2)[grep("age", names(cond2), ignore.case = T)] #alladmin_age
age2 = cond2 %>% select(alladmin_age)
names(cond3)[grep("age", names(cond3), ignore.case = T)] #stud_fu1_age
age3 = cond3 %>% select(stud_fu1_age)
ages = data.frame(age1, age2, age3) #no they are not consistent
ages %>% filter(stud_base_age == alladmin_age) %>% filter(stud_fu1_age == stud_base_age) %>% dim #1861 are consistent
dim(ages)[[1]] - 1861 #41113

#3
varnames[grep("BECE", varnames, ignore.case = T)]
expectations = dat %>% select(stud_base_bece_best, stud_base_bece_worst, stud_base_bece_likely)
inconsistent = expectations %>% filter(stud_base_bece_best < stud_base_bece_worst)
overconf = expectations %>% filter(stud_base_bece_best < stud_base_bece_likely)
underconf = expectations %>% filter(stud_base_bece_worst > stud_base_bece_likely)

###############
# [4] BALANCE #
###############
varnames[grep("treat", varnames, ignore.case = T)] #alladmin_treatgroup
dat$stud_base_GENDER = dat$stud_base_GENDER %>% recode("M"=0, "F"=1)
dat$stud_base_GENDER = as.integer(dat$stud_base_GENDER) 
dat$SHSregioncode = as.integer(dat$SHSregioncode)
dat %>% select(alladmin_treatgroup, stud_base_GENDER, stud_base_age, SHSregioncode) %>% group_by(alladmin_treatgroup) %>% summarise(mean_gender=mean(stud_base_GENDER, na.rm = T)) #not balanced
dat %>% select(alladmin_treatgroup, stud_base_GENDER, stud_base_age, SHSregioncode) %>% group_by(alladmin_treatgroup) %>% summarise(mean_age=mean(stud_base_age, na.rm = T)) #relatively balanced
dat %>% select(alladmin_treatgroup, stud_base_GENDER, stud_base_age, SHSregioncode) %>% group_by(alladmin_treatgroup) %>% summarise(mean_code=mean(SHSregioncode, na.rm = T)) #relatively balanced

##############################
# [5] RECODING AND HISTOGRAM #
##############################
#1
varnames[grep("educ|want", varnames, ignore.case = T)]
dat$stud_base_educ_want = as.factor(dat$stud_base_educ_want)
levels(dat$stud_base_educ_want) = c("Junior high school", "technical or vocational training", "senior high school", "nursing or teacher training", "polytechnic", "university")
#2
dat$stud_base_educ_want = as.numeric(dat$stud_base_educ_want)
tib1 = dat %>% filter(stud_base_GENDER == 1) %>% select(stud_base_educ_want) 
hist(tib1$stud_base_educ_want)
tib2 = dat %>% filter(stud_base_GENDER == 0) %>% select(stud_base_educ_want) 
hist(tib2$stud_base_educ_want)

#############################
# [6] MANIPULATING THE DATA #
#############################

#1
varnames[grep("mychoice", varnames, ignore.case = T)]
df = data.frame(names = names(table(dat$stud_base_mychoice_1_pgm)), pgm1 = as.vector(table(dat$stud_base_mychoice_1_pgm)), pgm2 = as.vector(table(dat$stud_base_mychoice_2_pgm)), pgm3 = as.vector(table(dat$stud_base_mychoice_3_pgm)), pgm4 = as.vector(table(dat$stud_base_mychoice_4_pgm))) %>% filter(names != "")

#2
varnames[grep("region", varnames, ignore.case = T)]
regionChoice = data.frame(names=names(table(dat$stud_base_mychoice_1_region)),reg1 = as.vector(table(dat$stud_base_mychoice_1_region)),reg2 = as.vector(table(dat$stud_base_mychoice_2_region)),reg3 = as.vector(table(dat$stud_base_mychoice_3_region)),reg4 = as.vector(table(dat$stud_base_mychoice_4_region))) %>% filter(names!="")

#3
varnames[grep("mychoice|bece", varnames, ignore.case = T)]
# [9] "stud_base_mychoice_1_bece"       
# [10] "stud_base_mychoice_2_bece"       
# [11] "stud_base_mychoice_3_bece"       
# [12] "stud_base_mychoice_4_bece"
pct25 = sapply(dat %>% select(stud_base_mychoice_1_bece, stud_base_mychoice_2_bece, stud_base_mychoice_3_bece, stud_base_mychoice_4_bece) %>% na.omit(), quantile)[2,1:4]
plot(c(1:4), pct25)
pct50 = sapply(dat %>% select(stud_base_mychoice_1_bece, stud_base_mychoice_2_bece, stud_base_mychoice_3_bece, stud_base_mychoice_4_bece) %>% na.omit(), quantile)[3,1:4]
plot(c(1:4), pct50)
pct75 = sapply(dat %>% select(stud_base_mychoice_1_bece, stud_base_mychoice_2_bece, stud_base_mychoice_3_bece, stud_base_mychoice_4_bece) %>% na.omit(), quantile)[4,1:4]
plot(c(1:4), pct7s5)
plot(c(1:4), as.vector(dat %>% select(stud_base_mychoice_1_bece, stud_base_mychoice_2_bece, stud_base_mychoice_3_bece, stud_base_mychoice_4_bece) %>% na.omit() %>% summarise_all(list(sd))))

require(tidyr)
availDat = dat %>% drop_na(stud_base_mychoice_1_bece, stud_base_mychoice_2_bece, stud_base_mychoice_3_bece, stud_base_mychoice_4_bece)
aggScore = availDat[sample.int(dim(availDat)[[1]],size=5),1:dim(availDat)[[2]]]
fivepeople = as.data.frame(cbind(matrix(c(aggScore$stud_base_mychoice_1_bece, aggScore$stud_base_mychoice_2_bece, aggScore$stud_base_mychoice_3_bece, aggScore$stud_base_mychoice_4_bece), nrow=4, byrow = T), c(1:4)))
lm(V6 ~., fivepeople) #in general, as rank of choice increases, expected score decreases

reverting = as.integer(availDat$stud_base_mychoice_1_bece >= availDat$stud_base_mychoice_2_bece)*as.integer(availDat$stud_base_mychoice_2_bece >= availDat$stud_base_mychoice_3_bece)*as.integer(availDat$stud_base_mychoice_3_bece >= availDat$stud_base_mychoice_4_bece)

##############
# [7] PROBIT #
##############
#Create inconsistent, overconfident, and underconfident variables
newDat = availDat %>% mutate(inconsistent = stud_base_bece_best < stud_base_bece_worst, overconf = stud_base_bece_best < stud_base_bece_likely, underconf = stud_base_bece_worst > stud_base_bece_likely)
probit=glm(reverting ~ stud_base_educ_want + stud_base_age + GENDER + inconsistent + overconf + underconf, family = binomial(link = "probit"), 
                data = newDat)

