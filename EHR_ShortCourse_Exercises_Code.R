## R code for Analysis of Big Healthcare Databases exercises

##-------------------------------------------------------------##
## Read in data
##-------------------------------------------------------------##

#### Encounter data
encounter = read.csv("https://raw.githubusercontent.com/rhubb/ASA_EHR_ShortCourse/master/data/encounter.csv", head=T)
 
#### Prescription medication data
meds = read.csv("https://raw.githubusercontent.com/rhubb/ASA_EHR_ShortCourse/master/data/meds.csv", head=T)

#### Measures data
measures = read.csv("https://raw.githubusercontent.com/rhubb/ASA_EHR_ShortCourse/master/data/measures.csv", head=T)

#### Validation data
validation = read.csv("https://raw.githubusercontent.com/rhubb/ASA_EHR_ShortCourse/master/data/validation.csv", head=T)

##-------------------------------------------------------------##
## Install R packages
##-------------------------------------------------------------##

install.packages("rpart")
install.packages("pROC")
install.packages("boot")
install.packages("gee")

library(rpart)
library(pROC)
library(boot)
library(gee)

##-------------------------------------------------------------##
## Exercises
##-------------------------------------------------------------##

##-------------------------------------------------------------
## 1. Data Quality Evaluation

## Use summary statistics and plots to investigate basic characteristics of the data

summary(measures)

# check types of measures available
unique(measures$measuretype)

# separate variables by measurement type
height <- measures[measures$measuretype == "height",-4]
names(height) <- c("patientid","servicedate","height")

weight <- measures[measures$measuretype == "weight",-4]
names(weight) <- c("patientid","servicedate","weight")

glucose <- measures[measures$measuretype == "glucose",-4]
names(glucose) <- c("patientid","servicedate","glucose")

hba1c <- measures[measures$measuretype == "hba1c",-4]
names(hba1c) <- c("patientid","servicedate","hba1c")

chol <- measures[measures$measuretype == "chol",-4]
names(chol) <- c("patientid","servicedate","chol")

# explore number of observations available per patient for each measurement type

summary(c(table(factor(height$patientid, levels = unique(encounter$patientid)))))
summary(c(table(factor(weight$patientid, levels = unique(encounter$patientid)))))
summary(c(table(factor(glucose$patientid, levels = unique(encounter$patientid)))))
summary(c(table(factor(hba1c$patientid, levels = unique(encounter$patientid)))))
summary(c(table(factor(chol$patientid, levels = unique(encounter$patientid)))))

# number of children with no measures available

sum(c(table(factor(height$patientid, levels = unique(encounter$patientid)))) == 0)
sum(c(table(factor(weight$patientid, levels = unique(encounter$patientid)))) == 0)
sum(c(table(factor(glucose$patientid, levels = unique(encounter$patientid)))) == 0)
sum(c(table(factor(hba1c$patientid, levels = unique(encounter$patientid)))) == 0)
sum(c(table(factor(chol$patientid, levels = unique(encounter$patientid)))) == 0)

# summarize distribution of variables across all patients, looking for values outside the plausible range

summary(height$height)
summary(weight$weight)
summary(glucose$glucose)
summary(hba1c$hba1c)
summary(chol$chol)

# values that are clearly outside the plausible range can be eliminated, those that seem unlikely should be
# noted for discussion with clinical collaborators

# remove negative weights as clearly lying outside the plausible range

weight$weight <- ifelse(weight$weight < 0, NA, weight$weight)

# identify patients with extreme heights and weights

extreme.heights <- height$patientid[height$height < 100] # flag patients with height < 1 m
extreme.weights <- weight$patientid[weight$weight > 200] # flag patients with weight > 200 kg

# implausible patterns in longitudinal measurements provide an additional means of identifying data errors

height.s <- split(data.frame(height$servicedate,height$height),height$patientid)
weight.s <- split(data.frame(as.Date(weight$servicedate),weight$weight),weight$patientid)
glucose.s <- split(data.frame(as.Date(glucose$servicedate),glucose$glucose),glucose$patientid)
hba1c.s <- split(data.frame(as.Date(hba1c$servicedate),hba1c$hba1c),hba1c$patientid)

# summarize rate of change and within-patient variability

# function to estimate rate of change and residual variability for each child's data
longrate <- function(x){
  days <- as.numeric(as.Date(x[,1]))
  measure <- x[,2]
  mod <- lm(measure ~ days)
  rate <- mod$coef[2]
  residsd <- summary(mod)$sigma
  return(c(rate,residsd))
}

height.lm <- t(sapply(height.s,longrate))

# take a look at a few patients with implausible trajectories

height.change.ind <- which(abs(height.lm[,1]) > 0.5)     
par(mfrow = c(2,3))
for (i in 1:6){
  plot(as.numeric(as.Date(height.s[[height.change.ind[i]]][,1])),height.s[[height.change.ind[i]]][,2], xlab = "date", ylab = "height", type = "l")
}

# a few of these measures look very suspicious, as if one measurement is about 2.5 times the other
# take a closer look at an example case

height[height$patientid == names(height.change.ind[4]),]

# generate BMI and look for implausible values

height$iddate <- paste(height$patientid,height$servicedate)
weight$iddate <- paste(weight$patientid,weight$servicedate)
bmi <- merge(height,weight,by = "iddate") # merge height and weight data 
bmi$bmi <- bmi$weight/(bmi$height/100)^2

par(mfrow = c(2,2))
hist(bmi$bmi)
plot(bmi$weight,bmi$bmi)
plot(bmi$height,bmi$bmi)

# unusual groupings in BMI plots suggest patients with wrong units for height or weight
# select a rule for eliminating these heights or weights

bmi$height <- ifelse(bmi$height < 100 & bmi$bmi > 100, bmi$height*2.54, bmi$height)
bmi$bmi <- bmi$weight/(bmi$height/100)^2

par(mfrow = c(2,2))
hist(bmi$bmi)
plot(bmi$weight,bmi$bmi)
plot(bmi$height,bmi$bmi)

##-------------------------------------------------------------
# 2. Phenotype Extraction

## Aggregate data to the patient level

# aggregate numeric measurements using the earliest, mean and maximum observed values

bmi$bmimean <- unsplit(sapply(split(bmi$bmi,bmi$patientid.x),mean,na.rm = T),bmi$patientid.x)
bmi$bmimax <- unsplit(sapply(split(bmi$bmi,bmi$patientid.x),max,na.rm = T),bmi$patientid.x)
bmi$firstbmi <- unsplit(sapply(split(bmi$bmi,bmi$patientid.x),function(x){x[1]}),bmi$patientid.x)

glucose$glucosemean <- unsplit(sapply(split(glucose$glucose,glucose$patientid),mean,na.rm = T),glucose$patientid)
glucose$glucosemax <- unsplit(sapply(split(glucose$glucose,glucose$patientid),max,na.rm = T),glucose$patientid)

hba1c$hba1cmean <- unsplit(sapply(split(hba1c$hba1c,hba1c$patientid),mean,na.rm = T),hba1c$patientid)
hba1c$hba1cmax <- unsplit(sapply(split(hba1c$hba1c,hba1c$patientid),max,na.rm = T),hba1c$patientid)

chol$cholmean <- unsplit(sapply(split(chol$chol,chol$patientid),mean,na.rm = T),chol$patientid)
chol$cholmax <- unsplit(sapply(split(chol$chol,chol$patientid),max,na.rm = T),chol$patientid)

encounter$agemean <- unsplit(sapply(split(encounter$age,encounter$patientid),mean,na.rm = T),encounter$patientid)
encounter$firstage <- unsplit(sapply(split(encounter$age,encounter$patientid),min,na.rm = T),encounter$patientid)

# look for any occurence of diabetes diagnosis codes, insulin, metformin,
# or visit to an endocrinologist within the period of interest
# T2DM ICD-9 = "250.00", T2DM ICD-10 = "E11.9", T1DM ICD-9 = "250.01", T1DM ICD-10 = "E10.9"
# Endocrinologist Medicare specialty code = 46

anycode <- function(x,code){
  code.present <- x %in% code
  return(sum(code.present)>0)
}

# Count number of occurences of code

sumcode <- function(x,code){
  code.present <- x %in% code
  return(sum(code.present))
}

# Determine whether metformin prescription precedes insulin prescription
# Returns 1 if only metformin prescribed or metformin prescribed before insulin
# otherwise returns 0

codeorder <- function(x){
  metdates <- as.Date(x$dates[x$drugs == "metformin"])
  insdates <- as.Date(x$dates[x$drugs == "insulin"])
  if (length(metdates) == 0) metfirst <- 0
  else if (length(metdates) > 0 & length(insdates) == 0) metfirst <- 1
  else if (length(metdates) == 0 & length(insdates) == 0) metfirst <- 0
  else metfirst <- suppressMessages(1*(min(metdates) < min(insdates)))
  return(metfirst)
}

# any T2DM code
encounter$T2DM <- unsplit(sapply(split(encounter$diag,encounter$patientid),anycode,code = c("250.00","E11.9")),encounter$patientid)

# number of T2DM codes
encounter$T2DMnum <- unsplit(sapply(split(encounter$diag,encounter$patientid),sumcode,code = c("250.00","E11.9")),encounter$patientid) # number of occurence of T2DM code

# any T1DM code
encounter$T1DM <- unsplit(sapply(split(encounter$diag,encounter$patientid),anycode,code = c("250.01","E10.9")),encounter$patientid)

# any visit to an endocrinologist
encounter$endo <- unsplit(sapply(split(encounter$prov,encounter$patientid),anycode,code = "46"),encounter$patientid)

# any depression diagnosis
encounter$dep <- unsplit(sapply(split(encounter$diag,encounter$patientid),anycode,code = c("296.2","296.9","296.3","300.4","F32.9","F41.8","F33.9")),encounter$patientid)

# any insulin prescription
meds$anyinsulin <- unsplit(sapply(split(meds$drug,meds$patientid),anycode,code = "insulin"),meds$patientid)

# any metformin prescription
meds$anymetformin <- unsplit(sapply(split(meds$drug,meds$patientid),anycode,code = "metformin"),meds$patientid)

# metformin prescription precedes insulin prescription
meds$metforminfirst <- unsplit(sapply(split(data.frame(dates=meds$presdate,drugs=meds$drug),
                                            meds$patientid),codeorder),meds$patientid)

## Create merged dataset with one observation per patient and aggregate variables

encounter1 <- encounter[!duplicated(encounter$patientid),c("patientid","agemean","firstage","race","gender","T2DM","T1DM","endo","T2DMnum","dep")]
bmi1 <- bmi[!duplicated(bmi$patientid.x),c("patientid.x","bmimean","bmimax","firstbmi")]
names(bmi1) <- c("patientid","bmimean","bmimax","firstbmi")
glucose1 <- glucose[!duplicated(glucose$patientid),c("patientid","glucosemean","glucosemax")]
hba1c1 <- hba1c[!duplicated(hba1c$patientid),c("patientid","hba1cmean","hba1cmax")]
chol1 <- chol[!duplicated(chol$patientid),c("patientid","cholmean","cholmax")]
meds1  <- meds[!duplicated(meds$patientid),c("patientid","anyinsulin","anymetformin","metforminfirst")]

data1 <- Reduce(function(x,y){merge(x,y, all = T)},list(encounter1,bmi1,glucose1,hba1c1,chol1,meds1,validation))

# create indicators for availability of any glucose or HbA1c measures

data1$anyglucose <- !is.na(data1$glucosemean)
data1$anyhba1c   <- !is.na(data1$hba1cmean)

# set insulin and metformin to false for patients with no medication data

data1$anyinsulin <- ifelse(is.na(data1$anyinsulin),FALSE,data1$anyinsulin)
data1$anymetformin <- ifelse(is.na(data1$anymetformin),FALSE,data1$anymetformin)

## Phenotyping models using gold standard labels from validation data set to construct prediction models for T2DM

## Logistic regression

mod.glm <- glm(T2DMv ~ T2DM + T1DM + bmimean + anyglucose + anyhba1c + anyinsulin + anymetformin, data = data1, family = "binomial")

# logistic regression-based phenotype
data1$T2DMglm <- predict(mod.glm, newdata = data1)

# evaluate performance of logistic regression phenotype
pred.glm <- na.omit(data.frame(pred = data1$T2DMglm,true = data1$T2DMv))
perf.glm <- roc(pred.glm$true, pred.glm$pred, auc = TRUE, print.auc = TRUE, show.thres = TRUE)

plot(perf.glm)

# Logistic regression AUC
perf.glm$auc

## CART

set.seed(20190805)
mod.cart <- rpart(T2DMv ~ T2DM + T1DM  + bmimean + anyglucose + anyhba1c, data = data1, method = "class")
mod.pruned<- prune(mod.cart, cp= mod.cart$cptable[which.min(mod.cart$cptable[,"xerror"]),"CP"])

par(xpd = NA) # prevent text labels from being cut off
plot(mod.pruned)
text(mod.pruned)

# predicted probabilities of T2DM based on CART
data1$T2DMcart <- predict(mod.pruned, newdata = data1, type = "prob")

# binary T2DM phenotype based on CART
data1$T2DMcart.class <- as.numeric(as.character(predict(mod.pruned, newdata = data1, type = "class")))

# evaluate performance of continuous CART phenotype
pred.cart <- na.omit(data.frame(pred = data1$T2DMcart[,2],true = data1$T2DMv))
perf.cart <- roc(pred.cart$true, pred.cart$pred, auc = TRUE, print.auc = TRUE, show.thres = TRUE)

par(xpd = FALSE)
plot(perf.cart)

# CART AUC
perf.cart$auc

## eMERGE T2DM rule

T2DM.rule <- function(x){
  if (x$T1DM ==1) T2DM <- 0
  else{
    if (x$T2DM ==1){
      if (x$anyinsulin == 1){
        if (x$anymetformin == 0){
          if (x$T2DMnum < 2){
            T2DM <- 0
          } else{
            T2DM <- 1
          }
        } else{
          if (x$metforminfirst == 0){
            T2DM <- 0
          } else{
            T2DM <- 1
          }
        }
      } else{
        if (x$anymetformin == 1){
          T2DM <- 1
        } else{
          if ((!is.na(x$glucosemax) & x$glucosemax > 200) | (!is.na(x$hba1cmax) & x$hba1cmax > 6.5)){
            T2DM <- 1
          } else{
            T2DM <- 0
          }
        }
      }
    } else{
      if (x$anymetformin == 0){
        T2DM <- 0
      } else{
        if ((!is.na(x$glucosemax) & x$glucosemax > 200) | (!is.na(x$hba1cmax) & x$hba1cmax > 6.5)){
          T2DM <- 1
        } else{
          T2DM <- 0
        }
      }
    }
  }
  return(T2DM)
}

data1$T2DMemerge <- unsplit(sapply(split(data1,data1$patientid),T2DM.rule),data1$patientid)

# eMERGE specificity
1-mean(data1$T2DMemerge[data1$T2DMv == 0 & !is.na(data1$T2DMv)])

# eMERGE sensitivity
mean(data1$T2DMemerge[data1$T2DMv == 1 & !is.na(data1$T2DMv)])

# eMERGE PPV
mean(data1$T2DMv[data1$T2DMemerge == 1],na.rm = T)

# eMERGE NPV
1 - mean(data1$T2DMv[data1$T2DMemerge == 0],na.rm = T)


##-------------------------------------------------------------
## 3. Missing Data

## For most real examples we would want to define our exposure (cholesterol) in a window around
## cohort entry. For this toy example we will just use all available data. 

# Percent missing cholesterol
mean(is.na(data1$cholmean))

# Number of encounters per patient
encounter$numvisit <- rep(c(table(encounter$patientid)), times = c(table(encounter$patientid)))
summary(c(table(encounter$patientid)))

# Merge number of encounters onto data set with one observation per patient
numvisit <- encounter[!duplicated(encounter$patientid),c("patientid","numvisit")]
data1 <- merge(data1,numvisit)

# Look for factors associated with missing cholesterol
data1$misschol <- is.na(data1$cholmean)
misschol.mod <- glm(misschol ~ firstage + race + gender + firstbmi, data = data1, family = binomial)
summary(misschol.mod)

# Generate probability of missingness from this model
data1$pmisschol1[!is.na(data1$firstbmi)] <- 1-predict(misschol.mod, type = "response", data = data1)

# Estimate probability of missingness using a two stage model
# first estimate probability of missingness conditional on making an endocrinologist visit
misschol.mod.2 <- glm(misschol ~ endo, data = data1, family = binomial)
summary(misschol.mod.2)

data1$pmisschol2 <- 1-predict(misschol.mod.2, type = "response")

# next estimate probability of missingness among those with and without endocrinologist visit
misschol.mod.20 <- glm(misschol ~ firstage + race + gender + firstbmi, data = data1[data1$endo == 0,], family = binomial)
summary(misschol.mod.20)

misschol.mod.21 <- glm(misschol ~ firstage + race + gender + firstbmi, data = data1[data1$endo == 1,], family = binomial)
summary(misschol.mod.21)

data1$pmisschol20[!is.na(data1$firstbmi)] <- 1-predict(misschol.mod.20, type = "response", newdata = data1[!is.na(data1$firstbmi),])
data1$pmisschol21[!is.na(data1$firstbmi)] <- 1-predict(misschol.mod.21, type = "response", newdata = data1[!is.na(data1$firstbmi),])

# create combined probability of having an observed cholesterol value given these two modules
data1$pmisschol.mod <- data1$pmisschol2*data1$pmisschol21+(1-data1$pmisschol2)*data1$pmisschol20

## Compare one module and two module probabilities of being observed
plot(data1$pmisschol1, data1$pmisschol.mod)

## Fit regression models using IPW to account for missingness in cholesterol

# Model using 1 step weights
data1$w1 <- 1/data1$pmisschol1
data1$w1 <- sum(!is.na(data1$w1))*data1$w1/sum(data1$w1,na.rm = T) # normalize weights to maintain sample size
chol.mod1 <- glm(T2DMcart.class~ firstage + factor(race) + gender + cholmean, data = data1, weights = w1, family = "binomial")
summary(chol.mod1)

# Model using 2 step weights
data1$w.mod <- 1/data1$pmisschol.mod
data1$w.mod <- sum(!is.na(data1$w.mod))*data1$w.mod/sum(data1$w.mod,na.rm = T) # normalize weights to maintain sample size
chol.mod2 <- glm(T2DMcart.class~ firstage + factor(race) + gender + cholmean, data = data1, weights = w.mod, family = "binomial")
summary(chol.mod2)

##-------------------------------------------------------------
## 4. Confounding by Utilization Intensity

## Analyze association between depression and T2DM with and without conditioning on visit intensity
dep.glm1 <- glm(T2DMcart.class ~ firstage + factor(race) + gender + dep, data = data1, family = "binomial")
summary(dep.glm1)

dep.glm2 <- glm(T2DMcart.class ~ firstage + factor(race) + gender + dep + numvisit, data = data1, family = "binomial")
summary(dep.glm2)

# compare odds ratios before and after adjustment
cbind(c(exp(dep.glm1$coef),NA),exp(dep.glm2$coef))

##-------------------------------------------------------------
## 5. Outcome Misclassification

## Analysis without additional adjustment variables

# first compute sensitivity and specificity using validation data
sens <- mean(as.numeric(as.character(data1$T2DMcart.class[data1$T2DMv == 1 & !is.na(data1$T2DMv)])))
sens

spec <- 1-mean(as.numeric(as.character(data1$T2DMcart.class[data1$T2DMv == 0 & !is.na(data1$T2DMv)])))
spec

# compute odds ratios based on 2x2 table
a <- sum(data1$T2DMcart.class == 1 & data1$dep == 1)
b <- sum(data1$T2DMcart.class == 0 & data1$dep == 1)
c <- sum(data1$T2DMcart.class == 1 & data1$dep == 0)
d <- sum(data1$T2DMcart.class == 0 & data1$dep == 0)

or.std <- a*d/(b*c) # standard odds ratio
or.mh <- (a/(a+b)-(1-spec))/(c/(c+d)-(1-spec))*(sens-c/(c+d))/(sens-a/(a+b)) # Magder and Hughes adjusted odds ratio

or.std
or.mh

## Adjusted analysis via logistic regression using EM algorithm

# posterior probability of Y
post.prob <- function(phat,S,sens,spec){
  post.probY <- ifelse(S== 1, sens*phat/(sens*phat+(1-spec)*(1-phat)),
                       (1-spec)*phat/((1-spec)*phat+sens*(1-phat)))
  return(post.probY)
}

# EM algorithm proposed by Magder and Hughes
mh.EM <- function(fmla, sens, spec, tol = 10^-4, maxit = 10){
  data1$Y <- data1$T2DMcart.class
  or1 <- glm(fmla, data = data1, family = "binomial")
  p0 <- predict(or1, type = "response")
  dif <- 1
  j <- 0
  while (dif > tol & j < maxit){
    w <- post.prob(p0,data1$T2DMcart.class,sens,spec)
    data2 <- rbind(data1, data1)
    data2$w <- c(w,1-w)
    data2$Y <- c(rep(1,nrow(data1)),rep(0,nrow(data1)))
    suppressWarnings(or2 <- glm(fmla, data = data2, family = "binomial", weights = w))
    p0 <- predict(or2, type = "response", newdata = data1)
    dif <- max(abs(or1$coef-or2$coef))
    or1 <- or2
    j <- j+1
  }
  if (dif > tol) return("Did not converge")
  else return(or2)
}

# fit model
fmla.dep <- formula("Y ~ firstage + factor(race) + gender + dep")
mod.MH <- mh.EM(fmla.dep, sens, spec, maxit = 100)
summary(mod.MH)

# naive model for comparison
mod.cart <- glm(T2DMcart.class ~ firstage + factor(race) + gender + dep, data = data1, family = "binomial")
summary(mod.cart)


# model based on validation data only
mod.valid <- glm(T2DMv ~ firstage + factor(race) + gender + dep, data = data1, family = "binomial")
summary(mod.valid)


# Compare odds ratios from all three models
data.frame(MH = exp(mod.MH$coef), Naive = exp(mod.cart$coef), Validation = exp(mod.valid$coef))

##-------------------------------------------------------------##
## 6. Using Probabilistic Phenotypes

# Function for bias correction with known values for mu0 and mu1
# link can take values "ident", "log", or "logistic"
bias.adjust.prob <- function(fmla,mu0,mu1,p0,link = "ident"){
  
  # regress probabilistic phenotype on predictors
  fitp = lm(fmla, data = data1)
  
  # make bias correctioon
  betastar = fitp$coef/(mu1 - mu0)
  
  if (link == "ident"){
    betastar = betastar
  } else if (link == "log"){
    betastar = betastar/p0
  } else if (link == "logit"){
    betastar <- betastar/(p0*(1-p0))
  } else return("unsupported link function")
  
  # return association parameters (drop intercept)
  return(betastar[-1])
}

# use validation data to compute mean phenotype probability among true cases and controls
data1$prob <- inv.logit(data1$T2DMglm)

mu0 <- mean(data1$prob[data1$T2DMv == 0 & !is.na(data1$T2DMv)],na.rm = T)
mu1 <- mean(data1$prob[data1$T2DMv == 1 & !is.na(data1$T2DMv)], na.rm = T)

# use mean of predicted probabilities to estimate prevalence
p0 <- mean(inv.logit(data1$T2DMglm),na.rm = T)

# fit model
fmla.prob <- formula("prob ~ firstage + factor(race) + gender + dep")

mod.prob <- bias.adjust.prob(fmla.prob, mu0, mu1, p0, link = "logit")

# compare with results using validation data
data.frame(Adj = exp(mod.prob), Validation = exp(mod.valid$coef)[-1])