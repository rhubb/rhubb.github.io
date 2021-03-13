## R code for ENAR EHR short course tutorials

##-------------------------------------------------------------##
## Install and load R packages
##-------------------------------------------------------------##

install.packages("rpart")
install.packages("pROC")
install.packages("MatchIt")
install.packages("survival")
install.packages("colorspace")

library(rpart)
library(pROC)
library(MatchIt)
library(survival)
library(colorspace)

##-------------------------------------------------------------##
## Tutorial 1 - Addressing phenotyping error in the study of diabetes
##-------------------------------------------------------------##

# Read in data
case1 <- read.csv("https://raw.githubusercontent.com/rhubb/ENAR_Short_Course/master/data/case1.csv", head=T)

# A. Phenotype Extraction
  
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

case1$T2DMemerge <- unsplit(sapply(split(case1,case1$patientid),T2DM.rule),case1$patientid)

# eMERGE specificity
1-mean(case1$T2DMemerge[case1$T2DMv == 0 & !is.na(case1$T2DMv)])

# eMERGE sensitivity
mean(case1$T2DMemerge[case1$T2DMv == 1 & !is.na(case1$T2DMv)])

# eMERGE PPV
mean(case1$T2DMv[case1$T2DMemerge == 1],na.rm = T)

# eMERGE NPV
1 - mean(case1$T2DMv[case1$T2DMemerge == 0],na.rm = T)

## Phenotyping models using gold standard labels from validated subset of observations to construct prediction models for T2DM

## Logistic regression

mod.glm <- glm(T2DMv ~ T2DM + T1DM + bmimean + anyglucose + anyhba1c + anyinsulin + anymetformin, 
               data = case1, family = "binomial")

# logistic regression-based phenotype
case1$T2DMglm <- predict(mod.glm, newdata = case1)

# evaluate performance of logistic regression phenotype
pred.glm <- na.omit(data.frame(pred = case1$T2DMglm,true = case1$T2DMv))
perf.glm <- roc(pred.glm$true, pred.glm$pred, auc = TRUE, print.auc = TRUE, show.thres = TRUE)

plot(perf.glm)

# Logistic regression AUC
perf.glm$auc

## CART

set.seed(20210303)
mod.cart <- rpart(T2DMv ~ T2DM + T1DM  + bmimean + anyglucose + anyhba1c + anyinsulin + anymetformin, data = case1, method = "class")
mod.pruned<- prune(mod.cart, cp= mod.cart$cptable[which.min(mod.cart$cptable[,"xerror"]),"CP"])

par(xpd = NA) # prevent text labels from being cut off
plot(mod.pruned)
text(mod.pruned)

# predicted probabilities of T2DM based on CART
case1$T2DMcart <- predict(mod.pruned, newdata = case1, type = "prob")

# binary T2DM phenotype based on CART
case1$T2DMcart.class <- as.numeric(as.character(predict(mod.pruned, newdata = case1, type = "class")))

# evaluate performance of continuous CART phenotype
pred.cart <- na.omit(data.frame(pred = case1$T2DMcart[,2],true = case1$T2DMv))
perf.cart <- roc(pred.cart$true, pred.cart$pred, auc = TRUE, print.auc = TRUE, show.thres = TRUE)

par(xpd = FALSE)
plot(perf.cart)

# CART AUC
perf.cart$auc

# sensitivity and specificity of the binary CART phenotype 
sens <- mean(as.numeric(as.character(case1$T2DMcart.class[case1$T2DMv == 1 & !is.na(case1$T2DMv)])))
sens

spec <- 1-mean(as.numeric(as.character(case1$T2DMcart.class[case1$T2DMv == 0 & !is.na(case1$T2DMv)])))
spec

# B. Outcome Misclassification

## Analysis without additional adjustment variables

# compute odds ratios based on 2x2 table
a <- sum(case1$T2DMcart.class == 1 & case1$rs7903146 == 1)
b <- sum(case1$T2DMcart.class == 0 & case1$rs7903146 == 1)
c <- sum(case1$T2DMcart.class == 1 & case1$rs7903146 == 0)
d <- sum(case1$T2DMcart.class == 0 & case1$rs7903146 == 0)

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
  case1$Y <- case1$T2DMcart.class
  or1 <- glm(fmla, data = case1, family = "binomial")
  p0 <- predict(or1, type = "response")
  dif <- 1
  j <- 0
  while (dif > tol & j < maxit){
    w <- post.prob(p0,case1$T2DMcart.class,sens,spec)
    data2 <- rbind(case1, case1)
    data2$w <- c(w,1-w)
    data2$Y <- c(rep(1,nrow(case1)),rep(0,nrow(case1)))
    suppressWarnings(or2 <- glm(fmla, data = data2, family = "binomial", weights = w))
    p0 <- predict(or2, type = "response", newdata = case1)
    dif <- max(abs(or1$coef-or2$coef))
    or1 <- or2
    j <- j+1
  }
  if (dif > tol) return("Did not converge")
  else return(or2)
}

# fit model
fmla.gen <- formula("Y ~ firstage + factor(race) + gender + rs7903146")
mod.MH <- mh.EM(fmla.gen, sens, spec, maxit = 100)
summary(mod.MH)

# naive model for comparison
mod.cart <- glm(T2DMcart.class ~ firstage + factor(race) + gender + rs7903146, data = case1, family = "binomial")
summary(mod.cart)

# model based on validation data only
mod.valid <- glm(T2DMv ~ firstage + factor(race) + gender + rs7903146, data = case1, family = "binomial")
summary(mod.valid)

# Compare odds ratios from all three models
data.frame(MH = exp(mod.MH$coef), Naive = exp(mod.cart$coef), Validation = exp(mod.valid$coef))

##-------------------------------------------------------------##
## Tutorial 2 - Clinical trial with EHR-derived external control arm
##-------------------------------------------------------------##

# Read in data
case2 <- read.csv("https://raw.githubusercontent.com/rhubb/ENAR_Short_Course/master/data/case2.csv", head=T)

#  A. Matching

# Mahalanobis matching
mah.match <- matchit(trial ~ age + race + chemo + et + cap + pr, data = case2, method = "nearest", distance = "mahalanobis")
mah.data <- match.data(mah.match)

# Propensity score matching
ps.match <- matchit(trial ~ age + race + chemo + et + cap + pr, data = case2, method = "nearest", distance = "glm", caliper = 0.1)
ps.data <- match.data(ps.match)

# check balance in patient characteristics before and after matching
summary(mah.match, addlvariables = ~ et*race) # generate balance statistics before and after Mahalanobis matching
plot(summary(mah.match), var.order = "unmatched")

summary(ps.match, addlvariables = ~ et*race, un = FALSE) # generate balance statistics after propensity score matching
plot(summary(ps.match), var.order = "unmatched")
plot(ps.match, type = "histogram") # historgam of propensity scores before and after matching

# B. Analysis of matched data  

## unadjusted comparison of trial and EHR patients prior to matching
COL = rainbow_hcl(2)
plot(survfit(Surv(time,Event)~ trial, data = case2), xlim = c(0,28), col = COL, lwd = 2, mark.time = T)
legend("topright",col = COL, legend = c("Trial","RWD"), lty = 1, lwd = 2)

# Analysis of data prior to matching
mod1 <- coxph(Surv(time,Event)~ trial + age + race + chemo + et + cap + pr, data = case2)
summary(mod1)

# Mahalanobis matching
mod2 <- coxph(Surv(time,Event)~ trial + age + race + chemo + et + cap + pr, data = mah.data)
summary(mod2)

# Propensity score matching
mod3 <- coxph(Surv(time,Event)~ trial + age + race + chemo + et + cap + pr, data = ps.data)
summary(mod3)

##-------------------------------------------------------------##
## Tutorial 3 - Distributed analysis
##-------------------------------------------------------------##

