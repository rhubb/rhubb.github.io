#####################################################################
## R commands (not intended to be a comprehensive list of commands!) 
#####################################################################

## To obtain help for a specific command --------------------------------------------- 
help(command_name)  

## Set up your working directory ----------------------------------------------------- 
# For PC
setwd("C:/your_working_directory/")  
# For Mac
setwd("/your_working_directory/")

## Installing packages and loading libraries ----------------------------------------- 
install.packages("package_name")
library("package_name")

## Reading an ASCII data set with the first row containing variable names and attributing
## to object called data -------------------------------------------------------------
data = read.table("your_dataset", header=T)  

## Reading a .csv dataset ------------------------------------------------------------
data = read.csv("your_dataset.csv", header=T)  

## Descriptive analysis: ------------------------------------------------------------- 
## summary statistics 
summary(data)  

## Summarize data by groups 
by(data, data$group, summary)  

## Histogram 
hist(data$varname)  

## Density plot 
plot(density(data$varname))  

## Box plot  
boxplot(data$varname) 

## Box plot by groups 
boxplot(data$response ~data$group)  

## Barplot (a little forced here...) 
barplot(table(data$varname))  

## Stem and leaf 
stem(data$varname)  

## Scatterplot 
plot(data$varname_x, data$varname_y)  

## Contingency table  
table(data$varname)  

## Contingency table 2 
table(data$varname_row, data$varname_column)

## Correlation
## Pearson
cor(data$varname1, data$varname2, method="pearson") 
## Spearman:
cor(data$varname1, data$varname2, method="spearman")  

## Inferential analysis:
## t-test for one sample
t.test(x=data$varname, y= NULL, mu=0) 

## t-test for two independent samples
t.test(data$response ~ data$groups, paired=FALSE)  
## alternative form
t.test(x=data$response[data$groups==0],y=data$response[data$groups==1],paired=FALSE) 

## Paired t-test 
t.test(x=data$var1, y=data$var2, paired=TRUE)  

## Proportion test 
prop.test(x=3, n=10, p=0.5) 

## binomial test (exact) 
binom.test(x=3, n=10, p=0.5)  

## chisquare test for one sample
chisq.test(x = c(4,3,4), p=c(1/3, 1/3, 1/3))  

## chisquare test 
chisq.test(data$var1, data$var2) 

## Fisher's exact test 
fisher.test(x=data$group, y=data$response)  

## MacNemar's test
mcnemar.test(data$group, data$response)  

## Wilcoxon rank sum test (ou Mann-Whitney)  
wilcox.test(data$response ~data$group)  

## Wilcoxon signed rank test 
wilcox.test(x=data$var1, y=data$var2, paired=TRUE)  

## Linear regression -----------------------------------------------------------------
fit = lm(response ~ predictor1 + predictor2, data=data) 
## now predictor2 is a categorical variable 
fit = lm(response ~ predictor1 + factor(predictor2), data=data)  

## model including interaction
fit = lm(response ~ predictor1 + predictor2 +predictor1*predictor2, data=data) 
## or alternatively 
fit = lm(response ~ predictor1*predictor2, data=data) 
#  (the above includes main effects!)

## summary of the model fit 
summary(fit)   

## estimated coefficients
coef(fit)

## fitted values 
fitted(fit)

## residuals
resid(fit)

## predictions
predict(fit)

## confidence intervals
confint(fit)  

## anova table 
anova(fit)  

## Scatterplot and regression line
plot(data$x, data$y, xlab="x", ylab="y") 
abline(fit, col=2, lwd=2)  

## Prediction intervalos for new observations with predictor values
## equal to 2, 3 e 4 
predict.lm(fit, newdata=data.frame(x=c(2,3,4)), interval="prediction")  

## residual analysis 
plot(fit, which=1:2)  

##  studentised rediduals
res = rstudent(fit)  

## Scatterplot of residuals
plot(res)    

## checking normality
qqnorm(res)  
qqline(res)  

## dfbeta 
dfb = dfbeta(fit) 
summary(dfb) 
plot(dfb[,2],ylab="dfbeta")   

## One-way analysis of variance ------------------------------------------------------
fit = oneway(data$y ~ data$x)

## summary of model fit
summary(fit)

##    ... similar output from lm command... 

## Logistic regression  -------------------------------------------------------------- 
fit = glm(response ~ predictor1 + predictor2, data=data, family=binomial)  

## summary of the model fit
summary(fit)   

## confidence interval for model parameters
confint(fit)  

## anova
anova(fit)  

## ppredictive probabilities 
pred = predict.glm(fit,type="response")  

## Pearson residuals 
resid.pearson = residuals(fit, type="pearson") 

## Deviance 
resid.dev = residuals(fit, type="deviance")  

## Risk difference regression  -------------------------------------------------------- 
library(gee)
fit = gee(response ~ predictor1 + predictor2, data=data, id = seq(1,nrow(data)))  

## Relative risk regression  ---------------------------------------------------------- 
fit = gee(response ~ predictor1 + predictor2, data=data, id = seq(1,nrow(data)),
	family = "poisson")  

## to quit an R session
q()