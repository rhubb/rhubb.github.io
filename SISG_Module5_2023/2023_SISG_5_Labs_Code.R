## SISG2023 R code for exercises

## Load data
cholesterol = read.csv("https://raw.githubusercontent.com/rhubb/SISG2023/master/data/SISG-Data-cholesterol.csv", head=T)
attach(cholesterol)

## Install and load R packages
install.packages("multcomp")
install.packages("lmtest")
install.packages("sandwich")
library(multcomp)
library(lmtest)
library(sandwich)

### Exercises

#1.
summary(TG)
summary(BMI)
group = 1*(BMI > 25)
group=factor(group,levels=c(0,1), labels=c("<=25",">25"))
table(group)
by(TG, group, mean)
by(TG, group, sd)
hist(TG)
hist(BMI)
boxplot(TG~group,ylab="Triglycerides (mg/dl)")
plot(TG ~ BMI, xlab = "BMI (kg/m2)", ylab = "Triglycerides (mg/dl)")

#2.
fit1 = lm(TG ~ BMI)
summary(fit1)
confint(fit1)
lines(BMI, fit1$fitted.values)

#3.
predict(fit1, newdata = data.frame(BMI = 23), interval = "confidence")
predict(fit1, newdata = data.frame(BMI = 23), interval = "prediction")

#4.
fit1 = lm(TG ~ BMI)
summary(fit1)$r.squared

#5
# Scatterplot of triglycerides vs BMI
plot(TG ~ BMI, xlab = "BMI (kg/m2)", ylab = "Triglycerides (mg/dl)")

# Identify observations with BMI <=37
bmi37 = which(BMI<=37)

# Consider again the regression of TG on BMI
fit1=lm(TG~BMI)
summary(fit1)

# excluding subjects with BMI > 37
fit2 = lm(TG[bmi37] ~ BMI[bmi37])
summary(fit2)

#6.
# Plot residuals vs fitted values
plot(fit1$fitted, fit1$residuals,xlab="Fitted values",ylab="Residuals")
abline(0,0)

# Quantile-quantile plot
qqnorm(fit1$residuals)
qqline(fit1$residuals)

# Deletion diagnostics
dfb=dfbeta(fit1)
index=order(abs(dfb[,2]),decreasing=T)
cbind(dfb[index[1:15],],BMI[index[1:15]],TG[index[1:15]])

# fit a linear regression model with robust standard errors
fit.robust = coeftest(fit1, vcov = sandwich)
fit.robust

#7
# Summarize the variable APOE
table_APOE=table(APOE)
table_APOE
prop.table(table_APOE)

# binary variable indicating presence of APOE4
APOE4 = ifelse(APOE %in% c(3,5,6), 1, 0)

## Linear regression analyses for association of APOE4 and BMI with TG
# multiple linear regression of triglycerides on BMI and APOE4
fit3=lm(TG~BMI+APOE4)
summary(fit3)

#8
# scatterplot with subjects stratified by APOE4
par(mfrow = c(1,1))
plot(BMI[APOE4 == 0], TG[APOE4 == 0], pch = 1, col=75,xlab = "BMI (kg/m2)", ylab = "Triglycerides (mg/dl)")
points(BMI[APOE4 == 1], TG[APOE4 == 1], pch = 1, col=34)

# multiple linear regression of triglycerides on BMI, APOE4, and interaction
fit4 = lm(TG ~ BMI*APOE4)
summary(fit4)

# Compare the models with and without interaction
anova(fit3,fit4)

# Compare with the model without APOE4
anova(fit1,fit4)

#9.
# descriptive statistics
summary(chol)
table(rs4775401)
hist(chol)

# graphical display: boxplot 
boxplot(chol ~ factor(rs4775401))

# numeric descriptives 
tapply(chol, factor(rs4775401), mean)
tapply(chol, factor(rs4775401), sd)

#10.
# ANOVA for cholesterol and rs4775401
fit1 = lm(chol ~ factor(rs4775401))
summary(fit1)
anova(fit1)

# One-way ANOVA (not assuming equal variances)
oneway.test(chol ~ factor(rs4775401))

# Using robust standard errors
coeftest(fit1, vcov = sandwich)

# Non-parametric ANOVA
kruskal.test(chol ~ factor(rs4775401))

#11.
# construct contrasts for all pairwise comparisons
M2 = contrMat(table(rs4775401), type="Tukey")
fit2 = lm(chol ~ -1 + factor(rs4775401))

# explore options to correct for multiple comparisons
mc2 = glht(fit2, linfct =M2)
summary(mc2, test=adjusted("none"))
summary(mc2, test=adjusted("bonferroni"))
summary(mc2, test=adjusted("hochberg"))
summary(mc2, test=adjusted("fdr"))

#12.
# exploratory data analysis
table(rs174548, APOE)
tapply(chol, list(factor(rs174548), factor(APOE)), mean)
tapply(chol, list(factor(rs174548), factor(APOE)), sd)

par(mfrow = c(1,1))
plot.design(chol ~ factor(rs174548) + factor(APOE))

# model with interaction
fit1 = lm(chol ~ factor(rs174548)*factor(APOE))
summary(fit1)

# model without interaction
fit2 = lm(chol ~ factor(rs174548) + factor(APOE))
summary(fit2)

# compare models with and without interaction
anova(fit2,fit1)

#13.
# Descriptive statistics for hypertension
table(HTN)
table(HTN,rs174548)
chisq.test(HTN,rs174548)
by(TG,HTN,mean)

# Logistic regression analysis for the association between rs174548 and hypertension
glm.mod1 <- glm(HTN ~ factor(rs174548), family = "binomial")
summary(glm.mod1)
exp(glm.mod1$coef)
exp(confint(glm.mod1))

#14.
# Logistic regression analysis for the association between triglycerides and hypertension
glm.mod2 <- glm(HTN ~ TG, family = "binomial")
summary(glm.mod2)
exp(glm.mod2$coef)
exp(confint(glm.mod2))

#15
# logistic regression analysis for the association between rs174548 and hypertension
# adjusting for triglycerides
glm.mod3 <- glm(HTN ~ TG+factor(rs174548), family = "binomial")
summary(glm.mod3)
exp(glm.mod3$coef)
exp(confint(glm.mod3))

lrtest(glm.mod2,glm.mod3)

#16.
# relative risk regression for the association between rs174548 and hypertension
# adjusting for triglycerides
glm.mod4 <- glm(HTN ~ TG+factor(rs174548), family = "poisson")
coeftest(glm.mod4, vcov = sandwich)
exp(glm.mod4$coef)

#17.
# risk difference regression for the association between rs174548 and hypertension
# adjusting for triglycerides
glm.mod5 <- lm(HTN ~ TG+factor(rs174548))
coeftest(glm.mod5, vcov = sandwich)
