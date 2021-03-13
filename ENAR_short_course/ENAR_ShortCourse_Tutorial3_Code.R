## R code for ENAR EHR short course tutorials

# set your working directory 
setwd("xxxx")

##-------------------------------------------------------------##
## Install and load R packages
##-------------------------------------------------------------##

install.packages("pda")
install.packages("survival")
install.packages('pscl')
library(pda)
library(survival)
library(pscl)

# # Or you can install pda via github:
# install.packages("devtools")
# library(devtools)
# devtools::install_github("penncil/pda")
# library(pda)


##-------------------------------------------------------------##
## ODAL
##-------------------------------------------------------------##
# ############################  read data ###############################
OUD_data = read.csv("https://raw.githubusercontent.com/rhubb/ENAR_Short_Course/master/data/case3_ODAL.csv")

# ############################  fit logistic reg using pooled data ###############################
fit.pool <- glm(diagnosis_OUD ~ age + sex + BMI + Depression + Alcohol_use_disorders, family = 'binomial', data = OUD_data)

# ############################  data for each site ###############################
site1_data = OUD_data[which(OUD_data$site == 'site1'),]
site2_data = OUD_data[which(OUD_data$site == 'site2'),]
site3_data = OUD_data[which(OUD_data$site == 'site3'),]


# ############################  STEP 0: Model specification using meta-data (in json file)  ###############################

## step 0a. specification of the model using json file (i.e., meta-data)
## the purpose of this is to ensure exactly the same model is implemented at all sites, by sharing the same control (json) file (i.e., meta-data)
control <- list(project_name = 'OUD study',
                step = 'initialize',
                sites = c('site1', 'site2', 'site3'),
                heterogeneity = FALSE,
                model = 'ODAL',
                family = 'binomial',
                outcome = "diagnosis_OUD",
                variables = c('age', 'sex', 'BMI', 'Depression','Alcohol_use_disorders'),
                optim_maxit = 100,
                lead_site = 'site1',
                upload_date = as.character(Sys.time()) )

## step 0b: the lead site informs all collaborating sites what model will be fitted, using the json file "control.json" 
## Specifically, at the lead site, we run the following command.
## enter "1" to allow transferring the control (json) file (i.e., meta-data)
## the file named "control.json" is automatically added to the cloud (or a local) folder, to be shared with all sites
pda(site_id = 'site1', control = control, dir = getwd())

## path of the folder:
getwd()



# ############################  STEP 1: initialize  ###############################

## at each collaborating site (i.e., site 2 and site 3), we 
## (1) find out what model to run, by retrieving the control file "control.json"
## (2) run this model, and share the regression coefficients and variance to the cloud (or a local) folder

## at site3: enter "1" to allow transferring of the regression coefficients and variance
## a json file is saved as "site3_initialize.json in the cloud folder
pda(site_id = 'site3', ipdata = site3_data, dir=getwd())

## ditto
## a json file is saved as "site2_initialize.json in the cloud folder
pda(site_id = 'site2', ipdata = site2_data, dir=getwd())


## at the lead site, update the control file ("control.json") by adding weighted-averaged regression coefficients (as the initial value for the ODAL)
## a json file is saved as "site1_initialize.json in the cloud folder
## control.json is automatically updated 
pda(site_id = 'site1', ipdata = site1_data, dir=getwd())


# ############################'  STEP 2: derivative  ###############################

## at each collaborating site (i.e., site 2 and site 3), we 
## (1) read the control file, which contains the weighted-averaged regression coefficients
## (2) plug-in the weighted-averaged regression coefficients into the log-likelihood and obtain the first and second derivative of the log-likelihood
## (3) share the intermediate results (i.e., first and second derivative of the log-likelihood) in the cloud folder

## at site3: enter "1" to allow transferring of the derivative and Hessian matrix
## a json file is saved as "site3_derive.json in the cloud folder
pda(site_id = 'site3', ipdata = site3_data, dir=getwd())

## ditto
## a json file is saved as "site2_derive.json in the cloud folder
pda(site_id = 'site2', ipdata = site2_data, dir=getwd())

## ditto
## a json file is saved as "site1_derive.json in the cloud folder
pda(site_id = 'site1', ipdata = site1_data, dir=getwd())


#' ############################'  STEP 3: estimate  ###############################

## at the lead site, we 
## (1) read the xxx_derive.json files and construct the surrogate likelihood
## (2) obtain the final estimates by optimizing the surrogate likelihood
pda(site_id = 'site1', ipdata = site1_data, dir=getwd())

##' the PDA ODAL is now completed!
##' All the sites can still run their own surrogate estimates and broadcast them.
##' compare the surrogate estimate with the pooled estimate
config <- getCloudConfig(site_id = 'site1', dir=getwd())
fit.odal <- pdaGet(name = 'site1_estimate', config = config)
control <- pdaGet(name = 'control', config)

# site 1: regression coefficients estimates by ODAL
fit.odal 

## for this tutorial, we can compare the ODAL estimates with the pooled estimate
cbind(b.pool=fit.pool$coef,
      b.odal=fit.odal$btilde,
      sd.pool=summary(fit.pool)$coef[,2],
      sd.odal=sqrt(diag(solve(fit.odal$Htilde)/nrow(OUD_data))))




##-------------------------------------------------------------##
## ODAC
##-------------------------------------------------------------##
# ############################  read data ###############################
ADRD_data = read.csv("https://raw.githubusercontent.com/rhubb/ENAR_Short_Course/master/data/case3_ODAC.csv")


# ############################  fit cox model using pooled data ###############################
fit.pool <- coxph(Surv(time, status) ~ age + sex + BMI + Depression + sleep_disorders, data = ADRD_data)

# ############################  data for each site ###############################
site1_data = ADRD_data[which(ADRD_data$site == 'site1'),]
site2_data = ADRD_data[which(ADRD_data$site == 'site2'),]
site3_data = ADRD_data[which(ADRD_data$site == 'site3'),]

# ############################  STEP 0: Model specification using meta-data (in json file)  ###############################

## step 0a. specification of the model using json file (i.e., meta-data)
## the purpose of this is to ensure exactly the same model is implemented at all sites, by sharing the same control (json) file (i.e., meta-data)
control <- list(project_name = 'ADRD study',
                step = 'initialize',
                sites = c('site1', 'site2', 'site3'),
                heterogeneity = FALSE,
                model = 'ODAC',
                family = 'cox',
                outcome = "Surv(time, status)",
                variables = c('age', 'sex', 'BMI', 'Depression', 'sleep_disorders'),
                optim_maxit = 100,
                lead_site = 'site1',
                upload_date = as.character(Sys.time()) )

## step 0b: the lead site informs all collaborating sites what model will be fitted, using the json file "control.json" 
## Specifically, at the lead site, we run the following command.
## enter "1" to allow transferring the control (json) file (i.e., meta-data)
## the file named "control.json" is automatically added to the cloud (or a local) folder, to be shared with all sites
pda(site_id = 'site1', control = control, dir = getwd())


# ############################  STEP 1: initialize  ###############################

## at each collaborating site (i.e., site 2 and site 3), we 
## (1) find out what model to run, by retrieving the control file "control.json"
## (2) run this model, and share the time points, regression coefficients, and variance to the cloud (or a local) folder

## at site3: enter "1" to allow transferring of the time points, regression coefficients, and variance
## a json file is saved as "site3_initialize.json in the cloud folder
pda(site_id = 'site3', ipdata = site3_data, dir=getwd())

## ditto
## a json file is saved as "site2_initialize.json in the cloud folder
pda(site_id = 'site2', ipdata = site2_data, dir=getwd())

## at the lead site, update the control file ("control.json") by adding weighted-averaged regression coefficients (as the initial value for the ODAL)
## a json file is saved as "site1_initialize.json in the cloud folder
## control.json is automatically updated 
pda(site_id = 'site1', ipdata = site1_data, dir=getwd())


# ###################   STEP 2: summary stats at each time point  #################
## at each site, we 
## (1) read the control file, which contains the weighted-averaged regression coefficients
## (2) obtain the preliminary summary-level statistics at each time point
## (3) share these preliminary summary-level statistics in the cloud folder

## at site3: enter "1" to allow transferring of the summary statistics at each time point
## a json file is saved as "site3_derive.json in the cloud folder
pda(site_id = 'site3', ipdata = site3_data, dir=getwd())

## at site2: enter "1" to allow transferring of the summary statistics at each time point
## a json file is saved as "site2_derive.json in the cloud folder
pda(site_id = 'site2', ipdata = site2_data, dir=getwd())

## at site1: enter "1" to allow transferring of the summary statistics at each time point
## a json file is saved as "site1_derive.json in the cloud folder
pda(site_id = 'site1', ipdata = site1_data, dir=getwd())



# ############################  STEP 3: derivative  ###############################

## at each site, we 
## (1) read the summary statistics from the xxx_derive.json files, and calculate the first and second derivatives of the log-likelihood function
## (2) share the derivatives in the cloud folder

## at site3: enter "1" to allow transferring of the derivatives
## a json file is saved as "site3_derive_UWZ.json in the cloud folder
pda(site_id = 'site3', ipdata = site3_data, dir=getwd())

## at site2: enter "1" to allow transferring of the derivatives
## a json file is saved as "site2_derive_UWZ.json in the cloud folder
pda(site_id = 'site2', ipdata = site2_data, dir=getwd())

## at site1: enter "1" to allow transferring of the derivatives
## a json file is saved as "site1_derive_UWZ.json in the cloud folder
pda(site_id = 'site1', ipdata = site1_data, dir=getwd())



# ############################  STEP 4: estimate  ###############################
## at the lead site, we 
## (1) read the derivatives from the xxx_derive_UWZ.json files and construct the surrogate likelihood
## (2) obtain the final estimates by optimizing the surrogate likelihood
pda(site_id = 'site1', ipdata = site1_data, dir=getwd())

## the PDA ODAC is now completed!
## All the sites can still run their own surrogate estimates and broadcast them.
config <- getCloudConfig(site_id = 'site1', dir=getwd())
fit.pda <- pdaGet(name = 'site1_estimate', config = config)
control <- pdaGet('control', config)

# site 1: regression coefficients estimates by ODAC
fit.pda 

## compare the surrogate estimate with the pooled estimate
cbind(b.pool=fit.pool$coef,
      b.pda=fit.pda$btilde )


##-------------------------------------------------------------##
## ODAH
##-------------------------------------------------------------##
# ############################  read data ###############################
SAE_data = read.csv("https://raw.githubusercontent.com/rhubb/ENAR_Short_Course/master/data/case3_ODAH.csv")

# ############################  fit hurdle regression using pooled data ###############################
fit.pool <- hurdle(SAE_freq ~ age+sex+CCI, data = SAE_data)

# ############################  data for each site ###############################
site1_data = SAE_data[which(SAE_data$site == 'site1'),]
site2_data = SAE_data[which(SAE_data$site == 'site2'),]
site3_data = SAE_data[which(SAE_data$site == 'site3'),]


# ############################  STEP 0: Model specification using meta-data (in json file)  ###############################

## step 0a. specification of the model using json file (i.e., meta-data)
## the purpose of this is to ensure exactly the same model is implemented at all sites, by sharing the same control (json) file (i.e., meta-data)
control <- list(project_name = 'Colorectal cancer study',
                step = 'initialize',
                sites = c('site1', 'site2', 'site3'),
                heterogeneity = FALSE,
                model = 'ODAH',
                family = 'hurdle',
                outcome = "SAE_freq",
                variables_hurdle_count = c("age", "sex", "CCI"),
                variables_hurdle_zero = c("age", "sex", "CCI"),
                offset = NULL,
                optim_maxit = 100,
                lead_site = 'site1',
                upload_date = as.character(Sys.time()) )

## step 0b: the lead site informs all collaborating sites what model will be fitted, using the json file "control.json" 
## Specifically, at the lead site, we run the following command.
## enter "1" to allow transferring the control (json) file (i.e., meta-data)
## the file named "control.json" is automatically added to the cloud (or a local) folder, to be shared with all sites
pda(site_id = 'site1', control = control, dir = getwd())



# ############################  STEP 1: initialize  ###############################

## at each collaborating site (i.e., site 2 and site 3), we 
## (1) find out what model to run, by retrieving the control file "control.json"
## (2) run this model, and share the regression coefficients and variance to the cloud (or a local) folder

## at site3: enter "1" to allow transferring of the regression coefficients and variance
## a json file is saved as "site3_initialize.json in the cloud folder
pda(site_id = 'site3', ipdata = site3_data, dir=getwd())

## ditto
## a json file is saved as "site2_initialize.json in the cloud folder
pda(site_id = 'site2', ipdata = site2_data, dir=getwd())


## at the lead site, update the control file ("control.json") by adding weighted-averaged regression coefficients (as the initial value for the ODAH)
## a json file is saved as "site1_initialize.json in the cloud folder
## control.json is automatically updated 
pda(site_id = 'site1', ipdata = site1_data, dir=getwd())


# ############################'  STEP 2: derivative  ###############################

## at each collaborating site (i.e., site 2 and site 3), we 
## (1) read the control file, which contains the weighted-averaged regression coefficients
## (2) plug-in the weighted-averaged regression coefficients into the log-likelihood and obtain the first and second derivatives of the log-likelihood
## (3) share the intermediate results (i.e., first and second derivatives of the log-likelihood) in the cloud folder

## at site3: enter "1" to allow transferring of the derivatives
## a json file is saved as "site3_derive.json in the cloud folder
pda(site_id = 'site3', ipdata = site3_data, dir=getwd())

## ditto
## a json file is saved as "site2_derive.json in the cloud folder
pda(site_id = 'site2', ipdata = site2_data, dir=getwd())

## ditto
## a json file is saved as "site1_derive.json in the cloud folder
pda(site_id = 'site1', ipdata = site1_data, dir=getwd())


#' ############################'  STEP 3: estimate  ###############################

## at the lead site, we 
## (1) read the xxx_derive.json files and construct the surrogate likelihood
## (2) obtain the final estimates by optimizing the surrogate likelihood
pda(site_id = 'site1', ipdata = site1_data, dir=getwd())

##' the PDA ODAH is now completed!
##' All the sites can still run their own surrogate estimates and broadcast them.
##' compare the surrogate estimate with the pooled estimate
config <- getCloudConfig(site_id = 'site1', dir=getwd())
fit.odah <- pdaGet(name = 'site1_estimate', config = config)
control <- pdaGet(name = 'control', config = config)

# site 1: regression coefficients estimates by ODAH
fit.odah 

## for this tutorial, we can compare the ODAH estimates with the pooled estimate
cbind(b.count.pool=fit.pool$coef$count,
      b.count.odah=fit.odah$btilde_count,
      b.zero.pool=fit.pool$coef$zero,
      b.zero.odah=fit.odah$btilde_zero)
