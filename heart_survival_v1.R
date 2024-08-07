# Survival Analysis

# Heart Survival

# Participants:
# Amit Agarwal
# Lorena Romeo
# Nikolai Len
# Quentin Camilleri

# Data set: heart from survival package.
# Survival of patients on the waiting list for the Stanford heart transplant 
# program.

# Objective:
# Assess the effect on survival of transplantation, treating the patient 
# population as homogeneous, while in the influence of a number of covariates 
# was investigated through pairwise correlations and explore techniques to see 
# the simultaneous effect of several covariates and see for what values of these 
# covariates, if any, transplantation is likely to prolong survival.


# Three dataset are available.

# 1. jasa, the original data with 103 observations with the following variables:
# birth.dt, accept.dt, tx.date, fu.date, fustat, surgery, age, futime, 
# wait.time,transplant, mismatch, hla.a2, mscore, reject

# jasa1 with 170 observations containing the following variables:
# id, start, stop, event, transplant, age, year, surgery

# 2. heart (the main data set) with 172 observations with the following 
# variables:
# start, stop, event, age, year, surgery, transplant, id
# J Crowley and M Hu (1977), Covariance analysis of heart transplant survival data. Journal of the
# American Statistical Association, 72, 27–36.

# 3. stanford2 (Stanford Heart Transplant data in a different format) with184
# observations with the following variables:
# id, time, status, age, t5
# LA Escobar and WQ Meeker Jr (1992), Assessing influence in regression analysis with censored
# data. Biometrics 48, 507–528. Page 519.

# jasa, jasa1, heart, all have 103 subjects while stanford2 has 184 subjects.


# Libraries seen in class or intended to be used.
library(dplyr)
library(tidyr)
library(ggplot2)
library(magrittr)
library(lubridate)
library(tidyverse)
library(broom)
library(xtable)
library(BSDA) # to install this library run.  devtools::install_github('alanarnholt/BSDA')
library(webr)  #.   devtools::install_github("cardiomoon/webr")
library(gtsummary)
library(modelsummary)
library(epiDisplay)
library(mgcv)
library(survival)
library(survivalROC)
library(timeROC)
library(survAUC)
library(ggfortify)
library(gridExtra)
library(survminer) 
library(epiR)
library(swimplot)
library(muhaz)
library(asaur)
library(maxLik)
library(plyr)
library(glmnet)
library(tidycmprsk)
library(mstate)
library(cmprsk)
library(tidycmprsk)
library(VSURF)
library(Hmisc)
library(pec)
library(riskRegression)

data1 = heart
data2 = stanford2
data3 = jasa
data4 = jasa1

write.csv(data1, file = "heart.csv", row.names = FALSE)
write.csv(data2, file = "stanford2.csv", row.names = FALSE)
write.csv(data3, file = "jasa.csv", row.names = FALSE)
write.csv(data4, file = "jasa1.csv", row.names = FALSE)

data(heart,package="survival")
?heart
str(heart)
head(heart)
summary(heart)

?stanford2
str(stanford2)
head(stanford2)
summary(stanford2)

?jasa
str(jasa)
head(jasa)
summary(jasa)

?jasa1
str(jasa1)
head(jasa1)
summary(jasa1)

# Predictions
fit <- survreg(Surv(time,status) ~ age + I(age^2), data=stanford2, 
               dist='lognormal') 
with(stanford2, plot(age, time, xlab='Age', ylab='Days', 
                     xlim=c(0,65), ylim=c(.1, 10^5), log='y', type='n'))
with(stanford2, points(age, time, pch=c(2,4)[status+1], cex=.7))
pred <- predict(fit, newdata=list(age=1:65), type='quantile', 
                p=c(.1, .5, .9)) 
matlines(1:65, pred, lty=c(2,1,2), col=1) 

# Predicted Weibull survival curve for heart transplant subject with
#  average age of 48.
lfit <- survreg(Surv(time, status) ~ age + I(age^2), data=stanford2)
pct <- 1:98/100   # The 100th percentile of predicted survival is at +infinity
ptime <- predict(lfit, newdata=data.frame(age=48), type='quantile',
                 p=pct, se=TRUE)
matplot(cbind(ptime$fit, ptime$fit + 2*ptime$se.fit,
              ptime$fit - 2*ptime$se.fit)/30.5, 1-pct,
        xlab="Months", ylab="Survival", type='l', lty=c(1,2,2), col=1)


# Fit an exponential model: the two fits are the same
survreg(Surv(time, status) ~ age + I(age^2), data=stanford2, dist='weibull',
        scale=1)
survreg(Surv(time, status) ~ age + I(age^2), data=stanford2,
        dist="exponential")

# A model with different baseline survival shapes for two groups, i.e.,
#   two different scale parameters
survreg(Surv(time, status) ~ age + I(age^2), data=stanford2)

fit <- survreg(Surv(time, status) ~ age + I(age^2), data=stanford2)
summary(fit)   # age and sex are both important

rr  <- residuals(fit, type='matrix')

### To reliagn code for our data set
# sum(rr[,1]) - with(mgus2, sum(log(futime[death==1]))) # loglik
# plot(mgus2$age, rr[,2], col= (1+mgus2$death)) # ldresp

# There are multiple ways to parameterize a Weibull distribution.
# The survreg function embeds it in a general location-scale family, which is a 
# different parameterization than the rweibull function, and often leads
# to confusion.
#   survreg's scale  =    1/(rweibull shape)
#   survreg's intercept = log(rweibull scale)
#   For the log-likelihood all parameterizations lead to the same value.
# y <- rweibull(1000, shape=2, scale=5)
# survreg(Surv(y)~1, dist="weibull")

# Economists fit a model called `tobit regression', which is a standard
# linear regression with Gaussian errors, and left censored data.
# tobinfit <- survreg(Surv(durable, durable>0, type='left') ~ age + quant,
#                     data=tobin, dist='gaussian')

# Using jasa data set to create transplant as a time dependent covariate
# The data set contains a few anomalies
# 1. A subject died on the day of entry. We cannot use (0,0) with coxph routine.
# So we assume the subject died at day 0.5. Alternatively, we can add a day to
# everyone's follow-up (The Statistical Analysis of Failure Time Data - 
# Kalbfleisch and Prentice). The result remains the same in both.
# 2. A subject transplanted on day 10 is considered to have been on medical 
# treatment for days 1-10 and as transplanted starting on day 11 but for 
# subject 38 who died on the same day as their procedure.
# This should be treated as a transplant death by advancing it by 0.5 days.
# 3. The treatment coefficients in table 6.1 (regression coefficients for model 
# fitted to the heart transplant data in The Statistical Analysis of Failure 
# Time Data - # Kalbfleisch and Prentice) can only be obtained if covariates 
# are defined in precisely the same way, since their models include 
# interactions (Table 5.2 of the book).
# For age this is (age in days)/ 365.25 - 48 years, and for year of enrollment 
# it is the number of years since the start of the study:
# (entry date - 1967/10/1)/365.25.


# Since time is in days the fractional time of 0.5 can be any value between 0 
# and 1.
jasa$subject <- 1:nrow(jasa) #we need an identifier variable
tdata <- with(jasa, data.frame(subject = subject,
futime= pmax(.5, fu.date - accept.dt),
txtime= ifelse(tx.date== fu.date,
               (tx.date -accept.dt) -.5,
               (tx.date - accept.dt)),
fustat = fustat
))

xdata <- tmerge(jasa, tdata, id=subject,
                death = event(futime, fustat),
                transplant = tdc(txtime),options= list(idname="subject"))

sdata <- tmerge(jasa, tdata, id=subject,death = event(futime, fustat),
                trt = tdc(txtime),options= list(idname="subject"))

attr(sdata, "tcount")

sdata$age <- sdata$age -48
sdata$year <- as.numeric(sdata$accept.dt - as.Date("1967-10-01"))/365.25

# model 6 of the table in The Statistical Analysis of Failure Time Data - 
# Kalbfleisch and Prentice
coxph(Surv(tstart, tstop, death) ~ age*trt + surgery + year,
      data= sdata, ties="breslow")

# This example is a special case for the tmerge function that is quite common:
# if the first created variable is an event then the time range for each 
# subject is inferred to be from 0 to that event time: explicit tstop and 
# tstart arguments are not required.
# It also makes use of a two argument form of event.
# Each of the event and cumevent functions may have a second argument, which if 
# present will be used as the value for the event code.
# If this second argument is not present a value of 1 is used.
# If an event variable is already a part of data1, then the updates make 
# changes to that variable, possibly adding more events.
# This feature allows for the infection indicator to be build up incrementally.
# The tdc and cumtdc arguments can have 1, 2 or three arguments.
# The first is always the time point, the second, if present, is the value to 
# be inserted, and an optional third argument is the initial value.
# If the tdc call has a single argument the result is always a 0/1 variable, 0 
# before the time point and 1 after.
# For the 2 or three argument form, the starting value before the first
# definition of the new variable (before the first time point) will be the 
# initial value.
# The default for the initial value is NA, the value of the tdcstart option.
# If the tdc variable being created is already a part of data1, the old value 
# is replaced wholesale, it is not updated.
# This differs from the behavior for events.
# Although there is a use case for updating an existing time-dependent value,
# say from a new data source, the much more common case was accidental reuse of
# an existing name, resulting in a mishmash between an existing baseline 
# covariate and the additions, and creating a column which was valid for neither.
# Equally, some examples showed that we can not reliably update a tdc, 
# i.e., cases where the correct answer is unclear.
# The tcount table for the above fit shows all the deaths at the trailing edge 
# of their interval, which is expected since the time of death or last follow-up
# was used to define each subject's interval of risk.
# Two of the transplants happened on day 0 and are listed as occurring on the
# leading edge of the first follow-up interval for the subject.
# The other 67 transplants were strictly within the (0, last follow up) interval
# of each subject.


# 1.1. Survival Analysis


# Calculating survival times


# Creating survival objects and curves



# Estimating x-year survival (surviving beyond a certain number of years)



# stimating median survival time


# Comparing survival times between groups


# 1.2. Hazard and Cumulative Hazard



# 1.3. Survival Function



# 1.4. Survival curves



# Surv()


# survfit()


# ggsurvplot()


# 1.5. Kaplan-Meier Curve



# 1.6. Censoring



# 1.7. Testing Proportional Hazards Assumption



# 1.8. Testing Influential Observations



# 1.9. Testing Non-Linearity



# 1.10. Cox Proportional Hazards Model


# log-rank tests


# 1.11. Parametric Survival Models


# 1.12. Binary Logistic Regression

