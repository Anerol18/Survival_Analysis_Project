# Survival Analysis

# Participants:
# Amit Agarwal
# Lorena Romeo
# Nikolai Len
# Quentin Camilleri

# Data set:  Survival of patients on the waiting list for the Stanford heart  
# transplant program present in the survival package.

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
# id and subject are the same.

# actual_age = accept.dt - birth.dt
# futime = fu.date - accept.dt
# wait.time = tx.date - accept.dt

# stop = fu.date - accept.dt + 1 = end of followup - acceptance into program + 1
# wait.time = tx.date - accept.dt + 1 = transplant date - acceptance into program + 1
# (stop - start) time between the events.

#------------------------------------------------------------------------------
# Libraries seen in class or intended to be used.
library(dplyr)
library(tidyr)
library(ggplot2)
library(magrittr)
library(lubridate)
library(tidyverse)
library(broom)
library(xtable)
library(BSDA)  # devtools::install_github('alanarnholt/BSDA')
library(webr)  # devtools::install_github("cardiomoon/webr")
library(gtsummary)
library(modelsummary)
library(epiDisplay)
library(mgcv)
library(survival)
library(ggfortify)
library(gridExtra)
library(survminer) 
library(epiR)
library(swimplot)
library(muhaz)
library(asaur)
library(maxLik)
library(survivalROC)
library(plyr)
library(glmnet)
library(tidycmprsk)
library(mstate)
library(cmprsk)
library(timeROC)
library(survAUC)
library(tidycmprsk)
library(openxlsx)

library(VSURF)
library(Hmisc)
library(pec)
library(riskRegression)
library(car)


#------------------------------------------------------------------------------
# Load the datasets
heart_data <- heart
stanford2_data = stanford2
jasa_data = jasa
jasa1_data = jasa1
#------------------------------------------------------------------------------
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
#------------------------------------------------------------------------------
# Check for missing values
missing_values <- heart_data %>% summarise_all(~sum(is.na(.)))
print(missing_values) # No missing values

# Add jasa columns to heart_data
heart_data <- heart_data %>%
  mutate(
    actual_age = jasa_data$age[match(id, 1:nrow(jasa_data))],
    accept_dt = jasa_data$accept.dt[match(id, 1:nrow(jasa_data))],
    tx_date = jasa_data$tx.date[match(id, 1:nrow(jasa_data))],
    fu_date = jasa_data$fu.date[match(id, 1:nrow(jasa_data))],
    reject  = jasa_data$reject[match(id, 1:nrow(jasa_data))],
    fustat = jasa_data$fustat[match(id, 1:nrow(jasa_data))]
  )

# Duration in years
heart_data <- heart_data %>%
  mutate(duration = (stop - start)/365.25)

# Convert reject and fustat to factors
# heart_data <- heart_data %>%
#   mutate(
#     reject = as.factor(reject),
#     fustat = as.factor(fustat)
#   )

# Age groups
heart_data$age_group <- cut(heart_data$actual_age, breaks = seq(0, 80, 10),
                        labels = c("0-10", "10-20", "20-30", "30-40", "40-50",
                                "50-60", "60-70", "70-80"))

# Descriptive statistics
summary(heart_data)

str(heart_data)

#------------------------------------------------------------------------------
wb <- createWorkbook()
addWorksheet(wb, "heart")
writeData(wb, "heart", heart_data)
addWorksheet(wb, "stanford2")
writeData(wb, "stanford2", stanford2_data)
addWorksheet(wb, "jasa")
writeData(wb, "jasa", jasa_data)
addWorksheet(wb, "jasa1")
writeData(wb, "jasa1", jasa1_data)
saveWorkbook(wb, "surviving_heart_transplant.xlsx", overwrite = TRUE)

write.csv(heart_data, file = "heart.csv", row.names = FALSE)
write.csv(stanford2_data, file = "stanford2.csv", row.names = FALSE)
write.csv(jasa_data, file = "jasa.csv", row.names = FALSE)
write.csv(jasa1_data, file = "jasa1.csv", row.names = FALSE)
#------------------------------------------------------------------------------
# Create survival object
# surv_object <- Surv(time = heart_data$stop, event = heart_data$event)
surv_object <- Surv(time = heart_data$stop, event = heart_data$event)

# Check for censoring in the dataset
heart_data$censored <- ifelse(heart_data$event == 0, TRUE, FALSE)
print(table(heart_data$censored))
data_uncensored <- heart_data %>% filter(event == 1)

# 75 observations are not censored, i.e. death has occurred.
# 97 observations are censored, i.e. the event has not occurred by the end of 
# the study period.

# Accrual and Follow-up Plot
ggplot(heart_data, aes(x = start, xend = stop, y = id, yend = id)) +
  geom_segment(size = 1) +
  geom_point(aes(x = start), size = 3) +
  geom_point(aes(x = stop, shape = factor(event)), linewidth = 3) +
  scale_shape_manual(values = c(16, 21)) +  # Filled circle for censored, open circle for event
  labs(x = "Year of entry - calendar time", y = "", 
       title = "Accrual and Follow-up", shape = "Censored") +
  theme_minimal()

# Survival Time Plot
ggplot(heart_data, aes(x = 0, xend = duration, y = id, yend = id)) +
  geom_segment(size = 1) +
  geom_point(aes(x = 0), size = 3) +
  geom_point(aes(x = duration, shape = factor(event)), size = 3) +
  scale_shape_manual(values = c(16, 4)) +  # Filled circle for censored, X for event
  labs(x = "Survival time in years - patient time", y = "", 
       title = "Survival Time") +
  theme_minimal()

#------------------------------------------------------------------------------
# The summary of fitted models provide the following information:
# time
# n.risk: number of subjects still under observation and at risk
# n.event: number of events (that occurred at each specific point in time
# survival: probability of surviving beyond a given time point.
# std.err: uncertainty in the survival estimate at each point in time
# lower 95% CI: lower 95% confidence intervals for the survival estimates
# upper 95% CI: upper 95% confidence intervals for the survival estimates
#------------------------------------------------------------------------------
# Kaplan-Meier estimator
# fitdata <- 
km_fit <- survfit(surv_object ~ 1, data = heart_data)
summary(km_fit)

# Over time, the probability of survival decreases.
# By 1387 (~3.8 years), only 18.1% of the individuals are expected to survive 
# beyond this.

# Plot Kaplan-Meier survival curve with censoring
ggsurvplot(km_fit, conf.int = TRUE, 
           ggtheme = theme_minimal(), 
           title = "Kaplan-Meier Survival Curve with Censoring", 
           censor.shape = '|', 
           censor.size = 4)

# Estimating x-year survival (surviving beyond a certain number of years)
x_year_survival <- summary(km_fit, times = 365.25 * 5)  # 5-year survival
print(paste("x-year survival: ", x_year_survival))

# Estimating median survival time
median_survival_time <- summary(km_fit)$table["median"]
print(paste("Median Survival Time: ", median_survival_time))

# Hazard function with complementary log-log transformation
ggsurvplot(km_fit, fun = "cloglog",
           ggtheme = theme_minimal(),
           title = "Complementary Log-Log Transformation")

# Cumulative Hazard function
ggsurvplot(km_fit, fun = "cumhaz",
           ggtheme = theme_minimal(),
           title = "Cumulative Hazard Function")

#------------------------------------------------------------------------------
# Kaplan-Meier estimator stratified by transplant status
km_fit_transplant <- survfit(surv_object ~ transplant, data = heart_data)
summary(km_fit_transplant)

# Over time, the probability of survival decreases.
# By 1387 (~3.8 years), only 19% of the individuals who had a transplant are 
# expected to survive  beyond this.
# Less than 23.6% of subjects who didnt have a transplant will survive over a 
# year.

# Plot Kaplan-Meier survival curve with censoring by transplant status
ggsurvplot(km_fit_transplant, conf.int = TRUE,
           ggtheme = theme_minimal(),
           title = "Survival Curve by Transplant Status with Censoring", 
           censor.shape = '|', censor.size = 4)

# Estimating x-year survival (surviving beyond a certain number of years)
x_year_survival <- summary(km_fit_transplant, times = 365.25 * 5)  # 5-year survival
print(paste("x-year survival: ", x_year_survival))
      
# Estimating median survival time
median_survival_time <- summary(km_fit_transplant)$table["median"]
print(paste("Median Survival Time: ", median_survival_time))

# Hazard function
ggsurvplot(km_fit_transplant, fun = "cloglog",
           ggtheme = theme_minimal(),
           title = "Complementary Log-Log Transformation")

# Cumulative Hazard function
ggsurvplot(km_fit_transplant, fun = "cumhaz",
           ggtheme = theme_minimal(),
           title = "Cumulative Hazard Function")

#------------------------------------------------------------------------------
# Kaplan-Meier estimator stratified by surgery status
km_fit_surgery <- survfit(surv_object ~ surgery, data = heart_data)
summary(km_fit_surgery)

# Over time, the probability of survival decreases.
# By 1387 (~3.8 years), only 15.5% of the individuals are expected to survive 
# beyond this who didnt have surgery and 34% are expected to survive beyond 
# 2.72 years among those who had a surgery.

# Plot Kaplan-Meier survival curve with censoring by surgery status
ggsurvplot(km_fit_surgery, conf.int = TRUE,
           ggtheme = theme_minimal(), 
           title = "Survival Curve by Surgery Status with Censoring",
           censor.shape = '|', censor.size = 4)

# Estimating x-year survival (surviving beyond a certain number of years)
x_year_survival <- summary(km_fit_surgery, times = 365.25 * 5)  # 5-year survival
print(paste("x-year survival: ", x_year_survival))
      
# Estimating median survival time
median_survival_time <- summary(km_fit_surgery)$table["median"]
print(paste("Median Survival Time: ", median_survival_time))

# Hazard function
ggsurvplot(km_fit_surgery, fun = "cloglog",
           ggtheme = theme_minimal(),
           title = "Complementary Log-Log Transformation")

# Cumulative Hazard function
ggsurvplot(km_fit_surgery, fun = "cumhaz",
           ggtheme = theme_minimal(),
           title = "Cumulative Hazard Function")

#------------------------------------------------------------------------------
# Kaplan-Meier estimator stratified by age groups
km_fit_age <- survfit(surv_object ~ age_group, data = heart_data)

summary(km_fit_age)

# Over time, the probability of survival decreases.
# See results for breakdown by age groups

# Plot Kaplan-Meier survival curve with censoring by age status
ggsurvplot(km_fit_age, conf.int = TRUE, 
           ggtheme = theme_minimal(), 
           title = "Survival Curve by Age Groups with Censoring",
           censor.shape = '|', censor.size = 4)

# Estimating x-year survival (surviving beyond a certain number of years)
x_year_survival <- summary(km_fit_age, times = 365.25 * 5)  # 5-year survival
print(paste("x-year survival: ", x_year_survival))
      
# Estimating median survival time
median_survival_time <- summary(km_fit_age)$table["median"]
print(paste("Median Survival Time: ", median_survival_time))

# Hazard function
ggsurvplot(km_fit_age, fun = "cloglog",
           ggtheme = theme_minimal(),
           title = "Complementary Log-Log Transformation")

# Cumulative Hazard function
ggsurvplot(km_fit_age, fun = "cumhaz",
           ggtheme = theme_minimal(),
           title = "Cumulative Hazard Function")

#------------------------------------------------------------------------------
# Cumulative hazard function using Nelson-Aalen estimator
na_fit <- survfit(Surv(stop, event) ~ 1, data = heart_data,
                  type = "fleming-harrington")

summary(na_fit)

# Over time, the probability of survival decreases.
# By 1387 (~3.8 years), only 19.1% of the individuals are expected to survive 
# beyond this.

# Plot cumulative hazard function using Nelson-Aalen estimator
ggsurvplot(na_fit, conf.int = TRUE,
           ggtheme = theme_minimal(),
           title = "Nelson-Aalen Cumulative Hazard Curve", fun = "cumhaz")

# Estimating x-year survival (surviving beyond a certain number of years)
x_year_survival <- summary(na_fit, times = 365.25 * 5)  # 5-year survival
print(paste("x-year survival: ", x_year_survival))
      
# Estimating median survival time
median_survival_time <- summary(na_fit)$table["median"]
print(paste("Median Survival Time: ", median_survival_time))

# Hazard function
ggsurvplot(na_fit, fun = "cloglog",
           ggtheme = theme_minimal(),
           title = "Complementary Log-Log Transformation")

# Cumulative Hazard function
ggsurvplot(na_fit, fun = "cumhaz",
           ggtheme = theme_minimal(),
           title = "Cumulative Hazard Function")

#------------------------------------------------------------------------------
# Cumulative hazard function using Nelson-Aalen estimator by transplant status
na_fit_transplant <- survfit(Surv(stop, event) ~ transplant,
                             data = heart_data, type = "fleming-harrington")

summary(na_fit_transplant)

# Over time, the probability of survival decreases.
# 26.7% of the subjects who didnt have a transplant are expected to survive 
# beyond a year.
# 20.2% are expected to survive beyond 1387 (~3.8 years) among those who had 
# a transplant.

# Plot cumulative hazard function using Nelson-Aalen estimator by transplant status
ggsurvplot(na_fit_transplant, conf.int = TRUE,
           ggtheme = theme_minimal(),
           title = "Nelson-Aalen Cumulative Hazard Curve by Transplant Status",
           fun = "cumhaz")

# Estimating x-year survival (surviving beyond a certain number of years)
x_year_survival <- summary(na_fit_transplant, times = 365.25 * 5)  # 5-year survival
print(paste("x-year survival: ", x_year_survival))
            
# Estimating median survival time
median_survival_time <- summary(na_fit_transplant)$table["median"]
print(paste("Median Survival Time: ", median_survival_time))

# Hazard function
ggsurvplot(na_fit_transplant, fun = "cloglog",
           ggtheme = theme_minimal(),
           title = "Complementary Log-Log Transformation")

# Cumulative Hazard function
ggsurvplot(na_fit_transplant, fun = "cumhaz",
           ggtheme = theme_minimal(),
           title = "Cumulative Hazard Function")

#------------------------------------------------------------------------------
# Cumulative hazard function using Nelson-Aalen estimator by surgery status
na_fit_surgery <- survfit(Surv(stop, event) ~ surgery, 
                          data = heart_data, 
                          type = "fleming-harrington")

summary(na_fit_surgery)

# Over time, the probability of survival decreases.
# 16.7% are expected to survive beyond 1387 (~3.8 years) among those who did 
# not have a surgery.
# Less than 37.2% of the subjects who had a surgery might survive for over a 
# year.

# Plot cumulative hazard function using Nelson-Aalen estimator by surgery Status
ggsurvplot(na_fit_surgery, conf.int = TRUE,
           ggtheme = theme_minimal(),
           title = "Nelson-Aalen Cumulative Hazard Curve by Surgery Status", 
           fun = "cumhaz")
            
# Estimating x-year survival (surviving beyond a certain number of years)
x_year_survival <- summary(na_fit_surgery, times = 365.25 * 5)  # 5-year survival
print(paste("x-year survival: ", x_year_survival))
      
# Estimating median survival time
median_survival_time <- summary(na_fit_surgery)$table["median"]
print(paste("Median Survival Time: ", median_survival_time))

# Hazard function
ggsurvplot(na_fit_surgery, fun = "cloglog",
           ggtheme = theme_minimal(),
           title = "Complementary Log-Log Transformation")

# Cumulative Hazard function
ggsurvplot(na_fit_surgery, fun = "cumhaz",
           ggtheme = theme_minimal(), 
           title = "Cumulative Hazard Function")

#------------------------------------------------------------------------------
# Cumulative hazard function using Nelson-Aalen estimator by age status
na_fit_age <- survfit(Surv(stop, event) ~ age_group, 
                      data = heart_data, type = "fleming-harrington")

summary(na_fit_age)
# Over time, the probability of survival decreases.
# See results for breakdown by age groups

# Plot cumulative hazard function using Nelson-Aalen estimator by age status
ggsurvplot(na_fit_age, conf.int = TRUE,
           ggtheme = theme_minimal(),
           title = "Nelson-Aalen Cumulative Hazard Curve by Age Groups", 
           fun = "cumhaz")

# Estimating x-year survival (surviving beyond a certain number of years)
x_year_survival <- summary(na_fit_age, times = 365.25 * 5)  # 5-year survival
print(paste("x-year survival: ", x_year_survival))
      
# Estimating median survival time
median_survival_time <- summary(na_fit_age)$table["median"]
print(paste("Median Survival Time: ", median_survival_time))

# Hazard function
ggsurvplot(na_fit_age, fun = "cloglog",
           ggtheme = theme_minimal(),
           title = "Complementary Log-Log Transformation")
      
# Cumulative Hazard function
ggsurvplot(na_fit_age, fun = "cumhaz",
           ggtheme = theme_minimal(),
           title = "Cumulative Hazard Function")
            
#------------------------------------------------------------------------------
# Comparing survival times between groups
# Log-rank tests for different groups

# Stratified log-rank test
stratified_logrank <- survdiff(surv_object ~ transplant + strata(surgery), 
                               data = heart_data)
print(paste("Log-rank test for transplant + surgery : ", stratified_logrank))


surv_obj_jasa <- Surv(time = jasa_data$futime, event = jasa_data$fustat)
stratified_logrank <- survdiff(surv_obj_jasa ~ transplant + strata(surgery), 
                               data = jasa_data)
print(paste("Log-rank test for transplant + surgery : ", stratified_logrank))

# issue here since it takes 103 subjects as without transplant and 69 subjects
# with transplant. Very different results with jasa: same for below
      
# Log-rank test for Transplant vs. no transplant
logrank_test_transplant <- survdiff(surv_object ~ transplant, 
                                    data = heart_data)
print(paste("Log-rank test for transplant: ", logrank_test_transplant))
      
# Log-rank test for Surgery vs. no surgery
logrank_test_surgery <- survdiff(surv_object ~ surgery, 
                                 data = heart_data)
print(paste("Log-rank test for surgery: ",logrank_test_surgery))
            
# Log-rank test for Age groups
logrank_test_age <- survdiff(surv_object ~ age_group,
                             data = heart_data)
print(paste("Log-rank test for age groups: ",logrank_test_age))

#------------------------------------------------------------------------------
# Testing Proportional Hazards Assumption using Cox Proportional Hazards Model
cox_model <- coxph(surv_object ~ age + year + surgery + transplant, 
                   data = heart_data)
summary(cox_model)

# cox_model <- coxph(surv_object ~ actual_age + year + surgery + transplant, 
#                    data = heart_data)
# No difference in results if actual_age is used instead of age.

# ToDO: explain results

cox.zph(cox_model)
ggcoxzph(cox.zph(cox_model))
ggforest(cox_model, data = heart_data)

# Testing Influential Observations
ggcoxdiagnostics(cox_model, type = "dfbeta",
                 linear.predictions = FALSE, ggtheme = theme_bw())

# Testing Non-Linearity
martingale_residuals <- residuals(cox_model, type = "martingale")
ggplot(heart_data, aes(x = age, y = martingale_residuals)) +
  geom_point() +
  geom_smooth() +
  labs(title = "Martingale Residuals vs Age",
       x = "Age", y = "Martingale Residuals")

# Evaluate model fit with residuals

# Martingale residuals
martingale_residuals <- residuals(cox_model, type = "martingale")
plot(heart_data$age, martingale_residuals, 
     main = "Martingale Residuals vs Age", 
     xlab = "Age", 
     ylab = "Martingale Residuals")
abline(h = 0, col = "red")

# Cox-Snell residuals
cox_snell_residuals <- heart_data$event - martingale_residuals
plot(cox_snell_residuals, 
     main = "Cox-Snell Residuals", 
     xlab = "Index", 
     ylab = "Cox-Snell Residuals")
abline(h = 0, col = "red")

# Deviance residuals
deviance_residuals <- residuals(cox_model, type = "deviance")
plot(deviance_residuals, 
     main = "Deviance Residuals", 
     xlab = "Index", 
     ylab = "Deviance Residuals")
abline(h = 0, col = "red")
                                          
# Schoenfeld residuals
schoenfeld_residuals <- residuals(cox_model, type = "schoenfeld")
plot(schoenfeld_residuals, 
     main = "Schoenfeld Residuals", 
     xlab = "Index", 
     ylab = "Schoenfeld Residuals")
abline(h = 0, col = "red")
      
# Summary of Results
summary_table <- tbl_regression(cox_model)
summary_table

#------------------------------------------------------------------------------
# Fit survival models

# Log-normal model
lognormal_fit <- survreg(surv_object ~ age +  + year + surgery + transplant,
                         data = heart_data, dist='lognormal')
summary(lognormal_fit)

with(heart_data, plot(age, duration, xlab='Age', ylab='Days', 
                      xlim=c(0,65), ylim=c(.1, 10^5), log='y', type='n'))
with(heart_data, points(age, duration, pch=c(2,4)[event+1], cex=.7))

rr  <- residuals(lognormal_fit, type='matrix')
plot(rr)

pred <- predict(lognormal_fit, 
                newdata = data.frame(
                  age = 1:65, 
                  year = mean(heart_data$year, na.rm = TRUE),
                  surgery = 0,
                  transplant = factor("0", levels = levels(heart_data$transplant))
                ), 
                type = 'quantile', 
                p = c(.1, .5, .9))

# Plot age vs. duration (survival time) on a log scale
with(heart_data, plot(age, duration, xlab = 'Age', ylab = 'Days', 
                      xlim = c(0, 65), ylim = c(.1, 10^5), log = 'y', type = 'n'))

# Add points for actual data
with(heart_data, points(age, duration, pch = c(2, 4)[event + 1], cex = .7))

# Add predicted quantiles to the plot
matlines(1:65, pred, lty = c(2, 1, 2), col = 1)

#------------------------------------------------------------------------------
# Weibull model
weibull_fit <- survreg(surv_object ~ age + year + surgery + transplant,
                       data = heart_data, dist = "weibull")
summary(weibull_fit)

with(heart_data, plot(age, duration, xlab='Age', ylab='Days', 
                      xlim=c(0,65), ylim=c(.1, 10^5), log='y', type='n'))
with(heart_data, points(age, duration, pch=c(2,4)[event+1], cex=.7))

rr  <- residuals(weibull_fit, type='matrix')
plot(rr)

pred <- predict(weibull_fit, 
                newdata = data.frame(
                  age = 1:65, 
                  year = mean(heart_data$year, na.rm = TRUE),
                  surgery = 0,
                  transplant = factor("0", levels = levels(heart_data$transplant))
                ), 
                type = 'quantile', 
                p = c(.1, .5, .9))

# Plot age vs. duration (survival time) on a log scale
with(heart_data, plot(age, duration, xlab = 'Age', ylab = 'Days', 
                      xlim = c(0, 65), ylim = c(.1, 10^5), log = 'y', type = 'n'))

# Add points for actual data
with(heart_data, points(age, duration, pch = c(2, 4)[event + 1], cex = .7))

# Add predicted quantiles to the plot
matlines(1:65, pred, lty = c(2, 1, 2), col = 1)

pct <- 1:98/100   # The 100th percentile of predicted survival is at +infinity
ptime <- predict(weibull_fit,
              newdata = data.frame(age = 48,
                year = mean(heart_data$year,na.rm = TRUE),surgery = 0,
                transplant = factor("0",levels = levels(heart_data$transplant))
), type='quantile',
                 p=pct, se=TRUE)
matplot(cbind(ptime$fit, ptime$fit + 2*ptime$se.fit,
              ptime$fit - 2*ptime$se.fit)/30.5, 1-pct,
        xlab="Months", ylab="Survival", type='l', lty=c(1,2,2), col=1)


#------------------------------------------------------------------------------
# Exponential distribution
exp_fit <- survreg(surv_object ~ age + year + surgery + transplant,
                   data = heart_data, dist = "exponential")
summary(exp_fit)

with(heart_data, plot(age, duration, xlab='Age', ylab='Days', 
                      xlim=c(0,65), ylim=c(.1, 10^5), log='y', type='n'))
with(heart_data, points(age, duration, pch=c(2,4)[event+1], cex=.7))

rr  <- residuals(exp_fit, type='matrix')
plot(rr)

pred <- predict(exp_fit, 
                newdata = data.frame(
                  age = 1:65, 
                  year = mean(heart_data$year, na.rm = TRUE),
                  surgery = 0,
                  transplant = factor("0", levels = levels(heart_data$transplant))
                ), 
                type = 'quantile', 
                p = c(.1, .5, .9))

# Plot age vs. duration (survival time) on a log scale
with(heart_data, plot(age, duration, xlab = 'Age', ylab = 'Days', 
                      xlim = c(0, 65), ylim = c(.1, 10^5), log = 'y', type = 'n'))

# Add points for actual data
with(heart_data, points(age, duration, pch = c(2, 4)[event + 1], cex = .7))

# Add predicted quantiles to the plot
matlines(1:65, pred, lty = c(2, 1, 2), col = 1)

#------------------------------------------------------------------------------
# Log-logistic model
loglog_fit <- survreg(surv_object ~ age + year + surgery + transplant,
                      data = heart_data, dist = "loglogistic")
summary(loglog_fit)
with(heart_data, plot(age, duration, xlab='Age', ylab='Days', 
                      xlim=c(0,65), ylim=c(.1, 10^5), log='y', type='n'))
with(heart_data, points(age, duration, pch=c(2,4)[event+1], cex=.7))

rr  <- residuals(loglog_fit, type='matrix')
plot(rr)

pred <- predict(loglog_fit, 
                newdata = data.frame(
                  age = 1:65, 
                  year = mean(heart_data$year, na.rm = TRUE),
                  surgery = 0,
                  transplant = factor("0", levels = levels(heart_data$transplant))
                ), 
                type = 'quantile', 
                p = c(.1, .5, .9))

# Plot age vs. duration (survival time) on a log scale
with(heart_data, plot(age, duration, xlab = 'Age', ylab = 'Days', 
                      xlim = c(0, 65), ylim = c(.1, 10^5), log = 'y', type = 'n'))

# Add points for actual data
with(heart_data, points(age, duration, pch = c(2, 4)[event + 1], cex = .7))

# Add predicted quantiles to the plot
matlines(1:65, pred, lty = c(2, 1, 2), col = 1)

#------------------------------------------------------------------------------
# Binary Logistic Regression

# Predict probabilities of transplant for ages 1 to 65
heart_data$surgery <- as.factor(heart_data$surgery)

logistic_fit <- glm(transplant ~ age + year + surgery,
                    data = heart_data, family = binomial)

summary(logistic_fit)
exp(coef(logistic_fit))
with(heart_data, plot(age, duration, xlab='Age', ylab='Days', 
                      xlim=c(0,65), ylim=c(.1, 10^5), log='y', type='n'))
with(heart_data, points(age, duration, pch=c(2,4)[event+1], cex=.7))

rr  <- residuals(logistic_fit, type='deviance')
plot(rr, main = "Deviance Residuals", ylab = "Residuals", xlab = "Index")
abline(h = 0, col = "red", lwd = 2)

# Predict probabilities of transplant for ages 1 to 65
pred <- predict(logistic_fit, 
                newdata = data.frame(
                  age = 1:65, 
                  year = mean(heart_data$year, na.rm = TRUE),
                  surgery = factor("0", levels = levels(heart_data$surgery))
                ), 
                type = 'response')

# Plot age vs. duration (survival time) on a log scale
with(heart_data, plot(age, duration, xlab = 'Age', ylab = 'Days', 
                      xlim = c(0, 65), ylim = c(.1, 10^5), log = 'y', type = 'n'))

# Add points for actual data
with(heart_data, points(age, duration, pch = c(2, 4)[event + 1], cex = .7))

# Add predicted probabilities to the plot
lines(1:65, pred, lty = 1, col = 1)


# Check for Multicollinearity
vif(logistic_fit)

#------------------------------------------------------------------------------


                     
#------------------------------------------------------------------------------
# miscellaneous code

#------------------------------------------------------------------------------


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

