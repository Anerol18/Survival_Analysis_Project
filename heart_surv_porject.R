#download the the survival library
library(survival)
#inspect the "heart" dataset
data(heart)
#heart is a dataset with 172 obs and 8 variables
head(heart)
str(heart)
?heart

#prepare the dataset for glm model
#prefer to work on the copy of the dataset not to modify the original one
data1 = (heart)

# Convert the categorical variables to factors, since this works better with glm ((event, surgery, transplant)
data1$event = as.factor(heart$event)
data1$surgery = as.factor(heart$surgery)
data1$transplant = as.factor(heart$transplant)

str(data1)

#verify the time unit since is not specified in the documentation, so we need to give our hypothesis giving the context
# the max(data1$year) is 6.47 that in days is ~2363. this is in line with survival studies that are generally long-term ones.
max(data1$year)
#the max(data1$stop) is 1800 day if we divide by years (365.25) gives ~5 years. that is still coherent with days time unit.
max(data1$stop)

#doing a double check based on the supposition that the time units is in  days
#calculation of follow-up time in different time units and check if seems coherent 
follow_up_days = data1$stop - data1$start
follow_up_months = follow_up_days / 30.44 #approx avg per month
follow_up_years = follow_up_days / 365.25 #approx avg per year

#summary on different time units
#the results is confirming that our hypotheis of days time unit is coherent
summary(follow_up_days)
summary(follow_up_months)
summary(follow_up_years)

#now that the time unit is validated we can start modelling the data

#glm model with binomial family
#set seed
set.seed(123)
heart_glm = glm(event ~ start + stop + age + year + surgery + transplant, data = data1, family = binomial)
#the summary of the glm model using all covariates is saying that a part from the intercept only stop, year and transplant1 are statistically significant
#the AIC is 185.3
summary(heart_glm)

#let's double-check what VSURF says about the variables importance
library('VSURF')
heart_vs=VSURF(event~.,data=data1)

#check variable importance result
heart_vs[[3]]
#plot it to visualize variable names. results: Transplant, stop, year
plot(heart_vs,step="pred",var.names=TRUE)

#rerun glm model with the selected variables confirmed by VSURF
heart_glm2 = glm(event ~ stop + year + transplant, data = data1, family = binomial)

#summary of second glm model
#the AIC is 180.25 so the model is better using these 3 variables
summary(heart_glm2)


#use Kaplan-Meier estimator to have survival curve
heart_KM = survfit(Surv(start, stop, event) ~ 1, data = data1)
summary(heart_KM)


# Plot Kaplan-Meier estimator
plot(heart_KM, mark.time = TRUE,
     main = "Kaplan-Meier Estimator",
     ylab = "Survival Probability",
     xlab = "Time (days)")

# used normal plot because ggplot not working :(  can you manage? 

library(tidyverse)

#other option: use Nelson-AAlen estimator to retrieve the cumulative hazard function
#ops, not working, giving error
#creating a new dataset because data1 has some factorizes variables and is giving error
data2 = heart
heart_NA = survfit(Surv(start, stop, event) ~ 1, data = data2, type = "fh")
summary (heart_NA)

# Plot Nelson-AAlen estimator
plot(heart_NA, mark.time = TRUE,
     main = "Nelson-AAlen Estimator",
     ylab = "Cumulative Hazard",
     xlab = "Time (days)")


#maybe good to add the logrank test having as two groups patients that had transplant yes or no
#to develop



#fit using COX hazard model
heart_cox = coxph(Surv(start, stop, event) ~ age + year + surgery + transplant, data = data2)
#based on this summary we can say that age, year (and surgery ? maybe) seems to be significant predictors of the event
#increasing of age increase the hazard and decreasing of year decrease the hazard ratio of the event
#concordance value is 0.363 with a standard error of 0.033 suggest we need to find a better model
summary(heart_cox)

#check Matringale residuals for cox model
data2$residuals = residuals(heart_cox, type = 'martingale')
summary(data2$residuals)

#plot cox model residuals
#####to review, plot is weird
with(data2, {
  plot(age, data2$residuals)
  lines(lowess(age, data2$residuals), lwd = 2)
  plot(year, data2$residuals)
  plot(surgery, data2$residuals)
})

#TO DO:
# comparing nested models (Likelihood)
# comparing non-nested models (AIC)
# assessing goodness of fit
# checking model assumptions
# Stratified Cox Regression
# Penalized Cox regression
# cross validation (lasso, ridge, elasticnet)
# predicting survival
# Interpretations
# embellishment using better plotting package and knit 



# not used code for the moment
# library(nnet)
# 
# library(gtsummary)
# tab1 = tbl_summary(heart)
# tab1
# 
# fit = lm(event ~ ., data=heart)
# summary(fit)
# 
# tbl_regression(fit)
# 
# #capture model summary as an object
# modelSummary = summary(fit)  
# 
# # model coefficients
# modelCoeffs = modelSummary$coefficients  

