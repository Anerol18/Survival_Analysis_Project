#PART 1 : LOAD DATASET, INSPECTION, DATA CLEANSING
#download the the survival library
library(survival)
#inspect the "heart" dataset
data(heart)
#heart is a dataset with 172 obs and 8 variables
head(heart)
str(heart)
#check documentation
?heart

#renaming dataset as best practice
data1 = heart

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
#the results is confirming that our hypothesis of days time unit is coherent
summary(follow_up_days)
summary(follow_up_months)
summary(follow_up_years)

###do we want to change time unit??

#PART 2 : CREATION OF SURVIVAL OBJECT AND FITTING SURVIVAL CURVE

#use Kaplan-Meier estimator to have survival curve
heart_KM = survfit(Surv(stop, event) ~ 1, data = data1)
#the summary is telling us that the survival rate at the end of the study (1800 days) is the ~18%
#at time 1987 with a standard error of 0.0567 and a CI ranges from 9.73% and 33.80%
summary(heart_KM)


# Plot Kaplan-Meier estimator
plot(heart_KM, mark.time = TRUE,
     main = "Kaplan-Meier Estimator",
     ylab = "Survival Probability",
     xlab = "Time (days)")


#plot nicer, not working
library(ggplot2)
library(survminer)
library(broom)
library(tidyverse)
ggplot(heart_KM, aes(x = time, y = surv) +
  geom_step() +
  geom_point(data = subset(heart_KM, n.censor == 1), aes(x = time, y = surv), shape = 3) +  # Adding censoring marks
  labs(title = "Kaplan-Meier Estimator", x = "Time (days)", y = "Survival Probability") +
  theme_minimal())
  #or theme_bw(base_size = 18)


#other option: use Nelson-AAlen estimator to retrieve the cumulative hazard function
#ops, not working, giving error
#creating a new dataset because data1 has some factorizes variables and is giving error
heart_NA = survfit(Surv(stop, event) ~ 1, data = data1, type = "fh")
summary (heart_NA)

# Plot Nelson-AAlen estimator
###shape is not correct
plot(heart_NA, mark.time = TRUE,
     main = "Nelson-AAlen Estimator",
     ylab = "Cumulative Hazard",
     xlab = "Time (days)")

###another try but is not cumulative, why???
heart.kphaz.fit1=kphaz.fit(data1$stop, data1$event, method = "nelson")

heart.kphaz.fit1 %>%
  as.data.frame() %>%
  ggplot(aes(x = time, y = haz)) +
  geom_line() +
  xlab("Time") +
  ylab("h(t)")

#maybe good to add the logrank test having as two groups patients that had transplant yes or no
#to develop




#PART 3 : ANALYSIS AND MODELLING
#prepare the dataset for glm model
#prefer to work on the copy of the dataset not to modify the original one
data2 = (heart)

# Convert the categorical variables to factors, since this works better with glm ((event, surgery, transplant)
data2$event = as.factor(heart$event)
data2$surgery = as.factor(heart$surgery)
data2$transplant = as.factor(heart$transplant)

#verify the new structure
str(data2)

#now that the time unit is validated we can start creating the survival object

#glm model with binomial family
#set seed
set.seed(123)
heart_glm = glm(event ~ start + stop + age + year + surgery + transplant, data = data2, family = binomial)
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
heart_glm2 = glm(event ~ stop + year + transplant, data = data2, family = binomial)

#summary of second glm model
#the AIC is 180.25 so the model is better using these 3 variables
summary(heart_glm2)



#fit using COX hazard model
#use data1 because the data2 is just for glm model (due to factorized variables)
#we will do the comparison using nested and non-nested models with cox function
#we start with the null model to have a baseline for the likelihood comparison
cox_null = coxph(Surv(stop, event) ~ 1, data = data1)
#baseline likelihood is -314.1482
summary(cox_null)
#we start with a non-nested model (aka full model) with all variables included
cox_full = coxph(Surv(stop, event) ~ age + year + surgery + transplant, data = data1)
#based on this summary we can say that age, year and transplant are significant predictors of the event
#increasing of age increase the hazard 3.2% and decreasing of year decrease the hazard ratio of the event of 15.8%
#transplant decrease the hazard ratio of the event of 46.9%, surgery decrease the hazard of 46.9% but p-value is >0.05 so it is not statistically significant
#concordance value is 0.679 so the model has good discriminatory power
#likelihood is 22.29 p=2e-04
summary(cox_full)

#we continue with nested models
cox_nested1 = coxph(Surv(stop, event) ~ age, data = data1)
cox_nested2 = coxph(Surv(stop, event) ~ year, data = data1)
cox_nested3 = coxph(Surv(stop, event) ~ age + year, data = data1)
cox_nested4 = coxph(Surv(stop, event) ~ age + year + transplant, data = data1)

#based on concordance and likelihood ratio, respectively of 0.665 and 18.88, p=3e-04 the best nest model is the n. 4
#however the best of all cox model seems still the cox_full
summary(cox_nested1)
summary(cox_nested2)
summary(cox_nested3)
summary(cox_nested4)

#we will then compare the full model with the nested one n.4
#comparing nested models: LRT
#based on obtained likelihood we can say that model1 (in which we add the variable "surgery") fits better the model but
#however the improvement is not statistically significant as indicated by the p-value that is >0.05
anova(cox_full, cox_nested4)

#we can also double check comparing using: AIC
fits = list(cox_full = cox_full, cox_nested1= cox_nested1, cox_nested2 = cox_nested2, cox_nested3 = cox_nested3, cox_nested4 = cox_nested4)
#best model is confirmed the cox_full
sapply(fits, AIC)

#we can  also apply the step function to the full model
#also in this case the best model is the full one with no variables removed with AIC of 614.01
step_AIC = step(cox_full)
step_AIC


###NEXT STEPS
## Split data into a training and a testing set
## Train candidate models
## Make predictions in the testing dataset
## Assess predictive performance
### 1. Check that the predictions go into the right direction
### 2. Extract the C-statistic from the 'summary' output
### 3. Flip directions as needed
# Measuring predictions performance - effect size
## Raw estimates
## Standardized predictions
## Discretized predictions - based on quantiles
# Model diagnostics

## Martingale residuals
## Case-deletion residuals
### Schoenfeld residuals
## Stratification
#. AUC/ROC. 

#TO DO:
# Stratified Cox Regression
# Penalized Cox regression
# cross validation (lasso, ridge, elasticnet)
# predicting survival
# Interpretations
# embellishment using better plotting package and knit


#NOT YET USED CODE


#check Matringale residuals for cox model
data1$residuals = residuals(heart_cox, type = 'martingale')
summary(data1$residuals)

#plot cox model residuals
#####to review, plot is weird
with(data1, {
  plot(age, data1$residuals)
  lines(lowess(age, data1$residuals), lwd = 2)
  plot(year, data1$residuals)
  plot(surgery, data1$residuals)
})

 



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

