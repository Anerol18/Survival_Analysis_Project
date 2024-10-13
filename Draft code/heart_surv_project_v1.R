# PART 1 : LOAD DATASET, INSPECTION, DATA CLEANSING
# download the the survival library
library(survival)
library(survivalROC)
library(timeROC)
library(ggplot2)
library(survminer)
library(broom)
library(tidyverse)
library(muhaz)
library(VSURF)
library(glmnet)
library(Hmisc)
library(survAUC)
library(pec)
library(riskRegression)
library(ggplot2)

# inspect the "heart" dataset
data(heart)
# heart is a dataset with 172 obs and 8 variables
head(heart)
str(heart)
#check documentation
?heart

#renaming dataset as best practice
data1 = heart

#verify the time unit since is not specified in the documentation, so we need 
# to give our hypothesis giving the context
# the max(data1$year) is 6.47 that in days is ~2363.
# This is in line with survival studies that are generally long-term ones.
max(data1$year)
# max(data1$stop) is 1800 day if we divide by years (365.25) gives ~5 years.
# that is still coherent with days time unit.
max(data1$stop)

# Double check based on the supposition that the time units is in days
# Calculation of follow-up time in different time units and check if seems coherent 
follow_up_days = data1$stop - data1$start
follow_up_months = follow_up_days / 30.44 # approx avg per month
follow_up_years = follow_up_days / 365.25 # approx avg per year

#summary on different time units
#the results is confirming that our hypothesis of days time unit is coherent
summary(follow_up_days)
summary(follow_up_months)
summary(follow_up_years)

### Do we want to change time unit??

# PART 2 : CREATION OF SURVIVAL OBJECT AND FITTING SURVIVAL CURVE

# Using Kaplan-Meier estimator to have survival curve
heart_KM = survfit(Surv(stop, event) ~ 1, data = data1)
#the summary is telling us that the survival rate at the end of the study 
# (1800 days) is ~18% at time 1987 with a standard error of 0.0567 and a CI 
# ranges from 9.73% and 33.80%
summary(heart_KM)

# Plot Kaplan-Meier estimator
plot(heart_KM, mark.time = TRUE,
     main = "Kaplan-Meier Estimator",
     ylab = "Survival Probability",
     xlab = "Time (days)")

# Plot nicer, not working

ggplot(heart_KM, aes(x = time, y = surv) +
  geom_step() +
  geom_point(data = subset(heart_KM, n.censor == 1),
             aes(x = time, y = surv), shape = 3) +  # Adding censoring marks
  labs(title = "Kaplan-Meier Estimator",
       x = "Time (days)",
       y = "Survival Probability") +
  theme_minimal())
  #or theme_bw(base_size = 18)

##################################################################################
# NOT WORKING
# Other option: use Nelson-AAlen estimator to retrieve the cumulative hazard
# function
#ops, not working, giving error

# Creating a new dataset because data1 has some factorizes variables and is 
# giving error
heart_NA = survfit(Surv(stop, event) ~ 1, data = data1, type = "fh")
summary (heart_NA)

# Plot Nelson-AAlen estimator
### shape is not correct
plot(heart_NA, mark.time = TRUE,
     main = "Nelson-AAlen Estimator",
     ylab = "Cumulative Hazard",
     xlab = "Time (days)")

### another try but is not cumulative, why???
heart.kphaz.fit1=kphaz.fit(data1$stop, data1$event, method = "nelson")

heart.kphaz.fit1 %>%
  as.data.frame() %>%
  ggplot(aes(x = time, y = haz)) +
  geom_line() +
  xlab("Time") +
  ylab("h(t)")

######################################################################################
# PART 3 : ANALYSIS AND MODELLING
# Prepare the data set for GLM model
# Prefer to work on the copy of the data set not to modify the original one
data2 = (heart)

# Convert the categorical variables to factors, since this works better with 
# glm ((event, surgery, transplant)
data2$event = as.factor(heart$event)
data2$surgery = as.factor(heart$surgery)
data2$transplant = as.factor(heart$transplant)

# Verify the new structure
str(data2)

# Now that the time unit is validated we can start creating the survival object

# GLM model with binomial family
set.seed(123)
heart_glm = glm(event ~ start + stop + age + year + surgery + 
                  transplant, data = data2, family = binomial)
# Summary of the glm model using all covariates says that a part from the 
# intercept only stop, year and transplant1 are statistically significant
# the AIC is 185.3
summary(heart_glm)

plot(heart_glm)

# let's double-check what VSURF says about the variables importance
heart_vs=VSURF(event~.,data=data1)

#check variable importance result
heart_vs[[3]]

# Plot it to visualize variable names. results: Transplant, stop, year
plot(heart_vs,step="pred",var.names=TRUE)

# Rerun glm model with the selected variables confirmed by VSURF
heart_glm2 = glm(event ~ stop + year + transplant, 
                 data = data2, family = binomial)

# Summary of second GLM model
# AIC is 180.25 so the model is better using these 3 variables
summary(heart_glm2)

plot(heart_glm2)

# Fit using COX hazard model
# Use data1 because the data2 is just for glm model (due to factorized variables)
# We will do the comparison using nested and non-nested models with cox function
# We start with the null model to have a baseline for the likelihood comparison
cox_null = coxph(Surv(stop, event) ~ 1, data = data1)

# Baseline likelihood is -314.1482
summary(cox_null)

# We start with a non-nested model (aka full model) with all variables included
cox_full = coxph(Surv(stop, event) ~ age + year + 
                   surgery + transplant, data = data1)

# Based on this summary we can say that age, year and transplant are significant 
# predictors of the event increasing of age increase the hazard 3.2% and 
# decreasing of year decrease the hazard ratio of the event of 15.8%
# Transplant decrease the hazard ratio of the event of 46.9%, surgery decrease 
# the hazard of 46.9% but p-value is >0.05 so it is not statistically 
# significant
# Concordance value is 0.679 so the model has good discriminatory power
# likelihood is 22.29 p=2e-04
summary(cox_full)

# We continue with nested models
cox_nested1 = coxph(Surv(stop, event) ~ age, data = data1)
cox_nested2 = coxph(Surv(stop, event) ~ year, data = data1)
cox_nested3 = coxph(Surv(stop, event) ~ age + year, data = data1)
cox_nested4 = coxph(Surv(stop, event) ~ age + year + transplant, data = data1)

# Based on concordance and likelihood ratio, respectively of 0.665 and 18.88, 
# p=3e-04 the best nest model is the n. 4
# However the best of all cox model seems still the cox_full
summary(cox_nested1)
summary(cox_nested2)
summary(cox_nested3)
summary(cox_nested4)

# We will then compare the full model with the nested one n.4
# Comparing nested models: LRT
# Based on obtained likelihood we can say that model1 (in which we add the 
# variable "surgery") fits better the model.
# However the improvement is not statistically significant as indicated by 
# the p-value that is > 0.05
anova(cox_full, cox_nested4)

# We can also double check comparing using: AIC
fits = list(cox_full = cox_full, 
            cox_nested1= cox_nested1, 
            cox_nested2 = cox_nested2, 
            cox_nested3 = cox_nested3, 
            cox_nested4 = cox_nested4)

# Best model is confirmed the cox_full
sapply(fits, AIC)

# We can  also apply the step function to the full model
# In this case the best model is the full one with no variables removed with 
# AIC of 614.01
step_AIC = step(cox_full)
step_AIC

# Next step is to verify if we can have a better model using stratification
# 1st test stratification applied on surgery
# Based on concordance and likelihood ratio, respectively of 0.647 and 16.43, 
# p=9e-04 the best model remains the cox_full
cox_strata1 = coxph(Surv(stop, event) ~ age + year + 
                      strata(surgery) + transplant, data = data1)
summary(cox_strata1)

# 2nd test stratification applied on transplant
#based on concordance and likelihood ratio, respectively of 0.645 and 17.38, 
# p=6e-04 the best model remains the cox_full
cox_strata2 = coxph(Surv(stop, event) ~ age + year + 
                      surgery + strata(transplant), data = data1)
summary(cox_strata2)

# As last test we can calculate the ROC curve
####NOT WORKING YET THE ROC
####################################################################

cox_ROC = coxph(Surv(stop, event) ~ age + year + 
                  surgery + transplant, data = data1)

#predict survival probability
predict_probs = predict(cox_ROC, type = "lp", se.fit = TRUE)

plot(predict_probs$fit)

#extract predicted probabilities and survival status
time = data1$stop
pred = predict_probs$fit
status = data1$event

#ROC curve calculation
roc = timeROC(status = status, marker = -pred, times = time)
# Error
######################################################################

#PART 4 : PREDICTION
# Split the data set in training and testing
set.seed(123)
i.training = sample.int(nrow(data1), size = floor(0.8 * nrow(data1)), replace = FALSE)
i.testing = setdiff(seq_len(nrow(data1)), i.training)
d_training = data1[i.training, ]
d_testing = data1[i.testing, ]

# Train the model using training data set
cox_full_tr = coxph(Surv(stop, event) ~ age + year + surgery + transplant, data = d_training)
cox_nested4_tr = coxph(Surv(stop, event) ~ age + year + transplant, data = d_training)

# We calculate prediction on testing data set
# Make predictions in the testing data set
d_testing$lp_A = predict(cox_full_tr, newdata = d_testing, type = "lp")
d_testing$lp_B = predict(cox_nested4_tr, newdata = d_testing, type = "lp")

d_testing

# Asses predictive performance
models = list(
  A = coxph(Surv(stop, event) ~ lp_A, data = d_testing),
  B = coxph(Surv(stop, event) ~ lp_B, data = d_testing)
)

summary(models$A)
summary(models$B)

plot(models$A$linear.predictors)

# Check that the predictions go into the right direction
sapply(models, coef)

# Extract the C-statistic from the 'summary' output
map_dbl(models, ~ summary(.)$concordance[1])

# Flip directions as needed
# Benchmark is selecting model A (that is our original cox_full)
benchmark =
  tibble(
    model = names(models),
    sign = sapply(models, coef) |> sign(),
    C_summary = map_dbl(models, ~ summary(.)$concordance[1]),
    C = ifelse(sign > 0, C_summary, 1 - C_summary)
  )
benchmark


# Measuring predictions performance - effect size
d1 = d_testing |> select(stop, event, lp_A, lp_B)
head(d1)

# Raw estimates
A = coxph(Surv(stop, event) ~ lp_A, data = d1) |> tidy()
B = coxph(Surv(stop, event) ~ lp_B, data = d1) |> tidy()
bind_rows(A, B)


#standardized predictions
d2 = mutate(d1,
             ZA = lp_A / sd(lp_A),
             ZB = lp_B / sd(lp_B))
head(d2)


# Discretized predictions - based on quantiles
d3 =
  mutate(d2,
         FA = factor((lp_A > median(lp_A)), levels = c(FALSE, TRUE), labels = c("low", "high")),
         FB = factor((lp_B > median(lp_B)), levels = c(FALSE, TRUE), labels = c("low", "high"))
  )

A = coxph(Surv(stop, event) ~ FA, data = d3) |> tidy()
B = coxph(Surv(stop, event) ~ FB, data = d3) |> tidy()
bind_rows(A, B) |>
  transmute(term, estimate, HR = exp(estimate), p.value)

fit.KMA = survfit(Surv(stop, event) ~ FA, data = d3)
fit.KMA
plot(fit.KMA, col = 1:2)

fit.KMB = survfit(Surv(stop, event) ~ FB, data = d3)
fit.KMB
plot(fit.KMB, col = 1:2)


# Cross-validation using Lasso model

# Prepare data for training dataset
x_train = as.matrix(d_training[, c("age", "year", "surgery", "transplant")]) 
y_train = as.matrix(Surv(d_training$stop, d_training$event))

# Prepare data for testing dataset
x_test = as.matrix(d_testing[, c("age", "year", "surgery", "transplant")]) 
y_test = as.matrix(Surv(d_testing$stop, d_testing$event))

# Set up cross-validation with Lasso
cvfit = cv.glmnet(x_train, y_train, alpha = 1, family = "cox")

# Retrieve optimal lambda
lambda_min = cvfit$lambda.min
lambda_min

# Fit with Lasso model
lasso_model = glmnet(x_train,y_train,alpha = 1,family="cox",lambda=lambda_min)

# Predict on test set
lasso_pred = predict(lasso_model,newx=x_test,s=lambda_min,type="link")

################################################################################
# Calculate C-index
### NOT WORKING YET
c_index <- survAUC::UnoC(Surv.rsp = y_test, Surv.rsp.new = y_test, lpnew = lasso_pred)
c_index

################################################################################

# Log-Rank Test
# The lasso model is statistically significant since p-value is very small
surv_diff = survdiff(Surv(d_testing$stop, d_testing$event) ~ lasso_pred)
p_value <- 1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1)
p_value

# The coefficients are all zero, this can happen with Lasso regression when the 
# regularization is strong
coef(cvfit, s = "lambda.1se")

################################################################################
# Create Calibration Plot
### NOT WORKING YET
pred_risk <- predictRisk(lasso_model, newdata = d_testing, 
                         times = seq(0, max(d_testing$stop), by = 1))
# Error

# Prepare the data for the calibration plot
cal_data <- data.frame(time = d_testing$stop, 
                       status = d_testing$event, pred_risk = pred_risk)

# Generate the calibration plot using ggplot2

ggplot(cal_data, aes(x = time, y = pred_risk)) +
  geom_point() +
  geom_smooth(method = "loess") +
  labs(title = "Calibration Plot", x = "Time (days)", y = "Predicted Risk") +
  theme_minimal()

################################################################################

## Martingale residuals??



# PART 5 : Conclusion
# fit.KMB might be considered better due to its higher median survival time 
# (FB=high) compared to fit.KMA (FA=high)?


# Interpretations
# embellishment using better plotting package and knit


# NOT YET USED CODE


# Check Matringale residuals for cox model
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

