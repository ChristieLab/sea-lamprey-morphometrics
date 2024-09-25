#---------------- Data Description -------------------------------------------------------------------------------------#
# script written by allison nalesnik on february 5, 2024; contact at al.nalesnik21@gmail.com
# R version 4.2.3
# data consists of multiple morphology measurements for fish preserved in ethanol 

#set working directory and load libraries
setwd("/Users/allison/Documents/Research/lamprey")
library(dplyr)
library(ggplot2)
library(ggthemes)
library(tidyr)
library(data.table)
library(nlme) # Non Linear Mixed Effect
library("lmerTest") # Linear Mixed Effect from package lme4

#----------------------------- Data loading and prep -------------------------------------------------------------#
dat_hbbs<- read.csv("dat_hbbs.csv", header=TRUE, sep=",", na.strings="?", dec=".", strip.white=TRUE)
dat_arl<- read.csv("dat_arl.csv", header=TRUE, sep=",", na.strings="?", dec=".", strip.white=TRUE)
dat_arl <- as.data.table(dat_arl)
dat_hbbs <- as.data.table(dat_hbbs)
dat_arl$mass_g <- as.numeric(dat_arl$mass_g)
dat_arl$condition <- as.numeric(dat_arl$condition)
dat_hbbs$site <- as.character("hbbs")
dat_arl$site <- as.character("arl")

# merged dataset of arl and hbbs data
all_dat <- rbind(dat_arl, dat_hbbs)

#--<>--<>--<>--<>--<>--<>-- Generate independent samples --<>--<>--<>--<>--<>--<>--<>
# randomly select individuals at random time points to model with independent samples
##  initiate empty data table
all_dat_sample <- data.table(unique_ID = character(0), init_length = numeric(0), init_mass = numeric(0),
                             init_condition = numeric(0), length_mm = numeric(0), mass_g = numeric(0),
                             condition = numeric(0), days_in_ethanol = numeric(0),
                             site= character(0))

## get list of unique IDs (one per sample)
uids_all <- unique(all_dat$unique_ID)

## randomly sample one row from each unique_ID to make all_dat_sample (datatable of independent samples)
for(i in uids_all){
  tall <- all_dat[unique_ID == i,]
  tall <- tall[sample(nrow(tall), 1), ] 
  all_dat_sample <- rbind(all_dat_sample, 
                          tall[,c("unique_ID","init_length", "init_mass", "init_condition", "length_mm", "mass_g",
                          "condition", "days_in_ethanol", "site")]) 
}

# Proceed with all_dat_sample (n=664)

###############################################################################################################
#---------------- EXPONENTIAL DECAY MODEL --------------------------------------------------------------------#

#---------------- Define the function ------------------------------------------------------------------------#
# A --> represents the initial size (e.g., length or mass) of the fish on day 0
# B --> represents the rate of decay. It determines how quickly the size decreases over time.
# A positive value of B indicates a decrease in size over time.
# The magnitude of B indicates how rapidly the size changes. A larger value of B means a faster decay.
# B is the Estimate when you look at summary(model)

#expDecayModel0 <- function(days_in_ethanol, init_length, B) {init_length * exp(-B * days_in_ethanol)} # calculates length given initial length and days in ethanol

expDecayModel1 <- function(days_in_ethanol, length_mm, B) {length_mm / (exp(-B * days_in_ethanol))} # calculates init_length

#---------------- Construction of Models Based on Predictor Variable -----------------------------------------#

# Model with all_dat
mod0 <- nls(init_length ~ expDecayModel1(days_in_ethanol, length_mm, B),  # predicting initial length
            fixed = site,
            start = c(B = 0.01), data = all_dat_sample)

summary(mod0)
# estimate = 4.483e-04
# std. error = 3.222e-05
# t-value = 13.91
# p-value = <2e-16 (8.750943e-39)
# Residual standard error: 3.81 on 663 degrees of freedom

#---------------- Leave-one-out Cross Evaluation -----------------------------------------#
##  # Initialize an empty vector to store predicted values
all_dat_sample$lpred <- numeric(nrow(all_dat_sample))

#  # Loop through each data point for cross-validation
for (i in 1:nrow(all_dat_sample)){
  
  train_data <- all_dat_sample[-i, ]                                           # Create a training dataset by removing the ith data point
  
  this_mod <- nls(init_length ~ expDecayModel1(days_in_ethanol, length_mm, B), # Fit the model using the training dataset
                  start = c(B = 0.01), data = train_data)
  
  prediction <- expDecayModel1(all_dat_sample$days_in_ethanol[i],              # Predict the initial length for the left-out data point
                               all_dat_sample$length_mm[i],
                               coef(this_mod)["B"])
  
  all_dat_sample$lpred[i] <- prediction # Store the prediction in the original dataset
}

# Plot
ggplot(all_dat_sample, aes(init_length, lpred)) + geom_point() + geom_abline()
plot(all_dat_sample$init_length, all_dat_sample$lpred) # *this looks good, better than the nls model*


#---------------- Difference between Predicted and Known values -----------------------------------------------#
# make pred_diff column. If -, then overpredicting. If +, then underpredicting
all_dat_sample$pred_diff <- (all_dat_sample$init_length - all_dat_sample$lpred)

# plot to see if this is related to days_in_ethanol
ggplot(all_dat_sample, aes(days_in_ethanol, pred_diff)) + geom_point()
range(all_dat_sample$pred_diff)
hist(all_dat_sample$pred_diff, breaks= 40)

# linear model between actual values and predicted values
mod1 <- lm(all_dat_sample$init_length ~ all_dat_sample$lpred)
summary(mod1)
summary(mod1)$coefficients[,4]

###############################################################################################################

#---------------- Make independent samples for modeling ----------------------------------------------------#
# HBBS RANDOM SAMPLES
dat_hbbs_sample <- data.table(unique_ID = character(0), init_length = numeric(0), init_mass = numeric(0),
                              init_condition = numeric(0), length_mm = numeric(0), mass_g = numeric(0),
                              condition = numeric(0), days_in_ethanol = numeric(0)) # initiate empty data table

uids_hbbs <- unique(dat_hbbs$unique_ID) # get list of unique IDs (one per sample)

# randomly sample one row from each unique_ID to make dat_arl_sample (dt of independent samples)
for(i in uids_hbbs){
  thbbs <- dat_hbbs[unique_ID == i,]
  thbbs <- thbbs[sample(nrow(thbbs), 1), ] 
  dat_hbbs_sample <- rbind(dat_hbbs_sample, 
                           thbbs[,c("unique_ID","init_length", "init_mass", "init_condition", "length_mm", "mass_g",
                                    "condition", "days_in_ethanol")]) 
}

# ARL RANDOM SAMPLES
dat_arl_sample <- data.table(unique_ID = character(0), init_length = numeric(0), init_mass = numeric(0),
                              init_condition = numeric(0), length_mm = numeric(0), mass_g = numeric(0),
                              condition = numeric(0), days_in_ethanol = numeric(0)) # initiate empty data table

uids_arl <- unique(dat_arl$unique_ID) # get list of unique IDs (one per sample)

# randomly sample one row from each unique_ID to make dat_arl_sample (dt of independent samples)
for(i in uids_arl){
  tarl <- dat_arl[unique_ID == i,]
  tarl <- tarl[sample(nrow(tarl), 1), ] 
  dat_arl_sample <- rbind(dat_arl_sample, 
                           tarl[,c("unique_ID","init_length", "init_mass", "init_condition", "length_mm", "mass_g",
                                    "condition", "days_in_ethanol")]) 
}

## Supplemental models

#Mass - NLS
expDecayModel2 <- function(days_in_ethanol, mass_g, B) {mass_g / (exp(-B * days_in_ethanol))} # calculates init_mass
mod0.1 <- nls(init_mass ~ expDecayModel2(days_in_ethanol, mass_g, B),
              start = c(B = 0.01),
              data = all_dat_sample)
summary(mod0.1)
summary(mod0.1)$coefficients[,4] # get unrounded p-value

#Mass - LM between known and predicted
all_dat_sample$mpred <- numeric(nrow(all_dat_sample))

#  # Loop through each data point for cross-validation
for (i in 1:nrow(all_dat_sample)){
  
  train_data <- all_dat_sample[-i, ]                                           # Create a training dataset by removing the ith data point
  
  this_mod <- nls(init_mass ~ expDecayModel2(days_in_ethanol, mass_g, B), # Fit the model using the training dataset
                  start = c(B = 0.01), data = train_data)
  
  prediction <- expDecayModel2(all_dat_sample$days_in_ethanol[i],              # Predict the initial length for the left-out data point
                               all_dat_sample$mass_g[i],
                               coef(this_mod)["B"])
  
  all_dat_sample$mpred[i] <- prediction # Store the prediction in the original dataset
}

mod1.1 <- lm(all_dat_sample$init_mass ~ all_dat_sample$mpred)
summary(mod1.1)
summary(mod1.1)$coefficients[,4] # get unrounded p-value


#CONDITION
expDecayModel3 <- function(days_in_ethanol, condition, B) {condition / (exp(-B * days_in_ethanol))} # calculates init_condition

mod_c1 <- nls(init_condition ~ expDecayModel3(days_in_ethanol, condition, B), # predicting initial condition
              start = c(B = 0.01),
              data = all_dat_sample)
summary(mod_c1)
summary(mod_c1)$coefficients[,4] # get unrounded p-value

#Condition - LM between known and predicted
all_dat_sample$cpred <- numeric(nrow(all_dat_sample))

#  # Loop through each data point for cross-validation
for (i in 1:nrow(all_dat_sample)){
  
  train_data <- all_dat_sample[-i, ]                                           # Create a training dataset by removing the ith data point
  
  this_mod <- nls(init_condition ~ expDecayModel3(days_in_ethanol, condition, B), # Fit the model using the training dataset
                  start = c(B = 0.01), data = train_data)
  
  prediction <- expDecayModel3(all_dat_sample$days_in_ethanol[i],              # Predict the initial length for the left-out data point
                               all_dat_sample$condition[i],
                               coef(this_mod)["B"])
  
  all_dat_sample$cpred[i] <- prediction # Store the prediction in the original dataset
}

mod_c2 <- lm(all_dat_sample$init_condition ~ all_dat_sample$cpred)
summary(mod_c2)
summary(mod_c2)$coefficients[,4]

