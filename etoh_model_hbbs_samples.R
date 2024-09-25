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
dat_hbbs <- as.data.table(dat_hbbs)
dat_hbbs$site <- as.character("hbbs")


#--<>--<>--<>--<>--<>--<>-- Generate independent samples --<>--<>--<>--<>--<>--<>--<>
# randomly select individuals at random time points to model with independent samples
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



###############################################################################################################
#---------------- EXPONENTIAL DECAY MODEL --------------------------------------------------------------------#

#---------------- Define the function ------------------------------------------------------------------------#
# A --> represents the initial size (e.g., length or mass) of the fish on day 0
# B --> represents the rate of decay. It determines how quickly the size decreases over time.
# A positive value of B indicates a decrease in size over time.
# The magnitude of B indicates how rapidly the size changes. A larger value of B means a faster decay.
# B is the Estimate when you look at summary(model)

expDecayModel1 <- function(days_in_ethanol, length_mm, B) {length_mm / (exp(-B * days_in_ethanol))} # calculates init_length

#---------------- Construction of Models Based on Predictor Variable -----------------------------------------#
mod0 <- nls(init_length ~ expDecayModel1(days_in_ethanol, length_mm, B),  # predicting initial length
            start = c(B = 0.01), data = dat_hbbs_sample)

summary(mod0)
summary(mod0)$coefficients[,4]

#---------------- Leave-one-out Cross Evaluation -----------------------------------------#
##  # Initialize an empty vector to store predicted values
dat_hbbs_sample$lpred <- numeric(nrow(dat_hbbs_sample))

#  # Loop through each data point for cross-validation
for (i in 1:nrow(dat_hbbs_sample)){
  
  train_data <- dat_hbbs_sample[-i, ]                                           # Create a training dataset by removing the ith data point
  
  this_mod <- nls(init_length ~ expDecayModel1(days_in_ethanol, length_mm, B), # Fit the model using the training dataset
                  start = c(B = 0.01), data = train_data)
  
  prediction <- expDecayModel1(dat_hbbs_sample$days_in_ethanol[i],              # Predict the initial length for the left-out data point
                               dat_hbbs_sample$length_mm[i],
                               coef(this_mod)["B"])
  
  dat_hbbs_sample$lpred[i] <- prediction # Store the prediction in the original dataset
}

# Plot
ggplot(dat_hbbs_sample, aes(init_length, lpred)) + geom_point() + geom_abline()
plot(dat_hbbs_sample$init_length, dat_hbbs_sample$lpred) # *this looks good, better than the nls model*


# linear model between actual values and predicted values
mod1 <- lm(dat_hbbs_sample$init_length ~ dat_hbbs_sample$lpred)
summary(mod1)
summary(mod1)$coefficients[,4]

###############################################################################################################

#Mass - NLS
expDecayModel2 <- function(days_in_ethanol, mass_g, B) {mass_g / (exp(-B * days_in_ethanol))} # calculates init_mass
mod0.1 <- nls(init_mass ~ expDecayModel2(days_in_ethanol, mass_g, B),
              start = c(B = 0.01),
              data = dat_hbbs_sample)
summary(mod0.1)
summary(mod0.1)$coefficients[,4] # get unrounded p-value

#Mass - LM between known and predicted
dat_hbbs_sample$mpred <- numeric(nrow(dat_hbbs_sample))

#  # Loop through each data point for cross-validation
for (i in 1:nrow(dat_hbbs_sample)){
  
  train_data <- dat_hbbs_sample[-i, ]                                           # Create a training dataset by removing the ith data point
  
  this_mod <- nls(init_mass ~ expDecayModel2(days_in_ethanol, mass_g, B), # Fit the model using the training dataset
                  start = c(B = 0.01), data = train_data)
  
  prediction <- expDecayModel2(dat_hbbs_sample$days_in_ethanol[i],              # Predict the initial length for the left-out data point
                               dat_hbbs_sample$mass_g[i],
                               coef(this_mod)["B"])
  
  dat_hbbs_sample$mpred[i] <- prediction # Store the prediction in the original dataset
}

mod1.1 <- lm(dat_hbbs_sample$init_mass ~ dat_hbbs_sample$mpred)
summary(mod1.1)
summary(mod1.1)$coefficients[,4]


#CONDITION
expDecayModel3 <- function(days_in_ethanol, condition, B) {condition / (exp(-B * days_in_ethanol))} # calculates init_condition

mod_c1 <- nls(init_condition ~ expDecayModel3(days_in_ethanol, condition, B), # predicting initial condition
              start = c(B = 0.01),
              data = dat_hbbs_sample)
summary(mod_c1)
summary(mod_c1)$coefficients[,4] # get unrounded p-value

#Condition - LM between known and predicted
dat_hbbs_sample$cpred <- numeric(nrow(dat_hbbs_sample))

#  # Loop through each data point for cross-validation
for (i in 1:nrow(dat_hbbs_sample)){
  
  train_data <- dat_hbbs_sample[-i, ]                                           # Create a training dataset by removing the ith data point
  
  this_mod <- nls(init_condition ~ expDecayModel3(days_in_ethanol, condition, B), # Fit the model using the training dataset
                  start = c(B = 0.01), data = train_data)
  
  prediction <- expDecayModel3(dat_hbbs_sample$days_in_ethanol[i],              # Predict the initial length for the left-out data point
                               dat_hbbs_sample$condition[i],
                               coef(this_mod)["B"])
  
  dat_hbbs_sample$cpred[i] <- prediction # Store the prediction in the original dataset
}

mod_c2 <- lm(dat_hbbs_sample$init_condition ~ dat_hbbs_sample$cpred)
summary(mod_c2)
summary(mod_c2)$coefficients[,4]

