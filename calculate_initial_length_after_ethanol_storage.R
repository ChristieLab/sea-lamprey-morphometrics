#---------------- Lamprey Shrinkage Model -------------------------------------------------------------------------#

# Written by Allison Nalesnik on February 5, 2024
# Contact at al.nalesnik21@gmail.com or analesni@purdue.edu
# R version 4.2.3

# Task: to calculate initial length of larval lamprey based on number of days stored in ethanol and total length

####################################################################################################################

setwd("/lamprey")

#---------------- Load data ---------------------------------------------------------------------------------------#
dat <- read.csv("dat.csv", header=TRUE, sep=",", na.strings="?", dec=".", strip.white=TRUE) # your data file goes here

#---------------- Define the function -----------------------------------------------------------------------------#
# length is the measured length of the fish after given days in ethanol
# days is the number of days that the fish has been stored in ethanol
# B is the rate of size decay

expDecayModel <- function(days, length, B) {length / (exp(-B * days))}

#---------------- Using the function ------------------------------------------------------------------------------#
# This will make a column called initial_length which will contain the calculated initial length values.
# 4.483e-04 is the rate of decay (calculated by Allison at Purdue, manuscript submitted in CJFAS 09/2024)

dat$initial_length <- expDecayModel(dat$days, dat$length, 4.483e-04)

################################################################################################################ end