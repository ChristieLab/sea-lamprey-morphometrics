# 8000+ larval lamprey were exposed to TFM over the course of 5 experiments. The below code works to write and analyze
# a model for this data.
  
setwd("/Users/allison/Documents/Research/lamprey") # Set directory
library(dplyr)
library(tidyr)
library(data.table)
library(lme4)
library(MuMIn)

# Add filtered data for making the model
dat<- read.csv("data/experiment_fish/all_experiment_fish.csv", header=TRUE, sep=",", na.strings="?", dec=".", strip.white=TRUE)
dat <- as.data.table(dat)
dat <- dat[,-1]
colnames(dat)[3] <- "bottle_id"
dat[,mass:=as.numeric(mass)]
dat[,condition:= as.numeric(condition)]

## MODEL A
mod <- lmerTest::lmer(hour ~ total_length*mass + (1 | experiment/tank/bag), data = dat, REML = FALSE)

## MODEL B
mod1 <- lmerTest::lmer(hour ~ total_length + mass + condition + (1 | experiment/tank/bag), data = dat, REML = FALSE)

## MODEL C
mod4 <- lmerTest::lmer(hour ~ condition + (1 | experiment/tank/bag), data = dat, REML = FALSE) 

                       
############################################################################################
# <>-<>-<>-<>-<>-<>-<>-<>-<>- Compare models with ANOVA <>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-

# ANOVA to test two models. Model with total_length*mass is best
anova(mod, mod1, mod4)

# Summary of the GLMM
summary(mod4)                                                                         # Wald Test is the estimates -- get p values
r.squaredGLMM(mod4)                                                                   # R^2 value
summary(mod4)$coefficients[,4]                                                        # unrounded p-value
conf <- confint(mod4, parm = "condition", level = 0.95)   # get confidence intervals
print(conf)

# <>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-


## MODEL D
lmod <- lmerTest::lmer(hour ~ total_length + (1 | experiment/tank/bag), data = dat, REML = FALSE)
summary(lmod)
r.squaredGLMM(lmod)                                            # condition r2 = ______
summary(lmod)$coefficients[,5]                                 # get unrounded p-value
confl <- confint(lmod, parm = "total_length", level = 0.95)    # confidence intervals
print(confl)

## MODEL E

mmod <- lmerTest::lmer(hour ~ mass + (1 | experiment/tank/bag), data = dat, REML = FALSE)
summary(mmod)
r.squaredGLMM(mmod)
summary(mmod)$coefficients[,5]                                 # get unrounded p-value
confm <- confint(mmod, parm = "mass", level = 0.95)     # confidence intervals
print(confm)

## MODEL F
lmmod <- lmerTest::lmer(hour ~ total_length + mass+ (1 | experiment/tank/bag), data = dat, REML = FALSE)
summary(lmmod)
r.squaredGLMM(lmmod)
summary(lmmod)$coefficients[,5]     
conflm <- confint(lmmod, parm = c("total_length","mass"), level = 0.95)
print(conflm)


# <>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-

# Predicted percentile survival times

## Calculate percentiles
percentiles <- c(0.10, 0.25, 0.50, 0.75, 0.90)
total_length_percentiles <- quantile(dat$total_length, percentiles)
mass_percentiles <- quantile(dat$mass, percentiles, na.rm = TRUE)

# Print the results
print("Total Length Percentiles:")
print(total_length_percentiles)
print("Mass Percentiles:")
print(mass_percentiles)

# Length -->   10%   25%   50%   75%   90% 
#             70.0  77.5  85.0  94.0 102.0 

# Mass -->     10%   25%   50%   75%   90% 
#            0.399 0.569 0.796 1.102 1.445 


mod <- lmerTest::lmer(hour ~ total_length*mass + (1 | experiment/tank/bag), data = dat, REML = FALSE)

# make dataframe of percentiles
pdata <- expand.grid(total_length = total_length_percentiles, mass = mass_percentiles)
pdata$predicted_hour <- predict(mod, newdata = pdata, re.form = NA)  # re.form=NA excludes random effects

# make matrix
library(tidyr)
wide_data <- pivot_wider(pdata, names_from = total_length, values_from = predicted_hour, id_cols = mass)

# <>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-

# Wilkie's data model

mw_fish <- read.csv("admin/Survivorship_manuscript/archive/materials_from_Wilke/mw_data.csv", header = TRUE)
mw_fish <- mw_fish[1:48,]
mw_fish$new_condition <- (mw_fish$mass * (10^5)) / (mw_fish$total_length ^ 2.6)


# linear model
mw_mod <- lm(hour ~ total_length*mass, data = mw_fish)                       # LENGTH, MASS, INTERACTION
mw_mod <- lm(hour ~ total_length + mass + total_length*mass, data = mw_fish) # LENGTH, MASS, INTERACTION -- report effect size, report in minutes (cm, g)
mw_mod <- lm(hour ~ total_length + mass, data = mw_fish)  # LENGTH AND MASS
summary(mw_mod)

# <>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-
#dat2 <- subset(dat, subset = entry != 47 | 155 | 1  | 176 | 55 | 52 | 9)
dat2 <- dat[8:8604,]

only_length <- lmerTest::lmer(hour ~ total_length + (1 | experiment/tank/bag), data = dat2, REML = FALSE)
only_mass <- lmerTest::lmer(hour ~ mass + (1 | experiment/tank/bag), data = dat2, REML = FALSE)
only_condition <- lmerTest::lmer(hour ~ condition + (1 | experiment/tank/bag), data = dat2, REML = FALSE)
length_mass <- lmerTest::lmer(hour ~ total_length + mass + (1 | experiment/tank/bag), data = dat2, REML = FALSE)
length_mass_int <- lmerTest::lmer(hour ~ total_length*mass + (1 | experiment/tank/bag), data = dat2, REML = FALSE)

anova(only_condition, only_length, only_mass, length_mass_int, length_mass)

# Summary of the GLMM
summary(mod)  

only_lengthmw <- lm(hour ~ total_length, data = mw_fish)
only_massmw <- lm(hour ~ mass, data = mw_fish)
only_conditionmw <- lm(hour ~ new_condition, data = mw_fish)
length_massmw <- lm(hour ~ total_length + mass, data = mw_fish)
length_mass_intmw <- lm(hour ~ total_length*mass, data = mw_fish)

anova(only_conditionmw, only_lengthmw, only_massmw, length_mass_intmw, length_massmw)

summary(length_mass_intmw)

## Final MW model reported in text
modmw <- lm(hour ~ total_length*mass, data = mw_fish)
summary(modmw)

## Interaction estimate is 0.05599





