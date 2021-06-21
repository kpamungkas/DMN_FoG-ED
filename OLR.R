#### Ordinal Logistic Regression ####
#
# Hypothesis 1: DMN connectivity would be predictive of 
# DBS-induced FoG changes in PD patients with STN-DBS
#
# Classification of DBS-induced freezing of gait changes
# based on VTA-based Average Resting-State Functional Connectivity
# of Default Mode Network
# Control region: Basal Ganglia
# Control region 2: Motor cortex
# Atlas: AAL3
#
# Effective Sample: Anne's Cohort (n=29) and Andy's Cohort (n = 22)
# Outcome: DBS-induced FoG changes (worsening, stay, and improved)
#
# Written by Kristia Pamungkas
# Schweinfurt, May 4th 2021
# Edited: June 9th 2021


#### DATA PREPARATION ####

## DMN ##

# install.packages("dplyr")
library(dplyr)

# Load clinical data from training set
setwd("...") # set working directory

# Clinical Data Cohort 1
clindata_c1 <- read.csv("clindata_c1_v2.csv", sep =";")
str(clindata_c1)
names(clindata_c1)
clindata_c1[1:2] <- lapply(clindata_c1[1:2], as.factor)
clindata_c1[8] <- lapply(clindata_c1[8], as.factor)
str(clindata_c1)
levels(clindata_c1$Outcome)
clindata_c1$Outcome <- factor(clindata_c1$Outcome, levels= c("Worsening","Stay","Improved"))
clindata_c1[clindata_c1 == -999] <- NA
clindata_c1[6] <- lapply(clindata_c1[6], as.factor)

# Clinical Data Cohort 2
clindata_c2<- read.csv("clindata_c2_v2.csv", sep=";")
str(clindata_c2)
names(clindata_c2)
clindata_c2[1:2] <- lapply(clindata_c2[1:2], as.factor)
clindata_c2[8] <- lapply(clindata_c2[8], as.factor)
str(clindata_c2)
levels(clindata_c2$Outcome)
clindata_c2$Outcome <- factor(clindata_c2$Outcome, levels= c("Worsening","Stay","Improved"))
clindata_c2[clindata_c2 == -999] <- NA
clindata_c2[6] <- lapply(clindata_c2[6], as.factor)

# Load resting state functional connectivity data
setwd("...") # set working directory to functional connectivity maps

# Cohort 1
rsfc_c1 <- read.csv("DMN_AvgRSFC_table_C1.csv")
str(rsfc_c1)
colnames(rsfc_c1) <- c("FSM","ACC", "PCC", "AngularGyrus", "Precuneus")
# Join clinical data with DMN resting-state functional connectivity
names(clindata_c1)
cohort1 <- select(clindata_c1, c(ID, Outcome))
names(cohort1)
cohort1[, 3:7] <- rsfc_c1

# Cohort 2
rsfc_c2 <- read.csv("DMN_AvgRSFC_table_C2.csv")
str(rsfc_c2)
colnames(rsfc_c2) <- c("FSM","ACC", "PCC", "AngularGyrus", "Precuneus")
# Join clinical data with DMN resting-state functional connectivity
names(clindata_c2)
cohort2 <- select(clindata_c2, c(ID, Outcome))
names(cohort2)
cohort2[, 3:7] <- rsfc_c2

# Join Cohort 1 and Cohort 2
trainingset <- rbind(cohort1, cohort2)


## Basal Ganglia ##

# Cohort 1
bg_c1 <- read.csv("BG_AvgRSFC_table_C1.csv")
str(bg_c1)
colnames(bg_c1) <- c("Caudate","Putamen", "Pallidum", "SN")
# Join clinical data with BGresting-state functional connectivity
names(clindata_c1)
bg_cohort1 <- dplyr::select(clindata_c1, c(ID, Outcome))
names(bg_cohort1)
bg_cohort1[, 3:6] <- bg_c1

# Cohort 2
bg_c2 <- read.csv("BG_AvgRSFC_table_C2_excl.csv")
str(bg_c2)
colnames(bg_c2) <- c("Caudate","Putamen", "Pallidum", "SN")
# Join clinical data with BG resting-state functional connectivity
names(clindata_c2)
bg_cohort2 <- dplyr::select(clindata_c2, c(ID, Outcome))
bg_cohort2[, 3:6] <- bg_c2

# Join Cohort 1 and Cohort 2
trainingset <- rbind(cohort1, cohort2)

# Join Cohort 1 and Cohort 2
controlset <- rbind(bg_cohort1, bg_cohort2)


## Motor ##

# Cohort 1
m_c1 <- read.csv("Motor_AvgRSFC_table_C1.csv")
str(m_c1)
colnames(m_c1) <- c("M1","SMA")
# Join clinical data with BGresting-state functional connectivity
names(clindata_c1)
m_cohort1 <- dplyr::select(clindata_c1, c(ID, Outcome))
names(m_cohort1)
m_cohort1[, 3:4] <- m_c1

# Cohort 2
m_c2 <- read.csv("Motor_AvgRSFC_table_C2_excl.csv")
str(m_c2)
colnames(m_c2) <- c("M1","SMA")
# Join clinical data with BGresting-state functional connectivity
names(clindata_c2)
m_cohort2 <- dplyr::select(clindata_c2, c(ID, Outcome))
names(m_cohort2)
m_cohort2[, 3:4] <- m_c2

# Join Cohort 1 and Cohort 2
controlset2 <- rbind(m_cohort1, m_cohort2)


# Add total average column
trainingset$DMN <- rowMeans(trainingset[3:7], na.rm = TRUE)
controlset$BG <- rowMeans(controlset[3:6], na.rm = TRUE)
controlset2$Motor <- rowMeans(controlset2[3:4], na.rm = TRUE)

# Standardization
trainingset_st <- trainingset
for(i in 1:nrow(trainingset)){
  for(j in 3:ncol(trainingset)){
    new  = (trainingset[i,j] - mean(trainingset[,j])) / sd(trainingset[,j])
    trainingset_st[i, j] = new
  }
}

controlset_st <- controlset
for(i in 1:nrow(controlset)){
  for(j in 3:ncol(controlset)){
    new  = (controlset[i,j] - mean(controlset[,j])) / sd(controlset[,j])
    controlset_st[i, j] = new
  }
}

controlset2_st <- controlset2
for(i in 1:nrow(controlset2)){
  for(j in 3:ncol(controlset2)){
    new  = (controlset2[i,j] - mean(controlset2[,j])) / sd(controlset2[,j])
    controlset2_st[i, j] = new
  }
}

#### ANALYSIS ####

# Checking the assumptions:
# 1. The dependent variable are ordered
levels(trainingset_st$Outcome) # checked
# > "Worsening" "Stay"      "Improved" 

# 2. One or more of the IV are either continuous, categorical or ordinal
str(trainingset_st) # checked
# > FSM, ACC, PCC, AG, Prec, DMN are continuous variables

# 3. No multicollinearity between IVs # will be checked during model specification

# 4. Proportional odds # will be checked during model specification


# Model Specification
# install.packages("MASS")
library(MASS)

# Model 1: Outcome ~ DMN
# Model 1.1: Outcome ~ FSM
# Model 1.2: Outcome ~ ACC
# Model 1.3: Outcome ~ PCC
# Model 1.4: Outcome ~ AG
# Model 1.5: Outcome ~ Prec
polr.dmn <- polr(Outcome ~ DMN, trainingset_st, Hess = TRUE)
polr.fsm <- polr(Outcome ~ FSM, trainingset_st, Hess = TRUE)
polr.acc <- polr(Outcome ~ ACC, trainingset_st, Hess = TRUE)
polr.pcc <- polr(Outcome ~ PCC, trainingset_st, Hess = TRUE)
polr.ag_ <- polr(Outcome ~ AngularGyrus, trainingset_st, Hess = TRUE)
polr.prec <- polr(Outcome ~ Precuneus, trainingset_st, Hess = TRUE)

# Model 2full: Outcome ~ FSM + ACC + PCC + AngularGyrus + Precuneus
polr.dmnhubsfull <- polr(Outcome ~ FSM + ACC + PCC + AngularGyrus + Precuneus,
                         trainingset_st, Hess = TRUE)
# Check multicollinearity
require(corrplot)
dmnhubsfull_corrmtx <- cor(trainingset_st[3:7], method= c("spearman"), use = "pairwise.complete.obs")
corrplot(dmnhubsfull_corrmtx, method="number")
# AngularGyrus is highly correlated with FSM
# Let's check the VIF to decide which IVs to eliminate

# Check VIF (Variation Inflation Factors)
require(car)
vif(polr.dmnhubsfull)
# FSM          ACC          PCC AngularGyrus    Precuneus 
# 9.637016     4.619440    16.590144     2.354605     2.213031 
# VIF > 5 is problematic; let's try eliminate PCC

# Model 2 minus AG: Outcome ~ FSM + ACC + PCC + Precuneus
polr.dmnhubsfull_2 <- polr(Outcome ~ FSM + ACC + PCC + Precuneus,
                           trainingset_st, Hess = TRUE)

# Check VIF (Variation Inflation Factors)
vif(polr.dmnhubsfull_2)
# FSM       ACC       PCC Precuneus 
# 9.275899  3.187420  2.242235  2.209390 
# Let's try removing ACC

# Model 2 minus AngularGyrus & ACC: Outcome ~ FSM + PCC + Precuneus
polr.dmnhubsfull_3 <- polr(Outcome ~ FSM + PCC + Precuneus,
                           trainingset_st, Hess = TRUE)

# Check VIF (Variation Inflation Factors)
vif(polr.dmnhubsfull_3)
# FSM       PCC Precuneus 
# 3.006146  1.387115  2.198018 
# All VIFs are below 3, thus we keep FSM, PCC, and Prec as DMN predictors


# Model 2: Outcome ~ FSM +  PCC + Precuneus
# Model 2.1: Outcome ~ FSM + PCC
# Model 2.2: Outcome ~ FSM + Precuneus
# Model 2.3: Outcome ~ PCC + Precuneus
polr.dmnhubs <- polr(Outcome ~ FSM + PCC + Precuneus,
                     trainingset_st, Hess = TRUE)
polr.dmnhubs_1 <- polr(Outcome ~ FSM + PCC,
                       trainingset_st, Hess = TRUE)
polr.dmnhubs_2 <- polr(Outcome ~ FSM + Precuneus,
                       trainingset_st, Hess = TRUE)
polr.dmnhubs_3 <- polr(Outcome ~ PCC + Precuneus,
                       trainingset_st, Hess = TRUE)

# Model 3: Outcome ~ FSM *  PCC * Precuneus
# Model 3.1: Outcome ~ FSM * PCC + Precuneus
# Model 3.2: Outcome ~ FSM * Precuneus + PCC
# Model 3.3: Outcome ~ FSM * PCC + Precuneus
polr.dmnhubsint <- polr(Outcome ~ FSM * PCC * Precuneus,
                        trainingset_st, Hess = TRUE)
polr.dmnhubsint_1 <- polr(Outcome ~ FSM * PCC + Precuneus,
                          trainingset_st, Hess = TRUE)
polr.dmnhubsint_2 <- polr(Outcome ~ FSM *Precuneus + PCC,
                          trainingset_st, Hess = TRUE)
polr.dmnhubsint_3 <- polr(Outcome ~ FSM * PCC + Precuneus,
                          trainingset_st, Hess = TRUE)

# Model Evaluation
anova(polr.dmn, polr.fsm)
anova(polr.dmn, polr.pcc)
anova(polr.dmn, polr.prec)
# > dmn is our simplest model

anova(polr.dmnhubs_1, polr.dmnhubs_2)
anova(polr.dmnhubs_1, polr.dmnhubs_3)
anova(polr.dmnhubs_3, polr.dmnhubs)

anova(polr.dmn, polr.dmnhubs_1)
anova(polr.dmn, polr.dmnhubs_2)
anova(polr.dmn, polr.dmnhubs_3)
anova(polr.dmn, polr.dmnhubs)
# > the additive models are not better than DMN model

anova(polr.dmn, polr.dmnhubsint_1)
anova(polr.dmn, polr.dmnhubsint_2)
anova(polr.dmn, polr.dmnhubsint_3)
# > the full additive models with 2 way interaction are not better than DMN model

anova(polr.dmn, polr.dmnhubsint)
# > the most complex model is better than the DMN model (not significant)
# > Thus, the dmn most complex model is the winning DMN model

anova(polr.dmn, polr.dmnhubs, polr.dmnhubsint)
#                    Model Resid. df Resid. Dev   Test    Df  LR stat.     Pr(Chi)
# 1                   DMN        48   86.06356
# 2 FSM + PCC + Precuneus        46   86.03317 1 vs 2     2  0.03038943 0.98492014
# 3 FSM * PCC * Precuneus        42   74.63462 2 vs 3     4 11.39855474 0.02243175

# DMN model summary
summary(polr.dmn)
(ctable <- coef(summary(polr.dmn)))
p <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
(ctable <- cbind(ctable, "p value" = p))
#                     Value Std. Error   t value      p value
# DMN            -0.2867457  0.2965420 -0.9669651 0.3335614920
# Worsening|Stay -1.2952713  0.3433773 -3.7721523 0.0001618454
# Stay|Improved  -0.6867531  0.2995059 -2.2929532 0.0218506994


# Winning model summary
summary(polr.dmnhubsint)
(ctable <- coef(summary(polr.dmnhubsint)))
p <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
(ctable <- cbind(ctable, "p value" = p))
#                        Value Std. Error    t value      p value
# FSM                0.04229756  0.7224629  0.05854634 0.9533134556
# PCC               -0.88722199  0.6892657 -1.28719889 0.1980249753
# Precuneus         -0.16719818  0.5187626 -0.32230190 0.7472239959
# FSM:PCC           -0.75244419  0.4929609 -1.52637708 0.1269159832
# FSM:Precuneus      0.39242444  0.7087300  0.55370090 0.5797835582
# PCC:Precuneus     -0.69253747  0.5799182 -1.19419854 0.2324003130
# FSM:PCC:Precuneus  1.06563317  0.5061903  2.10520253 0.0352736824 *
# Worsening|Stay    -2.05748011  0.5774271 -3.56318616 0.0003663807
# Stay|Improved     -1.30044890  0.5217651 -2.49240302 0.0126881975

# Get confidence intervals for the parameter estimates
ci.polr.dmnhubsint <- confint(polr.dmnhubsint)
# Convert the coefficients into OR
exp(coef(polr.dmnhubsint))
exp(cbind(OR = coef(polr.dmnhubsint), ci.polr.dmnhubsint))
#                          OR      2.5 %   97.5 %
# FSM               1.0432048 0.26948592 5.051756
# PCC               0.4117981 0.08678004 1.474181
# Precuneus         0.8460319 0.26870886 2.192163
# FSM:PCC           0.4712134 0.15956061 1.140697
# FSM:Precuneus     1.4805660 0.34754749 6.275834
# PCC:Precuneus     0.5003049 0.15547482 1.652832
# FSM:PCC:Precuneus 2.9026763 1.26222935 9.630575

# Effectplot
library(effects)
eff.dmnhubsint <- predictorEffects(polr.dmnhubsint, Outcome ~ FSM * PCC * Precuneus)
# Interaction effect plots
plot(eff.dmnhubsint$FSM, style = "stacked")
?plot
plot(eff.dmnhubsint$PCC, style = "stacked")
plot(eff.dmnhubsint$Precuneus, style = "stacked")


## BG ##

# Model Specification (Ordinal Logistic Regression)

# Checking the assumptions:
# 1. The dependent variable are ordered
levels(controlset_st$Outcome) # checked
# > "Worsening" "Stay"      "Improved" 

# 2. One or more of the IV are either continuous, categorical or ordinal
str(controlset_st) # checked
# > Caudate, Putamen, Pallidum, SN, and BG are continuous variables

# 3. No multicollinearity between IVs # will be checked during model specification

# 4. Proportional odds # will be checked during model specification


# Model Specification

# Model 1: Outcome ~ BG
# Model 1.1: Outcome ~ Caudate
# Model 1.2: Outcome ~ Putamen
# Model 1.3: Outcome ~ Pallidum
# Model 1.4: Outcome ~ SN
polr.bg <- polr(Outcome ~ BG, controlset_st, Hess = TRUE)
polr.caudate <- polr(Outcome ~ Caudate, controlset_st, Hess = TRUE)
polr.putamen <- polr(Outcome ~ Putamen, controlset_st, Hess = TRUE)
polr.pallidum <- polr(Outcome ~ Pallidum, controlset_st, Hess = TRUE)
polr.sn <- polr(Outcome ~ SN, controlset_st, Hess = TRUE)

# Model 2full: Outcome ~ Caudate + Putamen + Pallidum + SN
polr.bghubsfull <- polr(Outcome ~ Caudate + Putamen + Pallidum + SN,
                        controlset_st, Hess = TRUE)
# Check multicollinearity
bghubsfull_corrmtx <- cor(controlset_st[3:6], method= c("spearman"), use = "pairwise.complete.obs")
corrplot(bghubsfull_corrmtx, method="number")
# High correlation between Putamen and Pallidum and SN
# Let's check the VIF to decide which IVs to eliminate

# Check VIF (Variation Inflation Factors)
vif(polr.bghubsfull)
# Caudate   Putamen  Pallidum        SN 
# 31.617544 24.152923  5.849510  2.198293 
# Let's try eliminate Putamen

# Model 2 minus Putamen: Outcome ~ Caudate + Pallidum + SN
polr.bghubsfull_2 <- polr(Outcome ~ Caudate + Pallidum + SN,
                          controlset_st, Hess = TRUE)

# Check VIF (Variation Inflation Factors)
vif(polr.bghubsfull_2)
# Caudate Pallidum       SN 
# 4.391505 4.048452 2.199126 
# VIFs below 5 are acceptable
# > So we keep Caudate, Pallidum, and SN as BG predictors

# Model 2: Outcome ~ Caudate + Pallidum + SN
# Model 2.1: Outcome ~ Caudate + Pallidum
# Model 2.2: Outcome ~ Caudate + SN
# Model 2.3: Outcome ~ Pallidum + SN
polr.bghubs <- polr(Outcome ~ Caudate + Pallidum + SN,
                    controlset_st, Hess = TRUE)
polr.bghubs_1 <- polr(Outcome ~ Caudate + Pallidum,
                      controlset_st, Hess = TRUE)
polr.bghubs_2 <- polr(Outcome ~ Caudate + SN,
                      controlset_st, Hess = TRUE)
polr.bghubs_3 <- polr(Outcome ~ Pallidum + SN,
                      controlset_st, Hess = TRUE)

# Model 3: Outcome ~ Caudate * SN * Pallidum
# Model 3.1: Outcome ~ Caudate * Pallidum + SN
# Model 3.2: Outcome ~ Caudate * SN + Pallidum
# Model 3.3: Outcome ~ Caudate + SN * Pallidum
polr.bghubsint <- polr(Outcome ~ Caudate * SN * Pallidum,
                       controlset_st, Hess = TRUE)
polr.bghubsint_1 <- polr(Outcome ~ Caudate * Pallidum + SN,
                         controlset_st, Hess = TRUE)
polr.bghubsint_2 <- polr(Outcome ~ Caudate * SN + Pallidum,
                         controlset_st, Hess = TRUE)
polr.bghubsint_3 <- polr(Outcome ~ Caudate + SN * Pallidum,
                         controlset_st, Hess = TRUE)


# Model Evaluation
anova(polr.bg, polr.caudate)
anova(polr.bg, polr.pallidum)
anova(polr.bg, polr.sn) 
# > BG is our simplest model

anova(polr.bg, polr.bghubs_1)
anova(polr.bg, polr.bghubs_2)
anova(polr.bg, polr.bghubs_3)
anova(polr.bg, polr.bghubs)
# > the additive models are not better than BG model

anova(polr.bg, polr.bghubsint_1)
anova(polr.bg, polr.bghubsint_2)
anova(polr.bg, polr.bghubsint_3)
anova(polr.bg, polr.bghubsint)
# > the full additive models with 2 way interactions and the model with 3 way interactions
# > are not better than the BG model
# > Thus, our BG winning model is the BG model

anova(polr.bg, polr.bghubs, polr.bghubsint)
#                     Model Resid. df Resid. Dev   Test    Df LR stat.    Pr(Chi)
# 1                      BG        48   85.67222                                  
# 2 Caudate + Pallidum + SN        46   85.61549 1 vs 2     2 0.05673266 0.9720322
# 3 Caudate * SN * Pallidum        42   83.08373 2 vs 3     4 2.53175911 0.6389580


# Winning model summary
summary(polr.bg)
(ctable <- coef(summary(polr.bg)))
p <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
(ctable <- cbind(ctable, "p value" = p))
#                     Value Std. Error   t value      p value
# BG             -0.3411672  0.2978663 -1.145370 0.2520557345
# Worsening|Stay -1.3058099  0.3451253 -3.783582 0.0001545874
# Stay|Improved  -0.6932163  0.3009488 -2.303436 0.0212543216

# Get confidence intervals for the parameter estimates
ci.polr.bg <- confint(polr.bg)
# Convert the coefficients into OR
exp(coef(polr.bg))
exp(cbind(OR = coef(polr.bg), ci.polr.bg))
#            OR       2.5 %     97.5 %
# SN  0.71094   0.3887106  1.2725755

# Effectplot
eff.bg <- predictorEffects(polr.bg, Outcome ~ BG)
# Interaction effect plots
plot(eff.bg$BG, style = "stacked")


## MOTOR ##

# Model Specification (Ordinal Logistic Regression)

# Checking the assumptions:
# 1. The dependent variable are ordered
levels(controlset2_st$Outcome) # checked
# > "Worsening" "Stay"      "Improved" 

# 2. One or more of the IV are either continuous, categorical or ordinal
str(controlset2_st) # checked
# > M1, SMA, and M are continuous variables

# 3. No multicollinearity between IVs # will be checked during model specification

# 4. Proportional odds # will be checked during model specification


# Model Specification

# Model 1: Outcome ~ M
# Model 1.1: Outcome ~ M1
# Model 1.2: Outcome ~ SMA
polr.m <- polr(Outcome ~ Motor, controlset2_st, Hess = TRUE)
polr.m1 <- polr(Outcome ~ M1, controlset2_st, Hess = TRUE)
polr.sma <- polr(Outcome ~ SMA, controlset2_st, Hess = TRUE)

# Model 2full: Outcome ~ M1 + SMA
polr.mhubs <- polr(Outcome ~ M1 + SMA,
                   controlset2_st, Hess = TRUE)
# Check multicollinearity
mhubsfull_corrmtx <- cor(controlset2_st[3:4], method= c("spearman"), use = "pairwise.complete.obs")
corrplot(mhubsfull_corrmtx, method="number")
# Okay correlation between M1 and SMA

# Check VIF (Variation Inflation Factors)
vif(polr.mhubs)
# M1      SMA 
# 3.068608 2.220911 
# They are okay


# Model 3: Outcome ~ M1 * SMA
polr.mhubsint <- polr(Outcome ~ M1 * SMA,
                      controlset2_st, Hess = TRUE)


# Model Evaluation
anova(polr.m, polr.m1) # M1 better
anova(polr.m, polr.sma) # SMA is not better
# > let's keep M1 as our simplest model

anova(polr.m1, polr.mhubs)
# > the additive model is not better than M1 model

anova(polr.m1, polr.mhubsint)
# > the most complex model is not better than M1 model

anova(polr.m, polr.mhubs, polr.mhubsint)
#      Model Resid. df Resid. Dev   Test    Df   LR stat.      Pr(Chi)
# 1    Motor        48   86.95701                                  
# 2 M1 + SMA        47   86.94195 1 vs 2     1 0.01506112 0.9023259
# 3 M1 * SMA        46   84.49433 2 vs 3     1 2.44762620 0.1177027

anova(polr.m1, polr.mhubs, polr.mhubsint)
#      Model Resid. df Resid. Dev   Test    Df   LR stat.      Pr(Chi)
# 1       M1        48   86.94585                                   
# 2 M1 + SMA        47   86.94195 1 vs 2     1 0.003896288 0.9502282
# 3 M1 * SMA        46   84.49433 2 vs 3     1 2.447626201 0.1177027


# Winning model summary
summary(polr.m1)
(ctable <- coef(summary(polr.m1)))
p <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
(ctable <- cbind(ctable, "p value" = p))
#                   Value Std. Error    t value     p value
# M1              0.06932741  0.3024843  0.2291934 0.818718573
# Worsening|Stay -1.29125592  0.3406253 -3.7908395 0.000150139
# Stay|Improved  -0.69287411  0.2971911 -2.3314094 0.019731780

# Get confidence intervals for the parameter estimates
ci.polr.m1 <- confint(polr.m1)
# Convert the coefficients into OR
exp(coef(polr.m1))
exp(cbind(OR = coef(polr.m1), ci.polr.m1))
#              OR      2.5 %     97.5 %
# M1     1.071787 0.6006824 2.0048932

# Effectplot
eff.m1 <- predictorEffects(polr.m1, Outcome ~ M1)
# Interaction effect plots
plot(eff.m1, style = "stacked")


#### COMPARE ALL WINNING MODELS AGAINST EACH OTHER ####
anova(polr.bg, polr.m1, polr.dmnhubsint)
#                   Model Resid. df Resid. Dev   Test    Df  LR stat.   Pr(Chi)
# 1                    BG        48   85.67222                                  
# 2                    M1        48   86.94585 1 vs 2     0 -1.273624 1.00000000
# 3 FSM * PCC * Precuneus        42   74.63733 2 vs 3     6 12.308514 0.05542968

anova(polr.bg, polr.dmnhubsint)
#                   Model Resid. df Resid. Dev   Test    Df LR stat.    Pr(Chi)
# 1                    BG        48   85.67222                                 
# 2 FSM * PCC * Precuneus        42   74.63733 1 vs 2     6 11.03489 0.08730406

## DMN is non-significantly better than both SN and M1 models 


#### MODEL PERFORMANCE ####

# install.packages("MLmetrics")
library(MLmetrics)

## DMN ##

# Confusion Matrix
trainingset_2 <- dplyr::select(trainingset_st, -ACC, -AngularGyrus)
dmn.AP <- predict(polr.dmnhubsint, trainingset_2[, 3:5])
table(trainingset_2[, 2], dmn.AP)
#                  Predicted
#             Worsening Stay Improved
# Worsening         3    0        8
# Stay              1    0        5
# Improved          1    0       33

dmn.APprob <- predict(polr.dmnhubsint, trainingset_2[, 3:5], type = "probs")
MultiLogLoss(y_pred = dmn.APprob , y_true = trainingset_2$Outcome)
# 1.909761
Accuracy(dmn.AP, trainingset_2$Outcome)
# 71 %
F1_Score(dmn.AP, trainingset_2$Outcome)
# 38 %

# Balanced Accuracy
# ((3/11) + (0/6) + (33/34)/3) = 60 %


## BG ##

# Confusion Matrix
controlset_2 <- dplyr::select(controlset_st, -Putamen)
bg.AP <- predict(polr.bg, controlset_2[6])
table(controlset_2[, 2], bg.AP)
#                  Predicted
#            Worsening Stay Improved
# Worsening         0    0       11
# Stay              0    0        6
# Improved          0    0       34

bg.APprob <- predict(polr.bg, controlset_2[6], type = "probs")
MultiLogLoss(y_pred = bg.APprob , y_true = controlset_2$Outcome)
# 1.411606
Accuracy(bg.AP, controlset_2$Outcome)
# 67 %
F1_Score(bg.AP, controlset_2$Outcome)
# NaN %

# Balanced Accuracy
# ((0/11) + (0/6) + (34/34))/3 = 33 %


## Motor ##

# Confusion Matrix
m.AP <- predict(polr.m1, controlset2_st[3])
table(controlset2[, 2], m.AP)
#                  Predicted
#             Worsening Stay Improved
# Worsening         0    0       11
# Stay              0    0        6
# Improved          0    0       34

m.APprob <- predict(polr.m1, controlset2_st[3], type = "probs")
MultiLogLoss(y_pred = m.APprob , y_true = controlset2_st$Outcome)
# 1.363949
Accuracy(m.AP, controlset2_st$Outcome)
# 67 %
F1_Score(m.AP, controlset2_st$Outcome)
# NaN

# Balanced Accuracy
# ((0/11) + (0/6) + (34/34))/3 = 33 %


#### K-FOLD CROSSVALIDATION ####
# install.packages("caret")
library(caret)

## DMN ##

dmn.int <- model.matrix(lm(Outcome ~ -1 + FSM * PCC * Precuneus, data= trainingset_2)) 
colnames(dmn.int)

# Ordered Logistic Regression
set.seed(123)
control <- trainControl(method="cv",
                        number=5,
                        classProbs= TRUE,
                        summaryFunction = multiClassSummary,
                        savePredictions = TRUE)
metric <- "Accuracy"

set.seed(123)
mod.dmn_simplest <-train(Outcome ~ DMN,
                         trainingset_2,
                         method="polr",
                         metric=metric, 
                         trControl=control,
                         tuneGrid = expand.grid(method = "logistic"))
mod.dmn_simplest
summary(mod.dmn_simplest)
# Residual Deviance: 86.06422 
# AIC: 92.06422 

set.seed(123)
mod.dmn_additive <-train(Outcome ~ FSM + PCC + Precuneus,
                         trainingset_2,
                         method="polr",
                         metric=metric, 
                         trControl=control,
                         tuneGrid = expand.grid(method = "logistic"))
mod.dmn_additive
summary(mod.dmn_additive)
# Residual Deviance: 86.03354 
# AIC: 96.03354 

set.seed(123)
mod.dmn_additive2 <-train(Outcome ~ PCC + Precuneus,
                           trainingset_2,
                           method="polr",
                           metric=metric, 
                           trControl=control,
                           tuneGrid = expand.grid(method = "logistic"))
mod.dmn_additive2
summary(mod.dmn_additive2)
# Residual Deviance: 86.06558 
# AIC: 94.06558 

set.seed(123)
mod.dmn <- train(y = trainingset_2$Outcome,
                 x = dmn.int,
                 method="polr",
                 metric=metric, 
                 trControl=control,
                 tuneGrid = expand.grid(method = "logistic"))

mod.dmn
# logLoss    AUC        prAUC      Accuracy   Kappa      Mean_F1  Mean_Sensitivity  Mean_Specificity  Mean_Pos_Pred_Value
# 0.9802494  0.6531614  0.3748226  0.6682828  0.1445853  NaN      0.3825397         0.7111111         NaN     

# Mean_Neg_Pred_Value  Mean_Precision  Mean_Recall  Mean_Detection_Rate  Mean_Balanced_Accuracy
# 0.795202             NaN             0.3825397    0.2227609            0.5468254    

summary(mod.dmn)
# Residual Deviance: 74.63733 
# AIC: 92.63733 

(ctable <- coef(summary(mod.dmn)))
p <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
(ctable <- round(cbind(ctable, "p value" = p),3))
#                          Value Std. Error    t value      p value
# FSM                  0.04229756  0.7224629  0.05854634 0.9533134556
# PCC                 -0.88722199  0.6892657 -1.28719889 0.1980249753
# Precuneus           -0.16719818  0.5187626 -0.32230190 0.7472239959
# `FSM:PCC`           -0.75244419  0.4929609 -1.52637708 0.1269159832
# `FSM:Precuneus`      0.39242444  0.7087300  0.55370090 0.5797835582
# `PCC:Precuneus`     -0.69253747  0.5799182 -1.19419854 0.2324003130
# `FSM:PCC:Precuneus`  1.06563317  0.5061903  2.10520253 0.0352736824 *
# Worsening|Stay      -2.05748011  0.5774271 -3.56318616 0.0003663807
# Stay|Improved       -1.30044890  0.5217651 -2.49240302 0.0126881975

# Convert the coefficients into OR
pred <- coef(summary(mod.dmn))[, 1]
or <- exp(pred)
or

# FSM                 PCC           Precuneus           `FSM:PCC`     `FSM:Precuneus`     `PCC:Precuneus` 
# 1.0432048           0.4117981           0.8460319           0.4712134           1.4805660           0.5003049 

# `FSM:PCC:Precuneus`      Worsening|Stay       Stay|Improved 
# 2.9026763           0.1277755           0.2724095 

# Transform dataset
h1.dmn.kfcv <- dplyr::select(trainingset_2, Outcome)
h1.dmn.kfcv[, 2:8] <- dmn.int
names(h1.dmn.kfcv)
colnames(h1.dmn.kfcv) <- c("Outcome", "FSM", "PCC", "Precuneus",
                           "FSM:PCC","FSM:Precuneus","PCC:Precuneus",
                           "FSM:PCC:Precuneus")

# Confusion Matrix
dmn.kfcv.AP <- predict(mod.dmn, h1.dmn.kfcv[, 2:8])
dmn.kfcv.AP.prob <- predict(mod.dmn, h1.dmn.kfcv[, 2:8], type = "prob")
table(h1.dmn.kfcv[, 1], dmn.kfcv.AP)
#                  Predicted
#           Worsening Stay Improved
# Worsening         3    0        8
# Stay              1    0        5
# Improved          1    0       33

Accuracy(dmn.kfcv.AP, h1.dmn.kfcv$Outcome) # 71 %
MultiLogLoss(dmn.kfcv.AP.prob, h1.dmn.kfcv$Outcome) # 1.909761
confusionMatrix(dmn.kfcv.AP, h1.dmn.kfcv$Outcome)

# Features importance
dmn.imp <- varImp(mod.dmn, scale=FALSE)
print(dmn.imp)
plot(dmn.imp)


## BG ##

# Ordered Logistic Regression
set.seed(123)
control <- trainControl(method="cv",
                        number=5,
                        classProbs= TRUE,
                        summaryFunction = multiClassSummary,
                        savePredictions = TRUE)
metric <- "Accuracy"

set.seed(123)
mod.bg <- train(Outcome ~ BG,
                controlset_2,
                method="polr",
                metric=metric, 
                trControl=control,
                tuneGrid = expand.grid(method = "logistic"))
mod.bg
# logLoss    AUC        prAUC     Accuracy   Kappa        Mean_F1  Mean_Sensitivity  Mean_Specificity  Mean_Pos_Pred_Value
# 0.8751385  0.4080026  0.273078  0.6478788  -0.02857143  NaN      0.3238095         0.6583333         NaN           

# Mean_Neg_Pred_Value  Mean_Precision  Mean_Recall  Mean_Detection_Rate  Mean_Balanced_Accuracy
# 0.5592593            NaN             0.3238095    0.2159596            0.4910714   

summary(mod.bg)
# Residual Deviance: 85.67222 
# AIC: 91.67222 

(ctable <- coef(summary(mod.bg)))
p <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
(ctable <- cbind(ctable, "p value" = p))
#                     Value Std. Error   t value      p value
# BG             -0.3411672  0.2978663 -1.145370 0.2520557345
# Worsening|Stay -1.3058099  0.3451253 -3.783582 0.0001545874
# Stay|Improved  -0.6932163  0.3009488 -2.303436 0.0212543216

# Confusion Matrix
bg.kfcv.AP <- predict(mod.bg, controlset_2[6])
bg.kfcv.AP.prob <- predict(mod.bg, controlset_2[6], type = "prob")
table(controlset_2[, 2], bg.kfcv.AP)
#                  Predicted
#             Worsening Stay Improved
# Worsening         0    0       11
# Stay              0    0        6
# Improved          0    0       34

Accuracy(bg.kfcv.AP, controlset_2$Outcome) # 67 %
MultiLogLoss(bg.kfcv.AP.prob, controlset_2$Outcome) # 1.411606
confusionMatrix(bg.kfcv.AP, controlset_2$Outcome)

## Motor ##

# Ordered Logistic Regression
set.seed(123)
control <- trainControl(method="cv",
                        number=5,
                        classProbs= TRUE,
                        summaryFunction = multiClassSummary,
                        savePredictions = TRUE)
metric <- "Accuracy"

set.seed(123)
mod.m <- train(Outcome ~ M1,
               controlset2_st,
               method="polr",
               metric=metric, 
               trControl=control,
               tuneGrid = expand.grid(method = "logistic"))
mod.m
# logLoss    AUC        prAUC      Accuracy   Kappa  Mean_F1  Mean_Sensitivity  Mean_Specificity  Mean_Pos_Pred_Value
# 0.8884509  0.4556878  0.3009537  0.6678788  0      NaN      0.3333333         0.6666667         NaN     

# Mean_Neg_Pred_Value  Mean_Precision  Mean_Recall  Mean_Detection_Rate  Mean_Balanced_Accuracy
# NaN                  NaN             0.3333333    0.2226263            0.5                 

summary(mod.m)
# Residual Deviance: 86.94585 
# AIC: 92.94585 
(ctable <- coef(summary(mod.m)))
p <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
(ctable <- cbind(ctable, "p value" = p))
#                     Value Std. Error    t value      p value
# M1              0.06932741  0.3024843  0.2291934 0.818718573
# Worsening|Stay -1.29125592  0.3406253 -3.7908395 0.000150139
# Stay|Improved  -0.69287411  0.2971911 -2.3314094 0.019731780


# Confusion Matrix
m.kfcv.AP <- predict(mod.m, controlset2_st[3])
m.kfcv.AP.prob <- predict(mod.m, controlset2_st[3], type = "prob")
table(controlset2_st[, 2], m.kfcv.AP)
#                  Predicted
#             Worsening Stay Improved
# Worsening         0    0       11
# Stay              0    0        6
# Improved          0    0       34

Accuracy(m.kfcv.AP, controlset2_st$Outcome) # 67 %
MultiLogLoss(m.kfcv.AP.prob, controlset2_st$Outcome) # 1.363949
confusionMatrix(m.kfcv.AP, controlset2_st$Outcome)



############################################################
#### Ordinal Logistic Regression Final Version - Part 2 ####
############################################################
#
# Hypothesis 2: DMN connectivity associated with DBS-induced
# FoG changes would be predictive of 
# DBS-induced changes in executive dysfunction in PD patients with STN-DBS
#
# Classification of DBS-induced executive dysfunction changes
# based on VTA-based Average Resting-State Functional Connectivity
# of Default Mode Network
# Control region: Basal Ganglia
# Control region 2: Motor cortex
# Atlas: AAL3
#
# Sample: Frankenmolle's Cohort (n=10 @ 2 types of stimulation)
# Outcome: DBS-induced ED changes (worsening, stay, and improved)


#### DATA PREPARATION ####

setwd("...") # set working directory to clinical data of cohort 3

# Clinical Data Cohort 3
clindata_c3 <- read.csv("clindata_c3.csv", sep =";")
str(clindata_c3)
names(clindata_c3)
clindata_c3[1:2] <- lapply(clindata_c3[1:2], as.factor)
clindata_c3[8:9] <- lapply(clindata_c3[8:9], as.factor)
clindata_c3[13:14] <- lapply(clindata_c3[13:14], as.factor)
str(clindata_c3)
levels(clindata_c3$Single_2back_ClinStimEffect)
levels(clindata_c3$Single_2back_ClinStimEffect) <- c("Worsening", "Stay", "Improved")
levels(clindata_c3$Single_2back_ModStimEffect) <- c("Worsening", "Stay", "Improved")
levels(clindata_c3$Dual_2back_ClinStimEffect) <- c("Worsening", "Stay", "Improved")
levels(clindata_c3$Dual_2back_ModStimEffect) <- c("Worsening", "Stay", "Improved")
clindata_c3$DualDiff_Clin <- clindata_c3$Dual_2back_Clin - clindata_c3$Dual_2back_Off
clindata_c3$DualDiff_Mod <- clindata_c3$Dual_2back_Mod - clindata_c3$Dual_2back_Off


# Load resting state functional connectivity data
setwd("...") # set working directory to cohort 3 functional connectivity data
rsfc_c3clin <- read.csv("DMN_AvgRSFC_table_C3clin.csv")
str(rsfc_c3clin)
colnames(rsfc_c3clin) <- c("FSM","ACC", "PCC", "AngularGyrus", "Precuneus")

bg_c3clin <- read.csv("BG_AvgRSFC_table_C3clin.csv")
str(bg_c3clin)
colnames(bg_c3clin) <- c("Caudate","Putamen", "Pallidum", "SN")

motor_c3clin <- read.csv("Motor_AvgRSFC_table_C3clin.csv")
str(motor_c3clin)
colnames(motor_c3clin) <- c("M1","SMA")

# Join clinical data with DMN resting-state functional connectivity
names(clindata_c3)
cohort3clin <- dplyr::select(clindata_c3, c(ID, Dual_2back_ClinStimEffect, DualDiff_Clin))
names(cohort3clin)
cohort3clin[, 4:8] <- rsfc_c3clin
cohort3clin[, 9:12] <- bg_c3clin
cohort3clin[, 13:14] <- motor_c3clin
# Add 1 column of avg DMN and BG RSFC
names(cohort3clin)
cohort3clin$DMN <- rowMeans(cohort3clin[4:8], na.rm = TRUE)
cohort3clin$BG <- rowMeans(cohort3clin[9:12], na.rm = TRUE)
cohort3clin$Motor <- rowMeans(cohort3clin[13:14], na.rm = TRUE)
cohort3clin$Stim <- "clinical"
colnames(cohort3clin)[2:3] <- c("Dual_Outcome_1", "DualDiff") 


rsfc_c3mod <- read.csv("DMN_AvgRSFC_table_C3mod.csv")
str(rsfc_c3mod)
colnames(rsfc_c3mod) <- c("FSM","ACC", "PCC", "AngularGyrus", "Precuneus")

bg_c3mod <- read.csv("BG_AvgRSFC_table_C3mod.csv")
str(bg_c3mod)
colnames(bg_c3mod) <- c("Caudate","Putamen", "Pallidum", "SN")

motor_c3mod <- read.csv("Motor_AvgRSFC_table_C3mod.csv")
str(motor_c3mod)
colnames(motor_c3mod) <- c("M1","SMA")

# Join clinical data with DMN resting-state functional connectivity
names(clindata_c3)
cohort3mod <- dplyr::select(clindata_c3, c(ID, Dual_2back_ModStimEffect, DualDiff_Mod))
names(cohort3mod)
cohort3mod[, 4:8] <- rsfc_c3mod
cohort3mod[, 9:12] <- bg_c3mod
cohort3mod[, 13:14] <- motor_c3mod
names(cohort3mod)

# Add 1 column of avg DMN, BG, and Motor RSFC

cohort3mod$DMN <- rowMeans(cohort3mod[4:8], na.rm = TRUE)
cohort3mod$BG <- rowMeans(cohort3mod[9:12], na.rm = TRUE)
cohort3mod$Motor <- rowMeans(cohort3mod[13:14], na.rm = TRUE)
cohort3mod$Stim <- "model"
colnames(cohort3mod)[2:3] <- c("Dual_Outcome_1", "DualDiff") 

# Join cohort 3
cohort3 <- rbind(cohort3clin, cohort3mod)
cohort3 <- dplyr::select(cohort3, ID, Stim, Dual_Outcome_1, DualDiff,
                         FSM, ACC, PCC, AngularGyrus, Precuneus, DMN,
                         Caudate, Putamen, Pallidum, SN, BG,
                         M1, SMA, Motor)
str(cohort3)
cohort3$Stim <- as.factor(cohort3$Stim)

# Define Dual_Outcome_2 (No change group = performance difference <= +- 1 SD)
sd(cohort3$DualDiff)
names(cohort3)
cohort3$Dual_Outcome_2 <- cohort3$DualDiff - sd(cohort3$DualDiff)
for(i in 1:nrow(cohort3)){
  a = 0
  if((cohort3[i,4] <= 5.08148) & (cohort3[i, 4] >= -5.08148)){
    cohort3[i,19] = a
  } else if(cohort3[i,4] > 5.08148){
    a = 1
    cohort3[i, 19] = a
  } else{
    a = -1
    cohort3[i,19] = a
  }
}


# Subset to DMN, BG, Motor
names(cohort3)
cohort3.dmn <- dplyr::select(cohort3, ID, Stim, Dual_Outcome_1, Dual_Outcome_2,
                             FSM, ACC, PCC, AngularGyrus, Precuneus, DMN)
cohort3.bg <- dplyr::select(cohort3, ID, Stim, Dual_Outcome_1, Dual_Outcome_2,
                            Caudate, Putamen, Pallidum, SN, BG)
cohort3.m <- dplyr::select(cohort3, ID, Stim, Dual_Outcome_1, Dual_Outcome_2,
                           M1, SMA, Motor)


# Standardization
cohort3.dmn.st <- cohort3.dmn
for(i in 1:nrow(cohort3.dmn)){
  for(j in 5:ncol(cohort3.dmn)){
    new  = (cohort3.dmn[i,j] - mean(trainingset[,j-2])) / sd(trainingset[,j-2])
    cohort3.dmn.st[i, j] = new
  }
}

cohort3.bg.st <- cohort3.bg
for(i in 1:nrow(cohort3.bg)){
  for(j in 5:ncol(cohort3.bg)){
    new  = (cohort3.bg[i,j] - mean(controlset[,j-2])) / sd(controlset[,j-2])
    cohort3.bg.st[i, j] = new
  }
}

cohort3.m.st <- cohort3.m
for(i in 1:nrow(cohort3.m)){
  for(j in 5:ncol(cohort3.m)){
    new  = (cohort3.m[i,j] - mean(controlset2[,j-2])) / sd(controlset2[,j-2])
    cohort3.m.st[i, j] = new
  }
}


# Keep predictors of interest only
cohort3.dmn.st <- dplyr::select(cohort3.dmn.st, -ACC, -AngularGyrus, - DMN)
cohort3.bg.st <- dplyr::select(cohort3.bg.st, -Caudate, -Putamen, -Pallidum, -SN)
cohort3.m.st <- dplyr::select(cohort3.m.st, -SMA, -Motor)


#### DUAL TASK PREDICTION USING CROSSVALIDATED CLASSIFIERS ####
# install.packages("MLmetrics")
library("MLmetrics")

## Dual Outcome 1: +- 1% changes are grouped to stay group ##

# Create prediction matrices
dmn_pred.transf <- model.matrix(lm(Dual_Outcome_1~ -1 + FSM * PCC * Precuneus, data= cohort3.dmn.st)) 
colnames(dmn_pred.transf)
cohort3.dmn_transf <- dplyr::select(cohort3.dmn.st, Dual_Outcome_1)
cohort3.dmn_transf[,2:8] <- dmn_pred.transf
colnames(cohort3.dmn_transf)[2:8]<- c("FSM", "PCC", "Precuneus",
                                      "FSM:PCC", "FSM:Precuneus", "PCC:Precuneus",
                                      "FSM:PCC:Precuneus")

h2_1.dmn_pred <- predict(mod.dmn, cohort3.dmn_transf[, 2:8])
h2_1.dmn.prob <- predict(mod.dmn, cohort3.dmn_transf[, 2:8], type = "prob")

h2_1.bg_pred <- predict(mod.bg , cohort3.bg.st[5])
h2_1.bg.prob <- predict(mod.bg, cohort3.bg.st[5], type = "prob")

h2_1.m_pred <- predict(mod.m, cohort3.m.st[5])
h2_1.m.prob <- predict(mod.m, cohort3.m.st[5], type = "prob")

# DMN
# Confusion Matrix
table(cohort3.dmn_transf[, 1], h2_1.dmn_pred)
#                    dmn_pred
#           Worsening Stay Improved
# Worsening         0    0        6
# Stay              0    0        6
# Improved          0    0        8

MultiLogLoss(h2_1.dmn.prob , cohort3.dmn_transf$Dual_Outcome_1)
# 1.470226
Accuracy(h2_1.dmn_pred, cohort3.dmn_transf$Dual_Outcome_1)
# 40 %
F1_Score(h2_1.dmn_pred, cohort3.dmn_transf$Dual_Outcome_1)
# NaN %

## Balanced Accuracy
## (0/6 + 0/6 + 8/8)/3 = 33 %


# BG
# Confusion Matrix
table(cohort3.bg.st[, 3], h2_1.bg_pred)
#                    bg_pred
#             Worsening Stay Improved
# Worsening         0    0        6
# Stay              0    0        6
# Improved          0    0        8

MultiLogLoss(h2_1.bg.prob , cohort3.bg.st$Dual_Outcome_1)
# 1.529686
Accuracy(h2_1.bg_pred, cohort3.bg.st$Dual_Outcome_1)
# 40 %
F1_Score(h2_1.bg_pred, cohort3.bg.st$Dual_Outcome_1)
# NaN

## Balanced Accuracy
## (0/6 + 0/6 + 8/8)/3 = 33 %


# Motor
# Confusion Matrix
table(cohort3.m.st[, 3], h2_1.m_pred)
#                    m_pred
#             Worsening Stay Improved
# Worsening         0    0        6
# Stay              0    0        6
# Improved          0    0        8

MultiLogLoss(h2_1.m.prob , cohort3.m.st$Dual_Outcome_1)
# 1.42175
Accuracy(h2_1.m_pred, cohort3.m.st$Dual_Outcome_1)
# 40 %
F1_Score(h2_1.m_pred, cohort3.m.st$Dual_Outcome_1)
# NaN

## Balanced Accuracy
## (0/6 + 0/6 + 8/8)/3 = 33 %


## Dual Outcome 2: +- 1 SD (around 5%) changes are grouped to stay group ##

# install.packages("plyr")
library(plyr)
cohort3.dmn.st$Dual_Outcome_2 <- as.factor(cohort3.dmn.st$Dual_Outcome_2)
cohort3.dmn.st$Dual_Outcome_2 <- revalue(cohort3.dmn.st$Dual_Outcome_2, c("-1"="Worsening", "0"="Stay", "1" = "Improved"))

cohort3.bg.st$Dual_Outcome_2 <- as.factor(cohort3.bg.st$Dual_Outcome_2)
cohort3.bg.st$Dual_Outcome_2 <- revalue(cohort3.bg.st$Dual_Outcome_2, c("-1"="Worsening", "0"="Stay", "1" = "Improved"))

cohort3.m.st$Dual_Outcome_2 <- as.factor(cohort3.m.st$Dual_Outcome_2)
cohort3.m.st$Dual_Outcome_2 <- revalue(cohort3.m.st$Dual_Outcome_2, c("-1"="Worsening", "0"="Stay", "1" = "Improved"))


# Create prediction matrices
dmn_pred.transf2 <- model.matrix(lm(Dual_Outcome_2~ -1 + FSM * PCC * Precuneus, data= cohort3.dmn.st)) 
colnames(dmn_pred.transf2)
cohort3.dmn_transf2 <- dplyr::select(cohort3.dmn.st, Dual_Outcome_2)
cohort3.dmn_transf2[,2:8] <- dmn_pred.transf2
colnames(cohort3.dmn_transf2)[2:8]<- c("FSM", "PCC", "Precuneus",
                                       "FSM:PCC", "FSM:Precuneus", "PCC:Precuneus",
                                       "FSM:PCC:Precuneus")

h2_2.dmn_pred <- predict(mod.dmn, cohort3.dmn_transf2[, 2:8])
h2_2.dmn.prob <- predict(mod.dmn, cohort3.dmn_transf2[, 2:8], type = "prob")

h2_2.bg_pred <- predict(mod.bg , cohort3.bg.st[5])
h2_2.bg.prob <- predict(mod.bg, cohort3.bg.st[5], type = "prob")

h2_2.m_pred <- predict(mod.m, cohort3.m.st[5])
h2_2.m.prob <- predict(mod.m, cohort3.m.st[5], type = "prob")

# DMN
# Confusion Matrix
table(cohort3.dmn_transf2[, 1], h2_2.dmn_pred)
#                    dmn_pred
#           Worsening Stay Improved
# Worsening         0    0        1
# Stay              0    0       16
# Improved          0    0        3

MultiLogLoss(h2_2.dmn.prob , cohort3.dmn_transf2$Dual_Outcome_2)
# 1.984616
Accuracy(h2_2.dmn_pred, cohort3.dmn_transf2$Dual_Outcome_2)
# 15 %
F1_Score(h2_2.dmn_pred, cohort3.dmn_transf2$Dual_Outcome_2)
# NaN

## Balanced Accuracy
## (0/1 + 0/16 + 0/3)/3 = 33 %


# BG
# Confusion Matrix
table(cohort3.bg.st[, 4], h2_2.bg_pred)
#                    bg_pred
#             Worsening Stay Improved
# Worsening         0    0        1
# Stay              0    0       16
# Improved          0    0        3

MultiLogLoss(h2_2.bg.prob , cohort3.bg.st$Dual_Outcome_2)
# 2.119712
Accuracy(h2_2.bg_pred, cohort3.bg.st$Dual_Outcome_2)
# 15 %
F1_Score(h2_2.bg_pred, cohort3.bg.st$Dual_Outcome_2)
# NaN

## Balanced Accuracy
## (0/1 + 0/16 + 3/3)/3 = 33 %


# Motor
# Confusion Matrix
table(cohort3.m.st[, 4], h2_2.m_pred)
#                    m_pred
#             Worsening Stay Improved
# Worsening         0    0        1
# Stay              0    0       16
# Improved          0    0        3

MultiLogLoss(h2_2.m.prob , cohort3.m.st$Dual_Outcome_2)
# 2.00794
Accuracy(h2_2.m_pred, cohort3.m.st$Dual_Outcome_2)
# 15 %
F1_Score(h2_2.m_pred, cohort3.m.st$Dual_Outcome_2)
# NaN

## Balanced Accuracy
## (0/1 + 0/16 + 3/3)/3 = 33 %
