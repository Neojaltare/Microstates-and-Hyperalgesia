library(pbkrtest)
library(psych)
library(sjPlot)
library(readxl)
library(reghelper)
library(merTools)
library(lmerTest)
library(lme4)
library(ggeffects)
library(emmeans)
library(effects)
library(dbplyr)
library(readxl)
library(parameters)
library(ggplot2)
library(multcomp)
library(effectsize)
library(lattice)
library(plotly)
library(tidyverse)
library(papaja)
library(gt)
library(vtable)
library(plotrix)
library(car)
library(BayesFactor)
library(writexl)

Path = "Path_to_Microstate_results/StabilityAnalysis/"
runs = c("Together1","Together2","Together3","Together4","Together5")

thisrun = 1

fullpath = paste0(Path,runs[thisrun],"/Microstats.xlsx")
fullpath

# Load Single Clustering data ----
Microstats <- read_excel(fullpath)
# View(Microstats)
Microstats$Subject<- as.factor(Microstats$Subject)

# Run1 - A = 5, B = 2, C = 3, D = 1, E = 4
# Run2 - A = 5, B = 3, C = 2, D = 1, E = 4
# Run3 - A = 5, B = 2, C = 3, D = 1, E = 4
# Run4 - A = 3, B = 4, C = 2, D = 1, E = 5
# Run5 - A = 5, B = 1, C = 4, D = 2, E = 3


MS1 = Microstats[Microstats$MS == 'MS1',]
MS2 = Microstats[Microstats$MS == 'MS2',]
MS3 = Microstats[Microstats$MS == 'MS3',]
MS4 = Microstats[Microstats$MS == 'MS4',]
MS5 = Microstats[Microstats$MS == 'MS5',]

TotalGEVmod <- lm(TotalGEV ~ RS, data = MS1)
summary(TotalGEVmod)

mean(MS2$TotalGEV[MS2$RS == 'RS2'])
sd(MS2$TotalGEV[MS2$RS == 'RS2'])

# Create column and row names
col_names <- c("Duration", "Occurrence", "Coverage", "GEV", "TP1", "TP2", "TP3", "TP4", "TP5")
row_names <- c("MS1", "MS2", "MS3", "MS4", "MS5")
# Create empty data frame
ALL_Pvals <- data.frame(matrix(ncol = length(col_names), nrow = length(row_names)))
colnames(ALL_Pvals) <- col_names
rownames(ALL_Pvals) <- row_names

ALL_BFvals <- data.frame(matrix(ncol = length(col_names), nrow = length(row_names)))
colnames(ALL_BFvals) <- col_names
rownames(ALL_BFvals) <- row_names


# MS1
Durmod <- lmer(Duration ~ RS + (1|Subject), data = MS1,REML = FALSE)
summary(Durmod)

full_Durmod = lmBF(Duration ~  RS + Subject, data=MS1, whichRandom = 'Subject')
null_Durmod = lmBF(Duration ~ 1 + Subject, data = MS1, whichRandom = 'Subject')
BF_Duration = full_Durmod/null_Durmod

Occurmod <- lmer((Occurance) ~ RS + (1|Subject), data = MS1,REML = FALSE)
summary(Occurmod)

full_Occurmod = lmBF(Occurance ~  RS + Subject, data=MS1, whichRandom = 'Subject')
null_Occurmod = lmBF(Occurance ~ 1 + Subject, data = MS1, whichRandom = 'Subject')
BF_Occurance = full_Occurmod/null_Occurmod

Covermod <- lmer((Coverage) ~ RS + (1|Subject), data = MS1,REML = FALSE)
summary(Covermod)

full_Covermod = lmBF(Coverage ~  RS + Subject, data=MS1, whichRandom = 'Subject')
null_Covermod = lmBF(Coverage ~ 1 + Subject, data = MS1, whichRandom = 'Subject')
BF_Coverage = full_Covermod/null_Covermod

Gevmod <- lmer(GEV ~ RS + (1|Subject), data = MS1,REML = FALSE)
summary(Gevmod)

full_GEV = lmBF(GEV ~  RS + Subject, data=MS1, whichRandom = 'Subject')
null_GEV = lmBF(GEV ~ 1 + Subject, data = MS1, whichRandom = 'Subject')
BF_GEV = full_GEV/null_GEV

TP2mod <- lmer(TP2 ~ RS + (1|Subject), data = MS1,REML = FALSE)
summary(TP2mod)

full_TP2mod = lmBF(TP2 ~  RS + Subject, data=MS1, whichRandom = 'Subject')
null_TP2mod = lmBF(TP2 ~ 1 + Subject, data = MS1, whichRandom = 'Subject')
BF_TP2mod = full_TP2mod/null_TP2mod

TP3mod <- lmer(TP3 ~ RS + (1|Subject), data = MS1,REML = FALSE)
summary(TP3mod)

full_TP3mod = lmBF(TP3 ~  RS + Subject, data=MS1, whichRandom = 'Subject')
null_TP3mod = lmBF(TP3 ~ 1 + Subject, data = MS1, whichRandom = 'Subject')
BF_TP3mod = full_TP3mod/null_TP3mod

TP4mod <- lmer((TP4) ~ RS + (1|Subject), data = MS1,REML = FALSE)
summary(TP4mod)

full_TP4mod = lmBF(TP4 ~  RS + Subject, data=MS1, whichRandom = 'Subject')
null_TP4mod = lmBF(TP4 ~ 1 + Subject, data = MS1, whichRandom = 'Subject')
BF_TP4mod = full_TP4mod/null_TP4mod

TP5mod <- lmer((TP5) ~ RS + (1|Subject), data = MS1,REML = FALSE)
summary(TP5mod)

full_TP5mod = lmBF(TP5 ~  RS + Subject, data=MS1, whichRandom = 'Subject')
null_TP5mod = lmBF(TP5 ~ 1 + Subject, data = MS1, whichRandom = 'Subject')
BF_TP5mod = full_TP5mod/null_TP5mod

ALL_Pvals$Duration[1] = anova(Durmod)$`Pr(>F)`[1]
ALL_Pvals$Occurrence[1] = anova(Occurmod)$`Pr(>F)`[1]
ALL_Pvals$Coverage[1] = anova(Covermod)$`Pr(>F)`[1]
ALL_Pvals$GEV[1] = anova(Gevmod)$`Pr(>F)`[1]
ALL_Pvals$TP2[1] = anova(TP2mod)$`Pr(>F)`[1]
ALL_Pvals$TP3[1] = anova(TP3mod)$`Pr(>F)`[1]
ALL_Pvals$TP4[1] = anova(TP4mod)$`Pr(>F)`[1]
ALL_Pvals$TP5[1] = anova(TP5mod)$`Pr(>F)`[1]

ALL_BFvals$Duration[1] = extractBF(BF_Duration)$bf
ALL_BFvals$Occurrence[1] = extractBF(BF_Occurance)$bf
ALL_BFvals$Coverage[1] = extractBF(BF_Coverage)$bf
ALL_BFvals$GEV[1] = extractBF(BF_GEV)$bf
ALL_BFvals$TP2[1] = extractBF(BF_TP2mod)$bf
ALL_BFvals$TP3[1] = extractBF(BF_TP3mod)$bf
ALL_BFvals$TP4[1] = extractBF(BF_TP4mod)$bf
ALL_BFvals$TP5[1] = extractBF(BF_TP5mod)$bf

# MS2
Durmod <- lmer(Duration ~ RS + (1|Subject), data = MS2,REML = FALSE)
summary(Durmod)

Occurmod <- lmer(Occurance ~ RS + (1|Subject), data = MS2,REML = FALSE)
summary(Occurmod)

Covermod <- lmer(Coverage ~ RS + (1|Subject), data = MS2,REML = FALSE)
summary(Covermod)

Gevmod <- lmer(GEV ~ RS + (1|Subject), data = MS2,REML = FALSE)
summary(Gevmod)

TP1mod <- lmer(TP1 ~ RS + (1|Subject), data = MS2,REML = FALSE)
summary(TP1mod)

TP3mod <- lmer(TP3 ~ RS + (1|Subject), data = MS2,REML = FALSE)
summary(TP3mod)

TP4mod <- lmer(TP4 ~ RS + (1|Subject), data = MS2,REML = FALSE)
summary(TP4mod)

TP5mod <- lmer(TP5 ~ RS + (1|Subject), data = MS2,REML = FALSE)
summary(TP5mod)


# MS2 BF
full_Durmod = lmBF(Duration ~  RS + Subject, data=MS2, whichRandom = 'Subject')
null_Durmod = lmBF(Duration ~ 1 + Subject, data = MS2, whichRandom = 'Subject')
BF_Duration = full_Durmod/null_Durmod

full_Occurmod = lmBF(Occurance ~  RS + Subject, data=MS2, whichRandom = 'Subject')
null_Occurmod = lmBF(Occurance ~ 1 + Subject, data = MS2, whichRandom = 'Subject')
BF_Occurance = full_Occurmod/null_Occurmod

full_Covermod = lmBF(Coverage ~  RS + Subject, data=MS2, whichRandom = 'Subject')
null_Covermod = lmBF(Coverage ~ 1 + Subject, data = MS2, whichRandom = 'Subject')
BF_Coverage = full_Covermod/null_Covermod

full_GEV = lmBF(GEV ~  RS + Subject, data=MS2, whichRandom = 'Subject')
null_GEV = lmBF(GEV ~ 1 + Subject, data = MS2, whichRandom = 'Subject')
BF_GEV = full_GEV/null_GEV

full_TP1mod = lmBF(TP1 ~  RS + Subject, data=MS2, whichRandom = 'Subject')
null_TP1mod = lmBF(TP1 ~ 1 + Subject, data = MS2, whichRandom = 'Subject')
BF_TP1mod = full_TP1mod/null_TP1mod

full_TP3mod = lmBF(TP3 ~  RS + Subject, data=MS2, whichRandom = 'Subject')
null_TP3mod = lmBF(TP3 ~ 1 + Subject, data = MS2, whichRandom = 'Subject')
BF_TP3mod = full_TP3mod/null_TP3mod

full_TP4mod = lmBF(TP4 ~  RS + Subject, data=MS2, whichRandom = 'Subject')
null_TP4mod = lmBF(TP4 ~ 1 + Subject, data = MS2, whichRandom = 'Subject')
BF_TP4mod = full_TP4mod/null_TP4mod

full_TP5mod = lmBF(TP5 ~  RS + Subject, data=MS2, whichRandom = 'Subject')
null_TP5mod = lmBF(TP5 ~ 1 + Subject, data = MS2, whichRandom = 'Subject')
BF_TP5mod = full_TP5mod/null_TP5mod

ALL_Pvals$Duration[2] = anova(Durmod)$`Pr(>F)`[1]
ALL_Pvals$Occurrence[2] = anova(Occurmod)$`Pr(>F)`[1]
ALL_Pvals$Coverage[2] = anova(Covermod)$`Pr(>F)`[1]
ALL_Pvals$GEV[2] = anova(Gevmod)$`Pr(>F)`[1]
ALL_Pvals$TP1[2] = anova(TP1mod)$`Pr(>F)`[1]
ALL_Pvals$TP3[2] = anova(TP3mod)$`Pr(>F)`[1]
ALL_Pvals$TP4[2] = anova(TP4mod)$`Pr(>F)`[1]
ALL_Pvals$TP5[2] = anova(TP5mod)$`Pr(>F)`[1]

ALL_BFvals$Duration[2] = extractBF(BF_Duration)$bf
ALL_BFvals$Occurrence[2] = extractBF(BF_Occurance)$bf
ALL_BFvals$Coverage[2] = extractBF(BF_Coverage)$bf
ALL_BFvals$GEV[2] = extractBF(BF_GEV)$bf
ALL_BFvals$TP1[2] = extractBF(BF_TP1mod)$bf
ALL_BFvals$TP3[2] = extractBF(BF_TP3mod)$bf
ALL_BFvals$TP4[2] = extractBF(BF_TP4mod)$bf
ALL_BFvals$TP5[2] = extractBF(BF_TP5mod)$bf

# MS3
Durmod <- lmer(Duration ~ RS + (1|Subject), data = MS3,REML = FALSE)
summary(Durmod)

Occurmod <- lmer(Occurance ~ RS + (1|Subject), data = MS3,REML = FALSE)
summary(Occurmod)

Covermod <- lmer(Coverage ~ RS + (1|Subject), data = MS3,REML = FALSE)
summary(Covermod)

Gevmod <- lmer(GEV ~ RS + (1|Subject), data = MS3,REML = FALSE)
summary(Gevmod)

TP1mod <- lmer(TP1 ~ RS + (1|Subject), data = MS3,REML = FALSE)
summary(TP1mod)

TP2mod <- lmer(TP2 ~ RS + (1|Subject), data = MS3,REML = FALSE)
summary(TP2mod)

TP4mod <- lmer(TP4 ~ RS + (1|Subject), data = MS3,REML = FALSE)
summary(TP4mod)

TP5mod <- lmer(TP5 ~ RS + (1|Subject), data = MS3,REML = FALSE)
summary(TP5mod)

# MS3 BF
full_Durmod = lmBF(Duration ~  RS + Subject, data=MS3, whichRandom = 'Subject')
null_Durmod = lmBF(Duration ~ 1 + Subject, data = MS3, whichRandom = 'Subject')
BF_Duration = full_Durmod/null_Durmod

full_Occurmod = lmBF(Occurance ~  RS + Subject, data=MS3, whichRandom = 'Subject')
null_Occurmod = lmBF(Occurance ~ 1 + Subject, data = MS3, whichRandom = 'Subject')
BF_Occurance = full_Occurmod/null_Occurmod

full_Covermod = lmBF(Coverage ~  RS + Subject, data=MS3, whichRandom = 'Subject')
null_Covermod = lmBF(Coverage ~ 1 + Subject, data = MS3, whichRandom = 'Subject')
BF_Coverage = full_Covermod/null_Covermod

full_GEV = lmBF(GEV ~  RS + Subject, data=MS3, whichRandom = 'Subject')
null_GEV = lmBF(GEV ~ 1 + Subject, data = MS3, whichRandom = 'Subject')
BF_GEV = full_GEV/null_GEV

full_TP1mod = lmBF(TP1 ~  RS + Subject, data=MS3, whichRandom = 'Subject')
null_TP1mod = lmBF(TP1 ~ 1 + Subject, data = MS3, whichRandom = 'Subject')
BF_TP1mod = full_TP1mod/null_TP1mod

full_TP2mod = lmBF(TP2 ~  RS + Subject, data=MS3, whichRandom = 'Subject')
null_TP2mod = lmBF(TP2 ~ 1 + Subject, data = MS3, whichRandom = 'Subject')
BF_TP2mod = full_TP2mod/null_TP2mod

full_TP4mod = lmBF(TP4 ~  RS + Subject, data=MS3, whichRandom = 'Subject')
null_TP4mod = lmBF(TP4 ~ 1 + Subject, data = MS3, whichRandom = 'Subject')
BF_TP4mod = full_TP4mod/null_TP4mod

full_TP5mod = lmBF(TP5 ~  RS + Subject, data=MS3, whichRandom = 'Subject')
null_TP5mod = lmBF(TP5 ~ 1 + Subject, data = MS3, whichRandom = 'Subject')
BF_TP5mod = full_TP5mod/null_TP5mod

ALL_Pvals$Duration[3] = anova(Durmod)$`Pr(>F)`[1]
ALL_Pvals$Occurrence[3] = anova(Occurmod)$`Pr(>F)`[1]
ALL_Pvals$Coverage[3] = anova(Covermod)$`Pr(>F)`[1]
ALL_Pvals$GEV[3] = anova(Gevmod)$`Pr(>F)`[1]
ALL_Pvals$TP1[3] = anova(TP1mod)$`Pr(>F)`[1]
ALL_Pvals$TP2[3] = anova(TP2mod)$`Pr(>F)`[1]
ALL_Pvals$TP4[3] = anova(TP4mod)$`Pr(>F)`[1]
ALL_Pvals$TP5[3] = anova(TP5mod)$`Pr(>F)`[1]

ALL_BFvals$Duration[3] = extractBF(BF_Duration)$bf
ALL_BFvals$Occurrence[3] = extractBF(BF_Occurance)$bf
ALL_BFvals$Coverage[3] = extractBF(BF_Coverage)$bf
ALL_BFvals$GEV[3] = extractBF(BF_GEV)$bf
ALL_BFvals$TP1[3] = extractBF(BF_TP1mod)$bf
ALL_BFvals$TP2[3] = extractBF(BF_TP2mod)$bf
ALL_BFvals$TP4[3] = extractBF(BF_TP4mod)$bf
ALL_BFvals$TP5[3] = extractBF(BF_TP5mod)$bf

# MS4
Durmod <- lmer(Duration ~ RS + (1|Subject), data = MS4,REML = FALSE)
summary(Durmod)

Occurmod <- lmer(Occurance ~ RS + (1|Subject), data = MS4,REML = FALSE)
summary(Occurmod)

Covermod <- lmer(Coverage ~ RS + (1|Subject), data = MS4,REML = FALSE)
summary(Covermod)

hist(MS4$GEV)
Gevmod <- lmer(GEV ~ RS + (1|Subject), data = MS4,REML = FALSE)
summary(Gevmod)

TP1mod <- lmer(TP1 ~ RS + (1|Subject), data = MS4,REML = FALSE)
summary(TP1mod)

TP2mod <- lmer(TP2 ~RS + (1|Subject), data = MS4,REML = FALSE)
summary(TP2mod)

TP3mod <- lmer(TP3 ~ RS + (1|Subject), data = MS4,REML = FALSE)
summary(TP3mod)

TP5mod <- lmer(TP5 ~ RS + (1|Subject), data = MS4,REML = FALSE)
summary(TP5mod)

# MS4 BF
full_Durmod = lmBF(Duration ~  RS + Subject, data=MS4, whichRandom = 'Subject')
null_Durmod = lmBF(Duration ~ 1 + Subject, data = MS4, whichRandom = 'Subject')
BF_Duration = full_Durmod/null_Durmod

full_Occurmod = lmBF(Occurance ~  RS + Subject, data=MS4, whichRandom = 'Subject')
null_Occurmod = lmBF(Occurance ~ 1 + Subject, data = MS4, whichRandom = 'Subject')
BF_Occurance = full_Occurmod/null_Occurmod

full_Covermod = lmBF(Coverage ~  RS + Subject, data=MS4, whichRandom = 'Subject')
null_Covermod = lmBF(Coverage ~ 1 + Subject, data = MS4, whichRandom = 'Subject')
BF_Coverage = full_Covermod/null_Covermod

full_GEV = lmBF(GEV ~  RS + Subject, data=MS4, whichRandom = 'Subject')
null_GEV = lmBF(GEV ~ 1 + Subject, data = MS4, whichRandom = 'Subject')
BF_GEV = full_GEV/null_GEV

full_TP1mod = lmBF(TP1 ~  RS + Subject, data=MS4, whichRandom = 'Subject')
null_TP1mod = lmBF(TP1 ~ 1 + Subject, data = MS4, whichRandom = 'Subject')
BF_TP1mod = full_TP1mod/null_TP1mod

full_TP2mod = lmBF(TP2 ~  RS + Subject, data=MS4, whichRandom = 'Subject')
null_TP2mod = lmBF(TP2 ~ 1 + Subject, data = MS4, whichRandom = 'Subject')
BF_TP2mod = full_TP2mod/null_TP2mod

full_TP3mod = lmBF(TP3 ~  RS + Subject, data=MS4, whichRandom = 'Subject')
null_TP3mod = lmBF(TP3 ~ 1 + Subject, data = MS4, whichRandom = 'Subject')
BF_TP3mod = full_TP3mod/null_TP3mod

full_TP5mod = lmBF(TP5 ~  RS + Subject, data=MS4, whichRandom = 'Subject')
null_TP5mod = lmBF(TP5 ~ 1 + Subject, data = MS4, whichRandom = 'Subject')
BF_TP5mod = full_TP5mod/null_TP5mod

ALL_Pvals$Duration[4] = anova(Durmod)$`Pr(>F)`[1]
ALL_Pvals$Occurrence[4] = anova(Occurmod)$`Pr(>F)`[1]
ALL_Pvals$Coverage[4] = anova(Covermod)$`Pr(>F)`[1]
ALL_Pvals$GEV[4] = anova(Gevmod)$`Pr(>F)`[1]
ALL_Pvals$TP1[4] = anova(TP1mod)$`Pr(>F)`[1]
ALL_Pvals$TP2[4] = anova(TP2mod)$`Pr(>F)`[1]
ALL_Pvals$TP3[4] = anova(TP3mod)$`Pr(>F)`[1]
ALL_Pvals$TP5[4] = anova(TP5mod)$`Pr(>F)`[1]

ALL_BFvals$Duration[4] = extractBF(BF_Duration)$bf
ALL_BFvals$Occurrence[4] = extractBF(BF_Occurance)$bf
ALL_BFvals$Coverage[4] = extractBF(BF_Coverage)$bf
ALL_BFvals$GEV[4] = extractBF(BF_GEV)$bf
ALL_BFvals$TP1[4] = extractBF(BF_TP1mod)$bf
ALL_BFvals$TP2[4] = extractBF(BF_TP2mod)$bf
ALL_BFvals$TP3[4] = extractBF(BF_TP3mod)$bf
ALL_BFvals$TP5[4] = extractBF(BF_TP5mod)$bf

# MS5
Durmod <- lmer(Duration ~ RS + (1|Subject), data = MS5,REML = FALSE)
summary(Durmod)

Occurmod <- lmer(Occurance ~ RS + (1|Subject), data = MS5,REML = FALSE)
summary(Occurmod)

Covermod <- lmer(Coverage ~ RS + (1|Subject), data = MS5,REML = FALSE)
summary(Covermod)

Gevmod <- lmer(GEV ~ RS + (1|Subject), data = MS5,REML = FALSE)
summary(Gevmod)

TP1mod <- lmer(TP1 ~ RS + (1|Subject), data = MS5,REML = FALSE)
summary(TP1mod)

TP2mod <- lmer(TP2 ~ RS + (1|Subject), data = MS5,REML = FALSE)
summary(TP2mod)

TP3mod <- lmer(TP3 ~ RS + (1|Subject), data = MS5,REML = FALSE)
summary(TP3mod)

TP4mod <- lmer(TP4 ~RS + (1|Subject), data = MS5,REML = FALSE)
summary(TP4mod)

# MS5 BF
full_Durmod = lmBF(Duration ~  RS + Subject, data=MS5, whichRandom = 'Subject')
null_Durmod = lmBF(Duration ~ 1 + Subject, data = MS5, whichRandom = 'Subject')
BF_Duration = full_Durmod/null_Durmod

full_Occurmod = lmBF(Occurance ~  RS + Subject, data=MS5, whichRandom = 'Subject')
null_Occurmod = lmBF(Occurance ~ 1 + Subject, data = MS5, whichRandom = 'Subject')
BF_Occurance = full_Occurmod/null_Occurmod

full_Covermod = lmBF(Coverage ~  RS + Subject, data=MS5, whichRandom = 'Subject')
null_Covermod = lmBF(Coverage ~ 1 + Subject, data = MS5, whichRandom = 'Subject')
BF_Coverage = full_Covermod/null_Covermod

full_GEV = lmBF(GEV ~  RS + Subject, data=MS5, whichRandom = 'Subject')
null_GEV = lmBF(GEV ~ 1 + Subject, data = MS5, whichRandom = 'Subject')
BF_GEV = full_GEV/null_GEV

full_TP1mod = lmBF(TP1 ~  RS + Subject, data=MS5, whichRandom = 'Subject')
null_TP1mod = lmBF(TP1 ~ 1 + Subject, data = MS5, whichRandom = 'Subject')
BF_TP1mod = full_TP1mod/null_TP1mod

full_TP2mod = lmBF(TP2 ~  RS + Subject, data=MS5, whichRandom = 'Subject')
null_TP2mod = lmBF(TP2 ~ 1 + Subject, data = MS5, whichRandom = 'Subject')
BF_TP2mod = full_TP2mod/null_TP2mod

full_TP3mod = lmBF(TP3 ~  RS + Subject, data=MS5, whichRandom = 'Subject')
null_TP3mod = lmBF(TP3 ~ 1 + Subject, data = MS5, whichRandom = 'Subject')
BF_TP3mod = full_TP3mod/null_TP3mod

full_TP4mod = lmBF(TP4 ~  RS + Subject, data=MS5, whichRandom = 'Subject')
null_TP4mod = lmBF(TP4 ~ 1 + Subject, data = MS5, whichRandom = 'Subject')
BF_TP4mod = full_TP4mod/null_TP4mod


ALL_Pvals$Duration[5] = anova(Durmod)$`Pr(>F)`[1]
ALL_Pvals$Occurrence[5] = anova(Occurmod)$`Pr(>F)`[1]
ALL_Pvals$Coverage[5] = anova(Covermod)$`Pr(>F)`[1]
ALL_Pvals$GEV[5] = anova(Gevmod)$`Pr(>F)`[1]
ALL_Pvals$TP1[5] = anova(TP1mod)$`Pr(>F)`[1]
ALL_Pvals$TP2[5] = anova(TP2mod)$`Pr(>F)`[1]
ALL_Pvals$TP3[5] = anova(TP3mod)$`Pr(>F)`[1]
ALL_Pvals$TP4[5] = anova(TP4mod)$`Pr(>F)`[1]

ALL_BFvals$Duration[5] = extractBF(BF_Duration)$bf
ALL_BFvals$Occurrence[5] = extractBF(BF_Occurance)$bf
ALL_BFvals$Coverage[5] = extractBF(BF_Coverage)$bf
ALL_BFvals$GEV[5] = extractBF(BF_GEV)$bf
ALL_BFvals$TP1[5] = extractBF(BF_TP1mod)$bf
ALL_BFvals$TP2[5] = extractBF(BF_TP2mod)$bf
ALL_BFvals$TP3[5] = extractBF(BF_TP3mod)$bf
ALL_BFvals$TP4[5] = extractBF(BF_TP4mod)$bf

RS1vsRS2 = ALL_Pvals

# Run the FDR correction for all RS1 vs RS2 Pvals
ControlPvalues <- unlist(RS1vsRS2)
non_na_indices <- which(!is.na(ControlPvalues))
All_adj_p_values <- p.adjust(ControlPvalues[non_na_indices], method = "fdr")
ControlPvalues[non_na_indices] <- All_adj_p_values
SigControlled = ControlPvalues[ControlPvalues<.05]
ControlPvalues[ControlPvalues<.05]

ALL_BFvals

# Run1 - A = 5, B = 2, C = 3, D = 1, E = 4
# Run2 - A = 5, B = 3, C = 2, D = 1, E = 4
# Run3 - A = 5, B = 2, C = 3, D = 1, E = 4
# Run4 - A = 3, B = 4, C = 2, D = 1, E = 5
# Run5 - A = 5, B = 1, C = 4, D = 2, E = 3

# Save the results as dataframes
SigControlled = as.data.frame(SigControlled)
SigControlled$row_names <- rownames(SigControlled)
write_xlsx(SigControlled,paste0(Path,runs[thisrun],"/RS1vsRS2SigCorrected.xlsx"))


ControlPvalues = as.data.frame(ControlPvalues)
ControlPvalues$row_names <- rownames(ControlPvalues)
write_xlsx(ControlPvalues,paste0(Path,runs[thisrun],"/RS1vsRS2Corrected.xlsx"))


RS1vsRS2 = as.data.frame(RS1vsRS2)
RS1vsRS2$row_names <- rownames(RS1vsRS2)
write_xlsx(RS1vsRS2,paste0(Path,runs[thisrun],"/RS1vsRS2Uncorrected.xlsx"))


ALL_BFvals = as.data.frame(ALL_BFvals)
ALL_BFvals$row_names <- rownames(ALL_BFvals)
write_xlsx(ALL_BFvals,paste0(Path,runs[thisrun],"/RS1vsRS2BFvals.xlsx"))





## Microstates x Area ----
fullpath = paste0(Path,runs[thisrun],"/MicrostatsxArea.xlsx")
fullpath
MicrostatsxArea <- read_excel(fullpath)
MicrostatsxArea <- as.data.frame(MicrostatsxArea)
# MicrostatsxArea = MicrostatsxArea[MicrostatsxArea$Condition == 'N',]

MicrostatsxArea$GFP <- as.numeric(scale(MicrostatsxArea$GFP))
MicrostatsxArea$Duration <- as.numeric(scale(MicrostatsxArea$Duration))
MicrostatsxArea$Occurance <- as.numeric(scale(MicrostatsxArea$Occurance))
MicrostatsxArea$Coverage <- as.numeric((MicrostatsxArea$Coverage))

MS1RS1 = MicrostatsxArea[MicrostatsxArea$MS == 'MS1' & MicrostatsxArea$RS == 'RS1',]
MS2RS1 = MicrostatsxArea[MicrostatsxArea$MS == 'MS2' & MicrostatsxArea$RS == 'RS1',]
MS3RS1 = MicrostatsxArea[MicrostatsxArea$MS == 'MS3' & MicrostatsxArea$RS == 'RS1',]
MS4RS1 = MicrostatsxArea[MicrostatsxArea$MS == 'MS4' & MicrostatsxArea$RS == 'RS1',]
MS5RS1 = MicrostatsxArea[MicrostatsxArea$MS == 'MS5' & MicrostatsxArea$RS == 'RS1',]


# Function to compute p-values for individual predictors
compute_p_values_Area <- function(predictor_names, dependent_variable, data) {
  # Initialize a named vector to store p-values
  p_values <- numeric(length(predictor_names))
  p_values <- setNames(c(p_values), c(predictor_names))
  BF_values <- p_values
  
  # Loop over each predictor
  for (i in 1:length(predictor_names)) {
    # Create the formula for the model
    formula <- formula(paste(dependent_variable, "~", predictor_names[i]))
    formulaBF <- formula(paste(dependent_variable, "~", predictor_names[i]))
    
    # Fit the linear model
    model <- lm(formula, data = data)
    
    # Fit the BF model
    modelBF_full <- lmBF(formulaBF,data = data)
    
    print(summary(model))
    # Extract the p-value for the predictor and store it with the predictor name
    print(summary(model)$coefficients[predictor_names[i], "Pr(>|t|)"])
    p_values[predictor_names[i]] <- summary(model)$coefficients[2, "Pr(>|t|)"]
    BF_values[predictor_names[i]] <- extractBF(modelBF_full)$bf
  }
  # Return the named vector of p-values
  print(p_values)
  print(BF_values)
  return(list(p_values = p_values, BF_values = BF_values))
}


### Length ----

# Fit predictors separately to avoid the confusion of the other predictors
predictors <- c("Duration", "Occurance", "Coverage", "GEV", "TP2", "TP3", "TP4","TP5")
MS1_Length <- compute_p_values_Area(predictors, "Length", MS1RS1)

predictors <- c("Duration", "Occurance", "Coverage", "GEV", "TP1", "TP3", "TP4","TP5")
MS2_Length <- compute_p_values_Area(predictors, "Length", MS2RS1)

predictors <- c("Duration", "Occurance", "Coverage", "GEV", "TP2", "TP1", "TP4","TP5")
MS3_Length <- compute_p_values_Area(predictors, "Length", MS3RS1)

predictors <- c("Duration", "Occurance", "Coverage", "GEV", "TP2", "TP3", "TP1","TP5")
MS4_Length <- compute_p_values_Area(predictors, "Length", MS4RS1)

predictors <- c("Duration", "Occurance", "Coverage", "GEV", "TP1", "TP2", "TP3", "TP4")
MS5_Length <- compute_p_values_Area(predictors, "Length", MS5RS1)


## Multiple comparisons correction for Area Length
Pvalues <- unlist(c(MS1_Length$p_values,MS2_Length$p_values,MS3_Length$p_values,MS4_Length$p_values,MS5_Length$p_values))
CorrectedPvalues <- Pvalues
non_na_indices <- which(!is.na(CorrectedPvalues))
All_adj_p_values <- p.adjust(CorrectedPvalues[non_na_indices], method = "fdr")
CorrectedPvalues[non_na_indices] <- All_adj_p_values
SigCorrected <- CorrectedPvalues[CorrectedPvalues<.05]
CorrectedPvalues[CorrectedPvalues<.05]

BFvalues <- unlist(c(MS1_Length$BF_values,MS2_Length$BF_values,MS3_Length$BF_values,MS4_Length$BF_values,MS5_Length$BF_values))



# Save the data
Pvalues_df <- data.frame(
  Name = names(Pvalues),
  Value = Pvalues,
  row.names = NULL
)
write_xlsx(Pvalues_df,paste0(Path,runs[thisrun],"/AreaLengthUncorrected.xlsx"))


CorrectedPvalues_df <- data.frame(
  Name = names(CorrectedPvalues),
  Value = CorrectedPvalues,
  row.names = NULL
)
write_xlsx(CorrectedPvalues_df,paste0(Path,runs[thisrun],"/AreaLengthCorrected.xlsx"))

BFvalues_df <- data.frame(
  Name = names(BFvalues),
  Value = BFvalues,
  row.names = NULL
)
write_xlsx(BFvalues_df,paste0(Path,runs[thisrun],"/AreaLengthBFvals.xlsx"))

SigCorrected_df <- data.frame(
  Name = names(SigCorrected),
  Value = SigCorrected,
  row.names = NULL
)
write_xlsx(SigCorrected_df,paste0(Path,runs[thisrun],"/AreaLengthSigCorrected.xlsx"))


# Compute correlation among predictors
subset_cols <- c("GFP", "Duration", "Occurance", "Coverage", "GEV", "TP1", "TP2", "TP3", "TP4","TP5")
# Compute the correlation matrix for the subset of columns
cor_matrix <- cor(MS5RS1[, predictors])
# Print the correlation matrix
print(cor_matrix)



### Width ----
MS1RS1 <- MS1RS1[!is.na(MS1RS1$Width), ]
MS2RS1 <- MS2RS1[!is.na(MS2RS1$Width), ]
MS3RS1 <- MS3RS1[!is.na(MS3RS1$Width), ]
MS4RS1 <- MS4RS1[!is.na(MS4RS1$Width), ]
MS5RS1 <- MS5RS1[!is.na(MS5RS1$Width), ]

# Fit predictors separately to avoid the confusion of the other predictors
predictors <- c("Duration", "Occurance", "Coverage", "GEV", "TP2", "TP3", "TP4","TP5")
MS1_Width <- compute_p_values_Area(predictors, "Width", MS1RS1)

predictors <- c("Duration", "Occurance", "Coverage", "GEV", "TP1", "TP3", "TP4","TP5")
MS2_Width <- compute_p_values_Area(predictors, "Width", MS2RS1)

predictors <- c("Duration", "Occurance", "Coverage", "GEV", "TP2", "TP1", "TP4","TP5")
MS3_Width <- compute_p_values_Area(predictors, "Width", MS3RS1)

predictors <- c("Duration", "Occurance", "Coverage", "GEV", "TP2", "TP3", "TP1","TP5")
MS4_Width <- compute_p_values_Area(predictors, "Width", MS4RS1)

predictors <- c("Duration", "Occurance", "Coverage", "GEV", "TP1", "TP2", "TP3", "TP4")
MS5_Width <- compute_p_values_Area(predictors, "Width", MS5RS1)

# Multiple comparisons correction
Pvalues <- unlist(c(MS1_Width$p_values,MS2_Width$p_values,MS3_Width$p_values,MS4_Width$p_values,MS5_Width$p_values))
CorrectedPvalues <- Pvalues
non_na_indices <- which(!is.na(CorrectedPvalues))
All_adj_p_values <- p.adjust(CorrectedPvalues[non_na_indices], method = "fdr")
CorrectedPvalues[non_na_indices] <- All_adj_p_values
SigCorrected <- CorrectedPvalues[CorrectedPvalues<.05]
CorrectedPvalues[CorrectedPvalues<.05]

BFvalues <- unlist(c(MS1_Width$BF_values,MS2_Width$BF_values,MS3_Width$BF_values,MS4_Width$BF_values,MS5_Width$BF_values))


# Save the data
Pvalues_df <- data.frame(
  Name = names(Pvalues),
  Value = Pvalues,
  row.names = NULL
)
write_xlsx(Pvalues_df,paste0(Path,runs[thisrun],"/AreaWidthUncorrected.xlsx"))


CorrectedPvalues_df <- data.frame(
  Name = names(CorrectedPvalues),
  Value = CorrectedPvalues,
  row.names = NULL
)
write_xlsx(CorrectedPvalues_df,paste0(Path,runs[thisrun],"/AreaWidthCorrected.xlsx"))

BFvalues_df <- data.frame(
  Name = names(BFvalues),
  Value = BFvalues,
  row.names = NULL
)
write_xlsx(BFvalues_df,paste0(Path,runs[thisrun],"/AreaWidthBFvals.xlsx"))

SigCorrected_df <- data.frame(
  Name = names(SigCorrected),
  Value = SigCorrected,
  row.names = NULL
)
write_xlsx(SigCorrected_df,paste0(Path,runs[thisrun],"/AreawidthSigCorrected.xlsx"))



## Microstates x PPRatings ----
fullpath = paste0(Path,runs[thisrun],"/MicrostatsxPPratings.xlsx")
fullpath
MicrostatsxPP <- read_excel(fullpath)
# View(MicrostatsxPP)

# Convert the tibble to a data frame
MicrostatsxPP <- as.data.frame(MicrostatsxPP)

# Switch the Intensity ratings for participant 21
T0HFS21 = MicrostatsxPP$IntensityT0[MicrostatsxPP$Subject==21 & MicrostatsxPP$Arm == 'MFS']
T1HFS21 = MicrostatsxPP$IntensityRatings[MicrostatsxPP$Subject==21 & MicrostatsxPP$Arm == 'MFS']
T0Cont21 = MicrostatsxPP$IntensityT0[MicrostatsxPP$Subject==21 & MicrostatsxPP$Arm == 'Cont']
T1Cont21 =  MicrostatsxPP$IntensityRatings[MicrostatsxPP$Subject==21 & MicrostatsxPP$Arm == 'Cont']

MicrostatsxPP$IntensityT0[MicrostatsxPP$Subject==21 & MicrostatsxPP$Arm == 'MFS'] = T0Cont21
MicrostatsxPP$IntensityT0[MicrostatsxPP$Subject==21 & MicrostatsxPP$Arm == 'Cont'] = T0HFS21
MicrostatsxPP$IntensityRatings[MicrostatsxPP$Subject==21 & MicrostatsxPP$Arm == 'MFS'] = T1Cont21
MicrostatsxPP$IntensityRatings[MicrostatsxPP$Subject==21 & MicrostatsxPP$Arm == 'Cont'] = T1HFS21

# Switch the Unpleasantness ratings for participant 21
T0HFS21 = MicrostatsxPP$UnpleasantnessT0[MicrostatsxPP$Subject==21 & MicrostatsxPP$Arm == 'MFS']
T1HFS21 = MicrostatsxPP$UnpleasantnessRatings[MicrostatsxPP$Subject==21 & MicrostatsxPP$Arm == 'MFS']
T0Cont21 = MicrostatsxPP$UnpleasantnessT0[MicrostatsxPP$Subject==21 & MicrostatsxPP$Arm == 'Cont']
T1Cont21 =  MicrostatsxPP$UnpleasantnessRatings[MicrostatsxPP$Subject==21 & MicrostatsxPP$Arm == 'Cont']

MicrostatsxPP$UnpleasantnessT0[MicrostatsxPP$Subject==21 & MicrostatsxPP$Arm == 'MFS'] = T0Cont21
MicrostatsxPP$UnpleasantnessT0[MicrostatsxPP$Subject==21 & MicrostatsxPP$Arm == 'Cont'] = T0HFS21
MicrostatsxPP$UnpleasantnessRatings[MicrostatsxPP$Subject==21 & MicrostatsxPP$Arm == 'MFS'] = T1Cont21
MicrostatsxPP$UnpleasantnessRatings[MicrostatsxPP$Subject==21 & MicrostatsxPP$Arm == 'Cont'] = T1HFS21



MicrostatsxPP$GFP <- as.numeric(scale(MicrostatsxPP$GFP))
MicrostatsxPP$Duration <- as.numeric(scale(MicrostatsxPP$Duration))
MicrostatsxPP$Occurance <- as.numeric(scale(MicrostatsxPP$Occurance))
MicrostatsxPP$GEV <- as.numeric(scale(MicrostatsxPP$GEV))


MicrostatsxPP$Subject <- as.factor(MicrostatsxPP$Subject)
MS1RS1 = MicrostatsxPP[MicrostatsxPP$MS == 'MS1' & MicrostatsxPP$RS == 'RS1',]
MS2RS1 = MicrostatsxPP[MicrostatsxPP$MS == 'MS2' & MicrostatsxPP$RS == 'RS1',]
MS3RS1 = MicrostatsxPP[MicrostatsxPP$MS == 'MS3' & MicrostatsxPP$RS == 'RS1',]
MS4RS1 = MicrostatsxPP[MicrostatsxPP$MS == 'MS4' & MicrostatsxPP$RS == 'RS1',]
MS5RS1 = MicrostatsxPP[MicrostatsxPP$MS == 'MS5' & MicrostatsxPP$RS == 'RS1',]

# Function to compute p-values for individual predictors
Fit_Models_PP <- function(predictor_names, dependent_variable, data,MS) {
  # Initialize a named vector to store p-values
  p_values <- numeric(length(predictor_names))
  new_predictor_names <- paste0(MS, predictor_names)
  p_values <- setNames(c(p_values), c(new_predictor_names))
  BF_values <- p_values
  
  # Loop over each predictor
  for (i in 1:length(predictor_names)) {
    # Create the formula for the model
    formula <- formula(paste(dependent_variable, "~",'IntensityT0 + ', predictor_names[i],'*Arm + (1|Subject)'))
    formula_BFfull <- formula(paste(dependent_variable, "~",'IntensityT0 + ', predictor_names[i],'*Arm + Subject'))
    formula_BFnull <- formula(paste(dependent_variable, "~",'IntensityT0 + Arm + Subject'))
    
    # Fit the linear model
    model <- lmer(formula, data = data,REML = FALSE)
    BF_full <- lmBF(formula_BFfull, data = data,whichRandom = 'Subject')
    BF_null <- lmBF(formula_BFnull, data = data,whichRandom = 'Subject')
    BF <- BF_full/BF_null
    
    print(summary(model))
    # Extract the p-value for the predictor and store it with the predictor name
    print(summary(model)$coefficients[5, "Pr(>|t|)"])
    p_values[new_predictor_names[i]] <- summary(model)$coefficients[5, "Pr(>|t|)"]
    BF_values[new_predictor_names[i]] <- extractBF(BF)$bf
  }
  print(p_values)
  print(BF_values)
  return(list(p_values = p_values, BF_values = BF_values))
}


### Intensity ----

# Fit predictors individually

# Fit predictors separately to avoid the confusion of the other predictors
predictors <- c("Duration", "Occurance", "Coverage", "GEV", "TP2", "TP3", "TP4","TP5")
MS1_Intensity <- Fit_Models_PP(predictors, "IntensityRatings", MS1RS1,"MS1")

predictors <- c("Duration", "Occurance", "Coverage", "GEV", "TP1", "TP3", "TP4","TP5")
MS2_Intensity <- Fit_Models_PP(predictors, "IntensityRatings", MS2RS1,"MS2")

predictors <- c("Duration", "Occurance", "Coverage", "GEV", "TP2", "TP1", "TP4","TP5")
MS3_Intensity <- Fit_Models_PP(predictors, "IntensityRatings", MS3RS1,"MS3")

predictors <- c("Duration", "Occurance", "Coverage", "GEV", "TP2", "TP3", "TP1","TP5")
MS4_Intensity <- Fit_Models_PP(predictors, "IntensityRatings", MS4RS1,"MS4")

predictors <- c("Duration", "Occurance", "Coverage", "GEV", "TP1", "TP2", "TP3", "TP4")
MS5_Intensity <- Fit_Models_PP(predictors, "IntensityRatings", MS5RS1,"MS5")

## Multiple comparisons correction for Area Length
Pvalues <- unlist(c(MS1_Intensity$p_values,MS2_Intensity$p_values,MS3_Intensity$p_values,MS4_Intensity$p_values,MS5_Intensity$p_values))
CorrectedPvalues <- Pvalues
non_na_indices <- which(!is.na(CorrectedPvalues))
All_adj_p_values <- p.adjust(CorrectedPvalues[non_na_indices], method = "fdr")
CorrectedPvalues[non_na_indices] <- All_adj_p_values
SigCorrected <- CorrectedPvalues[CorrectedPvalues<.05]
CorrectedPvalues[CorrectedPvalues<.05]

BFvalues <- unlist(c(MS1_Intensity$BF_values,MS2_Intensity$BF_values,MS3_Intensity$BF_values,MS4_Intensity$BF_values,MS5_Intensity$BF_values))

# Save the data
Pvalues_df <- data.frame(
  Name = names(Pvalues),
  Value = Pvalues,
  row.names = NULL
)
write_xlsx(Pvalues_df,paste0(Path,runs[thisrun],"/PPIntensityUncorrected.xlsx"))


CorrectedPvalues_df <- data.frame(
  Name = names(CorrectedPvalues),
  Value = CorrectedPvalues,
  row.names = NULL
)
write_xlsx(CorrectedPvalues_df,paste0(Path,runs[thisrun],"/PPIntensityCorrected.xlsx"))

BFvalues_df <- data.frame(
  Name = names(BFvalues),
  Value = BFvalues,
  row.names = NULL
)
write_xlsx(BFvalues_df,paste0(Path,runs[thisrun],"/PPIntensityBFvals.xlsx"))

SigCorrected_df <- data.frame(
  Name = names(SigCorrected),
  Value = SigCorrected,
  row.names = NULL
)
write_xlsx(SigCorrected_df,paste0(Path,runs[thisrun],"/PPIntensitySigCorrected.xlsx"))


# Plotting for Intensity ----

dat = MS3RS1

axsize = 20
mytheme = theme(
  axis.title.x = element_text(size = axsize),
  axis.text.x = element_text(size = axsize-2),
  axis.title.y = element_text(size = axsize),
  axis.text.y = element_text(size = axsize-2),
  plot.title = element_text(size=25,hjust = 0.5),
  legend.text=element_text(size=15),
  legend.title=element_text(size=20))


# Duration
Fullmod <- lmer(IntensityRatings ~ IntensityT0 +  Duration*Arm + (1|Subject), data = dat,REML = FALSE)
summary(Fullmod)
Intplot <- ggpredict(Fullmod,terms = c('Duration','Arm'))

Intplot <- as.data.frame(Intplot)
ggplot()+
  geom_point(data = dat,
             aes(x = Duration, y = IntensityRatings, colour = Arm),
             position=position_jitter(h=0.05, w=.8), alpha = .5) + scale_colour_manual(values = c("black","red")) + 
  geom_line(data = Intplot,
            aes(x = x, y = predicted, colour = group),
            size = 1.2) + scale_fill_manual(values = c("black","red")) + 
  geom_ribbon(data = Intplot, aes(x = x, ymin = conf.low, ymax = conf.high, fill = group), alpha= 0.2, show.legend = FALSE) +
  theme_blank() + mytheme + theme(legend.position = 'top') + xlab("Duration (s)") + ylab("Intensity Ratings (NRS)")
 

# Occurance
Fullmod <- lmer(IntensityRatings ~ IntensityT0 +  Occurance*Arm + (1|Subject), data = dat,REML = FALSE)
summary(Fullmod)
Intplot <- ggpredict(Fullmod,terms = c('Occurance','Arm'))

Intplot <- as.data.frame(Intplot)
ggplot()+
  geom_point(data = dat,
             aes(x = Occurance, y = IntensityRatings, colour = Arm),
             position=position_jitter(h=0.05, w=.5), alpha = .5) + scale_colour_manual(values = c("black","red")) + 
  geom_line(data = Intplot,
            aes(x = x, y = predicted, colour = group),
            size = 1.2) + scale_fill_manual(values = c("black","red")) + 
  geom_ribbon(data = Intplot, aes(x = x, ymin = conf.low, ymax = conf.high, fill = group), alpha= 0.2, show.legend = FALSE) +
  theme_blank() + mytheme + theme(legend.position = 'top') + xlab("Occurrence/s") + ylab("Intensity Ratings (NRS)")


# Coverage
Fullmod <- lmer(IntensityRatings ~ IntensityT0 +  Coverage*Arm + (1|Subject), data = dat,REML = FALSE)
summary(Fullmod)
Intplot <- ggpredict(Fullmod,terms = c('Coverage','Arm'))

Intplot <- as.data.frame(Intplot)
ggplot()+
  geom_point(data = dat,
             aes(x = Coverage, y = IntensityRatings, colour = Arm),
             position=position_jitter(h=0.05, w=.03), alpha = .5) + scale_colour_manual(values = c("black","red")) + 
  geom_line(data = Intplot,
            aes(x = x, y = predicted, colour = group),
            size = 1.2) + scale_fill_manual(values = c("black","red")) + 
  geom_ribbon(data = Intplot, aes(x = x, ymin = conf.low, ymax = conf.high, fill = group), alpha= 0.2, show.legend = FALSE) +
  theme_blank() + mytheme + theme(legend.position = 'top') + xlab("Coverage (prop)") + ylab("Intensity Ratings (NRS)")

# GEV
Fullmod <- lmer(IntensityRatings ~ IntensityT0 +  GEV*Arm + (1|Subject), data = dat,REML = FALSE)
summary(Fullmod)
Intplot <- ggpredict(Fullmod,terms = c('GEV','Arm'))

Intplot <- as.data.frame(Intplot)
ggplot()+
  geom_point(data = dat,
             aes(x = GEV, y = IntensityRatings, colour = Arm),
             position=position_jitter(h=0.05, w=.5), alpha = .5) + scale_colour_manual(values = c("black","red")) + 
  geom_line(data = Intplot,
            aes(x = x, y = predicted, colour = group),
            size = 1.2) + scale_fill_manual(values = c("black","red")) + 
  geom_ribbon(data = Intplot, aes(x = x, ymin = conf.low, ymax = conf.high, fill = group), alpha= 0.2, show.legend = FALSE) +
  theme_blank() + mytheme + theme(legend.position = 'top') + xlab("GEV (prop)") + ylab("Intensity Ratings (NRS)")



# TP1
Fullmod <- lmer(IntensityRatings ~ IntensityT0 +  TP1*Arm + (1|Subject), data = dat,REML = FALSE)
summary(Fullmod)
Intplot <- ggpredict(Fullmod,terms = c('TP1','Arm'))

Intplot <- as.data.frame(Intplot)
ggplot()+
  geom_point(data = dat,
             aes(x = TP1, y = IntensityRatings, colour = Arm),
             position=position_jitter(h=0.05, w=.04), alpha = .5) + scale_colour_manual(values = c("black","red")) + 
  geom_line(data = Intplot,
            aes(x = x, y = predicted, colour = group),
            size = 1.2) + scale_fill_manual(values = c("black","red")) + 
  geom_ribbon(data = Intplot, aes(x = x, ymin = conf.low, ymax = conf.high, fill = group), alpha= 0.2, show.legend = FALSE) +
  theme_blank() + mytheme + theme(legend.position = 'top') + xlab("Transition Probability - MS1") + ylab("Intensity Ratings (NRS)")


# TP2
Fullmod <- lmer(IntensityRatings ~ IntensityT0 +  TP2*Arm + (1|Subject), data = dat,REML = FALSE)
summary(Fullmod)
Intplot <- ggpredict(Fullmod,terms = c('TP2','Arm'))

Intplot <- as.data.frame(Intplot)
ggplot()+
  geom_point(data = dat,
             aes(x = TP2, y = IntensityRatings, colour = Arm),
             position=position_jitter(h=0.05, w=.03), alpha = .5) + scale_colour_manual(values = c("black","red")) + 
  geom_line(data = Intplot,
            aes(x = x, y = predicted, colour = group),
            size = 1.2) + scale_fill_manual(values = c("black","red")) + 
  geom_ribbon(data = Intplot, aes(x = x, ymin = conf.low, ymax = conf.high, fill = group), alpha= 0.2, show.legend = FALSE) +
  theme_blank() + mytheme + theme(legend.position = 'top') + xlab("Transition Probability - MS2") + ylab("Intensity Ratings (NRS)")


# TP3
Fullmod <- lmer(IntensityRatings ~ IntensityT0 +  TP3*Arm + (1|Subject), data = dat,REML = FALSE)
summary(Fullmod)
Intplot <- ggpredict(Fullmod,terms = c('TP3','Arm'))

Intplot <- as.data.frame(Intplot)
ggplot()+
  geom_point(data = dat,
             aes(x = TP3, y = IntensityRatings, colour = Arm),
             position=position_jitter(h=0.05, w=.03), alpha = .5) + scale_colour_manual(values = c("black","red")) + 
  geom_line(data = Intplot,
            aes(x = x, y = predicted, colour = group),
            size = 1.2) + scale_fill_manual(values = c("black","red")) + 
  geom_ribbon(data = Intplot, aes(x = x, ymin = conf.low, ymax = conf.high, fill = group), alpha= 0.2, show.legend = FALSE) +
  theme_blank() + mytheme + theme(legend.position = 'top')+ xlab("Transition Probability - MS3") + ylab("Intensity Ratings (NRS)")

# TP4
Fullmod <- lmer(IntensityRatings ~ IntensityT0 +  TP4*Arm + (1|Subject), data = dat,REML = FALSE)
summary(Fullmod)
Intplot <- ggpredict(Fullmod,terms = c('TP4','Arm'))

Intplot <- as.data.frame(Intplot)
ggplot()+
  geom_point(data = dat,
             aes(x = TP4, y = IntensityRatings, colour = Arm),
             position=position_jitter(h=0.05, w=.03), alpha = .5) + scale_colour_manual(values = c("black","red")) + 
  geom_line(data = Intplot,
            aes(x = x, y = predicted, colour = group),
            size = 1.2) + scale_fill_manual(values = c("black","red")) + 
  geom_ribbon(data = Intplot, aes(x = x, ymin = conf.low, ymax = conf.high, fill = group), alpha= 0.2, show.legend = FALSE) +
  theme_blank() + mytheme + theme(legend.position = 'top')+ xlab("Transition Probability - MS4") + ylab("Intensity Ratings (NRS)")


# TP5
Fullmod <- lmer(IntensityRatings ~ IntensityT0 +  TP5*Arm + (1|Subject), data = dat,REML = FALSE)
summary(Fullmod)
Intplot <- ggpredict(Fullmod,terms = c('TP5','Arm'))

Intplot <- as.data.frame(Intplot)
ggplot()+
  geom_point(data = dat,
             aes(x = TP5, y = IntensityRatings, colour = Arm),
             position=position_jitter(h=0.05, w=.03), alpha = .5) + scale_colour_manual(values = c("black","red")) + 
  geom_line(data = Intplot,
            aes(x = x, y = predicted, colour = group),
            size = 1.2) + scale_fill_manual(values = c("black","red")) + 
  geom_ribbon(data = Intplot, aes(x = x, ymin = conf.low, ymax = conf.high, fill = group), alpha= 0.2, show.legend = FALSE) +
  theme_blank() + mytheme + theme(legend.position = 'top')+ xlab("Transition Probability - MS5") + ylab("Intensity Ratings (NRS)")




### Unpleasantness ----

# Function to compute p-values for individual predictors
Fit_Models_PP <- function(predictor_names, dependent_variable, data,MS) {
  # Initialize a named vector to store p-values
  p_values <- numeric(length(predictor_names))
  new_predictor_names <- paste0(MS, predictor_names)
  p_values <- setNames(c(p_values), c(new_predictor_names))
  BF_values <- p_values
  
  # Loop over each predictor
  for (i in 1:length(predictor_names)) {
    # Create the formula for the model
    formula <- formula(paste(dependent_variable, "~",'UnpleasantnessT0 + ', predictor_names[i],'*Arm + (1|Subject)'))
    formula_BFfull <- formula(paste(dependent_variable, "~",'UnpleasantnessT0 + ', predictor_names[i],'*Arm + Subject'))
    formula_BFnull <- formula(paste(dependent_variable, "~",'UnpleasantnessT0 + Arm + Subject'))
    
    # Fit the linear model
    model <- lmer(formula, data = data,REML = FALSE)
    BF_full <- lmBF(formula_BFfull, data = data,whichRandom = 'Subject')
    BF_null <- lmBF(formula_BFnull, data = data,whichRandom = 'Subject')
    BF <- BF_full/BF_null
    
    print(summary(model))
    # Extract the p-value for the predictor and store it with the predictor name
    print(summary(model)$coefficients[5, "Pr(>|t|)"])
    p_values[new_predictor_names[i]] <- summary(model)$coefficients[5, "Pr(>|t|)"]
    BF_values[new_predictor_names[i]] <- extractBF(BF)$bf
  }
  print(p_values)
  print(BF_values)
  return(list(p_values = p_values, BF_values = BF_values))
}

# Fit predictors individually
# Fit predictors separately to avoid the confusion of the other predictors
predictors <- c("Duration", "Occurance", "Coverage", "GEV", "TP2", "TP3", "TP4","TP5")
MS1_Unpleasantness <- Fit_Models_PP(predictors, "UnpleasantnessRatings", MS1RS1,'MS1')

predictors <- c("Duration", "Occurance", "Coverage", "GEV", "TP1", "TP3", "TP4","TP5")
MS2_Unpleasantness <- Fit_Models_PP(predictors, "UnpleasantnessRatings", MS2RS1,'MS2')

predictors <- c("Duration", "Occurance", "Coverage", "GEV", "TP2", "TP1", "TP4","TP5")
MS3_Unpleasantness <- Fit_Models_PP(predictors, "UnpleasantnessRatings", MS3RS1,'MS3')

predictors <- c("Duration", "Occurance", "Coverage", "GEV", "TP2", "TP3", "TP1","TP5")
MS4_Unpleasantness <- Fit_Models_PP(predictors, "UnpleasantnessRatings", MS4RS1,'MS4')

predictors <- c("Duration", "Occurance", "Coverage", "GEV", "TP1", "TP2", "TP3", "TP4")
MS5_Unpleasantness <- Fit_Models_PP(predictors, "UnpleasantnessRatings", MS5RS1,'MS5')


# Multiple comparisons correction
Pvalues <- unlist(c(MS1_Unpleasantness$p_values,MS2_Unpleasantness$p_values,MS3_Unpleasantness$p_values,MS4_Unpleasantness$p_values,MS5_Unpleasantness$p_values))
CorrectedPvalues <- Pvalues
non_na_indices <- which(!is.na(CorrectedPvalues))
All_adj_p_values <- p.adjust(CorrectedPvalues[non_na_indices], method = "fdr")
CorrectedPvalues[non_na_indices] <- All_adj_p_values
SigCorrected <- CorrectedPvalues[CorrectedPvalues<.05]
CorrectedPvalues[CorrectedPvalues<.05]

BFvalues <- unlist(c(MS1_Unpleasantness$BF_values,MS2_Unpleasantness$BF_values,MS3_Unpleasantness$BF_values,MS4_Unpleasantness$BF_values,MS5_Unpleasantness$BF_values))

# Save the data
Pvalues_df <- data.frame(
  Name = names(Pvalues),
  Value = Pvalues,
  row.names = NULL
)
write_xlsx(Pvalues_df,paste0(Path,runs[thisrun],"/PPUnpleasantnessUncorrected.xlsx"))


CorrectedPvalues_df <- data.frame(
  Name = names(CorrectedPvalues),
  Value = CorrectedPvalues,
  row.names = NULL
)
write_xlsx(CorrectedPvalues_df,paste0(Path,runs[thisrun],"/PPUnpleasantnessCorrected.xlsx"))

BFvalues_df <- data.frame(
  Name = names(BFvalues),
  Value = BFvalues,
  row.names = NULL
)
write_xlsx(BFvalues_df,paste0(Path,runs[thisrun],"/PPUnpleasantnessBFvals.xlsx"))

SigCorrected_df <- data.frame(
  Name = names(SigCorrected),
  Value = SigCorrected,
  row.names = NULL
)
write_xlsx(SigCorrected_df,paste0(Path,runs[thisrun],"/PPUnpleasantnessSigCorrected.xlsx"))



# Plotting for Unpleasantness ----
dat = MS1RS1
# Duration
Fullmod <- lmer(UnpleasantnessRatings ~ UnpleasantnessT0 +  Duration*Arm + (1|Subject), data = dat,REML = FALSE)
summary(Fullmod)
Unpplot <- ggpredict(Fullmod,terms = c('Duration','Arm'))

Unpplot <- as.data.frame(Unpplot)
ggplot()+
  geom_point(data = dat,
             aes(x = Duration, y = UnpleasantnessRatings, colour = Arm),
             position=position_jitter(h=0.05, w=.8), alpha = .5) + scale_colour_manual(values = c("black","red")) + 
  geom_line(data = Unpplot,
            aes(x = x, y = predicted, colour = group),
            size = 1.2) + scale_fill_manual(values = c("black","red")) + 
  geom_ribbon(data = Unpplot, aes(x = x, ymin = conf.low, ymax = conf.high, fill = group), alpha= 0.2, show.legend = FALSE) +
  theme_blank() + mytheme + theme(legend.position = 'top') + xlab("Duration (s)") + ylab("Unpleasantness Ratings (NRS)")


# Occurance
Fullmod <- lmer(UnpleasantnessRatings ~ UnpleasantnessT0 +  Occurance*Arm + (1|Subject), data = dat,REML = FALSE)
summary(Fullmod)
Unpplot <- ggpredict(Fullmod,terms = c('Occurance','Arm'))

Unpplot <- as.data.frame(Unpplot)
ggplot()+
  geom_point(data = dat,
             aes(x = Occurance, y = UnpleasantnessRatings, colour = Arm),
             position=position_jitter(h=0.05, w=.05), alpha = .5) + scale_colour_manual(values = c("black","red")) + 
  geom_line(data = Unpplot,
            aes(x = x, y = predicted, colour = group),
            size = 1.2) + scale_fill_manual(values = c("black","red")) + 
  geom_ribbon(data = Unpplot, aes(x = x, ymin = conf.low, ymax = conf.high, fill = group), alpha= 0.2, show.legend = FALSE) +
  theme_blank() + mytheme + theme(legend.position = 'top') + xlab("Occurance/s") + ylab("Unpleasantness Ratings (NRS)")


# Coverage
Fullmod <- lmer(UnpleasantnessRatings ~ UnpleasantnessT0 +  Coverage*Arm + (1|Subject), data = dat,REML = FALSE)
summary(Fullmod)
Unpplot <- ggpredict(Fullmod,terms = c('Coverage','Arm'))

Unpplot <- as.data.frame(Unpplot)
ggplot()+
  geom_point(data = dat,
             aes(x = Coverage, y = UnpleasantnessRatings, colour = Arm),
             position=position_jitter(h=0.05, w=.05), alpha = .5) + scale_colour_manual(values = c("black","red")) + 
  geom_line(data = Unpplot,
            aes(x = x, y = predicted, colour = group),
            size = 1.2) + scale_fill_manual(values = c("black","red")) + 
  geom_ribbon(data = Unpplot, aes(x = x, ymin = conf.low, ymax = conf.high, fill = group), alpha= 0.2, show.legend = FALSE) +
  theme_blank() + mytheme + theme(legend.position = 'top') + xlab("Coverage (prop)") + ylab("Unpleasantness Ratings (NRS)")

# GEV
Fullmod <- lmer(UnpleasantnessRatings ~ UnpleasantnessT0 +  GEV*Arm + (1|Subject), data = dat,REML = FALSE)
summary(Fullmod)
Unpplot <- ggpredict(Fullmod,terms = c('GEV','Arm'))

Unpplot <- as.data.frame(Unpplot)
ggplot()+
  geom_point(data = dat,
             aes(x = GEV, y = UnpleasantnessRatings, colour = Arm),
             position=position_jitter(h=0.05, w=.05), alpha = .5) + scale_colour_manual(values = c("black","red")) + 
  geom_line(data = Unpplot,
            aes(x = x, y = predicted, colour = group),
            size = 1.2) + scale_fill_manual(values = c("black","red")) + 
  geom_ribbon(data = Unpplot, aes(x = x, ymin = conf.low, ymax = conf.high, fill = group), alpha= 0.2, show.legend = FALSE) +
  theme_blank() + mytheme + theme(legend.position = 'top') + xlab("GEV (prop)") + ylab("Unpleasantness Ratings (NRS)")



# TP1
Fullmod <- lmer(UnpleasantnessRatings ~ UnpleasantnessT0 +  TP1*Arm + (1|Subject), data = dat,REML = FALSE)
summary(Fullmod)
Unpplot <- ggpredict(Fullmod,terms = c('TP1','Arm'))

Unpplot <- as.data.frame(Unpplot)
ggplot()+
  geom_point(data = dat,
             aes(x = TP1, y = UnpleasantnessRatings, colour = Arm),
             position=position_jitter(h=0.05, w=.03), alpha = .5) + scale_colour_manual(values = c("black","red")) + 
  geom_line(data = Unpplot,
            aes(x = x, y = predicted, colour = group),
            size = 1.2) + scale_fill_manual(values = c("black","red")) + 
  geom_ribbon(data = Unpplot, aes(x = x, ymin = conf.low, ymax = conf.high, fill = group), alpha= 0.2, show.legend = FALSE) +
  theme_blank() + mytheme + theme(legend.position = 'top') + xlab("Transition Probability - MS1") + ylab("Unpleasantness Ratings (NRS)")


# TP2
Fullmod <- lmer(UnpleasantnessRatings ~ UnpleasantnessT0 +  TP2*Arm + (1|Subject), data = dat,REML = FALSE)
summary(Fullmod)
Unpplot <- ggpredict(Fullmod,terms = c('TP2','Arm'))

Unpplot <- as.data.frame(Unpplot)
ggplot()+
  geom_point(data = dat,
             aes(x = TP2, y = UnpleasantnessRatings, colour = Arm),
             position=position_jitter(h=0.05, w=.03), alpha = .5) + scale_colour_manual(values = c("black","red")) + 
  geom_line(data = Unpplot,
            aes(x = x, y = predicted, colour = group),
            size = 1.2) + scale_fill_manual(values = c("black","red")) + 
  geom_ribbon(data = Unpplot, aes(x = x, ymin = conf.low, ymax = conf.high, fill = group), alpha= 0.2, show.legend = FALSE) +
  theme_blank() + mytheme + theme(legend.position = 'top') + xlab("Transition Probability - MS2") + ylab("Unpleasantness Ratings (NRS)")


# TP3
Fullmod <- lmer(UnpleasantnessRatings ~ UnpleasantnessT0 +  TP3*Arm + (1|Subject), data = dat,REML = FALSE)
summary(Fullmod)
Unpplot <- ggpredict(Fullmod,terms = c('TP3','Arm'))

Unpplot <- as.data.frame(Unpplot)
ggplot()+
  geom_point(data = dat,
             aes(x = TP3, y = UnpleasantnessRatings, colour = Arm),
             position=position_jitter(h=0.05, w=.03), alpha = .5) + scale_colour_manual(values = c("black","red")) + 
  geom_line(data = Unpplot,
            aes(x = x, y = predicted, colour = group),
            size = 1.2) + scale_fill_manual(values = c("black","red")) + 
  geom_ribbon(data = Unpplot, aes(x = x, ymin = conf.low, ymax = conf.high, fill = group), alpha= 0.2, show.legend = FALSE) +
  theme_blank() + mytheme + theme(legend.position = 'top')+ xlab("Transition Probability - MS3") + ylab("Unpleasantness Ratings (NRS)")

# TP4
Fullmod <- lmer(UnpleasantnessRatings ~ UnpleasantnessT0 +  TP4*Arm + (1|Subject), data = dat,REML = FALSE)
summary(Fullmod)
Unpplot <- ggpredict(Fullmod,terms = c('TP4','Arm'))

Unpplot <- as.data.frame(Unpplot)
ggplot()+
  geom_point(data = dat,
             aes(x = TP4, y = UnpleasantnessRatings, colour = Arm),
             position=position_jitter(h=0.05, w=.03), alpha = .5) + scale_colour_manual(values = c("black","red")) + 
  geom_line(data = Unpplot,
            aes(x = x, y = predicted, colour = group),
            size = 1.2) + scale_fill_manual(values = c("black","red")) + 
  geom_ribbon(data = Unpplot, aes(x = x, ymin = conf.low, ymax = conf.high, fill = group), alpha= 0.2, show.legend = FALSE) +
  theme_blank() + mytheme + theme(legend.position = 'top')+ xlab("Transition Probability - MS4") + ylab("Unpleasantness Ratings (NRS)")


# TP5
Fullmod <- lmer(UnpleasantnessRatings ~ UnpleasantnessT0 +  TP5*Arm + (1|Subject), data = dat,REML = FALSE)
summary(Fullmod)
Unpplot <- ggpredict(Fullmod,terms = c('TP5','Arm'))

Unpplot <- as.data.frame(Unpplot)
ggplot()+
  geom_point(data = dat,
             aes(x = TP5, y = UnpleasantnessRatings, colour = Arm),
             position=position_jitter(h=0.05, w=.03), alpha = .5) + scale_colour_manual(values = c("black","red")) + 
  geom_line(data = Unpplot,
            aes(x = x, y = predicted, colour = group),
            size = 1.2) + scale_fill_manual(values = c("black","red")) + 
  geom_ribbon(data = Unpplot, aes(x = x, ymin = conf.low, ymax = conf.high, fill = group), alpha= 0.2, show.legend = FALSE) +
  theme_blank() + mytheme + theme(legend.position = 'top')+ xlab("Transition Probability - MS5") + ylab("Unpleasantness Ratings (NRS)")



## Plotting for RS1 vs RS2 ----

Microstats <- read_excel("Path_to_Microstate_results/Microstats.xlsx")


Microstats$MS[Microstats$MS == 'MS1'] <- "D"
Microstats$MS[Microstats$MS == 'MS2'] <- "B"
Microstats$MS[Microstats$MS == 'MS3'] <- "E"
Microstats$MS[Microstats$MS == 'MS4'] <- "C"
Microstats$MS[Microstats$MS == 'MS5'] <- "A"

data = Microstats

axsize = 20
mytheme = theme(
  axis.title.x = element_text(size = axsize),
  axis.text.x = element_text(size = axsize-2),
  axis.title.y = element_text(size = axsize),
  axis.text.y = element_text(size = axsize-2),
  plot.title = element_text(size=25,hjust = 0.5),
  legend.text=element_text(size=15),
  legend.title=element_text(size=20))

plot <- ggplot(data, aes(x=MS, y=Coverage, fill=RS)) +
  geom_violin(scale="width", trim=FALSE,alpha = .6,position=position_dodge(.7),width=0.5) + scale_fill_manual(values = c("black","blue")) + 
  geom_boxplot(width=0.2, position=position_dodge(0.7),alpha = .1) +
  #geom_jitter(shape=21, position=position_dodge(0.7),alpha = .3) +
  labs(y="Coverage %") + ylim(c(0, .6)) +
  theme_blank() + mytheme
plot

plot <- ggplot(data, aes(x=MS, y=Duration, fill=RS)) +
  geom_violin(scale="width", trim=FALSE,alpha = .6,position=position_dodge(.7),width=0.5) + scale_fill_manual(values = c("black","blue")) + 
  geom_boxplot(width=0.2, position=position_dodge(0.7),alpha = .1) +
  #geom_jitter(shape=21, position=position_dodge(0.7),alpha = .3) +
  labs(y="Duration (ms)") + ylim(c(50, 130)) +
  theme_blank() + mytheme
plot

plot <- ggplot(data, aes(x=MS, y=Occurance, fill=RS)) +
  geom_violin(scale="width", trim=FALSE,alpha = .6,position=position_dodge(.7),width=0.5) + scale_fill_manual(values = c("black","blue")) + 
  geom_boxplot(width=0.2, position=position_dodge(0.7),alpha = .1) +
  #geom_jitter(shape=21, position=position_dodge(0.7),alpha = .3) +
  labs(y="Occurrence (ms/s)") + ylim(c(0, 6)) +
  theme_blank() + mytheme
plot

plot <- ggplot(data, aes(x=MS, y=GEV, fill=RS)) +
  geom_violin(scale="width", trim=FALSE,alpha = .6,position=position_dodge(.7),width=0.5) + scale_fill_manual(values = c("black","blue")) + 
  geom_boxplot(width=0.2, position=position_dodge(0.7),alpha = .1) +
  #geom_jitter(shape=21, position=position_dodge(0.7),alpha = .3) +
  labs(y="GEV (prop)") + ylim(c(0, .5)) +
  theme_blank() + mytheme
plot

plot <- ggplot(data, aes(x=RS, y=TotalGEV)) +
  geom_violin(scale="width", trim=FALSE,alpha = .6,position=position_dodge(.7),width=0.5,fill="dark green") + 
  geom_boxplot(width=0.2, position=position_dodge(0.7),alpha = .1) +
  #geom_jitter(shape=21, position=position_dodge(0.7),alpha = .3) +
  labs(y="Total GEV") + ylim(c(.5, 1)) +
  theme_blank() + mytheme
plot


