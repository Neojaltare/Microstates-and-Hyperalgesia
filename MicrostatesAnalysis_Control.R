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




# Load Single Clustering data (All) ----
# Here we do RS1 vs RS2 
Microstats <- read_excel("Path_to_Microstate_results/Microstats.xlsx")
View(Microstats)

Durmod <- lm(Duration ~ MS*RS, data = Microstats)
summary(Durmod)

Occurmod <- lm((Occurance) ~ MS*RS, data = Microstats)
summary(Occurmod)

hist((Microstats$Coverage))
Covermod <- lm((Coverage) ~ MS*RS, data = Microstats)
summary(Covermod)

hist((MS1$GEV))
Gevmod <- lm(GEV ~ MS*RS, data = Microstats)
summary(Gevmod)

hist((MS1$TP2))
TP2mod <- lm(TP2 ~ MS*RS, data = Microstats)
summary(TP2mod)

TP3mod <- lm(TP3 ~ MS*RS, data = Microstats)
summary(TP3mod)

TP4mod <- lm((TP4) ~ MS*RS, data = Microstats)
summary(TP4mod)

TP5mod <- lm((TP5) ~ MS*RS, data = Microstats)
summary(TP5mod)

# Split microstates
MS1 = Microstats[Microstats$MS == 'MS1',]
MS2 = Microstats[Microstats$MS == 'MS2',]
MS3 = Microstats[Microstats$MS == 'MS3',]
MS4 = Microstats[Microstats$MS == 'MS4',]
MS5 = Microstats[Microstats$MS == 'MS5',]

# Create column and row names
col_names <- c("Duration", "Occurrence", "Coverage", "GEV", "TP1", "TP2", "TP3", "TP4", "TP5")
row_names <- c("MS1", "MS2", "MS3", "MS4", "MS5")
# Create empty data frame
ALL_Pvals <- data.frame(matrix(ncol = length(col_names), nrow = length(row_names)))
colnames(ALL_Pvals) <- col_names
rownames(ALL_Pvals) <- row_names

# MS1
hist(log(MS1$Duration))
Durmod <- lmer(log(Duration) ~ Condition + RS + (1|Subject), data = MS1,REML = FALSE)
summary(Durmod)
ALL_Pvals$Duration[1] = anova(Durmod)$`Pr(>F)`[2]

hist(log10(MS1$Occurance))
Occurmod <- lmer(log(Occurance) ~ Condition + RS + (1|Subject), data = MS1,REML = FALSE)
summary(Occurmod)
ALL_Pvals$Occurrence[1] = anova(Occurmod)$`Pr(>F)`[2]

hist(log10(MS1$Coverage))
Covermod <- lmer(log(Coverage) ~ Condition + RS + (1|Subject), data = MS1,REML = FALSE)
summary(Covermod)
ALL_Pvals$Coverage[1] = anova(Covermod)$`Pr(>F)`[2]


hist(log(MS1$GEV))
Gevmod <- lmer(log(GEV) ~ Condition + RS + (1|Subject), data = MS1,REML = FALSE)
summary(Gevmod)
ALL_Pvals$GEV[1] = anova(Gevmod)$`Pr(>F)`[2]


hist(log10(MS1$TP2))
TP2mod <- lmer(log(TP2) ~ Condition + RS + (1|Subject), data = MS1,REML = FALSE)
summary(TP2mod)
ALL_Pvals$TP2[1] = anova(TP2mod)$`Pr(>F)`[2]


TP3mod <- lmer(log(TP3) ~ Condition + RS + (1|Subject), data = MS1,REML = FALSE)
summary(TP3mod)
ALL_Pvals$TP3[1] = anova(TP3mod)$`Pr(>F)`[2]


TP4mod <- lmer(log10(TP4) ~ Condition + RS + (1|Subject), data = MS1,REML = FALSE)
summary(TP4mod)
ALL_Pvals$TP4[1] = anova(TP4mod)$`Pr(>F)`[2]


TP5mod <- lmer(log(TP5) ~ Condition + RS + (1|Subject), data = MS1,REML = FALSE)
summary(TP5mod)
ALL_Pvals$TP5[1] = anova(TP5mod)$`Pr(>F)`[2]


# MS2
Durmod <- lmer(log(Duration) ~ Condition + RS + (1|Subject), data = MS2,REML = FALSE)
summary(Durmod)

Occurmod <- lmer(log(Occurance) ~ Condition + RS + (1|Subject), data = MS2,REML = FALSE)
summary(Occurmod)

Covermod <- lmer(log(Coverage) ~ Condition + RS + (1|Subject), data = MS2,REML = FALSE)
summary(Covermod)

Gevmod <- lmer(log(GEV) ~ Condition + RS + (1|Subject), data = MS2,REML = FALSE)
summary(Gevmod)

TP1mod <- lmer(log(TP1) ~ Condition + RS + (1|Subject), data = MS2,REML = FALSE)
summary(TP1mod)

TP3mod <- lmer(log(TP3) ~ Condition + RS + (1|Subject), data = MS2,REML = FALSE)
summary(TP3mod)

TP4mod <- lmer(log(TP4) ~ Condition + RS + (1|Subject), data = MS2,REML = FALSE)
summary(TP4mod)

TP5mod <- lmer(log(TP5) ~ Condition + RS + (1|Subject), data = MS2,REML = FALSE)
summary(TP5mod)

ALL_Pvals$Duration[2] = anova(Durmod)$`Pr(>F)`[2]
ALL_Pvals$Occurrence[2] = anova(Occurmod)$`Pr(>F)`[2]
ALL_Pvals$Coverage[2] = anova(Covermod)$`Pr(>F)`[2]
ALL_Pvals$GEV[2] = anova(Gevmod)$`Pr(>F)`[2]
ALL_Pvals$TP1[2] = anova(TP1mod)$`Pr(>F)`[2]
ALL_Pvals$TP3[2] = anova(TP3mod)$`Pr(>F)`[2]
ALL_Pvals$TP4[2] = anova(TP4mod)$`Pr(>F)`[2]
ALL_Pvals$TP5[2] = anova(TP5mod)$`Pr(>F)`[2]

# MS3
Durmod <- lmer(log(Duration) ~ Condition + RS + (1|Subject), data = MS3,REML = FALSE)
summary(Durmod)

Occurmod <- lmer(log(Occurance) ~ Condition + RS + (1|Subject), data = MS3,REML = FALSE)
summary(Occurmod)

Covermod <- lmer(log(Coverage) ~ Condition + RS + (1|Subject), data = MS3,REML = FALSE)
summary(Covermod)

Gevmod <- lmer(log(GEV) ~ Condition + RS + (1|Subject), data = MS3,REML = FALSE)
summary(Gevmod)

TP1mod <- lmer(log(TP1) ~ Condition + RS + (1|Subject), data = MS3,REML = FALSE)
summary(TP1mod)

TP2mod <- lmer(log(TP2) ~ Condition + RS + (1|Subject), data = MS3,REML = FALSE)
summary(TP2mod)

TP4mod <- lmer(log(TP4) ~ Condition + RS + (1|Subject), data = MS3,REML = FALSE)
summary(TP4mod)

TP5mod <- lmer(log(TP5) ~ Condition + RS + (1|Subject), data = MS3,REML = FALSE)
summary(TP5mod)


ALL_Pvals$Duration[3] = anova(Durmod)$`Pr(>F)`[2]
ALL_Pvals$Occurrence[3] = anova(Occurmod)$`Pr(>F)`[2]
ALL_Pvals$Coverage[3] = anova(Covermod)$`Pr(>F)`[2]
ALL_Pvals$GEV[3] = anova(Gevmod)$`Pr(>F)`[2]
ALL_Pvals$TP1[3] = anova(TP1mod)$`Pr(>F)`[2]
ALL_Pvals$TP2[3] = anova(TP2mod)$`Pr(>F)`[2]
ALL_Pvals$TP4[3] = anova(TP4mod)$`Pr(>F)`[2]
ALL_Pvals$TP5[3] = anova(TP5mod)$`Pr(>F)`[2]

# MS4
Durmod <- lmer(log(Duration) ~ Condition + RS + (1|Subject), data = MS4,REML = FALSE)
summary(Durmod)

Occurmod <- lmer(log(Occurance) ~ Condition + RS + (1|Subject), data = MS4,REML = FALSE)
summary(Occurmod)

Covermod <- lmer(log(Coverage) ~ Condition + RS + (1|Subject), data = MS4,REML = FALSE)
summary(Covermod)

hist(MS4$GEV)
Gevmod <- lmer(log(GEV) ~ Condition + RS + (1|Subject), data = MS4,REML = FALSE)
summary(Gevmod)

TP1mod <- lmer(log(TP1) ~ Condition + RS + (1|Subject), data = MS4,REML = FALSE)
summary(TP1mod)

TP2mod <- lmer(log(TP2) ~ Condition + RS + (1|Subject), data = MS4,REML = FALSE)
summary(TP2mod)

TP3mod <- lmer(log(TP3) ~ Condition + RS + (1|Subject), data = MS4,REML = FALSE)
summary(TP3mod)

TP5mod <- lmer(log(TP5) ~ Condition + RS + (1|Subject), data = MS4,REML = FALSE)
summary(TP5mod)

ALL_Pvals$Duration[4] = anova(Durmod)$`Pr(>F)`[2]
ALL_Pvals$Occurrence[4] = anova(Occurmod)$`Pr(>F)`[2]
ALL_Pvals$Coverage[4] = anova(Covermod)$`Pr(>F)`[2]
ALL_Pvals$GEV[4] = anova(Gevmod)$`Pr(>F)`[2]
ALL_Pvals$TP1[4] = anova(TP1mod)$`Pr(>F)`[2]
ALL_Pvals$TP2[4] = anova(TP2mod)$`Pr(>F)`[2]
ALL_Pvals$TP3[4] = anova(TP3mod)$`Pr(>F)`[2]
ALL_Pvals$TP5[4] = anova(TP5mod)$`Pr(>F)`[2]

# MS5
Durmod <- lmer(log(Duration) ~ Condition + RS + (1|Subject), data = MS5,REML = FALSE)
summary(Durmod)

Occurmod <- lmer(log(Occurance) ~ Condition + RS + (1|Subject), data = MS5,REML = FALSE)
summary(Occurmod)

Covermod <- lmer(log(Coverage) ~ Condition + RS + (1|Subject), data = MS5,REML = FALSE)
summary(Covermod)

Gevmod <- lmer(log(GEV) ~ Condition + RS + (1|Subject), data = MS5,REML = FALSE)
summary(Gevmod)

TP1mod <- lmer(log(TP1) ~ Condition + RS + (1|Subject), data = MS5,REML = FALSE)
summary(TP1mod)

TP2mod <- lmer(log(TP2) ~ Condition + RS + (1|Subject), data = MS5,REML = FALSE)
summary(TP2mod)

TP3mod <- lmer(log(TP3) ~ Condition + RS + (1|Subject), data = MS5,REML = FALSE)
summary(TP3mod)

TP4mod <- lmer(log(TP4) ~ Condition + RS + (1|Subject), data = MS5,REML = FALSE)
summary(TP4mod)

ALL_Pvals$Duration[5] = anova(Durmod)$`Pr(>F)`[2]
ALL_Pvals$Occurrence[5] = anova(Occurmod)$`Pr(>F)`[2]
ALL_Pvals$Coverage[5] = anova(Covermod)$`Pr(>F)`[2]
ALL_Pvals$GEV[5] = anova(Gevmod)$`Pr(>F)`[2]
ALL_Pvals$TP1[5] = anova(TP1mod)$`Pr(>F)`[2]
ALL_Pvals$TP2[5] = anova(TP2mod)$`Pr(>F)`[2]
ALL_Pvals$TP3[5] = anova(TP3mod)$`Pr(>F)`[2]
ALL_Pvals$TP4[5] = anova(TP4mod)$`Pr(>F)`[2]


Pvalues <- unlist(ALL_Pvals)
non_na_indices <- which(!is.na(Pvalues))
All_adj_p_values <- p.adjust(Pvalues[non_na_indices], method = "fdr")
Pvalues[non_na_indices] <- All_adj_p_values
Pvalues[Pvalues<.05]


## Microstates x Area ----
# Do microstates predict area
MicrostatsxArea <- read_excel("Path_to_Microstate_Area_data/MicrostatsxArea.xlsx")
View(MicrostatsxArea)

MS1RS1 = MicrostatsxArea[MicrostatsxArea$MS == 'MS1' & MicrostatsxArea$RS == 'RS1',]
MS2RS1 = MicrostatsxArea[MicrostatsxArea$MS == 'MS2' & MicrostatsxArea$RS == 'RS1',]
MS3RS1 = MicrostatsxArea[MicrostatsxArea$MS == 'MS3' & MicrostatsxArea$RS == 'RS1',]
MS4RS1 = MicrostatsxArea[MicrostatsxArea$MS == 'MS4' & MicrostatsxArea$RS == 'RS1',]
MS5RS1 = MicrostatsxArea[MicrostatsxArea$MS == 'MS5' & MicrostatsxArea$RS == 'RS1',]

MS1RS2 = MicrostatsxArea[MicrostatsxArea$MS == 'MS1' & MicrostatsxArea$RS == 'RS2',]
MS2RS2 = MicrostatsxArea[MicrostatsxArea$MS == 'MS2' & MicrostatsxArea$RS == 'RS2',]
MS3RS2 = MicrostatsxArea[MicrostatsxArea$MS == 'MS3' & MicrostatsxArea$RS == 'RS2',]
MS4RS2 = MicrostatsxArea[MicrostatsxArea$MS == 'MS4' & MicrostatsxArea$RS == 'RS2',]
MS5RS2 = MicrostatsxArea[MicrostatsxArea$MS == 'MS5' & MicrostatsxArea$RS == 'RS2',]

MS1diff = MS1RS1
MS1diff[,c(6:16)] = MS1RS1[,c(6:16)] - MS1RS2[,c(6:16)]
MS2diff = MS2RS1
MS2diff[,c(6:16)] = MS2RS1[,c(6:16)] - MS2RS2[,c(6:16)]
MS3diff = MS3RS1
MS3diff[,c(6:16)] = MS3RS1[,c(6:16)] - MS3RS2[,c(6:16)]
MS4diff = MS4RS1
MS4diff[,c(6:16)] = MS4RS1[,c(6:16)] - MS4RS2[,c(6:16)]
MS5diff = MS5RS1
MS5diff[,c(6:16)] = MS5RS1[,c(6:16)] - MS5RS2[,c(6:16)]

### Length ----
# Initialize variables

# Create column and row names
col_names <- c("RS1","RS2","Diff")
row_names <- c("MS1", "MS2", "MS3", "MS4", "MS5")
# Create empty data frames
Duration_Pvals <- data.frame(matrix(ncol = length(col_names), nrow = length(row_names)))
colnames(Duration_Pvals) <- col_names
rownames(Duration_Pvals) <- row_names
Occurance_Pvals = Duration_Pvals
Coverage_Pvals = Duration_Pvals
GEV_Pvals = Duration_Pvals
TP1_Pvals = Duration_Pvals
TP2_Pvals = Duration_Pvals
TP3_Pvals = Duration_Pvals
TP4_Pvals = Duration_Pvals
TP5_Pvals = Duration_Pvals


## MS1

# RS1
FullmodMS1 <- lm(Length ~ Condition + Duration + Occurance + Coverage + GEV + TP2 + TP3 + TP4 + TP5, data = MS1RS1)
summary(FullmodMS1)
anova(FullmodMS1)
Duration_Pvals["MS1","RS1"]  = anova(FullmodMS1)['Duration','Pr(>F)']
Occurance_Pvals["MS1","RS1"] = anova(FullmodMS1)['Occurance','Pr(>F)']
Coverage_Pvals["MS1","RS1"]  = anova(FullmodMS1)['Coverage','Pr(>F)']
GEV_Pvals["MS1","RS1"]       = anova(FullmodMS1)['GEV','Pr(>F)']
TP1_Pvals["MS1","RS1"]       = anova(FullmodMS1)['TP1','Pr(>F)']
TP2_Pvals["MS1","RS1"]       = anova(FullmodMS1)['TP2','Pr(>F)']
TP3_Pvals["MS1","RS1"]       = anova(FullmodMS1)['TP3','Pr(>F)']
TP4_Pvals["MS1","RS1"]       = anova(FullmodMS1)['TP4','Pr(>F)']
TP5_Pvals["MS1","RS1"]       = anova(FullmodMS1)['TP5','Pr(>F)']

# RS2
FullmodMS1 <- lm(Length ~ Condition + Duration + Occurance + Coverage + GEV + TP2 + TP3 + TP4 + TP5, data = MS1RS2)
summary(FullmodMS1)
anova(FullmodMS1)
Duration_Pvals["MS1","RS2"]  = anova(FullmodMS1)['Duration','Pr(>F)']
Occurance_Pvals["MS1","RS2"] = anova(FullmodMS1)['Occurance','Pr(>F)']
Coverage_Pvals["MS1","RS2"]  = anova(FullmodMS1)['Coverage','Pr(>F)']
GEV_Pvals["MS1","RS2"]       = anova(FullmodMS1)['GEV','Pr(>F)']
TP1_Pvals["MS1","RS2"]       = anova(FullmodMS1)['TP1','Pr(>F)']
TP2_Pvals["MS1","RS2"]       = anova(FullmodMS1)['TP2','Pr(>F)']
TP3_Pvals["MS1","RS2"]       = anova(FullmodMS1)['TP3','Pr(>F)']
TP4_Pvals["MS1","RS2"]       = anova(FullmodMS1)['TP4','Pr(>F)']
TP5_Pvals["MS1","RS2"]       = anova(FullmodMS1)['TP5','Pr(>F)']

# diff
DiffmodMS1 <- lm(Length ~ Condition + Duration + Occurance + Coverage + GEV + TP2 + TP3 + TP4 + TP5, data = MS1diff)
summary(DiffmodMS1)
anova(DiffmodMS1)
Duration_Pvals["MS1","Diff"]  = anova(DiffmodMS1)['Duration','Pr(>F)']
Occurance_Pvals["MS1","Diff"] = anova(DiffmodMS1)['Occurance','Pr(>F)']
Coverage_Pvals["MS1","Diff"]  = anova(DiffmodMS1)['Coverage','Pr(>F)']
GEV_Pvals["MS1","Diff"]       = anova(DiffmodMS1)['GEV','Pr(>F)']
TP1_Pvals["MS1","Diff"]       = anova(DiffmodMS1)['TP1','Pr(>F)']
TP2_Pvals["MS1","Diff"]       = anova(DiffmodMS1)['TP2','Pr(>F)']
TP3_Pvals["MS1","Diff"]       = anova(DiffmodMS1)['TP3','Pr(>F)']
TP4_Pvals["MS1","Diff"]       = anova(DiffmodMS1)['TP4','Pr(>F)']
TP5_Pvals["MS1","Diff"]       = anova(DiffmodMS1)['TP5','Pr(>F)']


## MS2
# MS2RS1
Fullmod <- lm(Length ~ Condition + Duration + Occurance + Coverage + GEV + TP1 + TP3 + TP4 + TP5, data = MS2RS1)
summary(Fullmod)
anova(Fullmod)
Duration_Pvals["MS2","RS1"]  = anova(Fullmod)['Duration','Pr(>F)']
Occurance_Pvals["MS2","RS1"] = anova(Fullmod)['Occurance','Pr(>F)']
Coverage_Pvals["MS2","RS1"]  = anova(Fullmod)['Coverage','Pr(>F)']
GEV_Pvals["MS2","RS1"]       = anova(Fullmod)['GEV','Pr(>F)']
TP1_Pvals["MS2","RS1"]       = anova(Fullmod)['TP1','Pr(>F)']
TP2_Pvals["MS2","RS1"]       = anova(Fullmod)['TP2','Pr(>F)']
TP3_Pvals["MS2","RS1"]       = anova(Fullmod)['TP3','Pr(>F)']
TP4_Pvals["MS2","RS1"]       = anova(Fullmod)['TP4','Pr(>F)']
TP5_Pvals["MS2","RS1"]       = anova(Fullmod)['TP5','Pr(>F)']

# MS2RS2
Fullmod <- lm(Length ~ Condition + Duration + Occurance + Coverage + GEV + TP1 + TP3 + TP4 + TP5, data = MS2RS2)
summary(Fullmod)
anova(Fullmod)
Duration_Pvals["MS2","RS2"]  = anova(Fullmod)['Duration','Pr(>F)']
Occurance_Pvals["MS2","RS2"] = anova(Fullmod)['Occurance','Pr(>F)']
Coverage_Pvals["MS2","RS2"]  = anova(Fullmod)['Coverage','Pr(>F)']
GEV_Pvals["MS2","RS2"]       = anova(Fullmod)['GEV','Pr(>F)']
TP1_Pvals["MS2","RS2"]       = anova(Fullmod)['TP1','Pr(>F)']
TP2_Pvals["MS2","RS2"]       = anova(Fullmod)['TP2','Pr(>F)']
TP3_Pvals["MS2","RS2"]       = anova(Fullmod)['TP3','Pr(>F)']
TP4_Pvals["MS2","RS2"]       = anova(Fullmod)['TP4','Pr(>F)']
TP5_Pvals["MS2","RS2"]       = anova(Fullmod)['TP5','Pr(>F)']

# diff
Diffmod <- lm(Length ~ Condition + Duration + Occurance + Coverage + GEV + TP1 + TP3 + TP4 + TP5, data = MS2diff)
summary(Diffmod)
anova(Diffmod)
Duration_Pvals["MS2","Diff"]  = anova(Diffmod)['Duration','Pr(>F)']
Occurance_Pvals["MS2","Diff"] = anova(Diffmod)['Occurance','Pr(>F)']
Coverage_Pvals["MS2","Diff"]  = anova(Diffmod)['Coverage','Pr(>F)']
GEV_Pvals["MS2","Diff"]       = anova(Diffmod)['GEV','Pr(>F)']
TP1_Pvals["MS2","Diff"]       = anova(Diffmod)['TP1','Pr(>F)']
TP2_Pvals["MS2","Diff"]       = anova(Diffmod)['TP2','Pr(>F)']
TP3_Pvals["MS2","Diff"]       = anova(Diffmod)['TP3','Pr(>F)']
TP4_Pvals["MS2","Diff"]       = anova(Diffmod)['TP4','Pr(>F)']
TP5_Pvals["MS2","Diff"]       = anova(Diffmod)['TP5','Pr(>F)']


## MS3
# MS3RS1
Fullmod <- lm(Length ~ Condition + Duration + Occurance + Coverage + GEV + TP1 + TP2 + TP4 + TP5, data = MS3RS1)
summary(Fullmod)
anova(Fullmod)
Duration_Pvals["MS3","RS1"]  = anova(Fullmod)['Duration','Pr(>F)']
Occurance_Pvals["MS3","RS1"] = anova(Fullmod)['Occurance','Pr(>F)']
Coverage_Pvals["MS3","RS1"]  = anova(Fullmod)['Coverage','Pr(>F)']
GEV_Pvals["MS3","RS1"]       = anova(Fullmod)['GEV','Pr(>F)']
TP1_Pvals["MS3","RS1"]       = anova(Fullmod)['TP1','Pr(>F)']
TP2_Pvals["MS3","RS1"]       = anova(Fullmod)['TP2','Pr(>F)']
TP3_Pvals["MS3","RS1"]       = anova(Fullmod)['TP3','Pr(>F)']
TP4_Pvals["MS3","RS1"]       = anova(Fullmod)['TP4','Pr(>F)']
TP5_Pvals["MS3","RS1"]       = anova(Fullmod)['TP5','Pr(>F)']

# MS3RS2
Fullmod <- lm(Length ~ Condition + Duration + Occurance + Coverage + GEV + TP1 + TP2 + TP4 + TP5, data = MS3RS2)
summary(Fullmod)
anova(Fullmod)
Duration_Pvals["MS3","RS2"]  = anova(Fullmod)['Duration','Pr(>F)']
Occurance_Pvals["MS3","RS2"] = anova(Fullmod)['Occurance','Pr(>F)']
Coverage_Pvals["MS3","RS2"]  = anova(Fullmod)['Coverage','Pr(>F)']
GEV_Pvals["MS3","RS2"]       = anova(Fullmod)['GEV','Pr(>F)']
TP1_Pvals["MS3","RS2"]       = anova(Fullmod)['TP1','Pr(>F)']
TP2_Pvals["MS3","RS2"]       = anova(Fullmod)['TP2','Pr(>F)']
TP3_Pvals["MS3","RS2"]       = anova(Fullmod)['TP3','Pr(>F)']
TP4_Pvals["MS3","RS2"]       = anova(Fullmod)['TP4','Pr(>F)']
TP5_Pvals["MS3","RS2"]       = anova(Fullmod)['TP5','Pr(>F)']

# diff
Diffmod <- lm(Length ~ Condition + Duration + Occurance + Coverage + GEV + TP1 + TP2 + TP4 + TP5, data = MS3diff)
summary(Diffmod)
anova(Diffmod)
Duration_Pvals["MS3","Diff"]  = anova(Diffmod)['Duration','Pr(>F)']
Occurance_Pvals["MS3","Diff"] = anova(Diffmod)['Occurance','Pr(>F)']
Coverage_Pvals["MS3","Diff"]  = anova(Diffmod)['Coverage','Pr(>F)']
GEV_Pvals["MS3","Diff"]       = anova(Diffmod)['GEV','Pr(>F)']
TP1_Pvals["MS3","Diff"]       = anova(Diffmod)['TP1','Pr(>F)']
TP2_Pvals["MS3","Diff"]       = anova(Diffmod)['TP2','Pr(>F)']
TP3_Pvals["MS3","Diff"]       = anova(Diffmod)['TP3','Pr(>F)']
TP4_Pvals["MS3","Diff"]       = anova(Diffmod)['TP4','Pr(>F)']
TP5_Pvals["MS3","Diff"]       = anova(Diffmod)['TP5','Pr(>F)']

## MS4
# MS4RS1
Fullmod <- lm(Length ~ Condition + Duration + Occurance + Coverage + GEV + TP1 + TP2 + TP3 + TP5, data = MS4RS1)
summary(Fullmod)
anova(Fullmod)
Duration_Pvals["MS4","RS1"]  = anova(Fullmod)['Duration','Pr(>F)']
Occurance_Pvals["MS4","RS1"] = anova(Fullmod)['Occurance','Pr(>F)']
Coverage_Pvals["MS4","RS1"]  = anova(Fullmod)['Coverage','Pr(>F)']
GEV_Pvals["MS4","RS1"]       = anova(Fullmod)['GEV','Pr(>F)']
TP1_Pvals["MS4","RS1"]       = anova(Fullmod)['TP1','Pr(>F)']
TP2_Pvals["MS4","RS1"]       = anova(Fullmod)['TP2','Pr(>F)']
TP3_Pvals["MS4","RS1"]       = anova(Fullmod)['TP3','Pr(>F)']
TP4_Pvals["MS4","RS1"]       = anova(Fullmod)['TP4','Pr(>F)']
TP5_Pvals["MS4","RS1"]       = anova(Fullmod)['TP5','Pr(>F)']

# MS4RS2
Fullmod <- lm(Length ~ Condition + Duration + Occurance + Coverage + GEV + TP1 + TP2 + TP3 + TP5, data = MS4RS2)
summary(Fullmod)
anova(Fullmod)
Duration_Pvals["MS4","RS2"]  = anova(Fullmod)['Duration','Pr(>F)']
Occurance_Pvals["MS4","RS2"] = anova(Fullmod)['Occurance','Pr(>F)']
Coverage_Pvals["MS4","RS2"]  = anova(Fullmod)['Coverage','Pr(>F)']
GEV_Pvals["MS4","RS2"]       = anova(Fullmod)['GEV','Pr(>F)']
TP1_Pvals["MS4","RS2"]       = anova(Fullmod)['TP1','Pr(>F)']
TP2_Pvals["MS4","RS2"]       = anova(Fullmod)['TP2','Pr(>F)']
TP3_Pvals["MS4","RS2"]       = anova(Fullmod)['TP3','Pr(>F)']
TP4_Pvals["MS4","RS2"]       = anova(Fullmod)['TP4','Pr(>F)']
TP5_Pvals["MS4","RS2"]       = anova(Fullmod)['TP5','Pr(>F)']

# diff
Diffmod <- lm(Length ~ Condition + Duration + Occurance + Coverage + GEV + TP1 + TP2 + TP3 + TP5, data = MS4diff)
summary(Diffmod)
anova(Diffmod)
Duration_Pvals["MS4","Diff"]  = anova(Diffmod)['Duration','Pr(>F)']
Occurance_Pvals["MS4","Diff"] = anova(Diffmod)['Occurance','Pr(>F)']
Coverage_Pvals["MS4","Diff"]  = anova(Diffmod)['Coverage','Pr(>F)']
GEV_Pvals["MS4","Diff"]       = anova(Diffmod)['GEV','Pr(>F)']
TP1_Pvals["MS4","Diff"]       = anova(Diffmod)['TP1','Pr(>F)']
TP2_Pvals["MS4","Diff"]       = anova(Diffmod)['TP2','Pr(>F)']
TP3_Pvals["MS4","Diff"]       = anova(Diffmod)['TP3','Pr(>F)']
TP4_Pvals["MS4","Diff"]       = anova(Diffmod)['TP4','Pr(>F)']
TP5_Pvals["MS4","Diff"]       = anova(Diffmod)['TP5','Pr(>F)']

## MS5
# MS5RS1
Fullmod <- lm(Length ~ Condition + Duration + Occurance + Coverage + GEV + TP1 + TP2 + TP3 + TP4, data = MS5RS1)
summary(Fullmod)
anova(Fullmod)
Duration_Pvals["MS5","RS1"]  = anova(Fullmod)['Duration','Pr(>F)']
Occurance_Pvals["MS5","RS1"] = anova(Fullmod)['Occurance','Pr(>F)']
Coverage_Pvals["MS5","RS1"]  = anova(Fullmod)['Coverage','Pr(>F)']
GEV_Pvals["MS5","RS1"]       = anova(Fullmod)['GEV','Pr(>F)']
TP1_Pvals["MS5","RS1"]       = anova(Fullmod)['TP1','Pr(>F)']
TP2_Pvals["MS5","RS1"]       = anova(Fullmod)['TP2','Pr(>F)']
TP3_Pvals["MS5","RS1"]       = anova(Fullmod)['TP3','Pr(>F)']
TP4_Pvals["MS5","RS1"]       = anova(Fullmod)['TP4','Pr(>F)']
TP5_Pvals["MS5","RS1"]       = anova(Fullmod)['TP5','Pr(>F)']

# MS5RS2
Fullmod <- lm(Length ~ Condition + Duration + Occurance + Coverage + GEV + TP1 + TP2 + TP3 + TP4, data = MS5RS2)
summary(Fullmod)
anova(Fullmod)
Duration_Pvals["MS5","RS2"]  = anova(Fullmod)['Duration','Pr(>F)']
Occurance_Pvals["MS5","RS2"] = anova(Fullmod)['Occurance','Pr(>F)']
Coverage_Pvals["MS5","RS2"]  = anova(Fullmod)['Coverage','Pr(>F)']
GEV_Pvals["MS5","RS2"]       = anova(Fullmod)['GEV','Pr(>F)']
TP1_Pvals["MS5","RS2"]       = anova(Fullmod)['TP1','Pr(>F)']
TP2_Pvals["MS5","RS2"]       = anova(Fullmod)['TP2','Pr(>F)']
TP3_Pvals["MS5","RS2"]       = anova(Fullmod)['TP3','Pr(>F)']
TP4_Pvals["MS5","RS2"]       = anova(Fullmod)['TP4','Pr(>F)']
TP5_Pvals["MS5","RS2"]       = anova(Fullmod)['TP5','Pr(>F)']

# diff
Diffmod <- lm(Length ~ Condition + Duration + Occurance + Coverage + GEV + TP1 + TP2 + TP3 + TP4, data = MS5diff)
summary(Diffmod)
anova(Diffmod)
Duration_Pvals["MS5","Diff"]  = anova(Diffmod)['Duration','Pr(>F)']
Occurance_Pvals["MS5","Diff"] = anova(Diffmod)['Occurance','Pr(>F)']
Coverage_Pvals["MS5","Diff"]  = anova(Diffmod)['Coverage','Pr(>F)']
GEV_Pvals["MS5","Diff"]       = anova(Diffmod)['GEV','Pr(>F)']
TP1_Pvals["MS5","Diff"]       = anova(Diffmod)['TP1','Pr(>F)']
TP2_Pvals["MS5","Diff"]       = anova(Diffmod)['TP2','Pr(>F)']
TP3_Pvals["MS5","Diff"]       = anova(Diffmod)['TP3','Pr(>F)']
TP4_Pvals["MS5","Diff"]       = anova(Diffmod)['TP4','Pr(>F)']
TP5_Pvals["MS5","Diff"]       = anova(Diffmod)['TP5','Pr(>F)']

# Check for significance after FDR correction
Pvalues <- unlist(TP5_Pvals)
non_na_indices <- which(!is.na(Pvalues))
All_adj_p_values <- p.adjust(Pvalues[non_na_indices], method = "fdr")
Pvalues[non_na_indices] <- All_adj_p_values
Pvalues[Pvalues<.05]


### Width ----
# Initialize variables

# Create column and row names
col_names <- c("RS1","RS2","Diff")
row_names <- c("MS1", "MS2", "MS3", "MS4", "MS5")
# Create empty data frames
Duration_Pvals <- data.frame(matrix(ncol = length(col_names), nrow = length(row_names)))
colnames(Duration_Pvals) <- col_names
rownames(Duration_Pvals) <- row_names
Occurance_Pvals = Duration_Pvals
Coverage_Pvals = Duration_Pvals
GEV_Pvals = Duration_Pvals
TP1_Pvals = Duration_Pvals
TP2_Pvals = Duration_Pvals
TP3_Pvals = Duration_Pvals
TP4_Pvals = Duration_Pvals
TP5_Pvals = Duration_Pvals

## MS1

# RS1
FullmodMS1 <- lm(Width ~ Condition + Duration + Occurance + Coverage + GEV + TP2 + TP3 + TP4 + TP5, data = MS1RS1)
summary(FullmodMS1)
anova(FullmodMS1)
Duration_Pvals["MS1","RS1"]  = anova(FullmodMS1)['Duration','Pr(>F)']
Occurance_Pvals["MS1","RS1"] = anova(FullmodMS1)['Occurance','Pr(>F)']
Coverage_Pvals["MS1","RS1"]  = anova(FullmodMS1)['Coverage','Pr(>F)']
GEV_Pvals["MS1","RS1"]       = anova(FullmodMS1)['GEV','Pr(>F)']
TP1_Pvals["MS1","RS1"]       = anova(FullmodMS1)['TP1','Pr(>F)']
TP2_Pvals["MS1","RS1"]       = anova(FullmodMS1)['TP2','Pr(>F)']
TP3_Pvals["MS1","RS1"]       = anova(FullmodMS1)['TP3','Pr(>F)']
TP4_Pvals["MS1","RS1"]       = anova(FullmodMS1)['TP4','Pr(>F)']
TP5_Pvals["MS1","RS1"]       = anova(FullmodMS1)['TP5','Pr(>F)']

# RS2
FullmodMS1 <- lm(Width ~ Condition + Duration + Occurance + Coverage + GEV + TP2 + TP3 + TP4 + TP5, data = MS1RS2)
summary(FullmodMS1)
anova(FullmodMS1)
Duration_Pvals["MS1","RS2"]  = anova(FullmodMS1)['Duration','Pr(>F)']
Occurance_Pvals["MS1","RS2"] = anova(FullmodMS1)['Occurance','Pr(>F)']
Coverage_Pvals["MS1","RS2"]  = anova(FullmodMS1)['Coverage','Pr(>F)']
GEV_Pvals["MS1","RS2"]       = anova(FullmodMS1)['GEV','Pr(>F)']
TP1_Pvals["MS1","RS2"]       = anova(FullmodMS1)['TP1','Pr(>F)']
TP2_Pvals["MS1","RS2"]       = anova(FullmodMS1)['TP2','Pr(>F)']
TP3_Pvals["MS1","RS2"]       = anova(FullmodMS1)['TP3','Pr(>F)']
TP4_Pvals["MS1","RS2"]       = anova(FullmodMS1)['TP4','Pr(>F)']
TP5_Pvals["MS1","RS2"]       = anova(FullmodMS1)['TP5','Pr(>F)']

# diff
DiffmodMS1 <- lm(Width ~ Condition + Duration + Occurance + Coverage + GEV + TP2 + TP3 + TP4 + TP5, data = MS1diff)
summary(DiffmodMS1)
anova(DiffmodMS1)
Duration_Pvals["MS1","Diff"]  = anova(DiffmodMS1)['Duration','Pr(>F)']
Occurance_Pvals["MS1","Diff"] = anova(DiffmodMS1)['Occurance','Pr(>F)']
Coverage_Pvals["MS1","Diff"]  = anova(DiffmodMS1)['Coverage','Pr(>F)']
GEV_Pvals["MS1","Diff"]       = anova(DiffmodMS1)['GEV','Pr(>F)']
TP1_Pvals["MS1","Diff"]       = anova(DiffmodMS1)['TP1','Pr(>F)']
TP2_Pvals["MS1","Diff"]       = anova(DiffmodMS1)['TP2','Pr(>F)']
TP3_Pvals["MS1","Diff"]       = anova(DiffmodMS1)['TP3','Pr(>F)']
TP4_Pvals["MS1","Diff"]       = anova(DiffmodMS1)['TP4','Pr(>F)']
TP5_Pvals["MS1","Diff"]       = anova(DiffmodMS1)['TP5','Pr(>F)']

## MS2
# MS2RS1
Fullmod <- lm(Width ~ Condition + Duration + Occurance + Coverage + GEV + TP1 + TP3 + TP4 + TP5, data = MS2RS1)
summary(Fullmod)
anova(Fullmod)
Duration_Pvals["MS2","RS1"]  = anova(Fullmod)['Duration','Pr(>F)']
Occurance_Pvals["MS2","RS1"] = anova(Fullmod)['Occurance','Pr(>F)']
Coverage_Pvals["MS2","RS1"]  = anova(Fullmod)['Coverage','Pr(>F)']
GEV_Pvals["MS2","RS1"]       = anova(Fullmod)['GEV','Pr(>F)']
TP1_Pvals["MS2","RS1"]       = anova(Fullmod)['TP1','Pr(>F)']
TP2_Pvals["MS2","RS1"]       = anova(Fullmod)['TP2','Pr(>F)']
TP3_Pvals["MS2","RS1"]       = anova(Fullmod)['TP3','Pr(>F)']
TP4_Pvals["MS2","RS1"]       = anova(Fullmod)['TP4','Pr(>F)']
TP5_Pvals["MS2","RS1"]       = anova(Fullmod)['TP5','Pr(>F)']

# MS2RS2
Fullmod <- lm(Width ~ Condition + Duration + Occurance + Coverage + GEV + TP1 + TP3 + TP4 + TP5, data = MS2RS2)
summary(Fullmod)
anova(Fullmod)
Duration_Pvals["MS2","RS2"]  = anova(Fullmod)['Duration','Pr(>F)']
Occurance_Pvals["MS2","RS2"] = anova(Fullmod)['Occurance','Pr(>F)']
Coverage_Pvals["MS2","RS2"]  = anova(Fullmod)['Coverage','Pr(>F)']
GEV_Pvals["MS2","RS2"]       = anova(Fullmod)['GEV','Pr(>F)']
TP1_Pvals["MS2","RS2"]       = anova(Fullmod)['TP1','Pr(>F)']
TP2_Pvals["MS2","RS2"]       = anova(Fullmod)['TP2','Pr(>F)']
TP3_Pvals["MS2","RS2"]       = anova(Fullmod)['TP3','Pr(>F)']
TP4_Pvals["MS2","RS2"]       = anova(Fullmod)['TP4','Pr(>F)']
TP5_Pvals["MS2","RS2"]       = anova(Fullmod)['TP5','Pr(>F)']

# diff
Diffmod <- lm(Width ~ Condition + Duration + Occurance + Coverage + GEV + TP1 + TP3 + TP4 + TP5, data = MS2diff)
summary(Diffmod)
anova(Diffmod)
Duration_Pvals["MS2","Diff"]  = anova(Diffmod)['Duration','Pr(>F)']
Occurance_Pvals["MS2","Diff"] = anova(Diffmod)['Occurance','Pr(>F)']
Coverage_Pvals["MS2","Diff"]  = anova(Diffmod)['Coverage','Pr(>F)']
GEV_Pvals["MS2","Diff"]       = anova(Diffmod)['GEV','Pr(>F)']
TP1_Pvals["MS2","Diff"]       = anova(Diffmod)['TP1','Pr(>F)']
TP2_Pvals["MS2","Diff"]       = anova(Diffmod)['TP2','Pr(>F)']
TP3_Pvals["MS2","Diff"]       = anova(Diffmod)['TP3','Pr(>F)']
TP4_Pvals["MS2","Diff"]       = anova(Diffmod)['TP4','Pr(>F)']
TP5_Pvals["MS2","Diff"]       = anova(Diffmod)['TP5','Pr(>F)']

## MS3
# MS3RS1
Fullmod <- lm(Width ~ Condition + Duration + Occurance + Coverage + GEV + TP1 + TP2 + TP4 + TP5, data = MS3RS1)
summary(Fullmod)
anova(Fullmod)
Duration_Pvals["MS3","RS1"]  = anova(Fullmod)['Duration','Pr(>F)']
Occurance_Pvals["MS3","RS1"] = anova(Fullmod)['Occurance','Pr(>F)']
Coverage_Pvals["MS3","RS1"]  = anova(Fullmod)['Coverage','Pr(>F)']
GEV_Pvals["MS3","RS1"]       = anova(Fullmod)['GEV','Pr(>F)']
TP1_Pvals["MS3","RS1"]       = anova(Fullmod)['TP1','Pr(>F)']
TP2_Pvals["MS3","RS1"]       = anova(Fullmod)['TP2','Pr(>F)']
TP3_Pvals["MS3","RS1"]       = anova(Fullmod)['TP3','Pr(>F)']
TP4_Pvals["MS3","RS1"]       = anova(Fullmod)['TP4','Pr(>F)']
TP5_Pvals["MS3","RS1"]       = anova(Fullmod)['TP5','Pr(>F)']

# MS3RS2
Fullmod <- lm(Width ~ Condition + Duration + Occurance + Coverage + GEV + TP1 + TP2 + TP4 + TP5, data = MS3RS2)
summary(Fullmod)
anova(Fullmod)
Duration_Pvals["MS3","RS2"]  = anova(Fullmod)['Duration','Pr(>F)']
Occurance_Pvals["MS3","RS2"] = anova(Fullmod)['Occurance','Pr(>F)']
Coverage_Pvals["MS3","RS2"]  = anova(Fullmod)['Coverage','Pr(>F)']
GEV_Pvals["MS3","RS2"]       = anova(Fullmod)['GEV','Pr(>F)']
TP1_Pvals["MS3","RS2"]       = anova(Fullmod)['TP1','Pr(>F)']
TP2_Pvals["MS3","RS2"]       = anova(Fullmod)['TP2','Pr(>F)']
TP3_Pvals["MS3","RS2"]       = anova(Fullmod)['TP3','Pr(>F)']
TP4_Pvals["MS3","RS2"]       = anova(Fullmod)['TP4','Pr(>F)']
TP5_Pvals["MS3","RS2"]       = anova(Fullmod)['TP5','Pr(>F)']

# diff
Diffmod <- lm(Width ~ Condition + Duration + Occurance + Coverage + GEV + TP1 + TP2 + TP4 + TP5, data = MS3diff)
summary(Diffmod)
anova(Diffmod)
Duration_Pvals["MS3","Diff"]  = anova(Diffmod)['Duration','Pr(>F)']
Occurance_Pvals["MS3","Diff"] = anova(Diffmod)['Occurance','Pr(>F)']
Coverage_Pvals["MS3","Diff"]  = anova(Diffmod)['Coverage','Pr(>F)']
GEV_Pvals["MS3","Diff"]       = anova(Diffmod)['GEV','Pr(>F)']
TP1_Pvals["MS3","Diff"]       = anova(Diffmod)['TP1','Pr(>F)']
TP2_Pvals["MS3","Diff"]       = anova(Diffmod)['TP2','Pr(>F)']
TP3_Pvals["MS3","Diff"]       = anova(Diffmod)['TP3','Pr(>F)']
TP4_Pvals["MS3","Diff"]       = anova(Diffmod)['TP4','Pr(>F)']
TP5_Pvals["MS3","Diff"]       = anova(Diffmod)['TP5','Pr(>F)']

## MS4
# MS4RS1
Fullmod <- lm(Width ~ Condition + Duration + Occurance + Coverage + GEV + TP1 + TP2 + TP3 + TP5, data = MS4RS1)
summary(Fullmod)
anova(Fullmod)
Duration_Pvals["MS4","RS1"]  = anova(Fullmod)['Duration','Pr(>F)']
Occurance_Pvals["MS4","RS1"] = anova(Fullmod)['Occurance','Pr(>F)']
Coverage_Pvals["MS4","RS1"]  = anova(Fullmod)['Coverage','Pr(>F)']
GEV_Pvals["MS4","RS1"]       = anova(Fullmod)['GEV','Pr(>F)']
TP1_Pvals["MS4","RS1"]       = anova(Fullmod)['TP1','Pr(>F)']
TP2_Pvals["MS4","RS1"]       = anova(Fullmod)['TP2','Pr(>F)']
TP3_Pvals["MS4","RS1"]       = anova(Fullmod)['TP3','Pr(>F)']
TP4_Pvals["MS4","RS1"]       = anova(Fullmod)['TP4','Pr(>F)']
TP5_Pvals["MS4","RS1"]       = anova(Fullmod)['TP5','Pr(>F)']

# MS4RS2
Fullmod <- lm(Width ~ Condition + Duration + Occurance + Coverage + GEV + TP1 + TP2 + TP3 + TP5, data = MS4RS2)
summary(Fullmod)
anova(Fullmod)
Duration_Pvals["MS4","RS2"]  = anova(Fullmod)['Duration','Pr(>F)']
Occurance_Pvals["MS4","RS2"] = anova(Fullmod)['Occurance','Pr(>F)']
Coverage_Pvals["MS4","RS2"]  = anova(Fullmod)['Coverage','Pr(>F)']
GEV_Pvals["MS4","RS2"]       = anova(Fullmod)['GEV','Pr(>F)']
TP1_Pvals["MS4","RS2"]       = anova(Fullmod)['TP1','Pr(>F)']
TP2_Pvals["MS4","RS2"]       = anova(Fullmod)['TP2','Pr(>F)']
TP3_Pvals["MS4","RS2"]       = anova(Fullmod)['TP3','Pr(>F)']
TP4_Pvals["MS4","RS2"]       = anova(Fullmod)['TP4','Pr(>F)']
TP5_Pvals["MS4","RS2"]       = anova(Fullmod)['TP5','Pr(>F)']

# diff
Diffmod <- lm(Width ~ Condition + Duration + Occurance + Coverage + GEV + TP1 + TP2 + TP3 + TP5, data = MS4diff)
summary(Diffmod)
anova(Diffmod)
Duration_Pvals["MS4","Diff"]  = anova(Diffmod)['Duration','Pr(>F)']
Occurance_Pvals["MS4","Diff"] = anova(Diffmod)['Occurance','Pr(>F)']
Coverage_Pvals["MS4","Diff"]  = anova(Diffmod)['Coverage','Pr(>F)']
GEV_Pvals["MS4","Diff"]       = anova(Diffmod)['GEV','Pr(>F)']
TP1_Pvals["MS4","Diff"]       = anova(Diffmod)['TP1','Pr(>F)']
TP2_Pvals["MS4","Diff"]       = anova(Diffmod)['TP2','Pr(>F)']
TP3_Pvals["MS4","Diff"]       = anova(Diffmod)['TP3','Pr(>F)']
TP4_Pvals["MS4","Diff"]       = anova(Diffmod)['TP4','Pr(>F)']
TP5_Pvals["MS4","Diff"]       = anova(Diffmod)['TP5','Pr(>F)']

## MS5
# MS5RS1
Fullmod <- lm(Width ~ Condition + Duration + Occurance + Coverage + GEV + TP1 + TP2 + TP3 + TP4, data = MS5RS1)
summary(Fullmod)
anova(Fullmod)
Duration_Pvals["MS5","RS1"]  = anova(Fullmod)['Duration','Pr(>F)']
Occurance_Pvals["MS5","RS1"] = anova(Fullmod)['Occurance','Pr(>F)']
Coverage_Pvals["MS5","RS1"]  = anova(Fullmod)['Coverage','Pr(>F)']
GEV_Pvals["MS5","RS1"]       = anova(Fullmod)['GEV','Pr(>F)']
TP1_Pvals["MS5","RS1"]       = anova(Fullmod)['TP1','Pr(>F)']
TP2_Pvals["MS5","RS1"]       = anova(Fullmod)['TP2','Pr(>F)']
TP3_Pvals["MS5","RS1"]       = anova(Fullmod)['TP3','Pr(>F)']
TP4_Pvals["MS5","RS1"]       = anova(Fullmod)['TP4','Pr(>F)']
TP5_Pvals["MS5","RS1"]       = anova(Fullmod)['TP5','Pr(>F)']

# MS5RS2
Fullmod <- lm(Width ~ Condition + Duration + Occurance + Coverage + GEV + TP1 + TP2 + TP3 + TP4, data = MS5RS2)
summary(Fullmod)
anova(Fullmod)
Duration_Pvals["MS5","RS2"]  = anova(Fullmod)['Duration','Pr(>F)']
Occurance_Pvals["MS5","RS2"] = anova(Fullmod)['Occurance','Pr(>F)']
Coverage_Pvals["MS5","RS2"]  = anova(Fullmod)['Coverage','Pr(>F)']
GEV_Pvals["MS5","RS2"]       = anova(Fullmod)['GEV','Pr(>F)']
TP1_Pvals["MS5","RS2"]       = anova(Fullmod)['TP1','Pr(>F)']
TP2_Pvals["MS5","RS2"]       = anova(Fullmod)['TP2','Pr(>F)']
TP3_Pvals["MS5","RS2"]       = anova(Fullmod)['TP3','Pr(>F)']
TP4_Pvals["MS5","RS2"]       = anova(Fullmod)['TP4','Pr(>F)']
TP5_Pvals["MS5","RS2"]       = anova(Fullmod)['TP5','Pr(>F)']

# diff
Diffmod <- lm(Width ~ Condition + Duration + Occurance + Coverage + GEV + TP1 + TP2 + TP3 + TP4, data = MS5diff)
summary(Diffmod)
anova(Diffmod)
Duration_Pvals["MS5","Diff"]  = anova(Diffmod)['Duration','Pr(>F)']
Occurance_Pvals["MS5","Diff"] = anova(Diffmod)['Occurance','Pr(>F)']
Coverage_Pvals["MS5","Diff"]  = anova(Diffmod)['Coverage','Pr(>F)']
GEV_Pvals["MS5","Diff"]       = anova(Diffmod)['GEV','Pr(>F)']
TP1_Pvals["MS5","Diff"]       = anova(Diffmod)['TP1','Pr(>F)']
TP2_Pvals["MS5","Diff"]       = anova(Diffmod)['TP2','Pr(>F)']
TP3_Pvals["MS5","Diff"]       = anova(Diffmod)['TP3','Pr(>F)']
TP4_Pvals["MS5","Diff"]       = anova(Diffmod)['TP4','Pr(>F)']
TP5_Pvals["MS5","Diff"]       = anova(Diffmod)['TP5','Pr(>F)']

# Check for significance after FDR correction
Pvalues <- unlist(TP5_Pvals)
non_na_indices <- which(!is.na(Pvalues))
All_adj_p_values <- p.adjust(Pvalues[non_na_indices], method = "fdr")
Pvalues[non_na_indices] <- All_adj_p_values
Pvalues[Pvalues<.05]


## Microstates x PPRatings ----
MicrostatsxPP <- read_excel("Path_to_Microstate_PPRatings_data/MicrostatsxPPratings.xlsx")
View(MicrostatsxPP)
MicrostatsxPP = as.data.frame(MicrostatsxPP)

MicrostatsxPP$IntensityRatings[MicrostatsxPP$Subject == 21] = NA
MicrostatsxPP$IntensityRatings[MicrostatsxPP$Subject == 30] = NA
MicrostatsxPP$IntensityRatings[MicrostatsxPP$Subject == 40] = NA
MicrostatsxPP$IntensityRatings[MicrostatsxPP$Subject == 46] = NA
MicrostatsxPP$IntensityRatings[MicrostatsxPP$Subject == 57] = NA

MS1RS1 = MicrostatsxPP[MicrostatsxPP$MS == 'MS1' & MicrostatsxPP$RS == 'RS1',]
MS2RS1 = MicrostatsxPP[MicrostatsxPP$MS == 'MS2' & MicrostatsxPP$RS == 'RS1',]
MS3RS1 = MicrostatsxPP[MicrostatsxPP$MS == 'MS3' & MicrostatsxPP$RS == 'RS1',]
MS4RS1 = MicrostatsxPP[MicrostatsxPP$MS == 'MS4' & MicrostatsxPP$RS == 'RS1',]
MS5RS1 = MicrostatsxPP[MicrostatsxPP$MS == 'MS5' & MicrostatsxPP$RS == 'RS1',]

MS1RS2 = MicrostatsxPP[MicrostatsxPP$MS == 'MS1' & MicrostatsxPP$RS == 'RS2',]
MS2RS2 = MicrostatsxPP[MicrostatsxPP$MS == 'MS2' & MicrostatsxPP$RS == 'RS2',]
MS3RS2 = MicrostatsxPP[MicrostatsxPP$MS == 'MS3' & MicrostatsxPP$RS == 'RS2',]
MS4RS2 = MicrostatsxPP[MicrostatsxPP$MS == 'MS4' & MicrostatsxPP$RS == 'RS2',]
MS5RS2 = MicrostatsxPP[MicrostatsxPP$MS == 'MS5' & MicrostatsxPP$RS == 'RS2',]

MS1diff = MS1RS1
MS1diff[,c(6:16)] = MS1RS1[,c(6:16)] - MS1RS2[,c(6:16)]
MS2diff = MS2RS1
MS2diff[,c(6:16)] = MS2RS1[,c(6:16)] - MS2RS2[,c(6:16)]
MS3diff = MS3RS1
MS3diff[,c(6:16)] = MS3RS1[,c(6:16)] - MS3RS2[,c(6:16)]
MS4diff = MS4RS1
MS4diff[,c(6:16)] = MS4RS1[,c(6:16)] - MS4RS2[,c(6:16)]
MS5diff = MS5RS1
MS5diff[,c(6:16)] = MS5RS1[,c(6:16)] - MS5RS2[,c(6:16)]


### Intensity ----
# Initialize variables

# Create column and row names
col_names <- c("RS1","RS2","Diff")
row_names <- c("MS1", "MS2", "MS3", "MS4", "MS5")
# Create empty data frames
Duration_Pvals <- data.frame(matrix(ncol = length(col_names), nrow = length(row_names)))
colnames(Duration_Pvals) <- col_names
rownames(Duration_Pvals) <- row_names
Occurance_Pvals = Duration_Pvals
Coverage_Pvals = Duration_Pvals
GEV_Pvals = Duration_Pvals
TP1_Pvals = Duration_Pvals
TP2_Pvals = Duration_Pvals
TP3_Pvals = Duration_Pvals
TP4_Pvals = Duration_Pvals
TP5_Pvals = Duration_Pvals

## MS1
# RS1
FullmodMS1 <- lmer(IntensityRatings ~ IntensityT0 + Condition + GFP:Arm + Duration:Arm + Occurance:Arm + Coverage:Arm + GEV:Arm + TP2:Arm + TP3:Arm + TP4:Arm + TP5:Arm + (1|Subject), data = MS1RS1,REML = FALSE)
summary(FullmodMS1)
anova(FullmodMS1)
Duration_Pvals["MS1","RS1"]  = anova(FullmodMS1)['Arm:Duration','Pr(>F)']
Occurance_Pvals["MS1","RS1"] = anova(FullmodMS1)['Arm:Occurance','Pr(>F)']
Coverage_Pvals["MS1","RS1"]  = anova(FullmodMS1)['Arm:Coverage','Pr(>F)']
GEV_Pvals["MS1","RS1"]       = anova(FullmodMS1)['Arm:GEV','Pr(>F)']
TP1_Pvals["MS1","RS1"]       = anova(FullmodMS1)['Arm:TP1','Pr(>F)']
TP2_Pvals["MS1","RS1"]       = anova(FullmodMS1)['Arm:TP2','Pr(>F)']
TP3_Pvals["MS1","RS1"]       = anova(FullmodMS1)['Arm:TP3','Pr(>F)']
TP4_Pvals["MS1","RS1"]       = anova(FullmodMS1)['Arm:TP4','Pr(>F)']
TP5_Pvals["MS1","RS1"]       = anova(FullmodMS1)['Arm:TP5','Pr(>F)']

# RS2
FullmodMS1 <- lmer(IntensityRatings ~ IntensityT0 + Condition + GFP:Arm + Duration:Arm + Occurance:Arm + Coverage:Arm + GEV:Arm + TP2:Arm + TP3:Arm + TP4:Arm + TP5:Arm + (1|Subject), data = MS1RS2,REML = FALSE)
summary(FullmodMS1)
anova(FullmodMS1)
Duration_Pvals["MS1","RS2"]  = anova(FullmodMS1)['Arm:Duration','Pr(>F)']
Occurance_Pvals["MS1","RS2"] = anova(FullmodMS1)['Arm:Occurance','Pr(>F)']
Coverage_Pvals["MS1","RS2"]  = anova(FullmodMS1)['Arm:Coverage','Pr(>F)']
GEV_Pvals["MS1","RS2"]       = anova(FullmodMS1)['Arm:GEV','Pr(>F)']
TP1_Pvals["MS1","RS2"]       = anova(FullmodMS1)['Arm:TP1','Pr(>F)']
TP2_Pvals["MS1","RS2"]       = anova(FullmodMS1)['Arm:TP2','Pr(>F)']
TP3_Pvals["MS1","RS2"]       = anova(FullmodMS1)['Arm:TP3','Pr(>F)']
TP4_Pvals["MS1","RS2"]       = anova(FullmodMS1)['Arm:TP4','Pr(>F)']
TP5_Pvals["MS1","RS2"]       = anova(FullmodMS1)['Arm:TP5','Pr(>F)']

# Diff
FullmodMS1 <- lmer(IntensityRatings ~ IntensityT0 + Condition + GFP:Arm + Duration:Arm + Occurance:Arm + Coverage:Arm + GEV:Arm + TP2:Arm + TP3:Arm + TP4:Arm + TP5:Arm + (1|Subject), data = MS1diff,REML = FALSE)
summary(FullmodMS1)
anova(FullmodMS1)
Duration_Pvals["MS1","Diff"]  = anova(FullmodMS1)['Arm:Duration','Pr(>F)']
Occurance_Pvals["MS1","Diff"] = anova(FullmodMS1)['Arm:Occurance','Pr(>F)']
Coverage_Pvals["MS1","Diff"]  = anova(FullmodMS1)['Arm:Coverage','Pr(>F)']
GEV_Pvals["MS1","Diff"]       = anova(FullmodMS1)['Arm:GEV','Pr(>F)']
TP1_Pvals["MS1","Diff"]       = anova(FullmodMS1)['Arm:TP1','Pr(>F)']
TP2_Pvals["MS1","Diff"]       = anova(FullmodMS1)['Arm:TP2','Pr(>F)']
TP3_Pvals["MS1","Diff"]       = anova(FullmodMS1)['Arm:TP3','Pr(>F)']
TP4_Pvals["MS1","Diff"]       = anova(FullmodMS1)['Arm:TP4','Pr(>F)']
TP5_Pvals["MS1","Diff"]       = anova(FullmodMS1)['Arm:TP5','Pr(>F)']


## MS2
# RS1
FullmodMS2 <- lmer(IntensityRatings ~ IntensityT0 + Condition + GFP:Arm + Duration:Arm + Occurance:Arm + Coverage:Arm + GEV:Arm + TP1:Arm + TP3:Arm + TP4:Arm + TP5:Arm + (1|Subject), data = MS2RS1,REML = FALSE)
summary(FullmodMS2)
anova(FullmodMS2)
Duration_Pvals["MS2","RS1"]  = anova(FullmodMS2)['Arm:Duration','Pr(>F)']
Occurance_Pvals["MS2","RS1"] = anova(FullmodMS2)['Arm:Occurance','Pr(>F)']
Coverage_Pvals["MS2","RS1"]  = anova(FullmodMS2)['Arm:Coverage','Pr(>F)']
GEV_Pvals["MS2","RS1"]       = anova(FullmodMS2)['Arm:GEV','Pr(>F)']
TP1_Pvals["MS2","RS1"]       = anova(FullmodMS2)['Arm:TP1','Pr(>F)']
TP2_Pvals["MS2","RS1"]       = anova(FullmodMS2)['Arm:TP2','Pr(>F)']
TP3_Pvals["MS2","RS1"]       = anova(FullmodMS2)['Arm:TP3','Pr(>F)']
TP4_Pvals["MS2","RS1"]       = anova(FullmodMS2)['Arm:TP4','Pr(>F)']
TP5_Pvals["MS2","RS1"]       = anova(FullmodMS2)['Arm:TP5','Pr(>F)']

# RS2
FullmodMS2 <- lmer(IntensityRatings ~ IntensityT0 + Condition + GFP:Arm + Duration:Arm + Occurance:Arm + Coverage:Arm + GEV:Arm + TP1:Arm + TP3:Arm + TP4:Arm + TP5:Arm + (1|Subject), data = MS2RS2,REML = FALSE)
summary(FullmodMS2)
anova(FullmodMS2)
Duration_Pvals["MS2","RS2"]  = anova(FullmodMS2)['Arm:Duration','Pr(>F)']
Occurance_Pvals["MS2","RS2"] = anova(FullmodMS2)['Arm:Occurance','Pr(>F)']
Coverage_Pvals["MS2","RS2"]  = anova(FullmodMS2)['Arm:Coverage','Pr(>F)']
GEV_Pvals["MS2","RS2"]       = anova(FullmodMS2)['Arm:GEV','Pr(>F)']
TP1_Pvals["MS2","RS2"]       = anova(FullmodMS2)['Arm:TP1','Pr(>F)']
TP2_Pvals["MS2","RS2"]       = anova(FullmodMS2)['Arm:TP2','Pr(>F)']
TP3_Pvals["MS2","RS2"]       = anova(FullmodMS2)['Arm:TP3','Pr(>F)']
TP4_Pvals["MS2","RS2"]       = anova(FullmodMS2)['Arm:TP4','Pr(>F)']
TP5_Pvals["MS2","RS2"]       = anova(FullmodMS2)['Arm:TP5','Pr(>F)']

# Diff
FullmodMS2 <- lmer(IntensityRatings ~ IntensityT0 + Condition + GFP:Arm + Duration:Arm + Occurance:Arm + Coverage:Arm + GEV:Arm + TP1:Arm + TP3:Arm + TP4:Arm + TP5:Arm + (1|Subject), data = MS2diff,REML = FALSE)
summary(FullmodMS2)
anova(FullmodMS2)
Duration_Pvals["MS2","Diff"]  = anova(FullmodMS2)['Arm:Duration','Pr(>F)']
Occurance_Pvals["MS2","Diff"] = anova(FullmodMS2)['Arm:Occurance','Pr(>F)']
Coverage_Pvals["MS2","Diff"]  = anova(FullmodMS2)['Arm:Coverage','Pr(>F)']
GEV_Pvals["MS2","Diff"]       = anova(FullmodMS2)['Arm:GEV','Pr(>F)']
TP1_Pvals["MS2","Diff"]       = anova(FullmodMS2)['Arm:TP1','Pr(>F)']
TP2_Pvals["MS2","Diff"]       = anova(FullmodMS2)['Arm:TP2','Pr(>F)']
TP3_Pvals["MS2","Diff"]       = anova(FullmodMS2)['Arm:TP3','Pr(>F)']
TP4_Pvals["MS2","Diff"]       = anova(FullmodMS2)['Arm:TP4','Pr(>F)']
TP5_Pvals["MS2","Diff"]       = anova(FullmodMS2)['Arm:TP5','Pr(>F)']

## MS3
# RS1
FullmodMS3 <- lmer(IntensityRatings ~ IntensityT0 + Condition + GFP:Arm + Duration:Arm + Occurance:Arm + Coverage:Arm + GEV:Arm + TP1:Arm + TP2:Arm + TP4:Arm + TP5:Arm + (1|Subject), data = MS3RS1,REML = FALSE)
summary(FullmodMS3)
anova(FullmodMS3)
Duration_Pvals["MS3","RS1"]  = anova(FullmodMS3)['Arm:Duration','Pr(>F)']
Occurance_Pvals["MS3","RS1"] = anova(FullmodMS3)['Arm:Occurance','Pr(>F)']
Coverage_Pvals["MS3","RS1"]  = anova(FullmodMS3)['Arm:Coverage','Pr(>F)']
GEV_Pvals["MS3","RS1"]       = anova(FullmodMS3)['Arm:GEV','Pr(>F)']
TP1_Pvals["MS3","RS1"]       = anova(FullmodMS3)['Arm:TP1','Pr(>F)']
TP2_Pvals["MS3","RS1"]       = anova(FullmodMS3)['Arm:TP2','Pr(>F)']
TP3_Pvals["MS3","RS1"]       = anova(FullmodMS3)['Arm:TP3','Pr(>F)']
TP4_Pvals["MS3","RS1"]       = anova(FullmodMS3)['Arm:TP4','Pr(>F)']
TP5_Pvals["MS3","RS1"]       = anova(FullmodMS3)['Arm:TP5','Pr(>F)']


FullmodMS3 <- lmer(IntensityRatings ~ IntensityT0 + GEV*Arm + (1|Subject), data = MS4RS1,REML = FALSE)
summary(FullmodMS3)

# RS2
FullmodMS3 <- lmer(IntensityRatings ~ IntensityT0 + Condition + GFP:Arm + Duration:Arm + Occurance:Arm + Coverage:Arm + GEV:Arm + TP1:Arm + TP2:Arm + TP4:Arm + TP5:Arm + (1|Subject), data = MS3RS2,REML = FALSE)
summary(FullmodMS3)
anova(FullmodMS3)
Duration_Pvals["MS3","RS2"]  = anova(FullmodMS3)['Arm:Duration','Pr(>F)']
Occurance_Pvals["MS3","RS2"] = anova(FullmodMS3)['Arm:Occurance','Pr(>F)']
Coverage_Pvals["MS3","RS2"]  = anova(FullmodMS3)['Arm:Coverage','Pr(>F)']
GEV_Pvals["MS3","RS2"]       = anova(FullmodMS3)['Arm:GEV','Pr(>F)']
TP1_Pvals["MS3","RS2"]       = anova(FullmodMS3)['Arm:TP1','Pr(>F)']
TP2_Pvals["MS3","RS2"]       = anova(FullmodMS3)['Arm:TP2','Pr(>F)']
TP3_Pvals["MS3","RS2"]       = anova(FullmodMS3)['Arm:TP3','Pr(>F)']
TP4_Pvals["MS3","RS2"]       = anova(FullmodMS3)['Arm:TP4','Pr(>F)']
TP5_Pvals["MS3","RS2"]       = anova(FullmodMS3)['Arm:TP5','Pr(>F)']

# Diff
FullmodMS3 <- lmer(IntensityRatings ~ IntensityT0 + Condition + GFP:Arm + Duration:Arm + Occurance:Arm + Coverage:Arm + GEV:Arm + TP1:Arm + TP2:Arm + TP4:Arm + TP5:Arm + (1|Subject), data = MS3diff,REML = FALSE)
summary(FullmodMS3)
anova(FullmodMS3)
Duration_Pvals["MS3","Diff"]  = anova(FullmodMS3)['Arm:Duration','Pr(>F)']
Occurance_Pvals["MS3","Diff"] = anova(FullmodMS3)['Arm:Occurance','Pr(>F)']
Coverage_Pvals["MS3","Diff"]  = anova(FullmodMS3)['Arm:Coverage','Pr(>F)']
GEV_Pvals["MS3","Diff"]       = anova(FullmodMS3)['Arm:GEV','Pr(>F)']
TP1_Pvals["MS3","Diff"]       = anova(FullmodMS3)['Arm:TP1','Pr(>F)']
TP2_Pvals["MS3","Diff"]       = anova(FullmodMS3)['Arm:TP2','Pr(>F)']
TP3_Pvals["MS3","Diff"]       = anova(FullmodMS3)['Arm:TP3','Pr(>F)']
TP4_Pvals["MS3","Diff"]       = anova(FullmodMS3)['Arm:TP4','Pr(>F)']
TP5_Pvals["MS3","Diff"]       = anova(FullmodMS3)['Arm:TP5','Pr(>F)']


## MS4
# RS1
FullmodMS4 <- lmer(IntensityRatings ~ IntensityT0 + Condition + GFP:Arm + Duration:Arm + Occurance:Arm + Coverage:Arm + GEV:Arm + TP2:Arm + TP3:Arm + TP1:Arm + TP5:Arm + (1|Subject), data = MS4RS1,REML = FALSE)
summary(FullmodMS4)
anova(FullmodMS4)
Duration_Pvals["MS4","RS1"]  = anova(FullmodMS4)['Arm:Duration','Pr(>F)']
Occurance_Pvals["MS4","RS1"] = anova(FullmodMS4)['Arm:Occurance','Pr(>F)']
Coverage_Pvals["MS4","RS1"]  = anova(FullmodMS4)['Arm:Coverage','Pr(>F)']
GEV_Pvals["MS4","RS1"]       = anova(FullmodMS4)['Arm:GEV','Pr(>F)']
TP1_Pvals["MS4","RS1"]       = anova(FullmodMS4)['Arm:TP1','Pr(>F)']
TP2_Pvals["MS4","RS1"]       = anova(FullmodMS4)['Arm:TP2','Pr(>F)']
TP3_Pvals["MS4","RS1"]       = anova(FullmodMS4)['Arm:TP3','Pr(>F)']
TP4_Pvals["MS4","RS1"]       = anova(FullmodMS4)['Arm:TP4','Pr(>F)']
TP5_Pvals["MS4","RS1"]       = anova(FullmodMS4)['Arm:TP5','Pr(>F)']

# RS2
FullmodMS4 <- lmer(IntensityRatings ~ IntensityT0 + Condition + GFP:Arm + Duration:Arm + Occurance:Arm + Coverage:Arm + GEV:Arm + TP2:Arm + TP3:Arm + TP1:Arm + TP5:Arm + (1|Subject), data = MS4RS2,REML = FALSE)
summary(FullmodMS4)
anova(FullmodMS4)
Duration_Pvals["MS4","RS2"]  = anova(FullmodMS4)['Arm:Duration','Pr(>F)']
Occurance_Pvals["MS4","RS2"] = anova(FullmodMS4)['Arm:Occurance','Pr(>F)']
Coverage_Pvals["MS4","RS2"]  = anova(FullmodMS4)['Arm:Coverage','Pr(>F)']
GEV_Pvals["MS4","RS2"]       = anova(FullmodMS4)['Arm:GEV','Pr(>F)']
TP1_Pvals["MS4","RS2"]       = anova(FullmodMS4)['Arm:TP1','Pr(>F)']
TP2_Pvals["MS4","RS2"]       = anova(FullmodMS4)['Arm:TP2','Pr(>F)']
TP3_Pvals["MS4","RS2"]       = anova(FullmodMS4)['Arm:TP3','Pr(>F)']
TP4_Pvals["MS4","RS2"]       = anova(FullmodMS4)['Arm:TP4','Pr(>F)']
TP5_Pvals["MS4","RS2"]       = anova(FullmodMS4)['Arm:TP5','Pr(>F)']

# Diff
FullmodMS4 <- lmer(IntensityRatings ~ IntensityT0 + Condition + GFP:Arm + Duration:Arm + Occurance:Arm + Coverage:Arm + GEV:Arm + TP2:Arm + TP3:Arm + TP1:Arm + TP5:Arm + (1|Subject), data = MS4diff,REML = FALSE)
summary(FullmodMS4)
anova(FullmodMS4)
Duration_Pvals["MS4","Diff"]  = anova(FullmodMS4)['Arm:Duration','Pr(>F)']
Occurance_Pvals["MS4","Diff"] = anova(FullmodMS4)['Arm:Occurance','Pr(>F)']
Coverage_Pvals["MS4","Diff"]  = anova(FullmodMS4)['Arm:Coverage','Pr(>F)']
GEV_Pvals["MS4","Diff"]       = anova(FullmodMS4)['Arm:GEV','Pr(>F)']
TP1_Pvals["MS4","Diff"]       = anova(FullmodMS4)['Arm:TP1','Pr(>F)']
TP2_Pvals["MS4","Diff"]       = anova(FullmodMS4)['Arm:TP2','Pr(>F)']
TP3_Pvals["MS4","Diff"]       = anova(FullmodMS4)['Arm:TP3','Pr(>F)']
TP4_Pvals["MS4","Diff"]       = anova(FullmodMS4)['Arm:TP4','Pr(>F)']
TP5_Pvals["MS4","Diff"]       = anova(FullmodMS4)['Arm:TP5','Pr(>F)']

## MS5
# RS1
FullmodMS5 <- lmer(IntensityRatings ~ IntensityT0 + Condition + GFP:Arm + Duration:Arm + Occurance:Arm + Coverage:Arm + GEV:Arm + TP2:Arm + TP3:Arm + TP4:Arm + TP1:Arm + (1|Subject), data = MS5RS1,REML = FALSE)
summary(FullmodMS5)
anova(FullmodMS5)
Duration_Pvals["MS5","RS1"]  = anova(FullmodMS5)['Arm:Duration','Pr(>F)']
Occurance_Pvals["MS5","RS1"] = anova(FullmodMS5)['Arm:Occurance','Pr(>F)']
Coverage_Pvals["MS5","RS1"]  = anova(FullmodMS5)['Arm:Coverage','Pr(>F)']
GEV_Pvals["MS5","RS1"]       = anova(FullmodMS5)['Arm:GEV','Pr(>F)']
TP1_Pvals["MS5","RS1"]       = anova(FullmodMS5)['Arm:TP1','Pr(>F)']
TP2_Pvals["MS5","RS1"]       = anova(FullmodMS5)['Arm:TP2','Pr(>F)']
TP3_Pvals["MS5","RS1"]       = anova(FullmodMS5)['Arm:TP3','Pr(>F)']
TP4_Pvals["MS5","RS1"]       = anova(FullmodMS5)['Arm:TP4','Pr(>F)']
TP5_Pvals["MS5","RS1"]       = anova(FullmodMS5)['Arm:TP5','Pr(>F)']

# RS2
FullmodMS5 <- lmer(IntensityRatings ~ IntensityT0 + Condition + GFP:Arm + Duration:Arm + Occurance:Arm + Coverage:Arm + GEV:Arm + TP2:Arm + TP3:Arm + TP4:Arm + TP1:Arm + (1|Subject), data = MS5RS2,REML = FALSE)
summary(FullmodMS5)
anova(FullmodMS5)
Duration_Pvals["MS5","RS2"]  = anova(FullmodMS5)['Arm:Duration','Pr(>F)']
Occurance_Pvals["MS5","RS2"] = anova(FullmodMS5)['Arm:Occurance','Pr(>F)']
Coverage_Pvals["MS5","RS2"]  = anova(FullmodMS5)['Arm:Coverage','Pr(>F)']
GEV_Pvals["MS5","RS2"]       = anova(FullmodMS5)['Arm:GEV','Pr(>F)']
TP1_Pvals["MS5","RS2"]       = anova(FullmodMS5)['Arm:TP1','Pr(>F)']
TP2_Pvals["MS5","RS2"]       = anova(FullmodMS5)['Arm:TP2','Pr(>F)']
TP3_Pvals["MS5","RS2"]       = anova(FullmodMS5)['Arm:TP3','Pr(>F)']
TP4_Pvals["MS5","RS2"]       = anova(FullmodMS5)['Arm:TP4','Pr(>F)']
TP5_Pvals["MS5","RS2"]       = anova(FullmodMS5)['Arm:TP5','Pr(>F)']

# Diff
FullmodMS5 <- lmer(IntensityRatings ~ IntensityT0 + Condition + GFP:Arm + Duration:Arm + Occurance:Arm + Coverage:Arm + GEV:Arm + TP2:Arm + TP3:Arm + TP4:Arm + TP1:Arm + (1|Subject), data = MS5diff,REML = FALSE)
summary(FullmodMS5)
anova(FullmodMS5)
Duration_Pvals["MS5","Diff"]  = anova(FullmodMS5)['Arm:Duration','Pr(>F)']
Occurance_Pvals["MS5","Diff"] = anova(FullmodMS5)['Arm:Occurance','Pr(>F)']
Coverage_Pvals["MS5","Diff"]  = anova(FullmodMS5)['Arm:Coverage','Pr(>F)']
GEV_Pvals["MS5","Diff"]       = anova(FullmodMS5)['Arm:GEV','Pr(>F)']
TP1_Pvals["MS5","Diff"]       = anova(FullmodMS5)['Arm:TP1','Pr(>F)']
TP2_Pvals["MS5","Diff"]       = anova(FullmodMS5)['Arm:TP2','Pr(>F)']
TP3_Pvals["MS5","Diff"]       = anova(FullmodMS5)['Arm:TP3','Pr(>F)']
TP4_Pvals["MS5","Diff"]       = anova(FullmodMS5)['Arm:TP4','Pr(>F)']
TP5_Pvals["MS5","Diff"]       = anova(FullmodMS5)['Arm:TP5','Pr(>F)']


# Check for significance after FDR correction
Pvalues <- unlist(Coverage_Pvals)
non_na_indices <- which(!is.na(Pvalues))
All_adj_p_values <- p.adjust(Pvalues[non_na_indices], method = "fdr")
Pvalues[non_na_indices] <- All_adj_p_values
Pvalues[Pvalues<.05]


### Unpleasantness ----
# Initialize variables

# Create column and row names
col_names <- c("RS1","RS2","Diff")
row_names <- c("MS1", "MS2", "MS3", "MS4", "MS5")
# Create empty data frames
Duration_Pvals <- data.frame(matrix(ncol = length(col_names), nrow = length(row_names)))
colnames(Duration_Pvals) <- col_names
rownames(Duration_Pvals) <- row_names
Occurance_Pvals = Duration_Pvals
Coverage_Pvals = Duration_Pvals
GEV_Pvals = Duration_Pvals
TP1_Pvals = Duration_Pvals
TP2_Pvals = Duration_Pvals
TP3_Pvals = Duration_Pvals
TP4_Pvals = Duration_Pvals
TP5_Pvals = Duration_Pvals

## MS1
# RS1
FullmodMS1 <- lmer(UnpleasantnessRatings ~ UnpleasantnessT0 + Condition + GFP:Arm + Duration:Arm + Occurance:Arm + Coverage:Arm + GEV:Arm + TP2:Arm + TP3:Arm + TP4:Arm + TP5:Arm + (1|Subject), data = MS1RS1,REML = FALSE,contrasts = list(Condition = "contr.sum", Arm = "contr.sum"))
summary(FullmodMS1)
anova(FullmodMS1)
Duration_Pvals["MS1","RS1"]  = anova(FullmodMS1)['Arm:Duration','Pr(>F)']
Occurance_Pvals["MS1","RS1"] = anova(FullmodMS1)['Arm:Occurance','Pr(>F)']
Coverage_Pvals["MS1","RS1"]  = anova(FullmodMS1)['Arm:Coverage','Pr(>F)']
GEV_Pvals["MS1","RS1"]       = anova(FullmodMS1)['Arm:GEV','Pr(>F)']
TP1_Pvals["MS1","RS1"]       = anova(FullmodMS1)['Arm:TP1','Pr(>F)']
TP2_Pvals["MS1","RS1"]       = anova(FullmodMS1)['Arm:TP2','Pr(>F)']
TP3_Pvals["MS1","RS1"]       = anova(FullmodMS1)['Arm:TP3','Pr(>F)']
TP4_Pvals["MS1","RS1"]       = anova(FullmodMS1)['Arm:TP4','Pr(>F)']
TP5_Pvals["MS1","RS1"]       = anova(FullmodMS1)['Arm:TP5','Pr(>F)']

# RS2
FullmodMS1 <- lmer(UnpleasantnessRatings ~ UnpleasantnessT0 + Condition + GFP:Arm + Duration:Arm + Occurance:Arm + Coverage:Arm + GEV:Arm + TP2:Arm + TP3:Arm + TP4:Arm + TP5:Arm + (1|Subject), data = MS1RS2,REML = FALSE)
summary(FullmodMS1)
anova(FullmodMS1)
Duration_Pvals["MS1","RS2"]  = anova(FullmodMS1)['Arm:Duration','Pr(>F)']
Occurance_Pvals["MS1","RS2"] = anova(FullmodMS1)['Arm:Occurance','Pr(>F)']
Coverage_Pvals["MS1","RS2"]  = anova(FullmodMS1)['Arm:Coverage','Pr(>F)']
GEV_Pvals["MS1","RS2"]       = anova(FullmodMS1)['Arm:GEV','Pr(>F)']
TP1_Pvals["MS1","RS2"]       = anova(FullmodMS1)['Arm:TP1','Pr(>F)']
TP2_Pvals["MS1","RS2"]       = anova(FullmodMS1)['Arm:TP2','Pr(>F)']
TP3_Pvals["MS1","RS2"]       = anova(FullmodMS1)['Arm:TP3','Pr(>F)']
TP4_Pvals["MS1","RS2"]       = anova(FullmodMS1)['Arm:TP4','Pr(>F)']
TP5_Pvals["MS1","RS2"]       = anova(FullmodMS1)['Arm:TP5','Pr(>F)']

# Diff
FullmodMS1 <- lmer(UnpleasantnessRatings ~ UnpleasantnessT0 + Condition + GFP:Arm + Duration:Arm + Occurance:Arm + Coverage:Arm + GEV:Arm + TP2:Arm + TP3:Arm + TP4:Arm + TP5:Arm + (1|Subject), data = MS1diff,REML = FALSE)
summary(FullmodMS1)
anova(FullmodMS1)
Duration_Pvals["MS1","Diff"]  = anova(FullmodMS1)['Arm:Duration','Pr(>F)']
Occurance_Pvals["MS1","Diff"] = anova(FullmodMS1)['Arm:Occurance','Pr(>F)']
Coverage_Pvals["MS1","Diff"]  = anova(FullmodMS1)['Arm:Coverage','Pr(>F)']
GEV_Pvals["MS1","Diff"]       = anova(FullmodMS1)['Arm:GEV','Pr(>F)']
TP1_Pvals["MS1","Diff"]       = anova(FullmodMS1)['Arm:TP1','Pr(>F)']
TP2_Pvals["MS1","Diff"]       = anova(FullmodMS1)['Arm:TP2','Pr(>F)']
TP3_Pvals["MS1","Diff"]       = anova(FullmodMS1)['Arm:TP3','Pr(>F)']
TP4_Pvals["MS1","Diff"]       = anova(FullmodMS1)['Arm:TP4','Pr(>F)']
TP5_Pvals["MS1","Diff"]       = anova(FullmodMS1)['Arm:TP5','Pr(>F)']

## MS2
# RS1
FullmodMS2 <- lmer(UnpleasantnessRatings ~ UnpleasantnessT0 + Condition + GFP:Arm + Duration:Arm + Occurance:Arm + Coverage:Arm + GEV:Arm + TP1:Arm + TP3:Arm + TP4:Arm + TP5:Arm + (1|Subject), data = MS2RS1,REML = FALSE)
summary(FullmodMS2)
anova(FullmodMS2)
Duration_Pvals["MS2","RS1"]  = anova(FullmodMS2)['Arm:Duration','Pr(>F)']
Occurance_Pvals["MS2","RS1"] = anova(FullmodMS2)['Arm:Occurance','Pr(>F)']
Coverage_Pvals["MS2","RS1"]  = anova(FullmodMS2)['Arm:Coverage','Pr(>F)']
GEV_Pvals["MS2","RS1"]       = anova(FullmodMS2)['Arm:GEV','Pr(>F)']
TP1_Pvals["MS2","RS1"]       = anova(FullmodMS2)['Arm:TP1','Pr(>F)']
TP2_Pvals["MS2","RS1"]       = anova(FullmodMS2)['Arm:TP2','Pr(>F)']
TP3_Pvals["MS2","RS1"]       = anova(FullmodMS2)['Arm:TP3','Pr(>F)']
TP4_Pvals["MS2","RS1"]       = anova(FullmodMS2)['Arm:TP4','Pr(>F)']
TP5_Pvals["MS2","RS1"]       = anova(FullmodMS2)['Arm:TP5','Pr(>F)']


# RS2
FullmodMS2 <- lmer(UnpleasantnessRatings ~ UnpleasantnessT0 + Condition + GFP:Arm + Duration:Arm + Occurance:Arm + Coverage:Arm + GEV:Arm + TP1:Arm + TP3:Arm + TP4:Arm + TP5:Arm + (1|Subject), data = MS2RS2,REML = FALSE)
summary(FullmodMS2)
anova(FullmodMS2)
Duration_Pvals["MS2","RS2"]  = anova(FullmodMS2)['Arm:Duration','Pr(>F)']
Occurance_Pvals["MS2","RS2"] = anova(FullmodMS2)['Arm:Occurance','Pr(>F)']
Coverage_Pvals["MS2","RS2"]  = anova(FullmodMS2)['Arm:Coverage','Pr(>F)']
GEV_Pvals["MS2","RS2"]       = anova(FullmodMS2)['Arm:GEV','Pr(>F)']
TP1_Pvals["MS2","RS2"]       = anova(FullmodMS2)['Arm:TP1','Pr(>F)']
TP2_Pvals["MS2","RS2"]       = anova(FullmodMS2)['Arm:TP2','Pr(>F)']
TP3_Pvals["MS2","RS2"]       = anova(FullmodMS2)['Arm:TP3','Pr(>F)']
TP4_Pvals["MS2","RS2"]       = anova(FullmodMS2)['Arm:TP4','Pr(>F)']
TP5_Pvals["MS2","RS2"]       = anova(FullmodMS2)['Arm:TP5','Pr(>F)']

# Diff
FullmodMS2 <- lmer(UnpleasantnessRatings ~ UnpleasantnessT0 + Condition + GFP:Arm + Duration:Arm + Occurance:Arm + Coverage:Arm + GEV:Arm + TP1:Arm + TP3:Arm + TP4:Arm + TP5:Arm + (1|Subject), data = MS2diff,REML = FALSE)
summary(FullmodMS2)
anova(FullmodMS2)
Duration_Pvals["MS2","Diff"]  = anova(FullmodMS2)['Arm:Duration','Pr(>F)']
Occurance_Pvals["MS2","Diff"] = anova(FullmodMS2)['Arm:Occurance','Pr(>F)']
Coverage_Pvals["MS2","Diff"]  = anova(FullmodMS2)['Arm:Coverage','Pr(>F)']
GEV_Pvals["MS2","Diff"]       = anova(FullmodMS2)['Arm:GEV','Pr(>F)']
TP1_Pvals["MS2","Diff"]       = anova(FullmodMS2)['Arm:TP1','Pr(>F)']
TP2_Pvals["MS2","Diff"]       = anova(FullmodMS2)['Arm:TP2','Pr(>F)']
TP3_Pvals["MS2","Diff"]       = anova(FullmodMS2)['Arm:TP3','Pr(>F)']
TP4_Pvals["MS2","Diff"]       = anova(FullmodMS2)['Arm:TP4','Pr(>F)']
TP5_Pvals["MS2","Diff"]       = anova(FullmodMS2)['Arm:TP5','Pr(>F)']


## MS3
# RS1
FullmodMS3 <- lmer(UnpleasantnessRatings ~ UnpleasantnessT0 + Condition + GFP:Arm + Duration:Arm + Occurance:Arm + Coverage:Arm + GEV:Arm + TP2:Arm + TP1:Arm + TP4:Arm + TP5:Arm + (1|Subject), data = MS3RS1,REML = FALSE)
summary(FullmodMS3)
anova(FullmodMS3)
Duration_Pvals["MS3","RS1"]  = anova(FullmodMS3)['Arm:Duration','Pr(>F)']
Occurance_Pvals["MS3","RS1"] = anova(FullmodMS3)['Arm:Occurance','Pr(>F)']
Coverage_Pvals["MS3","RS1"]  = anova(FullmodMS3)['Arm:Coverage','Pr(>F)']
GEV_Pvals["MS3","RS1"]       = anova(FullmodMS3)['Arm:GEV','Pr(>F)']
TP1_Pvals["MS3","RS1"]       = anova(FullmodMS3)['Arm:TP1','Pr(>F)']
TP2_Pvals["MS3","RS1"]       = anova(FullmodMS3)['Arm:TP2','Pr(>F)']
TP3_Pvals["MS3","RS1"]       = anova(FullmodMS3)['Arm:TP3','Pr(>F)']
TP4_Pvals["MS3","RS1"]       = anova(FullmodMS3)['Arm:TP4','Pr(>F)']
TP5_Pvals["MS3","RS1"]       = anova(FullmodMS3)['Arm:TP5','Pr(>F)']

# RS2
FullmodMS3 <- lmer(UnpleasantnessRatings ~ UnpleasantnessT0 + Condition + GFP:Arm + Duration:Arm + Occurance:Arm + Coverage:Arm + GEV:Arm + TP2:Arm + TP1:Arm + TP4:Arm + TP5:Arm + (1|Subject), data = MS3RS2,REML = FALSE)
summary(FullmodMS3)
anova(FullmodMS3)
Duration_Pvals["MS3","RS2"]  = anova(FullmodMS3)['Arm:Duration','Pr(>F)']
Occurance_Pvals["MS3","RS2"] = anova(FullmodMS3)['Arm:Occurance','Pr(>F)']
Coverage_Pvals["MS3","RS2"]  = anova(FullmodMS3)['Arm:Coverage','Pr(>F)']
GEV_Pvals["MS3","RS2"]       = anova(FullmodMS3)['Arm:GEV','Pr(>F)']
TP1_Pvals["MS3","RS2"]       = anova(FullmodMS3)['Arm:TP1','Pr(>F)']
TP2_Pvals["MS3","RS2"]       = anova(FullmodMS3)['Arm:TP2','Pr(>F)']
TP3_Pvals["MS3","RS2"]       = anova(FullmodMS3)['Arm:TP3','Pr(>F)']
TP4_Pvals["MS3","RS2"]       = anova(FullmodMS3)['Arm:TP4','Pr(>F)']
TP5_Pvals["MS3","RS2"]       = anova(FullmodMS3)['Arm:TP5','Pr(>F)']

# Diff
FullmodMS3 <- lmer(UnpleasantnessRatings ~ UnpleasantnessT0 + Condition + GFP:Arm + Duration:Arm + Occurance:Arm + Coverage:Arm + GEV:Arm + TP2:Arm + TP1:Arm + TP4:Arm + TP5:Arm + (1|Subject), data = MS3diff,REML = FALSE)
summary(FullmodMS3)
anova(FullmodMS3)
Duration_Pvals["MS3","Diff"]  = anova(FullmodMS3)['Arm:Duration','Pr(>F)']
Occurance_Pvals["MS3","Diff"] = anova(FullmodMS3)['Arm:Occurance','Pr(>F)']
Coverage_Pvals["MS3","Diff"]  = anova(FullmodMS3)['Arm:Coverage','Pr(>F)']
GEV_Pvals["MS3","Diff"]       = anova(FullmodMS3)['Arm:GEV','Pr(>F)']
TP1_Pvals["MS3","Diff"]       = anova(FullmodMS3)['Arm:TP1','Pr(>F)']
TP2_Pvals["MS3","Diff"]       = anova(FullmodMS3)['Arm:TP2','Pr(>F)']
TP3_Pvals["MS3","Diff"]       = anova(FullmodMS3)['Arm:TP3','Pr(>F)']
TP4_Pvals["MS3","Diff"]       = anova(FullmodMS3)['Arm:TP4','Pr(>F)']
TP5_Pvals["MS3","Diff"]       = anova(FullmodMS3)['Arm:TP5','Pr(>F)']


## MS4
# RS1
FullmodMS4 <- lmer(UnpleasantnessRatings ~ UnpleasantnessT0 + Condition + GFP:Arm + Duration:Arm + Occurance:Arm + Coverage:Arm + GEV:Arm + TP2:Arm + TP3:Arm + TP1:Arm + TP5:Arm + (1|Subject), data = MS4RS1,REML = FALSE)
summary(FullmodMS4)
anova(FullmodMS4)
Duration_Pvals["MS4","RS1"]  = anova(FullmodMS4)['Arm:Duration','Pr(>F)']
Occurance_Pvals["MS4","RS1"] = anova(FullmodMS4)['Arm:Occurance','Pr(>F)']
Coverage_Pvals["MS4","RS1"]  = anova(FullmodMS4)['Arm:Coverage','Pr(>F)']
GEV_Pvals["MS4","RS1"]       = anova(FullmodMS4)['Arm:GEV','Pr(>F)']
TP1_Pvals["MS4","RS1"]       = anova(FullmodMS4)['Arm:TP1','Pr(>F)']
TP2_Pvals["MS4","RS1"]       = anova(FullmodMS4)['Arm:TP2','Pr(>F)']
TP3_Pvals["MS4","RS1"]       = anova(FullmodMS4)['Arm:TP3','Pr(>F)']
TP4_Pvals["MS4","RS1"]       = anova(FullmodMS4)['Arm:TP4','Pr(>F)']
TP5_Pvals["MS4","RS1"]       = anova(FullmodMS4)['Arm:TP5','Pr(>F)']

# RS2
FullmodMS4 <- lmer(UnpleasantnessRatings ~ UnpleasantnessT0 + Condition + GFP:Arm + Duration:Arm + Occurance:Arm + Coverage:Arm + GEV:Arm + TP2:Arm + TP3:Arm + TP1:Arm + TP5:Arm + (1|Subject), data = MS4RS2,REML = FALSE)
summary(FullmodMS4)
anova(FullmodMS4)
Duration_Pvals["MS4","RS2"]  = anova(FullmodMS4)['Arm:Duration','Pr(>F)']
Occurance_Pvals["MS4","RS2"] = anova(FullmodMS4)['Arm:Occurance','Pr(>F)']
Coverage_Pvals["MS4","RS2"]  = anova(FullmodMS4)['Arm:Coverage','Pr(>F)']
GEV_Pvals["MS4","RS2"]       = anova(FullmodMS4)['Arm:GEV','Pr(>F)']
TP1_Pvals["MS4","RS2"]       = anova(FullmodMS4)['Arm:TP1','Pr(>F)']
TP2_Pvals["MS4","RS2"]       = anova(FullmodMS4)['Arm:TP2','Pr(>F)']
TP3_Pvals["MS4","RS2"]       = anova(FullmodMS4)['Arm:TP3','Pr(>F)']
TP4_Pvals["MS4","RS2"]       = anova(FullmodMS4)['Arm:TP4','Pr(>F)']
TP5_Pvals["MS4","RS2"]       = anova(FullmodMS4)['Arm:TP5','Pr(>F)']

# Diff
FullmodMS4 <- lmer(UnpleasantnessRatings ~ UnpleasantnessT0 + Condition + GFP:Arm + Duration:Arm + Occurance:Arm + Coverage:Arm + GEV:Arm + TP2:Arm + TP3:Arm + TP1:Arm + TP5:Arm + (1|Subject), data = MS4diff,REML = FALSE)
summary(FullmodMS4)
anova(FullmodMS4)
Duration_Pvals["MS4","Diff"]  = anova(FullmodMS4)['Arm:Duration','Pr(>F)']
Occurance_Pvals["MS4","Diff"] = anova(FullmodMS4)['Arm:Occurance','Pr(>F)']
Coverage_Pvals["MS4","Diff"]  = anova(FullmodMS4)['Arm:Coverage','Pr(>F)']
GEV_Pvals["MS4","Diff"]       = anova(FullmodMS4)['Arm:GEV','Pr(>F)']
TP1_Pvals["MS4","Diff"]       = anova(FullmodMS4)['Arm:TP1','Pr(>F)']
TP2_Pvals["MS4","Diff"]       = anova(FullmodMS4)['Arm:TP2','Pr(>F)']
TP3_Pvals["MS4","Diff"]       = anova(FullmodMS4)['Arm:TP3','Pr(>F)']
TP4_Pvals["MS4","Diff"]       = anova(FullmodMS4)['Arm:TP4','Pr(>F)']
TP5_Pvals["MS4","Diff"]       = anova(FullmodMS4)['Arm:TP5','Pr(>F)']


## MS5
# RS1
FullmodMS5 <- lmer(UnpleasantnessRatings ~ UnpleasantnessT0 + Condition + GFP:Arm + Duration:Arm + Occurance:Arm + Coverage:Arm + GEV:Arm + TP2:Arm + TP3:Arm + TP4:Arm + TP1:Arm + (1|Subject), data = MS5RS1,REML = FALSE)
summary(FullmodMS5)
anova(FullmodMS5)
Duration_Pvals["MS5","RS1"]  = anova(FullmodMS5)['Arm:Duration','Pr(>F)']
Occurance_Pvals["MS5","RS1"] = anova(FullmodMS5)['Arm:Occurance','Pr(>F)']
Coverage_Pvals["MS5","RS1"]  = anova(FullmodMS5)['Arm:Coverage','Pr(>F)']
GEV_Pvals["MS5","RS1"]       = anova(FullmodMS5)['Arm:GEV','Pr(>F)']
TP1_Pvals["MS5","RS1"]       = anova(FullmodMS5)['Arm:TP1','Pr(>F)']
TP2_Pvals["MS5","RS1"]       = anova(FullmodMS5)['Arm:TP2','Pr(>F)']
TP3_Pvals["MS5","RS1"]       = anova(FullmodMS5)['Arm:TP3','Pr(>F)']
TP4_Pvals["MS5","RS1"]       = anova(FullmodMS5)['Arm:TP4','Pr(>F)']
TP5_Pvals["MS5","RS1"]       = anova(FullmodMS5)['Arm:TP5','Pr(>F)']


# RS2
FullmodMS5 <- lmer(UnpleasantnessRatings ~ UnpleasantnessT0 + Condition + GFP:Arm + Duration:Arm + Occurance:Arm + Coverage:Arm + GEV:Arm + TP2:Arm + TP3:Arm + TP4:Arm + TP1:Arm + (1|Subject), data = MS5RS2,REML = FALSE)
summary(FullmodMS5)
anova(FullmodMS5)
Duration_Pvals["MS5","RS2"]  = anova(FullmodMS5)['Arm:Duration','Pr(>F)']
Occurance_Pvals["MS5","RS2"] = anova(FullmodMS5)['Arm:Occurance','Pr(>F)']
Coverage_Pvals["MS5","RS2"]  = anova(FullmodMS5)['Arm:Coverage','Pr(>F)']
GEV_Pvals["MS5","RS2"]       = anova(FullmodMS5)['Arm:GEV','Pr(>F)']
TP1_Pvals["MS5","RS2"]       = anova(FullmodMS5)['Arm:TP1','Pr(>F)']
TP2_Pvals["MS5","RS2"]       = anova(FullmodMS5)['Arm:TP2','Pr(>F)']
TP3_Pvals["MS5","RS2"]       = anova(FullmodMS5)['Arm:TP3','Pr(>F)']
TP4_Pvals["MS5","RS2"]       = anova(FullmodMS5)['Arm:TP4','Pr(>F)']
TP5_Pvals["MS5","RS2"]       = anova(FullmodMS5)['Arm:TP5','Pr(>F)']

# Diff
FullmodMS5 <- lmer(UnpleasantnessRatings ~ UnpleasantnessT0 + Condition + GFP:Arm + Duration:Arm + Occurance:Arm + Coverage:Arm + GEV:Arm + TP2:Arm + TP3:Arm + TP4:Arm + TP1:Arm + (1|Subject), data = MS5diff,REML = FALSE)
summary(FullmodMS5)
anova(FullmodMS5)
Duration_Pvals["MS5","Diff"]  = anova(FullmodMS5)['Arm:Duration','Pr(>F)']
Occurance_Pvals["MS5","Diff"] = anova(FullmodMS5)['Arm:Occurance','Pr(>F)']
Coverage_Pvals["MS5","Diff"]  = anova(FullmodMS5)['Arm:Coverage','Pr(>F)']
GEV_Pvals["MS5","Diff"]       = anova(FullmodMS5)['Arm:GEV','Pr(>F)']
TP1_Pvals["MS5","Diff"]       = anova(FullmodMS5)['Arm:TP1','Pr(>F)']
TP2_Pvals["MS5","Diff"]       = anova(FullmodMS5)['Arm:TP2','Pr(>F)']
TP3_Pvals["MS5","Diff"]       = anova(FullmodMS5)['Arm:TP3','Pr(>F)']
TP4_Pvals["MS5","Diff"]       = anova(FullmodMS5)['Arm:TP4','Pr(>F)']
TP5_Pvals["MS5","Diff"]       = anova(FullmodMS5)['Arm:TP5','Pr(>F)']

# Check for significance after FDR correction
Pvalues <- unlist(Duration_Pvals)
non_na_indices <- which(!is.na(Pvalues))
All_adj_p_values <- p.adjust(Pvalues[non_na_indices], method = "fdr")
Pvalues[non_na_indices] <- All_adj_p_values
Pvalues[Pvalues<.05]

## Double Subtracted ----
## Microstates x PPRatings ----
MicrostatsxPP <- read_excel("Path_to_Microstate_PPRatings_data/MicrostatsxPPratings.xlsx")
View(MicrostatsxPP)

# Create the double subtracted data
MicrostatsxPP$Intsubtr = MicrostatsxPP$IntensityRatings - MicrostatsxPP$IntensityT0
MicrostatsxPP$Unpsubtr = MicrostatsxPP$UnpleasantnessRatings - MicrostatsxPP$UnpleasantnessT0

HFS  = MicrostatsxPP[MicrostatsxPP$Arm == "MFS",]
Cont = MicrostatsxPP[MicrostatsxPP$Arm == "Cont",]

HFS$DSint = HFS$Intsubtr - Cont$Intsubtr
HFS$DSUnp = HFS$Unpsubtr - Cont$Unpsubtr

MicrostatsxPP = HFS


MS1RS1 = MicrostatsxPP[MicrostatsxPP$MS == 'MS1' & MicrostatsxPP$RS == 'RS1',]
MS2RS1 = MicrostatsxPP[MicrostatsxPP$MS == 'MS2' & MicrostatsxPP$RS == 'RS1',]
MS3RS1 = MicrostatsxPP[MicrostatsxPP$MS == 'MS3' & MicrostatsxPP$RS == 'RS1',]
MS4RS1 = MicrostatsxPP[MicrostatsxPP$MS == 'MS4' & MicrostatsxPP$RS == 'RS1',]
MS5RS1 = MicrostatsxPP[MicrostatsxPP$MS == 'MS5' & MicrostatsxPP$RS == 'RS1',]

MS1RS2 = MicrostatsxPP[MicrostatsxPP$MS == 'MS1' & MicrostatsxPP$RS == 'RS2',]
MS2RS2 = MicrostatsxPP[MicrostatsxPP$MS == 'MS2' & MicrostatsxPP$RS == 'RS2',]
MS3RS2 = MicrostatsxPP[MicrostatsxPP$MS == 'MS3' & MicrostatsxPP$RS == 'RS2',]
MS4RS2 = MicrostatsxPP[MicrostatsxPP$MS == 'MS4' & MicrostatsxPP$RS == 'RS2',]
MS5RS2 = MicrostatsxPP[MicrostatsxPP$MS == 'MS5' & MicrostatsxPP$RS == 'RS2',]

MS1diff = MS1RS1
MS1diff[,c(6:16)] = MS1RS1[,c(6:16)] - MS1RS2[,c(6:16)]
MS2diff = MS2RS1
MS2diff[,c(6:16)] = MS2RS1[,c(6:16)] - MS2RS2[,c(6:16)]
MS3diff = MS3RS1
MS3diff[,c(6:16)] = MS3RS1[,c(6:16)] - MS3RS2[,c(6:16)]
MS4diff = MS4RS1
MS4diff[,c(6:16)] = MS4RS1[,c(6:16)] - MS4RS2[,c(6:16)]
MS5diff = MS5RS1
MS5diff[,c(6:16)] = MS5RS1[,c(6:16)] - MS5RS2[,c(6:16)]

### Intensity ----
# Initialize variables

# Create column and row names
col_names <- c("RS1","RS2","Diff")
row_names <- c("MS1", "MS2", "MS3", "MS4", "MS5")
# Create empty data frames
Duration_Pvals <- data.frame(matrix(ncol = length(col_names), nrow = length(row_names)))
colnames(Duration_Pvals) <- col_names
rownames(Duration_Pvals) <- row_names
Occurance_Pvals = Duration_Pvals
Coverage_Pvals = Duration_Pvals
GEV_Pvals = Duration_Pvals
TP1_Pvals = Duration_Pvals
TP2_Pvals = Duration_Pvals
TP3_Pvals = Duration_Pvals
TP4_Pvals = Duration_Pvals
TP5_Pvals = Duration_Pvals

## MS1
# RS1
Fullmod <- lmer(DSint ~ GFP  + Duration  + Occurance  + Coverage  + GEV  + TP2  + TP3  + TP4  + TP5  + (1|Subject), data = MS1RS1,REML = FALSE)
anova(Fullmod)
Duration_Pvals["MS1","RS1"]  = anova(Fullmod)['Duration','Pr(>F)']
Occurance_Pvals["MS1","RS1"] = anova(Fullmod)['Occurance','Pr(>F)']
Coverage_Pvals["MS1","RS1"]  = anova(Fullmod)['Coverage','Pr(>F)']
GEV_Pvals["MS1","RS1"]       = anova(Fullmod)['GEV','Pr(>F)']
TP1_Pvals["MS1","RS1"]       = anova(Fullmod)['TP1','Pr(>F)']
TP2_Pvals["MS1","RS1"]       = anova(Fullmod)['TP2','Pr(>F)']
TP3_Pvals["MS1","RS1"]       = anova(Fullmod)['TP3','Pr(>F)']
TP4_Pvals["MS1","RS1"]       = anova(Fullmod)['TP4','Pr(>F)']
TP5_Pvals["MS1","RS1"]       = anova(Fullmod)['TP5','Pr(>F)']

# Intplot <- ggpredict(Fullmod,terms = c('Occurance','Arm'))
# plot(Intplot,add.data = TRUE)

# RS2
Fullmod <- lmer(DSint ~ GFP  + Duration  + Occurance  + Coverage  + GEV  + TP2  + TP3  + TP4  + TP5  + (1|Subject), data = MS1RS2,REML = FALSE)
anova(Fullmod)
Duration_Pvals["MS1","RS2"]  = anova(Fullmod)['Duration','Pr(>F)']
Occurance_Pvals["MS1","RS2"] = anova(Fullmod)['Occurance','Pr(>F)']
Coverage_Pvals["MS1","RS2"]  = anova(Fullmod)['Coverage','Pr(>F)']
GEV_Pvals["MS1","RS2"]       = anova(Fullmod)['GEV','Pr(>F)']
TP1_Pvals["MS1","RS2"]       = anova(Fullmod)['TP1','Pr(>F)']
TP2_Pvals["MS1","RS2"]       = anova(Fullmod)['TP2','Pr(>F)']
TP3_Pvals["MS1","RS2"]       = anova(Fullmod)['TP3','Pr(>F)']
TP4_Pvals["MS1","RS2"]       = anova(Fullmod)['TP4','Pr(>F)']
TP5_Pvals["MS1","RS2"]       = anova(Fullmod)['TP5','Pr(>F)']
# Diff
Fullmod <- lmer(DSint ~ GFP  + Duration  + Occurance  + Coverage  + GEV  + TP2  + TP3  + TP4  + TP5  + (1|Subject), data = MS1diff,REML = FALSE)
anova(Fullmod)
Duration_Pvals["MS1","Diff"]  = anova(Fullmod)['Duration','Pr(>F)']
Occurance_Pvals["MS1","Diff"] = anova(Fullmod)['Occurance','Pr(>F)']
Coverage_Pvals["MS1","Diff"]  = anova(Fullmod)['Coverage','Pr(>F)']
GEV_Pvals["MS1","Diff"]       = anova(Fullmod)['GEV','Pr(>F)']
TP1_Pvals["MS1","Diff"]       = anova(Fullmod)['TP1','Pr(>F)']
TP2_Pvals["MS1","Diff"]       = anova(Fullmod)['TP2','Pr(>F)']
TP3_Pvals["MS1","Diff"]       = anova(Fullmod)['TP3','Pr(>F)']
TP4_Pvals["MS1","Diff"]       = anova(Fullmod)['TP4','Pr(>F)']
TP5_Pvals["MS1","Diff"]       = anova(Fullmod)['TP5','Pr(>F)']


## MS2
# RS1
Fullmod <- lmer(DSint ~ GFP  + Duration  + Occurance  + Coverage  + GEV  + TP1  + TP3  + TP4  + TP5  + (1|Subject), data = MS2RS1,REML = FALSE)
anova(Fullmod)
Duration_Pvals["MS2","RS1"]  = anova(Fullmod)['Duration','Pr(>F)']
Occurance_Pvals["MS2","RS1"] = anova(Fullmod)['Occurance','Pr(>F)']
Coverage_Pvals["MS2","RS1"]  = anova(Fullmod)['Coverage','Pr(>F)']
GEV_Pvals["MS2","RS1"]       = anova(Fullmod)['GEV','Pr(>F)']
TP1_Pvals["MS2","RS1"]       = anova(Fullmod)['TP1','Pr(>F)']
TP2_Pvals["MS2","RS1"]       = anova(Fullmod)['TP2','Pr(>F)']
TP3_Pvals["MS2","RS1"]       = anova(Fullmod)['TP3','Pr(>F)']
TP4_Pvals["MS2","RS1"]       = anova(Fullmod)['TP4','Pr(>F)']
TP5_Pvals["MS2","RS1"]       = anova(Fullmod)['TP5','Pr(>F)']

# RS2
Fullmod <- lmer(DSint ~ GFP  + Duration  + Occurance  + Coverage  + GEV  + TP1  + TP3  + TP4  + TP5  + (1|Subject), data = MS2RS2,REML = FALSE)
anova(Fullmod)
Duration_Pvals["MS2","RS2"]  = anova(Fullmod)['Duration','Pr(>F)']
Occurance_Pvals["MS2","RS2"] = anova(Fullmod)['Occurance','Pr(>F)']
Coverage_Pvals["MS2","RS2"]  = anova(Fullmod)['Coverage','Pr(>F)']
GEV_Pvals["MS2","RS2"]       = anova(Fullmod)['GEV','Pr(>F)']
TP1_Pvals["MS2","RS2"]       = anova(Fullmod)['TP1','Pr(>F)']
TP2_Pvals["MS2","RS2"]       = anova(Fullmod)['TP2','Pr(>F)']
TP3_Pvals["MS2","RS2"]       = anova(Fullmod)['TP3','Pr(>F)']
TP4_Pvals["MS2","RS2"]       = anova(Fullmod)['TP4','Pr(>F)']
TP5_Pvals["MS2","RS2"]       = anova(Fullmod)['TP5','Pr(>F)']

# Diff
Fullmod <- lmer(DSint ~ GFP  + Duration  + Occurance  + Coverage  + GEV  + TP1  + TP3  + TP4  + TP5  + (1|Subject), data = MS2diff,REML = FALSE)
anova(Fullmod)
Duration_Pvals["MS2","Diff"]  = anova(Fullmod)['Duration','Pr(>F)']
Occurance_Pvals["MS2","Diff"] = anova(Fullmod)['Occurance','Pr(>F)']
Coverage_Pvals["MS2","Diff"]  = anova(Fullmod)['Coverage','Pr(>F)']
GEV_Pvals["MS2","Diff"]       = anova(Fullmod)['GEV','Pr(>F)']
TP1_Pvals["MS2","Diff"]       = anova(Fullmod)['TP1','Pr(>F)']
TP2_Pvals["MS2","Diff"]       = anova(Fullmod)['TP2','Pr(>F)']
TP3_Pvals["MS2","Diff"]       = anova(Fullmod)['TP3','Pr(>F)']
TP4_Pvals["MS2","Diff"]       = anova(Fullmod)['TP4','Pr(>F)']
TP5_Pvals["MS2","Diff"]       = anova(Fullmod)['TP5','Pr(>F)']


## MS3
# RS1
Fullmod <- lmer(DSint ~ GFP  + Duration  + Occurance  + Coverage  + GEV  + TP2  + TP1  + TP4  + TP5  + (1|Subject), data = MS3RS1,REML = FALSE)
anova(Fullmod)
Duration_Pvals["MS3","RS1"]  = anova(Fullmod)['Duration','Pr(>F)']
Occurance_Pvals["MS3","RS1"] = anova(Fullmod)['Occurance','Pr(>F)']
Coverage_Pvals["MS3","RS1"]  = anova(Fullmod)['Coverage','Pr(>F)']
GEV_Pvals["MS3","RS1"]       = anova(Fullmod)['GEV','Pr(>F)']
TP1_Pvals["MS3","RS1"]       = anova(Fullmod)['TP1','Pr(>F)']
TP2_Pvals["MS3","RS1"]       = anova(Fullmod)['TP2','Pr(>F)']
TP3_Pvals["MS3","RS1"]       = anova(Fullmod)['TP3','Pr(>F)']
TP4_Pvals["MS3","RS1"]       = anova(Fullmod)['TP4','Pr(>F)']
TP5_Pvals["MS3","RS1"]       = anova(Fullmod)['TP5','Pr(>F)']

# RS2
Fullmod <- lmer(DSint ~ GFP  + Duration  + Occurance  + Coverage  + GEV  + TP2  + TP1  + TP4  + TP5  + (1|Subject), data = MS3RS2,REML = FALSE)
anova(Fullmod)
Duration_Pvals["MS3","RS2"]  = anova(Fullmod)['Duration','Pr(>F)']
Occurance_Pvals["MS3","RS2"] = anova(Fullmod)['Occurance','Pr(>F)']
Coverage_Pvals["MS3","RS2"]  = anova(Fullmod)['Coverage','Pr(>F)']
GEV_Pvals["MS3","RS2"]       = anova(Fullmod)['GEV','Pr(>F)']
TP1_Pvals["MS3","RS2"]       = anova(Fullmod)['TP1','Pr(>F)']
TP2_Pvals["MS3","RS2"]       = anova(Fullmod)['TP2','Pr(>F)']
TP3_Pvals["MS3","RS2"]       = anova(Fullmod)['TP3','Pr(>F)']
TP4_Pvals["MS3","RS2"]       = anova(Fullmod)['TP4','Pr(>F)']
TP5_Pvals["MS3","RS2"]       = anova(Fullmod)['TP5','Pr(>F)']

# Diff
Fullmod <- lmer(DSint ~ GFP  + Duration  + Occurance  + Coverage  + GEV  + TP2  + TP1  + TP4  + TP5  + (1|Subject), data = MS3diff,REML = FALSE)
anova(Fullmod)
Duration_Pvals["MS3","Diff"]  = anova(Fullmod)['Duration','Pr(>F)']
Occurance_Pvals["MS3","Diff"] = anova(Fullmod)['Occurance','Pr(>F)']
Coverage_Pvals["MS3","Diff"]  = anova(Fullmod)['Coverage','Pr(>F)']
GEV_Pvals["MS3","Diff"]       = anova(Fullmod)['GEV','Pr(>F)']
TP1_Pvals["MS3","Diff"]       = anova(Fullmod)['TP1','Pr(>F)']
TP2_Pvals["MS3","Diff"]       = anova(Fullmod)['TP2','Pr(>F)']
TP3_Pvals["MS3","Diff"]       = anova(Fullmod)['TP3','Pr(>F)']
TP4_Pvals["MS3","Diff"]       = anova(Fullmod)['TP4','Pr(>F)']
TP5_Pvals["MS3","Diff"]       = anova(Fullmod)['TP5','Pr(>F)']


## MS4
# RS1
Fullmod <- lmer(DSint ~ GFP  + Duration  + Occurance  + Coverage  + GEV  + TP2  + TP3  + TP1  + TP5  + (1|Subject), data = MS4RS1,REML = FALSE)
anova(Fullmod)
Duration_Pvals["MS4","RS1"]  = anova(Fullmod)['Duration','Pr(>F)']
Occurance_Pvals["MS4","RS1"] = anova(Fullmod)['Occurance','Pr(>F)']
Coverage_Pvals["MS4","RS1"]  = anova(Fullmod)['Coverage','Pr(>F)']
GEV_Pvals["MS4","RS1"]       = anova(Fullmod)['GEV','Pr(>F)']
TP1_Pvals["MS4","RS1"]       = anova(Fullmod)['TP1','Pr(>F)']
TP2_Pvals["MS4","RS1"]       = anova(Fullmod)['TP2','Pr(>F)']
TP3_Pvals["MS4","RS1"]       = anova(Fullmod)['TP3','Pr(>F)']
TP4_Pvals["MS4","RS1"]       = anova(Fullmod)['TP4','Pr(>F)']
TP5_Pvals["MS4","RS1"]       = anova(Fullmod)['TP5','Pr(>F)']

# RS2
Fullmod <- lmer(DSint ~ GFP  + Duration  + Occurance  + Coverage  + GEV  + TP2  + TP3  + TP1  + TP5  + (1|Subject), data = MS4RS2,REML = FALSE)
anova(Fullmod)
Duration_Pvals["MS4","RS2"]  = anova(Fullmod)['Duration','Pr(>F)']
Occurance_Pvals["MS4","RS2"] = anova(Fullmod)['Occurance','Pr(>F)']
Coverage_Pvals["MS4","RS2"]  = anova(Fullmod)['Coverage','Pr(>F)']
GEV_Pvals["MS4","RS2"]       = anova(Fullmod)['GEV','Pr(>F)']
TP1_Pvals["MS4","RS2"]       = anova(Fullmod)['TP1','Pr(>F)']
TP2_Pvals["MS4","RS2"]       = anova(Fullmod)['TP2','Pr(>F)']
TP3_Pvals["MS4","RS2"]       = anova(Fullmod)['TP3','Pr(>F)']
TP4_Pvals["MS4","RS2"]       = anova(Fullmod)['TP4','Pr(>F)']
TP5_Pvals["MS4","RS2"]       = anova(Fullmod)['TP5','Pr(>F)']

# Diff
Fullmod <- lmer(DSint ~ GFP  + Duration  + Occurance  + Coverage  + GEV  + TP2  + TP3  + TP1  + TP5  + (1|Subject), data = MS4diff,REML = FALSE)
anova(Fullmod)
Duration_Pvals["MS4","Diff"]  = anova(Fullmod)['Duration','Pr(>F)']
Occurance_Pvals["MS4","Diff"] = anova(Fullmod)['Occurance','Pr(>F)']
Coverage_Pvals["MS4","Diff"]  = anova(Fullmod)['Coverage','Pr(>F)']
GEV_Pvals["MS4","Diff"]       = anova(Fullmod)['GEV','Pr(>F)']
TP1_Pvals["MS4","Diff"]       = anova(Fullmod)['TP1','Pr(>F)']
TP2_Pvals["MS4","Diff"]       = anova(Fullmod)['TP2','Pr(>F)']
TP3_Pvals["MS4","Diff"]       = anova(Fullmod)['TP3','Pr(>F)']
TP4_Pvals["MS4","Diff"]       = anova(Fullmod)['TP4','Pr(>F)']
TP5_Pvals["MS4","Diff"]       = anova(Fullmod)['TP5','Pr(>F)']


## MS5
# RS1
Fullmod <- lmer(DSint ~ GFP  + Duration  + Occurance  + Coverage  + GEV  + TP2  + TP3  + TP4  + TP1  + (1|Subject), data = MS5RS1,REML = FALSE)
anova(Fullmod)
anova(Fullmod)
Duration_Pvals["MS5","RS1"]  = anova(Fullmod)['Duration','Pr(>F)']
Occurance_Pvals["MS5","RS1"] = anova(Fullmod)['Occurance','Pr(>F)']
Coverage_Pvals["MS5","RS1"]  = anova(Fullmod)['Coverage','Pr(>F)']
GEV_Pvals["MS5","RS1"]       = anova(Fullmod)['GEV','Pr(>F)']
TP1_Pvals["MS5","RS1"]       = anova(Fullmod)['TP1','Pr(>F)']
TP2_Pvals["MS5","RS1"]       = anova(Fullmod)['TP2','Pr(>F)']
TP3_Pvals["MS5","RS1"]       = anova(Fullmod)['TP3','Pr(>F)']
TP4_Pvals["MS5","RS1"]       = anova(Fullmod)['TP4','Pr(>F)']
TP5_Pvals["MS5","RS1"]       = anova(Fullmod)['TP5','Pr(>F)']

# RS2
Fullmod <- lmer(DSint ~ GFP  + Duration  + Occurance  + Coverage  + GEV  + TP2  + TP3  + TP4  + TP1  + (1|Subject), data = MS5RS2,REML = FALSE)
anova(Fullmod)
Duration_Pvals["MS5","RS2"]  = anova(Fullmod)['Duration','Pr(>F)']
Occurance_Pvals["MS5","RS2"] = anova(Fullmod)['Occurance','Pr(>F)']
Coverage_Pvals["MS5","RS2"]  = anova(Fullmod)['Coverage','Pr(>F)']
GEV_Pvals["MS5","RS2"]       = anova(Fullmod)['GEV','Pr(>F)']
TP1_Pvals["MS5","RS2"]       = anova(Fullmod)['TP1','Pr(>F)']
TP2_Pvals["MS5","RS2"]       = anova(Fullmod)['TP2','Pr(>F)']
TP3_Pvals["MS5","RS2"]       = anova(Fullmod)['TP3','Pr(>F)']
TP4_Pvals["MS5","RS2"]       = anova(Fullmod)['TP4','Pr(>F)']
TP5_Pvals["MS5","RS2"]       = anova(Fullmod)['TP5','Pr(>F)']

# Diff
Fullmod <- lmer(DSint ~ GFP  + Duration  + Occurance  + Coverage  + GEV  + TP2  + TP3  + TP4  + TP1  + (1|Subject), data = MS5diff,REML = FALSE)
anova(Fullmod)
Duration_Pvals["MS5","Diff"]  = anova(Fullmod)['Duration','Pr(>F)']
Occurance_Pvals["MS5","Diff"] = anova(Fullmod)['Occurance','Pr(>F)']
Coverage_Pvals["MS5","Diff"]  = anova(Fullmod)['Coverage','Pr(>F)']
GEV_Pvals["MS5","Diff"]       = anova(Fullmod)['GEV','Pr(>F)']
TP1_Pvals["MS5","Diff"]       = anova(Fullmod)['TP1','Pr(>F)']
TP2_Pvals["MS5","Diff"]       = anova(Fullmod)['TP2','Pr(>F)']
TP3_Pvals["MS5","Diff"]       = anova(Fullmod)['TP3','Pr(>F)']
TP4_Pvals["MS5","Diff"]       = anova(Fullmod)['TP4','Pr(>F)']
TP5_Pvals["MS5","Diff"]       = anova(Fullmod)['TP5','Pr(>F)']

# Check for significance after FDR correction
Pvalues <- unlist(TP5_Pvals)
non_na_indices <- which(!is.na(Pvalues))
All_adj_p_values <- p.adjust(Pvalues[non_na_indices], method = "fdr")
Pvalues[non_na_indices] <- All_adj_p_values
Pvalues[Pvalues<.05]

### Unpleasantness ----
# Initialize variables

# Create column and row names
col_names <- c("RS1","RS2","Diff")
row_names <- c("MS1", "MS2", "MS3", "MS4", "MS5")
# Create empty data frames
Duration_Pvals <- data.frame(matrix(ncol = length(col_names), nrow = length(row_names)))
colnames(Duration_Pvals) <- col_names
rownames(Duration_Pvals) <- row_names
Occurance_Pvals = Duration_Pvals
Coverage_Pvals = Duration_Pvals
GEV_Pvals = Duration_Pvals
TP1_Pvals = Duration_Pvals
TP2_Pvals = Duration_Pvals
TP3_Pvals = Duration_Pvals
TP4_Pvals = Duration_Pvals
TP5_Pvals = Duration_Pvals

## MS1
# RS1
Fullmod <- lmer(DSUnp ~ GFP  + Duration  + Occurance  + Coverage  + GEV  + TP2  + TP3  + TP4  + TP5  + (1|Subject), data = MS1RS1,REML = FALSE)
anova(Fullmod)
Duration_Pvals["MS1","RS1"]  = anova(Fullmod)['Duration','Pr(>F)']
Occurance_Pvals["MS1","RS1"] = anova(Fullmod)['Occurance','Pr(>F)']
Coverage_Pvals["MS1","RS1"]  = anova(Fullmod)['Coverage','Pr(>F)']
GEV_Pvals["MS1","RS1"]       = anova(Fullmod)['GEV','Pr(>F)']
TP1_Pvals["MS1","RS1"]       = anova(Fullmod)['TP1','Pr(>F)']
TP2_Pvals["MS1","RS1"]       = anova(Fullmod)['TP2','Pr(>F)']
TP3_Pvals["MS1","RS1"]       = anova(Fullmod)['TP3','Pr(>F)']
TP4_Pvals["MS1","RS1"]       = anova(Fullmod)['TP4','Pr(>F)']
TP5_Pvals["MS1","RS1"]       = anova(Fullmod)['TP5','Pr(>F)']

# RS2
Fullmod <- lmer(DSUnp ~ GFP  + Duration  + Occurance  + Coverage  + GEV  + TP2  + TP3  + TP4  + TP5  + (1|Subject), data = MS1RS2,REML = FALSE)
anova(Fullmod)
Duration_Pvals["MS1","RS2"]  = anova(Fullmod)['Duration','Pr(>F)']
Occurance_Pvals["MS1","RS2"] = anova(Fullmod)['Occurance','Pr(>F)']
Coverage_Pvals["MS1","RS2"]  = anova(Fullmod)['Coverage','Pr(>F)']
GEV_Pvals["MS1","RS2"]       = anova(Fullmod)['GEV','Pr(>F)']
TP1_Pvals["MS1","RS2"]       = anova(Fullmod)['TP1','Pr(>F)']
TP2_Pvals["MS1","RS2"]       = anova(Fullmod)['TP2','Pr(>F)']
TP3_Pvals["MS1","RS2"]       = anova(Fullmod)['TP3','Pr(>F)']
TP4_Pvals["MS1","RS2"]       = anova(Fullmod)['TP4','Pr(>F)']
TP5_Pvals["MS1","RS2"]       = anova(Fullmod)['TP5','Pr(>F)']

# Diff
Fullmod <- lmer(DSUnp ~ GFP  + Duration  + Occurance  + Coverage  + GEV  + TP2  + TP3  + TP4  + TP5  + (1|Subject), data = MS1diff,REML = FALSE)
anova(Fullmod)
Duration_Pvals["MS1","Diff"]  = anova(Fullmod)['Duration','Pr(>F)']
Occurance_Pvals["MS1","Diff"] = anova(Fullmod)['Occurance','Pr(>F)']
Coverage_Pvals["MS1","Diff"]  = anova(Fullmod)['Coverage','Pr(>F)']
GEV_Pvals["MS1","Diff"]       = anova(Fullmod)['GEV','Pr(>F)']
TP1_Pvals["MS1","Diff"]       = anova(Fullmod)['TP1','Pr(>F)']
TP2_Pvals["MS1","Diff"]       = anova(Fullmod)['TP2','Pr(>F)']
TP3_Pvals["MS1","Diff"]       = anova(Fullmod)['TP3','Pr(>F)']
TP4_Pvals["MS1","Diff"]       = anova(Fullmod)['TP4','Pr(>F)']
TP5_Pvals["MS1","Diff"]       = anova(Fullmod)['TP5','Pr(>F)']


## MS2
# RS1
Fullmod <- lmer(DSUnp ~ GFP  + Duration  + Occurance  + Coverage  + GEV  + TP1  + TP3  + TP4  + TP5  + (1|Subject), data = MS2RS1,REML = FALSE)
anova(Fullmod)
Duration_Pvals["MS2","RS1"]  = anova(Fullmod)['Duration','Pr(>F)']
Occurance_Pvals["MS2","RS1"] = anova(Fullmod)['Occurance','Pr(>F)']
Coverage_Pvals["MS2","RS1"]  = anova(Fullmod)['Coverage','Pr(>F)']
GEV_Pvals["MS2","RS1"]       = anova(Fullmod)['GEV','Pr(>F)']
TP1_Pvals["MS2","RS1"]       = anova(Fullmod)['TP1','Pr(>F)']
TP2_Pvals["MS2","RS1"]       = anova(Fullmod)['TP2','Pr(>F)']
TP3_Pvals["MS2","RS1"]       = anova(Fullmod)['TP3','Pr(>F)']
TP4_Pvals["MS2","RS1"]       = anova(Fullmod)['TP4','Pr(>F)']
TP5_Pvals["MS2","RS1"]       = anova(Fullmod)['TP5','Pr(>F)']

# RS2
Fullmod <- lmer(DSUnp ~ GFP  + Duration  + Occurance  + Coverage  + GEV  + TP1  + TP3  + TP4  + TP5  + (1|Subject), data = MS2RS2,REML = FALSE)
anova(Fullmod)
Duration_Pvals["MS2","RS2"]  = anova(Fullmod)['Duration','Pr(>F)']
Occurance_Pvals["MS2","RS2"] = anova(Fullmod)['Occurance','Pr(>F)']
Coverage_Pvals["MS2","RS2"]  = anova(Fullmod)['Coverage','Pr(>F)']
GEV_Pvals["MS2","RS2"]       = anova(Fullmod)['GEV','Pr(>F)']
TP1_Pvals["MS2","RS2"]       = anova(Fullmod)['TP1','Pr(>F)']
TP2_Pvals["MS2","RS2"]       = anova(Fullmod)['TP2','Pr(>F)']
TP3_Pvals["MS2","RS2"]       = anova(Fullmod)['TP3','Pr(>F)']
TP4_Pvals["MS2","RS2"]       = anova(Fullmod)['TP4','Pr(>F)']
TP5_Pvals["MS2","RS2"]       = anova(Fullmod)['TP5','Pr(>F)']

# Diff
Fullmod <- lmer(DSUnp ~ GFP  + Duration  + Occurance  + Coverage  + GEV  + TP1  + TP3  + TP4  + TP5  + (1|Subject), data = MS2diff,REML = FALSE)
anova(Fullmod)
Duration_Pvals["MS2","Diff"]  = anova(Fullmod)['Duration','Pr(>F)']
Occurance_Pvals["MS2","Diff"] = anova(Fullmod)['Occurance','Pr(>F)']
Coverage_Pvals["MS2","Diff"]  = anova(Fullmod)['Coverage','Pr(>F)']
GEV_Pvals["MS2","Diff"]       = anova(Fullmod)['GEV','Pr(>F)']
TP1_Pvals["MS2","Diff"]       = anova(Fullmod)['TP1','Pr(>F)']
TP2_Pvals["MS2","Diff"]       = anova(Fullmod)['TP2','Pr(>F)']
TP3_Pvals["MS2","Diff"]       = anova(Fullmod)['TP3','Pr(>F)']
TP4_Pvals["MS2","Diff"]       = anova(Fullmod)['TP4','Pr(>F)']
TP5_Pvals["MS2","Diff"]       = anova(Fullmod)['TP5','Pr(>F)']


## MS3
# RS1
Fullmod <- lmer(DSUnp ~ GFP  + Duration  + Occurance  + Coverage  + GEV  + TP2  + TP1  + TP4  + TP5  + (1|Subject), data = MS3RS1,REML = FALSE)
anova(Fullmod)
Duration_Pvals["MS3","RS1"]  = anova(Fullmod)['Duration','Pr(>F)']
Occurance_Pvals["MS3","RS1"] = anova(Fullmod)['Occurance','Pr(>F)']
Coverage_Pvals["MS3","RS1"]  = anova(Fullmod)['Coverage','Pr(>F)']
GEV_Pvals["MS3","RS1"]       = anova(Fullmod)['GEV','Pr(>F)']
TP1_Pvals["MS3","RS1"]       = anova(Fullmod)['TP1','Pr(>F)']
TP2_Pvals["MS3","RS1"]       = anova(Fullmod)['TP2','Pr(>F)']
TP3_Pvals["MS3","RS1"]       = anova(Fullmod)['TP3','Pr(>F)']
TP4_Pvals["MS3","RS1"]       = anova(Fullmod)['TP4','Pr(>F)']
TP5_Pvals["MS3","RS1"]       = anova(Fullmod)['TP5','Pr(>F)']

# RS2
Fullmod <- lmer(DSUnp ~ GFP  + Duration  + Occurance  + Coverage  + GEV  + TP2  + TP1  + TP4  + TP5  + (1|Subject), data = MS3RS2,REML = FALSE)
anova(Fullmod)
Duration_Pvals["MS3","RS2"]  = anova(Fullmod)['Duration','Pr(>F)']
Occurance_Pvals["MS3","RS2"] = anova(Fullmod)['Occurance','Pr(>F)']
Coverage_Pvals["MS3","RS2"]  = anova(Fullmod)['Coverage','Pr(>F)']
GEV_Pvals["MS3","RS2"]       = anova(Fullmod)['GEV','Pr(>F)']
TP1_Pvals["MS3","RS2"]       = anova(Fullmod)['TP1','Pr(>F)']
TP2_Pvals["MS3","RS2"]       = anova(Fullmod)['TP2','Pr(>F)']
TP3_Pvals["MS3","RS2"]       = anova(Fullmod)['TP3','Pr(>F)']
TP4_Pvals["MS3","RS2"]       = anova(Fullmod)['TP4','Pr(>F)']
TP5_Pvals["MS3","RS2"]       = anova(Fullmod)['TP5','Pr(>F)']

# Diff
Fullmod <- lmer(DSUnp ~ GFP  + Duration  + Occurance  + Coverage  + GEV  + TP2  + TP1  + TP4  + TP5  + (1|Subject), data = MS3diff,REML = FALSE)
anova(Fullmod)
Duration_Pvals["MS3","Diff"]  = anova(Fullmod)['Duration','Pr(>F)']
Occurance_Pvals["MS3","Diff"] = anova(Fullmod)['Occurance','Pr(>F)']
Coverage_Pvals["MS3","Diff"]  = anova(Fullmod)['Coverage','Pr(>F)']
GEV_Pvals["MS3","Diff"]       = anova(Fullmod)['GEV','Pr(>F)']
TP1_Pvals["MS3","Diff"]       = anova(Fullmod)['TP1','Pr(>F)']
TP2_Pvals["MS3","Diff"]       = anova(Fullmod)['TP2','Pr(>F)']
TP3_Pvals["MS3","Diff"]       = anova(Fullmod)['TP3','Pr(>F)']
TP4_Pvals["MS3","Diff"]       = anova(Fullmod)['TP4','Pr(>F)']
TP5_Pvals["MS3","Diff"]       = anova(Fullmod)['TP5','Pr(>F)']



## MS4
# RS1
Fullmod <- lmer(DSUnp ~ GFP  + Duration  + Occurance  + Coverage  + GEV  + TP2  + TP3  + TP1  + TP5  + (1|Subject), data = MS4RS1,REML = FALSE)
anova(Fullmod)
Duration_Pvals["MS4","RS1"]  = anova(Fullmod)['Duration','Pr(>F)']
Occurance_Pvals["MS4","RS1"] = anova(Fullmod)['Occurance','Pr(>F)']
Coverage_Pvals["MS4","RS1"]  = anova(Fullmod)['Coverage','Pr(>F)']
GEV_Pvals["MS4","RS1"]       = anova(Fullmod)['GEV','Pr(>F)']
TP1_Pvals["MS4","RS1"]       = anova(Fullmod)['TP1','Pr(>F)']
TP2_Pvals["MS4","RS1"]       = anova(Fullmod)['TP2','Pr(>F)']
TP3_Pvals["MS4","RS1"]       = anova(Fullmod)['TP3','Pr(>F)']
TP4_Pvals["MS4","RS1"]       = anova(Fullmod)['TP4','Pr(>F)']
TP5_Pvals["MS4","RS1"]       = anova(Fullmod)['TP5','Pr(>F)']

# RS2
Fullmod <- lmer(DSUnp ~ GFP  + Duration  + Occurance  + Coverage  + GEV  + TP2  + TP3  + TP1  + TP5  + (1|Subject), data = MS4RS2,REML = FALSE)
anova(Fullmod)
Duration_Pvals["MS4","RS2"]  = anova(Fullmod)['Duration','Pr(>F)']
Occurance_Pvals["MS4","RS2"] = anova(Fullmod)['Occurance','Pr(>F)']
Coverage_Pvals["MS4","RS2"]  = anova(Fullmod)['Coverage','Pr(>F)']
GEV_Pvals["MS4","RS2"]       = anova(Fullmod)['GEV','Pr(>F)']
TP1_Pvals["MS4","RS2"]       = anova(Fullmod)['TP1','Pr(>F)']
TP2_Pvals["MS4","RS2"]       = anova(Fullmod)['TP2','Pr(>F)']
TP3_Pvals["MS4","RS2"]       = anova(Fullmod)['TP3','Pr(>F)']
TP4_Pvals["MS4","RS2"]       = anova(Fullmod)['TP4','Pr(>F)']
TP5_Pvals["MS4","RS2"]       = anova(Fullmod)['TP5','Pr(>F)']

# Diff
Fullmod <- lmer(DSUnp ~ GFP  + Duration  + Occurance  + Coverage  + GEV  + TP2  + TP3  + TP1  + TP5  + (1|Subject), data = MS4diff,REML = FALSE)
anova(Fullmod)
Duration_Pvals["MS4","Diff"]  = anova(Fullmod)['Duration','Pr(>F)']
Occurance_Pvals["MS4","Diff"] = anova(Fullmod)['Occurance','Pr(>F)']
Coverage_Pvals["MS4","Diff"]  = anova(Fullmod)['Coverage','Pr(>F)']
GEV_Pvals["MS4","Diff"]       = anova(Fullmod)['GEV','Pr(>F)']
TP1_Pvals["MS4","Diff"]       = anova(Fullmod)['TP1','Pr(>F)']
TP2_Pvals["MS4","Diff"]       = anova(Fullmod)['TP2','Pr(>F)']
TP3_Pvals["MS4","Diff"]       = anova(Fullmod)['TP3','Pr(>F)']
TP4_Pvals["MS4","Diff"]       = anova(Fullmod)['TP4','Pr(>F)']
TP5_Pvals["MS4","Diff"]       = anova(Fullmod)['TP5','Pr(>F)']

## MS5
# RS1
Fullmod <- lmer(DSUnp ~ GFP  + Duration  + Occurance  + Coverage  + GEV  + TP2  + TP3  + TP4  + TP1  + (1|Subject), data = MS5RS1,REML = FALSE)
anova(Fullmod)
Duration_Pvals["MS5","RS1"]  = anova(Fullmod)['Duration','Pr(>F)']
Occurance_Pvals["MS5","RS1"] = anova(Fullmod)['Occurance','Pr(>F)']
Coverage_Pvals["MS5","RS1"]  = anova(Fullmod)['Coverage','Pr(>F)']
GEV_Pvals["MS5","RS1"]       = anova(Fullmod)['GEV','Pr(>F)']
TP1_Pvals["MS5","RS1"]       = anova(Fullmod)['TP1','Pr(>F)']
TP2_Pvals["MS5","RS1"]       = anova(Fullmod)['TP2','Pr(>F)']
TP3_Pvals["MS5","RS1"]       = anova(Fullmod)['TP3','Pr(>F)']
TP4_Pvals["MS5","RS1"]       = anova(Fullmod)['TP4','Pr(>F)']
TP5_Pvals["MS5","RS1"]       = anova(Fullmod)['TP5','Pr(>F)']

# RS2
Fullmod <- lmer(DSUnp ~ GFP  + Duration  + Occurance  + Coverage  + GEV  + TP2  + TP3  + TP4  + TP1  + (1|Subject), data = MS5RS2,REML = FALSE)
anova(Fullmod)
Duration_Pvals["MS5","RS2"]  = anova(Fullmod)['Duration','Pr(>F)']
Occurance_Pvals["MS5","RS2"] = anova(Fullmod)['Occurance','Pr(>F)']
Coverage_Pvals["MS5","RS2"]  = anova(Fullmod)['Coverage','Pr(>F)']
GEV_Pvals["MS5","RS2"]       = anova(Fullmod)['GEV','Pr(>F)']
TP1_Pvals["MS5","RS2"]       = anova(Fullmod)['TP1','Pr(>F)']
TP2_Pvals["MS5","RS2"]       = anova(Fullmod)['TP2','Pr(>F)']
TP3_Pvals["MS5","RS2"]       = anova(Fullmod)['TP3','Pr(>F)']
TP4_Pvals["MS5","RS2"]       = anova(Fullmod)['TP4','Pr(>F)']
TP5_Pvals["MS5","RS2"]       = anova(Fullmod)['TP5','Pr(>F)']

# Diff
Fullmod <- lmer(DSUnp ~ GFP  + Duration  + Occurance  + Coverage  + GEV  + TP2  + TP3  + TP4  + TP1  + (1|Subject), data = MS5diff,REML = FALSE)
anova(Fullmod)
Duration_Pvals["MS5","Diff"]  = anova(Fullmod)['Duration','Pr(>F)']
Occurance_Pvals["MS5","Diff"] = anova(Fullmod)['Occurance','Pr(>F)']
Coverage_Pvals["MS5","Diff"]  = anova(Fullmod)['Coverage','Pr(>F)']
GEV_Pvals["MS5","Diff"]       = anova(Fullmod)['GEV','Pr(>F)']
TP1_Pvals["MS5","Diff"]       = anova(Fullmod)['TP1','Pr(>F)']
TP2_Pvals["MS5","Diff"]       = anova(Fullmod)['TP2','Pr(>F)']
TP3_Pvals["MS5","Diff"]       = anova(Fullmod)['TP3','Pr(>F)']
TP4_Pvals["MS5","Diff"]       = anova(Fullmod)['TP4','Pr(>F)']
TP5_Pvals["MS5","Diff"]       = anova(Fullmod)['TP5','Pr(>F)']

# Check for significance after FDR correction
Pvalues <- unlist(TP5_Pvals)
non_na_indices <- which(!is.na(Pvalues))
All_adj_p_values <- p.adjust(Pvalues[non_na_indices], method = "fdr")
Pvalues[non_na_indices] <- All_adj_p_values
Pvalues[Pvalues<.05]



# Load Single Clustering data (Control) ----
Microstats <- read_excel("Path_to_Microstate_results/Microstats.xlsx")
View(Microstats)
Microstats$Subject<- as.factor(Microstats$Subject)


Durmod <- lm(Duration ~ MS*RS, data = Microstats)
summary(Durmod)

Occurmod <- lm((Occurance) ~ MS*RS, data = Microstats)
summary(Occurmod)

hist((MS1$Coverage))
Covermod <- lm((Coverage) ~ MS*RS, data = Microstats)
summary(Covermod)

hist((MS1$GEV))
Gevmod <- lm(GEV ~ MS*RS, data = Microstats)
summary(Gevmod)

hist((MS1$TP2))
TP2mod <- lm(TP2 ~ MS*RS, data = Microstats)
summary(TP2mod)

TP3mod <- lm(TP3 ~ MS*RS, data = Microstats)
summary(TP3mod)

TP4mod <- lm((TP4) ~ MS*RS, data = Microstats)
summary(TP4mod)

TP5mod <- lm((TP5) ~ MS*RS, data = Microstats)
summary(TP5mod)


MS1 = Microstats[Microstats$MS == 'MS1',]
MS2 = Microstats[Microstats$MS == 'MS2',]
MS3 = Microstats[Microstats$MS == 'MS3',]
MS4 = Microstats[Microstats$MS == 'MS4',]
MS5 = Microstats[Microstats$MS == 'MS5',]

TotalGEVmod <- lm(TotalGEV ~ RS, data = MS1)
summary(TotalGEVmod)

t.test(TotalGEV ~ RS, data = MS1)

mean(MS2$TotalGEV[MS2$RS == 'RS2'])
sd(MS2$TotalGEV[MS2$RS == 'RS2'])

# Create column and row names
col_names <- c("Duration", "Occurrence", "Coverage", "GEV", "TP1", "TP2", "TP3", "TP4", "TP5")
row_names <- c("MS1", "MS2", "MS3", "MS4", "MS5")
# Create empty data frame
ALL_Pvals <- data.frame(matrix(ncol = length(col_names), nrow = length(row_names)))
colnames(ALL_Pvals) <- col_names
rownames(ALL_Pvals) <- row_names

# Incorporate this into the function
# full_BF = lmBF(IntensityRatings ~  Baselineratings + Condition*Arm*Time + Participant ,data=Intensity2, whichRandom = 'Participant')
# null_BF = lmBF(IntensityRatings ~  Baselineratings + Condition + Arm*Time + Participant ,data=Intensity2, whichRandom = 'Participant')
# full_BF/null_BF

# MS1
hist((MS1$Duration))
Durmod <- lmer(Duration ~ RS + (1|Subject), data = MS1,REML = FALSE)
summary(Durmod)


full_Durmod = lmBF(Duration ~  RS + Subject, data=MS1, whichRandom = 'Subject')
null_Durmod = lmBF(Duration ~ 1 + Subject, data = MS1, whichRandom = 'Subject')
BF_Duration = full_Durmod/null_Durmod

hist((MS1$Occurance))
Occurmod <- lmer(Occurance ~ RS + (1|Subject), data = MS1,REML = FALSE)
summary(Occurmod)

t.test(Occurance ~ RS, data = MS1, paired = TRUE)


full_Occurmod = lmBF(Occurance ~  RS + Subject, data=MS1, whichRandom = 'Subject')
null_Occurmod = lmBF(Occurance ~ 1 + Subject, data = MS1, whichRandom = 'Subject')
BF_Occurance = full_Occurmod/null_Occurmod

hist((MS1$Coverage))
Covermod <- lmer((Coverage) ~ RS + (1|Subject), data = MS1,REML = FALSE)
summary(Covermod)

full_Covermod = lmBF(Coverage ~  RS + Subject, data=MS1, whichRandom = 'Subject')
null_Covermod = lmBF(Coverage ~ 1 + Subject, data = MS1, whichRandom = 'Subject')
BF_Coverage = full_Covermod/null_Covermod

hist((MS1$GEV))
Gevmod <- lmer(GEV ~ RS + (1|Subject), data = MS1,REML = FALSE)
summary(Gevmod)

full_GEV = lmBF(GEV ~  RS + Subject, data=MS1, whichRandom = 'Subject')
null_GEV = lmBF(GEV ~ 1 + Subject, data = MS1, whichRandom = 'Subject')
BF_GEV = full_GEV/null_GEV

hist((MS1$TP2))
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

# MS5
Durmod <- lmer(Duration ~ RS + (1|Subject), data = MS5,REML = FALSE)
summary(Durmod)
t_to_eta2(-2.67,20)

Occurmod <- lmer(Occurance ~ RS + (1|Subject), data = MS5,REML = FALSE)
summary(Occurmod)
t_to_eta2(-3.095,20)

Covermod <- lmer(Coverage ~ RS + (1|Subject), data = MS5,REML = FALSE)
summary(Covermod)
t_to_eta2(-3.301,20)


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

RS1vsRS2 = ALL_Pvals

# Run the FDR correction for all RS1 vs RS2 Pvals
ControlPvalues <- unlist(RS1vsRS2)
non_na_indices <- which(!is.na(ControlPvalues))
All_adj_p_values <- p.adjust(ControlPvalues[non_na_indices], method = "fdr")
ControlPvalues[non_na_indices] <- All_adj_p_values
ControlPvalues[ControlPvalues<.05]


# Control_CorrectedPs <- data.frame(
#   fdr_adj_duration = p.adjust(ALL_Pvals$Duration, method = "fdr"),
#   fdr_adj_coverage = p.adjust(ALL_Pvals$Coverage, method = "fdr"),
#   fdr_adj_occurrence = p.adjust(ALL_Pvals$Occurrence, method = "fdr"),
#   fdr_adj_GEV = p.adjust(ALL_Pvals$GEV, method = "fdr"),
#   fdr_adj_TP1 = p.adjust(ALL_Pvals$TP1, method = "fdr"),
#   fdr_adj_TP2 = p.adjust(ALL_Pvals$TP2, method = "fdr"),
#   fdr_adj_TP3 = p.adjust(ALL_Pvals$TP3, method = "fdr"),
#   fdr_adj_TP4 = p.adjust(ALL_Pvals$TP4, method = "fdr"),
#   fdr_adj_TP5 = p.adjust(ALL_Pvals$TP5, method = "fdr")
# )


## Microstates x Area ----
MicrostatsxArea <- read_excel("Path_to_Microstate_Area_data/MicrostatsxArea.xlsx")
View(MicrostatsxArea)
MicrostatsxArea <- as.data.frame(MicrostatsxArea)
MicrostatsxArea = MicrostatsxArea[MicrostatsxArea$Condition == 'C',]

MicrostatsxArea$GFP <- as.numeric(scale(MicrostatsxArea$GFP))
MicrostatsxArea$Duration <- as.numeric(scale(MicrostatsxArea$Duration))
MicrostatsxArea$Occurance <- as.numeric(scale(MicrostatsxArea$Occurance))
MicrostatsxArea$Coverage <- as.numeric((MicrostatsxArea$Coverage))

MS1RS1 = MicrostatsxArea[MicrostatsxArea$MS == 'MS1' & MicrostatsxArea$RS == 'RS1',]
MS2RS1 = MicrostatsxArea[MicrostatsxArea$MS == 'MS2' & MicrostatsxArea$RS == 'RS1',]
MS3RS1 = MicrostatsxArea[MicrostatsxArea$MS == 'MS3' & MicrostatsxArea$RS == 'RS1',]
MS4RS1 = MicrostatsxArea[MicrostatsxArea$MS == 'MS4' & MicrostatsxArea$RS == 'RS1',]
MS5RS1 = MicrostatsxArea[MicrostatsxArea$MS == 'MS5' & MicrostatsxArea$RS == 'RS1',]

# MS1RS2 = MicrostatsxArea[MicrostatsxArea$MS == 'MS1' & MicrostatsxArea$RS == 'RS2',]
# MS2RS2 = MicrostatsxArea[MicrostatsxArea$MS == 'MS2' & MicrostatsxArea$RS == 'RS2',]
# MS3RS2 = MicrostatsxArea[MicrostatsxArea$MS == 'MS3' & MicrostatsxArea$RS == 'RS2',]
# MS4RS2 = MicrostatsxArea[MicrostatsxArea$MS == 'MS4' & MicrostatsxArea$RS == 'RS2',]
# MS5RS2 = MicrostatsxArea[MicrostatsxArea$MS == 'MS5' & MicrostatsxArea$RS == 'RS2',]
# 
# MS1diff = MS1RS1
# MS1diff[,c(6:16)] = MS1RS1[,c(6:16)] - MS1RS2[,c(6:16)]
# MS2diff = MS2RS1
# MS2diff[,c(6:16)] = MS2RS1[,c(6:16)] - MS2RS2[,c(6:16)]
# MS3diff = MS3RS1
# MS3diff[,c(6:16)] = MS3RS1[,c(6:16)] - MS3RS2[,c(6:16)]
# MS4diff = MS4RS1
# MS4diff[,c(6:16)] = MS4RS1[,c(6:16)] - MS4RS2[,c(6:16)]
# MS5diff = MS5RS1
# MS5diff[,c(6:16)] = MS5RS1[,c(6:16)] - MS5RS2[,c(6:16)]

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
# Initialize variables

# Create column and row names
col_names <- c("RS1","RS2","Diff")
row_names <- c("MS1", "MS2", "MS3", "MS4", "MS5")
# Create empty data frames
Duration_Pvals <- data.frame(matrix(ncol = length(col_names), nrow = length(row_names)))
colnames(Duration_Pvals) <- col_names
rownames(Duration_Pvals) <- row_names
Occurance_Pvals = Duration_Pvals
Coverage_Pvals = Duration_Pvals
GEV_Pvals = Duration_Pvals
GFP_Pvals = Duration_Pvals
TP1_Pvals = Duration_Pvals
TP2_Pvals = Duration_Pvals
TP3_Pvals = Duration_Pvals
TP4_Pvals = Duration_Pvals
TP5_Pvals = Duration_Pvals

## MS1

# RS1
Fullmod <- lm(Length ~ Duration + Occurance + Coverage + GFP + GEV + TP2 + TP3 + TP4 + TP5, data = MS1RS1)
summary(Fullmod)
anova(Fullmod)
Duration_Pvals["MS1","RS1"]  = summary(Fullmod)$coefficients['Duration','Pr(>|t|)']
Occurance_Pvals["MS1","RS1"] = summary(Fullmod)$coefficients['Occurance','Pr(>|t|)']
Coverage_Pvals["MS1","RS1"]  = summary(Fullmod)$coefficients['Coverage','Pr(>|t|)']
GEV_Pvals["MS1","RS1"]       = summary(Fullmod)$coefficients['GEV','Pr(>|t|)']
GFP_Pvals["MS1","RS1"]       = summary(Fullmod)$coefficients['GFP','Pr(>|t|)']
TP2_Pvals["MS1","RS1"]       = summary(Fullmod)$coefficients['TP2','Pr(>|t|)']
TP3_Pvals["MS1","RS1"]       = summary(Fullmod)$coefficients['TP3','Pr(>|t|)']
TP4_Pvals["MS1","RS1"]       = summary(Fullmod)$coefficients['TP4','Pr(>|t|)']
TP5_Pvals["MS1","RS1"]       = summary(Fullmod)$coefficients['TP5','Pr(>|t|)']

# # RS2
# Fullmod <- lm(Length ~ Condition + Duration + Occurance + Coverage + GEV + TP2 + TP3 + TP4 + TP5, data = MS1RS2)
# anova(Fullmod)
# Duration_Pvals["MS1","RS2"]  = summary(Fullmod)$coefficients['Duration','Pr(>|t|)']
# Occurance_Pvals["MS1","RS2"] = summary(Fullmod)$coefficients['Occurance','Pr(>|t|)']
# Coverage_Pvals["MS1","RS2"]  = summary(Fullmod)$coefficients['Coverage','Pr(>|t|)']
# GEV_Pvals["MS1","RS2"]       = summary(Fullmod)$coefficients['GEV','Pr(>|t|)']
# TP1_Pvals["MS1","RS2"]       = summary(Fullmod)$coefficients['TP1','Pr(>|t|)']
# TP2_Pvals["MS1","RS2"]       = summary(Fullmod)$coefficients['TP2','Pr(>|t|)']
# TP3_Pvals["MS1","RS2"]       = summary(Fullmod)$coefficients['TP3','Pr(>|t|)']
# TP4_Pvals["MS1","RS2"]       = summary(Fullmod)$coefficients['TP4','Pr(>|t|)']
# TP5_Pvals["MS1","RS2"]       = summary(Fullmod)$coefficients['TP5','Pr(>|t|)']
# 
# # diff
# Diffmod <- lm(Length ~ Condition + Duration + Occurance + Coverage + GEV + TP2 + TP3 + TP4 + TP5, data = MS1diff)
# anova(Diffmod)
# Duration_Pvals["MS1","Diff"]  = anova(Diffmod)['Duration','Pr(>|t|)']
# Occurance_Pvals["MS1","Diff"] = anova(Diffmod)['Occurance','Pr(>|t|)']
# Coverage_Pvals["MS1","Diff"]  = anova(Diffmod)['Coverage','Pr(>|t|)']
# GEV_Pvals["MS1","Diff"]       = anova(Diffmod)['GEV','Pr(>|t|)']
# TP1_Pvals["MS1","Diff"]       = anova(Diffmod)['TP1','Pr(>|t|)']
# TP2_Pvals["MS1","Diff"]       = anova(Diffmod)['TP2','Pr(>|t|)']
# TP3_Pvals["MS1","Diff"]       = anova(Diffmod)['TP3','Pr(>|t|)']
# TP4_Pvals["MS1","Diff"]       = anova(Diffmod)['TP4','Pr(>|t|)']
# TP5_Pvals["MS1","Diff"]       = anova(Diffmod)['TP5','Pr(>|t|)']


## MS2
# MS2RS1
Fullmod <- lm(Length ~ Duration + Occurance + Coverage + GFP + GEV + TP1 + TP3 + TP4 + TP5, data = MS2RS1)
summary(Fullmod)
anova(Fullmod)
Duration_Pvals["MS2","RS1"]  = summary(Fullmod)$coefficients['Duration','Pr(>|t|)']
Occurance_Pvals["MS2","RS1"] = summary(Fullmod)$coefficients['Occurance','Pr(>|t|)']
Coverage_Pvals["MS2","RS1"]  = summary(Fullmod)$coefficients['Coverage','Pr(>|t|)']
GEV_Pvals["MS2","RS1"]       = summary(Fullmod)$coefficients['GEV','Pr(>|t|)']
GFP_Pvals["MS2","RS1"]       = summary(Fullmod)$coefficients['GFP','Pr(>|t|)']
TP1_Pvals["MS2","RS1"]       = summary(Fullmod)$coefficients['TP1','Pr(>|t|)']
TP3_Pvals["MS2","RS1"]       = summary(Fullmod)$coefficients['TP3','Pr(>|t|)']
TP4_Pvals["MS2","RS1"]       = summary(Fullmod)$coefficients['TP4','Pr(>|t|)']
TP5_Pvals["MS2","RS1"]       = summary(Fullmod)$coefficients['TP5','Pr(>|t|)']

# # MS2RS2
# Fullmod <- lm(Length ~ Condition + Duration + Occurance + Coverage + GEV + TP1 + TP3 + TP4 + TP5, data = MS2RS2)
# anova(Fullmod)
# Duration_Pvals["MS2","RS2"]  = summary(Fullmod)$coefficients['Duration','Pr(>|t|)']
# Occurance_Pvals["MS2","RS2"] = summary(Fullmod)$coefficients['Occurance','Pr(>|t|)']
# Coverage_Pvals["MS2","RS2"]  = summary(Fullmod)$coefficients['Coverage','Pr(>|t|)']
# GEV_Pvals["MS2","RS2"]       = summary(Fullmod)$coefficients['GEV','Pr(>|t|)']
# TP1_Pvals["MS2","RS2"]       = summary(Fullmod)$coefficients['TP1','Pr(>|t|)']
# TP2_Pvals["MS2","RS2"]       = summary(Fullmod)$coefficients['TP2','Pr(>|t|)']
# TP3_Pvals["MS2","RS2"]       = summary(Fullmod)$coefficients['TP3','Pr(>|t|)']
# TP4_Pvals["MS2","RS2"]       = summary(Fullmod)$coefficients['TP4','Pr(>|t|)']
# TP5_Pvals["MS2","RS2"]       = summary(Fullmod)$coefficients['TP5','Pr(>|t|)']
# 
# # diff
# Diffmod <- lm(Length ~ Condition + Duration + Occurance + Coverage + GEV + TP1 + TP3 + TP4 + TP5, data = MS2diff)
# anova(Diffmod)
# Duration_Pvals["MS2","Diff"]  = anova(Diffmod)['Duration','Pr(>|t|)']
# Occurance_Pvals["MS2","Diff"] = anova(Diffmod)['Occurance','Pr(>|t|)']
# Coverage_Pvals["MS2","Diff"]  = anova(Diffmod)['Coverage','Pr(>|t|)']
# GEV_Pvals["MS2","Diff"]       = anova(Diffmod)['GEV','Pr(>|t|)']
# TP1_Pvals["MS2","Diff"]       = anova(Diffmod)['TP1','Pr(>|t|)']
# TP2_Pvals["MS2","Diff"]       = anova(Diffmod)['TP2','Pr(>|t|)']
# TP3_Pvals["MS2","Diff"]       = anova(Diffmod)['TP3','Pr(>|t|)']
# TP4_Pvals["MS2","Diff"]       = anova(Diffmod)['TP4','Pr(>|t|)']
# TP5_Pvals["MS2","Diff"]       = anova(Diffmod)['TP5','Pr(>|t|)']


## MS3
# MS3RS1
Fullmod <- lm(Length ~ Duration + Occurance + Coverage + GFP + GEV + TP1 + TP2 + TP4 + TP5, data = MS3RS1)
summary(Fullmod)
anova(Fullmod)
Duration_Pvals["MS3","RS1"]  = summary(Fullmod)$coefficients['Duration','Pr(>|t|)']
Occurance_Pvals["MS3","RS1"] = summary(Fullmod)$coefficients['Occurance','Pr(>|t|)']
Coverage_Pvals["MS3","RS1"]  = summary(Fullmod)$coefficients['Coverage','Pr(>|t|)']
GEV_Pvals["MS3","RS1"]       = summary(Fullmod)$coefficients['GEV','Pr(>|t|)']
GFP_Pvals["MS3","RS1"]       = summary(Fullmod)$coefficients['GFP','Pr(>|t|)']
TP1_Pvals["MS3","RS1"]       = summary(Fullmod)$coefficients['TP1','Pr(>|t|)']
TP2_Pvals["MS3","RS1"]       = summary(Fullmod)$coefficients['TP2','Pr(>|t|)']
TP4_Pvals["MS3","RS1"]       = summary(Fullmod)$coefficients['TP4','Pr(>|t|)']
TP5_Pvals["MS3","RS1"]       = summary(Fullmod)$coefficients['TP5','Pr(>|t|)']

# # MS3RS2
# Fullmod <- lm(Length ~ Condition + Duration + Occurance + Coverage + GEV + TP1 + TP2 + TP4 + TP5, data = MS3RS2)
# anova(Fullmod)
# Duration_Pvals["MS3","RS2"]  = summary(Fullmod)$coefficients['Duration','Pr(>|t|)']
# Occurance_Pvals["MS3","RS2"] = summary(Fullmod)$coefficients['Occurance','Pr(>|t|)']
# Coverage_Pvals["MS3","RS2"]  = summary(Fullmod)$coefficients['Coverage','Pr(>|t|)']
# GEV_Pvals["MS3","RS2"]       = summary(Fullmod)$coefficients['GEV','Pr(>|t|)']
# TP1_Pvals["MS3","RS2"]       = summary(Fullmod)$coefficients['TP1','Pr(>|t|)']
# TP2_Pvals["MS3","RS2"]       = summary(Fullmod)$coefficients['TP2','Pr(>|t|)']
# TP3_Pvals["MS3","RS2"]       = summary(Fullmod)$coefficients['TP3','Pr(>|t|)']
# TP4_Pvals["MS3","RS2"]       = summary(Fullmod)$coefficients['TP4','Pr(>|t|)']
# TP5_Pvals["MS3","RS2"]       = summary(Fullmod)$coefficients['TP5','Pr(>|t|)']
# 
# # diff
# Diffmod <- lm(Length ~ Condition + Duration + Occurance + Coverage + GEV + TP1 + TP2 + TP4 + TP5, data = MS3diff)
# anova(Diffmod)
# Duration_Pvals["MS3","Diff"]  = anova(Diffmod)['Duration','Pr(>|t|)']
# Occurance_Pvals["MS3","Diff"] = anova(Diffmod)['Occurance','Pr(>|t|)']
# Coverage_Pvals["MS3","Diff"]  = anova(Diffmod)['Coverage','Pr(>|t|)']
# GEV_Pvals["MS3","Diff"]       = anova(Diffmod)['GEV','Pr(>|t|)']
# TP1_Pvals["MS3","Diff"]       = anova(Diffmod)['TP1','Pr(>|t|)']
# TP2_Pvals["MS3","Diff"]       = anova(Diffmod)['TP2','Pr(>|t|)']
# TP3_Pvals["MS3","Diff"]       = anova(Diffmod)['TP3','Pr(>|t|)']
# TP4_Pvals["MS3","Diff"]       = anova(Diffmod)['TP4','Pr(>|t|)']
# TP5_Pvals["MS3","Diff"]       = anova(Diffmod)['TP5','Pr(>|t|)']


## MS4
# MS4RS1
Fullmod <- lm(Length ~ Duration + Occurance + Coverage + GFP + GEV + TP1 + TP2 + TP3 + TP5, data = MS4RS1)
summary(Fullmod)
anova(Fullmod)
Duration_Pvals["MS4","RS1"]  = summary(Fullmod)$coefficients['Duration','Pr(>|t|)']
Occurance_Pvals["MS4","RS1"] = summary(Fullmod)$coefficients['Occurance','Pr(>|t|)']
Coverage_Pvals["MS4","RS1"]  = summary(Fullmod)$coefficients['Coverage','Pr(>|t|)']
GEV_Pvals["MS4","RS1"]       = summary(Fullmod)$coefficients['GEV','Pr(>|t|)']
GFP_Pvals["MS4","RS1"]       = summary(Fullmod)$coefficients['GFP','Pr(>|t|)']
TP1_Pvals["MS4","RS1"]       = summary(Fullmod)$coefficients['TP1','Pr(>|t|)']
TP2_Pvals["MS4","RS1"]       = summary(Fullmod)$coefficients['TP2','Pr(>|t|)']
TP3_Pvals["MS4","RS1"]       = summary(Fullmod)$coefficients['TP3','Pr(>|t|)']
TP5_Pvals["MS4","RS1"]       = summary(Fullmod)$coefficients['TP5','Pr(>|t|)']

# # MS4RS2
# Fullmod <- lm(Length ~ Condition + Duration + Occurance + Coverage + GEV + TP1 + TP2 + TP3 + TP5, data = MS4RS2)
# anova(Fullmod)
# Duration_Pvals["MS4","RS2"]  = summary(Fullmod)$coefficients['Duration','Pr(>|t|)']
# Occurance_Pvals["MS4","RS2"] = summary(Fullmod)$coefficients['Occurance','Pr(>|t|)']
# Coverage_Pvals["MS4","RS2"]  = summary(Fullmod)$coefficients['Coverage','Pr(>|t|)']
# GEV_Pvals["MS4","RS2"]       = summary(Fullmod)$coefficients['GEV','Pr(>|t|)']
# TP1_Pvals["MS4","RS2"]       = summary(Fullmod)$coefficients['TP1','Pr(>|t|)']
# TP2_Pvals["MS4","RS2"]       = summary(Fullmod)$coefficients['TP2','Pr(>|t|)']
# TP3_Pvals["MS4","RS2"]       = summary(Fullmod)$coefficients['TP3','Pr(>|t|)']
# TP4_Pvals["MS4","RS2"]       = summary(Fullmod)$coefficients['TP4','Pr(>|t|)']
# TP5_Pvals["MS4","RS2"]       = summary(Fullmod)$coefficients['TP5','Pr(>|t|)']
# 
# # diff
# Diffmod <- lm(Length ~ Condition + Duration + Occurance + Coverage + GEV + TP2 + TP3 + TP1 + TP5, data = MS4diff)
# anova(Diffmod)
# Duration_Pvals["MS4","Diff"]  = anova(Diffmod)['Duration','Pr(>|t|)']
# Occurance_Pvals["MS4","Diff"] = anova(Diffmod)['Occurance','Pr(>|t|)']
# Coverage_Pvals["MS4","Diff"]  = anova(Diffmod)['Coverage','Pr(>|t|)']
# GEV_Pvals["MS4","Diff"]       = anova(Diffmod)['GEV','Pr(>|t|)']
# TP1_Pvals["MS4","Diff"]       = anova(Diffmod)['TP1','Pr(>|t|)']
# TP2_Pvals["MS4","Diff"]       = anova(Diffmod)['TP2','Pr(>|t|)']
# TP3_Pvals["MS4","Diff"]       = anova(Diffmod)['TP3','Pr(>|t|)']
# TP4_Pvals["MS4","Diff"]       = anova(Diffmod)['TP4','Pr(>|t|)']
# TP5_Pvals["MS4","Diff"]       = anova(Diffmod)['TP5','Pr(>|t|)']


## MS5
# MS5RS1 - Check if it is more reliable to fit separate models for each predictor
Fullmod <- lm(Length ~ Duration + Occurance + Coverage + GFP + GEV + TP1 + TP2 + TP3 + TP4, data = MS5RS1)
summary(Fullmod)
Duration_Pvals["MS5","RS1"]  = summary(Fullmod)$coefficients['Duration','Pr(>|t|)']
Occurance_Pvals["MS5","RS1"] = summary(Fullmod)$coefficients['Occurance','Pr(>|t|)']
Coverage_Pvals["MS5","RS1"]  = summary(Fullmod)$coefficients['Coverage','Pr(>|t|)']
GEV_Pvals["MS5","RS1"]       = summary(Fullmod)$coefficients['GEV','Pr(>|t|)']
GFP_Pvals["MS5","RS1"]       = summary(Fullmod)$coefficients['GFP','Pr(>|t|)']
TP1_Pvals["MS5","RS1"]       = summary(Fullmod)$coefficients['TP1','Pr(>|t|)']
TP2_Pvals["MS5","RS1"]       = summary(Fullmod)$coefficients['TP2','Pr(>|t|)']
TP3_Pvals["MS5","RS1"]       = summary(Fullmod)$coefficients['TP3','Pr(>|t|)']
TP4_Pvals["MS5","RS1"]       = summary(Fullmod)$coefficients['TP4','Pr(>|t|)']

Fullmod <- lm(Length ~ Coverage, data = MS5RS1)
summary(Fullmod)
Intplot <- ggpredict(Fullmod,terms = c('Coverage'))
plot(Intplot,add.data = TRUE)


# Now you have a list of models, each named after the predictor variable
# You can access the models like this: model_results$PredictorName


# # MS5RS2
# Fullmod <- lm(Length ~ Condition + Duration + Occurance + Coverage + GEV + TP1 + TP2 + TP3 + TP4, data = MS5RS2)
# anova(Fullmod)
# Duration_Pvals["MS5","RS2"]  = summary(Fullmod)$coefficients['Duration','Pr(>|t|)']
# Occurance_Pvals["MS5","RS2"] = summary(Fullmod)$coefficients['Occurance','Pr(>|t|)']
# Coverage_Pvals["MS5","RS2"]  = summary(Fullmod)$coefficients['Coverage','Pr(>|t|)']
# GEV_Pvals["MS5","RS2"]       = summary(Fullmod)$coefficients['GEV','Pr(>|t|)']
# TP1_Pvals["MS5","RS2"]       = summary(Fullmod)$coefficients['TP1','Pr(>|t|)']
# TP2_Pvals["MS5","RS2"]       = summary(Fullmod)$coefficients['TP2','Pr(>|t|)']
# TP3_Pvals["MS5","RS2"]       = summary(Fullmod)$coefficients['TP3','Pr(>|t|)']
# TP4_Pvals["MS5","RS2"]       = summary(Fullmod)$coefficients['TP4','Pr(>|t|)']
# TP5_Pvals["MS5","RS2"]       = summary(Fullmod)$coefficients['TP5','Pr(>|t|)']
# 
# # diff
# Diffmod <- lm(Length ~ Condition + Duration + Occurance + Coverage + GEV + TP2 + TP3 + TP4 + TP1, data = MS5diff)
# anova(Diffmod)
# Duration_Pvals["MS5","Diff"]  = anova(Diffmod)['Duration','Pr(>|t|)']
# Occurance_Pvals["MS5","Diff"] = anova(Diffmod)['Occurance','Pr(>|t|)']
# Coverage_Pvals["MS5","Diff"]  = anova(Diffmod)['Coverage','Pr(>|t|)']
# GEV_Pvals["MS5","Diff"]       = anova(Diffmod)['GEV','Pr(>|t|)']
# TP1_Pvals["MS5","Diff"]       = anova(Diffmod)['TP1','Pr(>|t|)']
# TP2_Pvals["MS5","Diff"]       = anova(Diffmod)['TP2','Pr(>|t|)']
# TP3_Pvals["MS5","Diff"]       = anova(Diffmod)['TP3','Pr(>|t|)']
# TP4_Pvals["MS5","Diff"]       = anova(Diffmod)['TP4','Pr(>|t|)']
# TP5_Pvals["MS5","Diff"]       = anova(Diffmod)['TP5','Pr(>|t|)']

# Check for significance after FDR correction
Pvalues <- unlist(TP5_Pvals)
non_na_indices <- which(!is.na(Pvalues))
All_adj_p_values <- p.adjust(Pvalues[non_na_indices], method = "fdr")
Pvalues[non_na_indices] <- All_adj_p_values
Pvalues[Pvalues<.05]


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
non_na_indices <- which(!is.na(Pvalues))
All_adj_p_values <- p.adjust(Pvalues[non_na_indices], method = "fdr")
Pvalues[non_na_indices] <- All_adj_p_values
Pvalues[Pvalues<.05]


# Compute correlation among predictors
subset_cols <- c("GFP", "Duration", "Occurance", "Coverage", "GEV", "TP1", "TP2", "TP3", "TP4","TP5")
# Compute the correlation matrix for the subset of columns
cor_matrix <- cor(MS5RS1[, predictors])
# Print the correlation matrix
print(cor_matrix)



### Width ----
# Initialize variables

# Create column and row names
col_names <- c("RS1","RS2","Diff")
row_names <- c("MS1", "MS2", "MS3", "MS4", "MS5")
# Create empty data frames
Duration_Pvals <- data.frame(matrix(ncol = length(col_names), nrow = length(row_names)))
colnames(Duration_Pvals) <- col_names
rownames(Duration_Pvals) <- row_names
Occurance_Pvals = Duration_Pvals
Coverage_Pvals = Duration_Pvals
GEV_Pvals = Duration_Pvals
GFP_Pvals = Duration_Pvals
TP1_Pvals = Duration_Pvals
TP2_Pvals = Duration_Pvals
TP3_Pvals = Duration_Pvals
TP4_Pvals = Duration_Pvals
TP5_Pvals = Duration_Pvals

## MS1
# RS1
Fullmod <- lm(Width ~ Duration + Occurance + Coverage + GEV + GFP + TP2 + TP3 + TP4 + TP5, data = MS1RS1)
summary(Fullmod)
Duration_Pvals["MS1","RS1"]  = summary(Fullmod)$coefficients['Duration','Pr(>|t|)']
Occurance_Pvals["MS1","RS1"] = summary(Fullmod)$coefficients['Occurance','Pr(>|t|)']
Coverage_Pvals["MS1","RS1"]  = summary(Fullmod)$coefficients['Coverage','Pr(>|t|)']
GEV_Pvals["MS1","RS1"]       = summary(Fullmod)$coefficients['GEV','Pr(>|t|)']
GFP_Pvals["MS1","RS1"]       = summary(Fullmod)$coefficients['GFP','Pr(>|t|)']
TP2_Pvals["MS1","RS1"]       = summary(Fullmod)$coefficients['TP2','Pr(>|t|)']
TP3_Pvals["MS1","RS1"]       = summary(Fullmod)$coefficients['TP3','Pr(>|t|)']
TP4_Pvals["MS1","RS1"]       = summary(Fullmod)$coefficients['TP4','Pr(>|t|)']
TP5_Pvals["MS1","RS1"]       = summary(Fullmod)$coefficients['TP5','Pr(>|t|)']

# # RS2
# Fullmod <- lm(Width ~ Condition + Duration + Occurance + Coverage + GEV + TP2 + TP3 + TP4 + TP5, data = MS1RS2)
# anova(Fullmod)
# Duration_Pvals["MS1","RS2"]  = summary(Fullmod)$coefficients['Duration','Pr(>|t|)']
# Occurance_Pvals["MS1","RS2"] = summary(Fullmod)$coefficients['Occurance','Pr(>|t|)']
# Coverage_Pvals["MS1","RS2"]  = summary(Fullmod)$coefficients['Coverage','Pr(>|t|)']
# GEV_Pvals["MS1","RS2"]       = summary(Fullmod)$coefficients['GEV','Pr(>|t|)']
# TP1_Pvals["MS1","RS2"]       = summary(Fullmod)$coefficients['TP1','Pr(>|t|)']
# TP2_Pvals["MS1","RS2"]       = summary(Fullmod)$coefficients['TP2','Pr(>|t|)']
# TP3_Pvals["MS1","RS2"]       = summary(Fullmod)$coefficients['TP3','Pr(>|t|)']
# TP4_Pvals["MS1","RS2"]       = summary(Fullmod)$coefficients['TP4','Pr(>|t|)']
# TP5_Pvals["MS1","RS2"]       = summary(Fullmod)$coefficients['TP5','Pr(>|t|)']
# 
# # diff
# Diffmod <- lm(Width ~ Condition + Duration + Occurance + Coverage + GEV + TP2 + TP3 + TP4 + TP5, data = MS1diff)
# anova(Diffmod)
# Duration_Pvals["MS1","Diff"]  = anova(Diffmod)['Duration','Pr(>|t|)']
# Occurance_Pvals["MS1","Diff"] = anova(Diffmod)['Occurance','Pr(>|t|)']
# Coverage_Pvals["MS1","Diff"]  = anova(Diffmod)['Coverage','Pr(>|t|)']
# GEV_Pvals["MS1","Diff"]       = anova(Diffmod)['GEV','Pr(>|t|)']
# TP1_Pvals["MS1","Diff"]       = anova(Diffmod)['TP1','Pr(>|t|)']
# TP2_Pvals["MS1","Diff"]       = anova(Diffmod)['TP2','Pr(>|t|)']
# TP3_Pvals["MS1","Diff"]       = anova(Diffmod)['TP3','Pr(>|t|)']
# TP4_Pvals["MS1","Diff"]       = anova(Diffmod)['TP4','Pr(>|t|)']
# TP5_Pvals["MS1","Diff"]       = anova(Diffmod)['TP5','Pr(>|t|)']


## MS2
# MS2RS1
Fullmod <- lm(Width ~ Duration + Occurance + Coverage + GEV + GFP + TP1 + TP3 + TP4 + TP5, data = MS2RS1)
summary(Fullmod)
Duration_Pvals["MS2","RS1"]  = summary(Fullmod)$coefficients['Duration','Pr(>|t|)']
Occurance_Pvals["MS2","RS1"] = summary(Fullmod)$coefficients['Occurance','Pr(>|t|)']
Coverage_Pvals["MS2","RS1"]  = summary(Fullmod)$coefficients['Coverage','Pr(>|t|)']
GEV_Pvals["MS2","RS1"]       = summary(Fullmod)$coefficients['GEV','Pr(>|t|)']
GFP_Pvals["MS2","RS1"]       = summary(Fullmod)$coefficients['GFP','Pr(>|t|)']
TP1_Pvals["MS2","RS1"]       = summary(Fullmod)$coefficients['TP1','Pr(>|t|)']
TP3_Pvals["MS2","RS1"]       = summary(Fullmod)$coefficients['TP3','Pr(>|t|)']
TP4_Pvals["MS2","RS1"]       = summary(Fullmod)$coefficients['TP4','Pr(>|t|)']
TP5_Pvals["MS2","RS1"]       = summary(Fullmod)$coefficients['TP5','Pr(>|t|)']

# # MS2RS2
# Fullmod <- lm(Width ~ Condition + Duration + Occurance + Coverage + GEV + TP1 + TP3 + TP4 + TP5, data = MS2RS2)
# anova(Fullmod)
# Duration_Pvals["MS2","RS2"]  = summary(Fullmod)$coefficients['Duration','Pr(>|t|)']
# Occurance_Pvals["MS2","RS2"] = summary(Fullmod)$coefficients['Occurance','Pr(>|t|)']
# Coverage_Pvals["MS2","RS2"]  = summary(Fullmod)$coefficients['Coverage','Pr(>|t|)']
# GEV_Pvals["MS2","RS2"]       = summary(Fullmod)$coefficients['GEV','Pr(>|t|)']
# TP1_Pvals["MS2","RS2"]       = summary(Fullmod)$coefficients['TP1','Pr(>|t|)']
# TP2_Pvals["MS2","RS2"]       = summary(Fullmod)$coefficients['TP2','Pr(>|t|)']
# TP3_Pvals["MS2","RS2"]       = summary(Fullmod)$coefficients['TP3','Pr(>|t|)']
# TP4_Pvals["MS2","RS2"]       = summary(Fullmod)$coefficients['TP4','Pr(>|t|)']
# TP5_Pvals["MS2","RS2"]       = summary(Fullmod)$coefficients['TP5','Pr(>|t|)']
# 
# # diff
# Diffmod <- lm(Width ~ Condition + Duration + Occurance + Coverage + GEV + TP1 + TP3 + TP4 + TP5, data = MS2diff)
# anova(Diffmod)
# Duration_Pvals["MS2","Diff"]  = anova(Diffmod)['Duration','Pr(>|t|)']
# Occurance_Pvals["MS2","Diff"] = anova(Diffmod)['Occurance','Pr(>|t|)']
# Coverage_Pvals["MS2","Diff"]  = anova(Diffmod)['Coverage','Pr(>|t|)']
# GEV_Pvals["MS2","Diff"]       = anova(Diffmod)['GEV','Pr(>|t|)']
# TP1_Pvals["MS2","Diff"]       = anova(Diffmod)['TP1','Pr(>|t|)']
# TP2_Pvals["MS2","Diff"]       = anova(Diffmod)['TP2','Pr(>|t|)']
# TP3_Pvals["MS2","Diff"]       = anova(Diffmod)['TP3','Pr(>|t|)']
# TP4_Pvals["MS2","Diff"]       = anova(Diffmod)['TP4','Pr(>|t|)']
# TP5_Pvals["MS2","Diff"]       = anova(Diffmod)['TP5','Pr(>|t|)']


## MS3
# MS3RS1
Fullmod <- lm(Width ~ Duration + Occurance + Coverage + GEV + GFP + TP1 + TP2 + TP4 + TP5, data = MS3RS1)
summary(Fullmod)
Duration_Pvals["MS3","RS1"]  = summary(Fullmod)$coefficients['Duration','Pr(>|t|)']
Occurance_Pvals["MS3","RS1"] = summary(Fullmod)$coefficients['Occurance','Pr(>|t|)']
Coverage_Pvals["MS3","RS1"]  = summary(Fullmod)$coefficients['Coverage','Pr(>|t|)']
GEV_Pvals["MS3","RS1"]       = summary(Fullmod)$coefficients['GEV','Pr(>|t|)']
GFP_Pvals["MS3","RS1"]       = summary(Fullmod)$coefficients['GFP','Pr(>|t|)']
TP1_Pvals["MS3","RS1"]       = summary(Fullmod)$coefficients['TP1','Pr(>|t|)']
TP2_Pvals["MS3","RS1"]       = summary(Fullmod)$coefficients['TP2','Pr(>|t|)']
TP4_Pvals["MS3","RS1"]       = summary(Fullmod)$coefficients['TP4','Pr(>|t|)']
TP5_Pvals["MS3","RS1"]       = summary(Fullmod)$coefficients['TP5','Pr(>|t|)']

Intplot <- ggpredict(Fullmod,terms = c('GEV'))
plot(Intplot,add.data = TRUE)

# # MS3RS2
# Fullmod <- lm(Width ~ Condition + Duration + Occurance + Coverage + GEV + TP1 + TP2 + TP4 + TP5, data = MS3RS2)
# anova(Fullmod)
# Duration_Pvals["MS3","RS2"]  = summary(Fullmod)$coefficients['Duration','Pr(>|t|)']
# Occurance_Pvals["MS3","RS2"] = summary(Fullmod)$coefficients['Occurance','Pr(>|t|)']
# Coverage_Pvals["MS3","RS2"]  = summary(Fullmod)$coefficients['Coverage','Pr(>|t|)']
# GEV_Pvals["MS3","RS2"]       = summary(Fullmod)$coefficients['GEV','Pr(>|t|)']
# TP1_Pvals["MS3","RS2"]       = summary(Fullmod)$coefficients['TP1','Pr(>|t|)']
# TP2_Pvals["MS3","RS2"]       = summary(Fullmod)$coefficients['TP2','Pr(>|t|)']
# TP3_Pvals["MS3","RS2"]       = summary(Fullmod)$coefficients['TP3','Pr(>|t|)']
# TP4_Pvals["MS3","RS2"]       = summary(Fullmod)$coefficients['TP4','Pr(>|t|)']
# TP5_Pvals["MS3","RS2"]       = summary(Fullmod)$coefficients['TP5','Pr(>|t|)']
# 
# # diff
# Diffmod <- lm(Width ~ Condition + Duration + Occurance + Coverage + GEV + TP2 + TP1 + TP4 + TP5, data = MS3diff)
# anova(Diffmod)
# Duration_Pvals["MS3","Diff"]  = anova(Diffmod)['Duration','Pr(>|t|)']
# Occurance_Pvals["MS3","Diff"] = anova(Diffmod)['Occurance','Pr(>|t|)']
# Coverage_Pvals["MS3","Diff"]  = anova(Diffmod)['Coverage','Pr(>|t|)']
# GEV_Pvals["MS3","Diff"]       = anova(Diffmod)['GEV','Pr(>|t|)']
# TP1_Pvals["MS3","Diff"]       = anova(Diffmod)['TP1','Pr(>|t|)']
# TP2_Pvals["MS3","Diff"]       = anova(Diffmod)['TP2','Pr(>|t|)']
# TP3_Pvals["MS3","Diff"]       = anova(Diffmod)['TP3','Pr(>|t|)']
# TP4_Pvals["MS3","Diff"]       = anova(Diffmod)['TP4','Pr(>|t|)']
# TP5_Pvals["MS3","Diff"]       = anova(Diffmod)['TP5','Pr(>|t|)']


## MS4
# MS4RS1
Fullmod <- lm(Width ~ Duration + Occurance + Coverage + GEV + GFP + TP1 + TP2 + TP3 + TP5, data = MS4RS1)
summary(Fullmod)
Duration_Pvals["MS4","RS1"]  = summary(Fullmod)$coefficients['Duration','Pr(>|t|)']
Occurance_Pvals["MS4","RS1"] = summary(Fullmod)$coefficients['Occurance','Pr(>|t|)']
Coverage_Pvals["MS4","RS1"]  = summary(Fullmod)$coefficients['Coverage','Pr(>|t|)']
GEV_Pvals["MS4","RS1"]       = summary(Fullmod)$coefficients['GEV','Pr(>|t|)']
GFP_Pvals["MS4","RS1"]       = summary(Fullmod)$coefficients['GFP','Pr(>|t|)']
TP1_Pvals["MS4","RS1"]       = summary(Fullmod)$coefficients['TP1','Pr(>|t|)']
TP2_Pvals["MS4","RS1"]       = summary(Fullmod)$coefficients['TP2','Pr(>|t|)']
TP3_Pvals["MS4","RS1"]       = summary(Fullmod)$coefficients['TP3','Pr(>|t|)']
TP5_Pvals["MS4","RS1"]       = summary(Fullmod)$coefficients['TP5','Pr(>|t|)']


# # MS4RS2
# Fullmod <- lm(Width ~ Condition + Duration + Occurance + Coverage + GEV + TP1 + TP2 + TP3 + TP5, data = MS4RS2)
# anova(Fullmod)
# Duration_Pvals["MS4","RS2"]  = summary(Fullmod)$coefficients['Duration','Pr(>|t|)']
# Occurance_Pvals["MS4","RS2"] = summary(Fullmod)$coefficients['Occurance','Pr(>|t|)']
# Coverage_Pvals["MS4","RS2"]  = summary(Fullmod)$coefficients['Coverage','Pr(>|t|)']
# GEV_Pvals["MS4","RS2"]       = summary(Fullmod)$coefficients['GEV','Pr(>|t|)']
# TP1_Pvals["MS4","RS2"]       = summary(Fullmod)$coefficients['TP1','Pr(>|t|)']
# TP2_Pvals["MS4","RS2"]       = summary(Fullmod)$coefficients['TP2','Pr(>|t|)']
# TP3_Pvals["MS4","RS2"]       = summary(Fullmod)$coefficients['TP3','Pr(>|t|)']
# TP4_Pvals["MS4","RS2"]       = summary(Fullmod)$coefficients['TP4','Pr(>|t|)']
# TP5_Pvals["MS4","RS2"]       = summary(Fullmod)$coefficients['TP5','Pr(>|t|)']
# 
# # diff
# Diffmod <- lm(Width ~ Condition + Duration + Occurance + Coverage + GEV + TP2 + TP3 + TP1 + TP5, data = MS4diff)
# anova(Diffmod)
# Duration_Pvals["MS4","Diff"]  = anova(Diffmod)['Duration','Pr(>|t|)']
# Occurance_Pvals["MS4","Diff"] = anova(Diffmod)['Occurance','Pr(>|t|)']
# Coverage_Pvals["MS4","Diff"]  = anova(Diffmod)['Coverage','Pr(>|t|)']
# GEV_Pvals["MS4","Diff"]       = anova(Diffmod)['GEV','Pr(>|t|)']
# TP1_Pvals["MS4","Diff"]       = anova(Diffmod)['TP1','Pr(>|t|)']
# TP2_Pvals["MS4","Diff"]       = anova(Diffmod)['TP2','Pr(>|t|)']
# TP3_Pvals["MS4","Diff"]       = anova(Diffmod)['TP3','Pr(>|t|)']
# TP4_Pvals["MS4","Diff"]       = anova(Diffmod)['TP4','Pr(>|t|)']
# TP5_Pvals["MS4","Diff"]       = anova(Diffmod)['TP5','Pr(>|t|)']


## MS5
# MS5RS1
Fullmod <- lm(Width ~ Duration + Occurance + Coverage + GEV + GFP + TP1 + TP2 + TP3 + TP4, data = MS5RS1)
summary(Fullmod)
Duration_Pvals["MS5","RS1"]  = summary(Fullmod)$coefficients['Duration','Pr(>|t|)']
Occurance_Pvals["MS5","RS1"] = summary(Fullmod)$coefficients['Occurance','Pr(>|t|)']
Coverage_Pvals["MS5","RS1"]  = summary(Fullmod)$coefficients['Coverage','Pr(>|t|)']
GEV_Pvals["MS5","RS1"]       = summary(Fullmod)$coefficients['GEV','Pr(>|t|)']
GFP_Pvals["MS5","RS1"]       = summary(Fullmod)$coefficients['GFP','Pr(>|t|)']
TP1_Pvals["MS5","RS1"]       = summary(Fullmod)$coefficients['TP1','Pr(>|t|)']
TP2_Pvals["MS5","RS1"]       = summary(Fullmod)$coefficients['TP2','Pr(>|t|)']
TP3_Pvals["MS5","RS1"]       = summary(Fullmod)$coefficients['TP3','Pr(>|t|)']
TP4_Pvals["MS5","RS1"]       = summary(Fullmod)$coefficients['TP4','Pr(>|t|)']

# # MS5RS2
# Fullmod <- lm(Width ~ Condition + Duration + Occurance + Coverage + GEV + TP1 + TP2 + TP3 + TP4, data = MS5RS2)
# anova(Fullmod)
# Duration_Pvals["MS5","RS2"]  = summary(Fullmod)$coefficients['Duration','Pr(>|t|)']
# Occurance_Pvals["MS5","RS2"] = summary(Fullmod)$coefficients['Occurance','Pr(>|t|)']
# Coverage_Pvals["MS5","RS2"]  = summary(Fullmod)$coefficients['Coverage','Pr(>|t|)']
# GEV_Pvals["MS5","RS2"]       = summary(Fullmod)$coefficients['GEV','Pr(>|t|)']
# TP1_Pvals["MS5","RS2"]       = summary(Fullmod)$coefficients['TP1','Pr(>|t|)']
# TP2_Pvals["MS5","RS2"]       = summary(Fullmod)$coefficients['TP2','Pr(>|t|)']
# TP3_Pvals["MS5","RS2"]       = summary(Fullmod)$coefficients['TP3','Pr(>|t|)']
# TP4_Pvals["MS5","RS2"]       = summary(Fullmod)$coefficients['TP4','Pr(>|t|)']
# TP5_Pvals["MS5","RS2"]       = summary(Fullmod)$coefficients['TP5','Pr(>|t|)']
# 
# # diff
# Diffmod <- lm(Width ~ Condition + Duration + Occurance + Coverage + GEV + TP2 + TP3 + TP4 + TP1, data = MS5diff)
# anova(Diffmod)
# Duration_Pvals["MS5","Diff"]  = anova(Diffmod)['Duration','Pr(>|t|)']
# Occurance_Pvals["MS5","Diff"] = anova(Diffmod)['Occurance','Pr(>|t|)']
# Coverage_Pvals["MS5","Diff"]  = anova(Diffmod)['Coverage','Pr(>|t|)']
# GEV_Pvals["MS5","Diff"]       = anova(Diffmod)['GEV','Pr(>|t|)']
# TP1_Pvals["MS5","Diff"]       = anova(Diffmod)['TP1','Pr(>|t|)']
# TP2_Pvals["MS5","Diff"]       = anova(Diffmod)['TP2','Pr(>|t|)']
# TP3_Pvals["MS5","Diff"]       = anova(Diffmod)['TP3','Pr(>|t|)']
# TP4_Pvals["MS5","Diff"]       = anova(Diffmod)['TP4','Pr(>|t|)']
# TP5_Pvals["MS5","Diff"]       = anova(Diffmod)['TP5','Pr(>|t|)']

# Check for significance after FDR correction
Pvalues <- unlist(TP5_Pvals)
non_na_indices <- which(!is.na(Pvalues))
All_adj_p_values <- p.adjust(Pvalues[non_na_indices], method = "fdr")
Pvalues[non_na_indices] <- All_adj_p_values
Pvalues[Pvalues<.05]

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
non_na_indices <- which(!is.na(Pvalues))
All_adj_p_values <- p.adjust(Pvalues[non_na_indices], method = "fdr")
Pvalues[non_na_indices] <- All_adj_p_values
Pvalues[Pvalues<.05]


## Microstates x PPRatings ----
MicrostatsxPP <- read_excel("Path_to_Microstate_PPRatings_data/MicrostatsxPPratings.xlsx")
# View(MicrostatsxPP)


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


# MS1RS2 = MicrostatsxPP[MicrostatsxPP$MS == 'MS1' & MicrostatsxPP$RS == 'RS2',]
# MS2RS2 = MicrostatsxPP[MicrostatsxPP$MS == 'MS2' & MicrostatsxPP$RS == 'RS2',]
# MS3RS2 = MicrostatsxPP[MicrostatsxPP$MS == 'MS3' & MicrostatsxPP$RS == 'RS2',]
# MS4RS2 = MicrostatsxPP[MicrostatsxPP$MS == 'MS4' & MicrostatsxPP$RS == 'RS2',]
# MS5RS2 = MicrostatsxPP[MicrostatsxPP$MS == 'MS5' & MicrostatsxPP$RS == 'RS2',]
# 
# MS1diff = MS1RS1
# MS1diff[,c(6:16)] = MS1RS1[,c(6:16)] - MS1RS2[,c(6:16)]
# MS2diff = MS2RS1
# MS2diff[,c(6:16)] = MS2RS1[,c(6:16)] - MS2RS2[,c(6:16)]
# MS3diff = MS3RS1
# MS3diff[,c(6:16)] = MS3RS1[,c(6:16)] - MS3RS2[,c(6:16)]
# MS4diff = MS4RS1
# MS4diff[,c(6:16)] = MS4RS1[,c(6:16)] - MS4RS2[,c(6:16)]
# MS5diff = MS5RS1
# MS5diff[,c(6:16)] = MS5RS1[,c(6:16)] - MS5RS2[,c(6:16)]

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


Controlarm = MicrostatsxPP[MicrostatsxPP$Arm == 'Cont',]
hist(Controlarm$IntensityRatings - Controlarm$IntensityT0, main = "Intensity Ratings(T1 - T0): Control Arm")

HFSarm = MicrostatsxPP[MicrostatsxPP$Arm == 'MFS',]
hist(HFSarm$IntensityRatings - HFSarm$IntensityT0, main = "Intensity Ratings(T1 - T0): HFS Arm")

# Ensure the lengths match for both groups
mfs_values <- T1$IntensityRatings[T1$Arm == 'MFS']
control_values <- T1$IntensityRatings[T1$Arm == 'Cont']
min_length <- min(length(mfs_values), length(control_values))
differences <- mfs_values[1:min_length] - control_values[1:min_length]
differences <- na.omit(differences)
hist(differences, main = "Difference in Ratings at T1 (MFS-Control)", xlab = "Difference in Intensity")

### Intensity ----
# Initialize variables

# Create column and row names
col_names <- c("RS1","RS2","Diff")
row_names <- c("MS1", "MS2", "MS3", "MS4", "MS5")
# Create empty data frames
Duration_Pvals <- data.frame(matrix(ncol = length(col_names), nrow = length(row_names)))
colnames(Duration_Pvals) <- col_names
rownames(Duration_Pvals) <- row_names
Occurance_Pvals = Duration_Pvals
Coverage_Pvals = Duration_Pvals
GEV_Pvals = Duration_Pvals
GFP_Pvals = Duration_Pvals
TP1_Pvals = Duration_Pvals
TP2_Pvals = Duration_Pvals
TP3_Pvals = Duration_Pvals
TP4_Pvals = Duration_Pvals
TP5_Pvals = Duration_Pvals

## MS1
# RS1
Fullmod <- lmer(IntensityRatings ~ IntensityT0 + GFP*Arm + Duration*Arm + Occurance*Arm + Coverage*Arm + GEV*Arm + TP2*Arm + TP3*Arm + TP4*Arm + TP5*Arm + (1|Subject), data = MS1RS1,REML = FALSE)
# Fullmod <- lmer(IntensityRatings ~ IntensityT0 + GFP:Arm + Duration:Arm + Occurance:Arm + Coverage:Arm + GEV:Arm + TP2:Arm + TP3:Arm + TP4:Arm + TP5:Arm + (1|Subject), data = MS1RS1,REML = FALSE)
summary(Fullmod)
Duration_Pvals["MS1","RS1"]  = summary(Fullmod)$coefficients['ArmMFS:Duration','Pr(>|t|)']
Occurance_Pvals["MS1","RS1"] = summary(Fullmod)$coefficients['ArmMFS:Occurance','Pr(>|t|)']
Coverage_Pvals["MS1","RS1"]  = summary(Fullmod)$coefficients['ArmMFS:Coverage','Pr(>|t|)']
GEV_Pvals["MS1","RS1"]       = summary(Fullmod)$coefficients['ArmMFS:GEV','Pr(>|t|)']
GFP_Pvals["MS1","RS1"]       = summary(Fullmod)$coefficients['GFP:ArmMFS','Pr(>|t|)']
TP2_Pvals["MS1","RS1"]       = summary(Fullmod)$coefficients['ArmMFS:TP2','Pr(>|t|)']
TP3_Pvals["MS1","RS1"]       = summary(Fullmod)$coefficients['ArmMFS:TP3','Pr(>|t|)']
TP4_Pvals["MS1","RS1"]       = summary(Fullmod)$coefficients['ArmMFS:TP4','Pr(>|t|)']
TP5_Pvals["MS1","RS1"]       = summary(Fullmod)$coefficients['ArmMFS:TP5','Pr(>|t|)']

Intplot <- ggpredict(Fullmod,terms = c('Coverage','Arm'))
plot(Intplot,add.data = TRUE)

## MS2
# RS1
# Fullmod <- lmer(IntensityRatings ~ IntensityT0 + GFP:Arm + Duration:Arm + Occurance:Arm + Coverage:Arm + GEV:Arm + TP1:Arm + TP3:Arm + TP4:Arm + TP5:Arm + (1|Subject), data = MS2RS1,REML = FALSE)
Fullmod <- lmer(IntensityRatings ~ IntensityT0 + GFP*Arm + Duration*Arm + Occurance*Arm + Coverage*Arm + GEV*Arm + TP1*Arm + TP3*Arm + TP4*Arm + TP5*Arm + (1|Subject), data = MS2RS1,REML = FALSE)
summary(Fullmod)
Duration_Pvals["MS2","RS1"]  = summary(Fullmod)$coefficients['ArmMFS:Duration','Pr(>|t|)']
Occurance_Pvals["MS2","RS1"] = summary(Fullmod)$coefficients['ArmMFS:Occurance','Pr(>|t|)']
Coverage_Pvals["MS2","RS1"]  = summary(Fullmod)$coefficients['ArmMFS:Coverage','Pr(>|t|)']
GEV_Pvals["MS2","RS1"]       = summary(Fullmod)$coefficients['ArmMFS:GEV','Pr(>|t|)']
GFP_Pvals["MS2","RS1"]       = summary(Fullmod)$coefficients['GFP:ArmMFS','Pr(>|t|)']
TP1_Pvals["MS2","RS1"]       = summary(Fullmod)$coefficients['ArmMFS:TP1','Pr(>|t|)']
TP3_Pvals["MS2","RS1"]       = summary(Fullmod)$coefficients['ArmMFS:TP3','Pr(>|t|)']
TP4_Pvals["MS2","RS1"]       = summary(Fullmod)$coefficients['ArmMFS:TP4','Pr(>|t|)']
TP5_Pvals["MS2","RS1"]       = summary(Fullmod)$coefficients['ArmMFS:TP5','Pr(>|t|)']

Intplot <- ggpredict(Fullmod,terms = c('Occurance','Arm'))
plot(Intplot,add.data = TRUE)


## MS3
# RS1
# Fullmod <- lmer(IntensityRatings ~ IntensityT0 +  GFP:Arm + Duration:Arm + Occurance:Arm + Coverage:Arm + GEV:Arm + TP2:Arm + TP1:Arm + TP4:Arm + TP5:Arm + (1|Subject), data = MS3RS1,REML = FALSE)
Fullmod <- lmer(IntensityRatings ~ IntensityT0 +  GFP*Arm + Duration*Arm + Occurance*Arm + Coverage*Arm + GEV*Arm + TP2*Arm + TP1*Arm + TP4*Arm + TP5*Arm + (1|Subject), data = MS3RS1,REML = FALSE)
summary(Fullmod)
Duration_Pvals["MS3","RS1"]  = summary(Fullmod)$coefficients['ArmMFS:Duration','Pr(>|t|)']
Occurance_Pvals["MS3","RS1"] = summary(Fullmod)$coefficients['ArmMFS:Occurance','Pr(>|t|)']
Coverage_Pvals["MS3","RS1"]  = summary(Fullmod)$coefficients['ArmMFS:Coverage','Pr(>|t|)']
GEV_Pvals["MS3","RS1"]       = summary(Fullmod)$coefficients['ArmMFS:GEV','Pr(>|t|)']
GFP_Pvals["MS3","RS1"]       = summary(Fullmod)$coefficients['GFP:ArmMFS','Pr(>|t|)']
TP1_Pvals["MS3","RS1"]       = summary(Fullmod)$coefficients['ArmMFS:TP1','Pr(>|t|)']
TP2_Pvals["MS3","RS1"]       = summary(Fullmod)$coefficients['ArmMFS:TP2','Pr(>|t|)']
TP4_Pvals["MS3","RS1"]       = summary(Fullmod)$coefficients['ArmMFS:TP4','Pr(>|t|)']
TP5_Pvals["MS3","RS1"]       = summary(Fullmod)$coefficients['ArmMFS:TP5','Pr(>|t|)']

Fullmod <- lmer(IntensityRatings ~ IntensityT0 + GFP*Arm + (1|Subject), data = MS3RS1,REML = FALSE)
summary(Fullmod)

Intplot <- ggpredict(Fullmod,terms = c('Arm','GFP'))
plot(Intplot,add.data = TRUE)

## MS4
# RS1
# Fullmod <- lmer(IntensityRatings ~ IntensityT0 + GFP:Arm + Duration:Arm + Occurance:Arm + Coverage:Arm + GEV:Arm + TP2:Arm + TP3:Arm + TP1:Arm + TP5:Arm + (1|Subject), data = MS4RS1,REML = FALSE)
Fullmod <- lmer(IntensityRatings ~ IntensityT0 + GFP*Arm + Duration*Arm + Occurance*Arm + Coverage*Arm + GEV*Arm + TP2*Arm + TP3*Arm + TP1*Arm + TP5*Arm + (1|Subject), data = MS4RS1,REML = FALSE)
anova(Fullmod)
Duration_Pvals["MS4","RS1"]  = summary(Fullmod)$coefficients['ArmMFS:Duration','Pr(>|t|)']
Occurance_Pvals["MS4","RS1"] = summary(Fullmod)$coefficients['ArmMFS:Occurance','Pr(>|t|)']
Coverage_Pvals["MS4","RS1"]  = summary(Fullmod)$coefficients['ArmMFS:Coverage','Pr(>|t|)']
GEV_Pvals["MS4","RS1"]       = summary(Fullmod)$coefficients['ArmMFS:GEV','Pr(>|t|)']
GFP_Pvals["MS4","RS1"]       = summary(Fullmod)$coefficients['GFP:ArmMFS','Pr(>|t|)']
TP1_Pvals["MS4","RS1"]       = summary(Fullmod)$coefficients['ArmMFS:TP1','Pr(>|t|)']
TP2_Pvals["MS4","RS1"]       = summary(Fullmod)$coefficients['ArmMFS:TP2','Pr(>|t|)']
TP3_Pvals["MS4","RS1"]       = summary(Fullmod)$coefficients['ArmMFS:TP3','Pr(>|t|)']
TP5_Pvals["MS4","RS1"]       = summary(Fullmod)$coefficients['ArmMFS:TP5','Pr(>|t|)']

Intplot <- ggpredict(Fullmod,terms = c('Duration','Arm'))
plot(Intplot,add.data = TRUE)


## MS5
# RS1
# Fullmod <- lmer(IntensityRatings ~ IntensityT0 + GFP:Arm + Duration:Arm + Occurance:Arm + Coverage:Arm + GEV:Arm + TP2:Arm + TP3:Arm + TP4:Arm + TP1:Arm + (1|Subject), data = MS5RS1,REML = FALSE)
Fullmod <- lmer(IntensityRatings ~ IntensityT0 + GFP*Arm + Duration*Arm + Occurance*Arm + Coverage*Arm + GEV*Arm + TP2*Arm + TP3*Arm + TP4*Arm + TP1*Arm + (1|Subject), data = MS5RS1,REML = FALSE)
summary(Fullmod)
Duration_Pvals["MS5","RS1"]  = summary(Fullmod)$coefficients['ArmMFS:Duration','Pr(>|t|)']
Occurance_Pvals["MS5","RS1"] = summary(Fullmod)$coefficients['ArmMFS:Occurance','Pr(>|t|)']
Coverage_Pvals["MS5","RS1"]  = summary(Fullmod)$coefficients['ArmMFS:Coverage','Pr(>|t|)']
GEV_Pvals["MS5","RS1"]       = summary(Fullmod)$coefficients['ArmMFS:GEV','Pr(>|t|)']
GFP_Pvals["MS5","RS1"]       = summary(Fullmod)$coefficients['GFP:ArmMFS','Pr(>|t|)']
TP1_Pvals["MS5","RS1"]       = summary(Fullmod)$coefficients['ArmMFS:TP1','Pr(>|t|)']
TP2_Pvals["MS5","RS1"]       = summary(Fullmod)$coefficients['ArmMFS:TP2','Pr(>|t|)']
TP3_Pvals["MS5","RS1"]       = summary(Fullmod)$coefficients['ArmMFS:TP3','Pr(>|t|)']
TP4_Pvals["MS5","RS1"]       = summary(Fullmod)$coefficients['ArmMFS:TP4','Pr(>|t|)']

# Check for significance after FDR correction
Pvalues <- unlist(Coverage_Pvals)
non_na_indices <- which(!is.na(Pvalues))
All_adj_p_values <- p.adjust(Pvalues[non_na_indices], method = "fdr")
Pvalues[non_na_indices] <- All_adj_p_values
Pvalues[Pvalues<.05]

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
non_na_indices <- which(!is.na(Pvalues))
All_adj_p_values <- p.adjust(Pvalues[non_na_indices], method = "fdr")
Pvalues[non_na_indices] <- All_adj_p_values
Pvalues[Pvalues<.05]


Fullmod <- lmer(IntensityRatings ~ IntensityT0 + Duration*Arm + (1|Subject), data = MS5RS1,REML = FALSE)
summary(Fullmod)
plot(Fullmod)
t_to_eta2(4.450,379.37021)




# Plotting for Intensity ----
dat = MS5RS1
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
             position=position_jitter(h=0.05, w=.05), alpha = .5) + scale_colour_manual(values = c("black","red")) + 
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
             position=position_jitter(h=0.05, w=.005), alpha = .5) + scale_colour_manual(values = c("black","red")) + 
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
# Initialize variables

# Create column and row names
col_names <- c("RS1","RS2","Diff")
row_names <- c("MS1", "MS2", "MS3", "MS4", "MS5")
# Create empty data frames
Duration_Pvals <- data.frame(matrix(ncol = length(col_names), nrow = length(row_names)))
colnames(Duration_Pvals) <- col_names
rownames(Duration_Pvals) <- row_names
Occurance_Pvals = Duration_Pvals
Coverage_Pvals = Duration_Pvals
GEV_Pvals = Duration_Pvals
GFP_Pvals = Duration_Pvals
TP1_Pvals = Duration_Pvals
TP2_Pvals = Duration_Pvals
TP3_Pvals = Duration_Pvals
TP4_Pvals = Duration_Pvals
TP5_Pvals = Duration_Pvals

## MS1
# RS1
Fullmod <- lmer(UnpleasantnessRatings ~ UnpleasantnessT0 + GFP*Arm + Duration*Arm + Occurance*Arm + Coverage*Arm + GEV*Arm + TP2*Arm + TP3*Arm + TP4*Arm + TP5*Arm + (1|Subject), data = MS1RS1,REML = FALSE)
anova(Fullmod)
summary(Fullmod)
Duration_Pvals["MS1","RS1"]  = summary(Fullmod)$coefficients['ArmMFS:Duration','Pr(>|t|)']
Occurance_Pvals["MS1","RS1"] = summary(Fullmod)$coefficients['ArmMFS:Occurance','Pr(>|t|)']
Coverage_Pvals["MS1","RS1"]  = summary(Fullmod)$coefficients['ArmMFS:Coverage','Pr(>|t|)']
GEV_Pvals["MS1","RS1"]       = summary(Fullmod)$coefficients['ArmMFS:GEV','Pr(>|t|)']
GFP_Pvals["MS1","RS1"]       = summary(Fullmod)$coefficients['GFP:ArmMFS','Pr(>|t|)']
TP2_Pvals["MS1","RS1"]       = summary(Fullmod)$coefficients['ArmMFS:TP2','Pr(>|t|)']
TP3_Pvals["MS1","RS1"]       = summary(Fullmod)$coefficients['ArmMFS:TP3','Pr(>|t|)']
TP4_Pvals["MS1","RS1"]       = summary(Fullmod)$coefficients['ArmMFS:TP4','Pr(>|t|)']
TP5_Pvals["MS1","RS1"]       = summary(Fullmod)$coefficients['ArmMFS:TP5','Pr(>|t|)']

# # RS2
# Fullmod <- lmer(UnpleasantnessRatings ~ UnpleasantnessT0 + GFP:Arm + Duration:Arm + Occurance:Arm + Coverage:Arm + GEV:Arm + TP2:Arm + TP3:Arm + TP4:Arm + TP5:Arm + (1|Subject), data = MS1RS2,REML = FALSE)
# anova(Fullmod)
# Duration_Pvals["MS1","RS2"]  = summary(Fullmod)$coefficients['Arm:Duration','Pr(>|t|)']
# Occurance_Pvals["MS1","RS2"] = summary(Fullmod)$coefficients['Arm:Occurance','Pr(>|t|)']
# Coverage_Pvals["MS1","RS2"]  = summary(Fullmod)$coefficients['Arm:Coverage','Pr(>|t|)']
# GEV_Pvals["MS1","RS2"]       = summary(Fullmod)$coefficients['Arm:GEV','Pr(>|t|)']
# TP1_Pvals["MS1","RS2"]       = summary(Fullmod)$coefficients['Arm:TP1','Pr(>|t|)']
# TP2_Pvals["MS1","RS2"]       = summary(Fullmod)$coefficients['Arm:TP2','Pr(>|t|)']
# TP3_Pvals["MS1","RS2"]       = summary(Fullmod)$coefficients['Arm:TP3','Pr(>|t|)']
# TP4_Pvals["MS1","RS2"]       = summary(Fullmod)$coefficients['Arm:TP4','Pr(>|t|)']
# TP5_Pvals["MS1","RS2"]       = summary(Fullmod)$coefficients['Arm:TP5','Pr(>|t|)']
# 
# # Diff
# Fullmod <- lmer(UnpleasantnessRatings ~ UnpleasantnessT0 + GFP:Arm + Duration:Arm + Occurance:Arm + Coverage:Arm + GEV:Arm + TP2:Arm + TP3:Arm + TP4:Arm + TP5:Arm + (1|Subject), data = MS1diff,REML = FALSE)
# anova(Fullmod)
# Duration_Pvals["MS1","Diff"]  = summary(Fullmod)$coefficients['Arm:Duration','Pr(>|t|)']
# Occurance_Pvals["MS1","Diff"] = summary(Fullmod)$coefficients['Arm:Occurance','Pr(>|t|)']
# Coverage_Pvals["MS1","Diff"]  = summary(Fullmod)$coefficients['Arm:Coverage','Pr(>|t|)']
# GEV_Pvals["MS1","Diff"]       = summary(Fullmod)$coefficients['Arm:GEV','Pr(>|t|)']
# TP1_Pvals["MS1","Diff"]       = summary(Fullmod)$coefficients['Arm:TP1','Pr(>|t|)']
# TP2_Pvals["MS1","Diff"]       = summary(Fullmod)$coefficients['Arm:TP2','Pr(>|t|)']
# TP3_Pvals["MS1","Diff"]       = summary(Fullmod)$coefficients['Arm:TP3','Pr(>|t|)']
# TP4_Pvals["MS1","Diff"]       = summary(Fullmod)$coefficients['Arm:TP4','Pr(>|t|)']
# TP5_Pvals["MS1","Diff"]       = summary(Fullmod)$coefficients['Arm:TP5','Pr(>|t|)']


## MS2
# RS1
Fullmod <- lmer(UnpleasantnessRatings ~ UnpleasantnessT0 + GFP:Arm + Duration:Arm + Occurance:Arm + Coverage:Arm + GEV:Arm + TP1:Arm + TP3:Arm + TP4:Arm + TP5:Arm + (1|Subject), data = MS2RS1,REML = FALSE)
anova(Fullmod)
Duration_Pvals["MS2","RS1"]  = summary(Fullmod)$coefficients['ArmMFS:Duration','Pr(>|t|)']
Occurance_Pvals["MS2","RS1"] = summary(Fullmod)$coefficients['ArmMFS:Occurance','Pr(>|t|)']
Coverage_Pvals["MS2","RS1"]  = summary(Fullmod)$coefficients['ArmMFS:Coverage','Pr(>|t|)']
GEV_Pvals["MS2","RS1"]       = summary(Fullmod)$coefficients['ArmMFS:GEV','Pr(>|t|)']
GFP_Pvals["MS1","RS1"]       = summary(Fullmod)$coefficients['GFP:ArmMFS','Pr(>|t|)']
TP1_Pvals["MS2","RS1"]       = summary(Fullmod)$coefficients['ArmMFS:TP1','Pr(>|t|)']
TP3_Pvals["MS2","RS1"]       = summary(Fullmod)$coefficients['ArmMFS:TP3','Pr(>|t|)']
TP4_Pvals["MS2","RS1"]       = summary(Fullmod)$coefficients['ArmMFS:TP4','Pr(>|t|)']
TP5_Pvals["MS2","RS1"]       = summary(Fullmod)$coefficients['ArmMFS:TP5','Pr(>|t|)']

# # RS2
# Fullmod <- lmer(UnpleasantnessRatings ~ UnpleasantnessT0 + GFP:Arm + Duration:Arm + Occurance:Arm + Coverage:Arm + GEV:Arm + TP1:Arm + TP3:Arm + TP4:Arm + TP5:Arm + (1|Subject), data = MS2RS2,REML = FALSE)
# anova(Fullmod)
# Duration_Pvals["MS2","RS2"]  = summary(Fullmod)$coefficients['Arm:Duration','Pr(>|t|)']
# Occurance_Pvals["MS2","RS2"] = summary(Fullmod)$coefficients['Arm:Occurance','Pr(>|t|)']
# Coverage_Pvals["MS2","RS2"]  = summary(Fullmod)$coefficients['Arm:Coverage','Pr(>|t|)']
# GEV_Pvals["MS2","RS2"]       = summary(Fullmod)$coefficients['Arm:GEV','Pr(>|t|)']
# TP1_Pvals["MS2","RS2"]       = summary(Fullmod)$coefficients['Arm:TP1','Pr(>|t|)']
# TP2_Pvals["MS2","RS2"]       = summary(Fullmod)$coefficients['Arm:TP2','Pr(>|t|)']
# TP3_Pvals["MS2","RS2"]       = summary(Fullmod)$coefficients['Arm:TP3','Pr(>|t|)']
# TP4_Pvals["MS2","RS2"]       = summary(Fullmod)$coefficients['Arm:TP4','Pr(>|t|)']
# TP5_Pvals["MS2","RS2"]       = summary(Fullmod)$coefficients['Arm:TP5','Pr(>|t|)']
# 
# # Diff
# Fullmod <- lmer(UnpleasantnessRatings ~ UnpleasantnessT0 + GFP:Arm + Duration:Arm + Occurance:Arm + Coverage:Arm + GEV:Arm + TP1:Arm + TP3:Arm + TP4:Arm + TP5:Arm + (1|Subject), data = MS2diff,REML = FALSE)
# anova(Fullmod)
# Duration_Pvals["MS2","Diff"]  = summary(Fullmod)$coefficients['Arm:Duration','Pr(>|t|)']
# Occurance_Pvals["MS2","Diff"] = summary(Fullmod)$coefficients['Arm:Occurance','Pr(>|t|)']
# Coverage_Pvals["MS2","Diff"]  = summary(Fullmod)$coefficients['Arm:Coverage','Pr(>|t|)']
# GEV_Pvals["MS2","Diff"]       = summary(Fullmod)$coefficients['Arm:GEV','Pr(>|t|)']
# TP1_Pvals["MS2","Diff"]       = summary(Fullmod)$coefficients['Arm:TP1','Pr(>|t|)']
# TP2_Pvals["MS2","Diff"]       = summary(Fullmod)$coefficients['Arm:TP2','Pr(>|t|)']
# TP3_Pvals["MS2","Diff"]       = summary(Fullmod)$coefficients['Arm:TP3','Pr(>|t|)']
# TP4_Pvals["MS2","Diff"]       = summary(Fullmod)$coefficients['Arm:TP4','Pr(>|t|)']
# TP5_Pvals["MS2","Diff"]       = summary(Fullmod)$coefficients['Arm:TP5','Pr(>|t|)']


## MS3
# RS1
Fullmod <- lmer(UnpleasantnessRatings ~ UnpleasantnessT0 + GFP:Arm + Duration:Arm + Occurance:Arm + Coverage:Arm + GEV:Arm + TP2:Arm + TP1:Arm + TP4:Arm + TP5:Arm + (1|Subject), data = MS3RS1,REML = FALSE)
anova(Fullmod)
Duration_Pvals["MS3","RS1"]  = summary(Fullmod)$coefficients['ArmMFS:Duration','Pr(>|t|)']
Occurance_Pvals["MS3","RS1"] = summary(Fullmod)$coefficients['ArmMFS:Occurance','Pr(>|t|)']
Coverage_Pvals["MS3","RS1"]  = summary(Fullmod)$coefficients['ArmMFS:Coverage','Pr(>|t|)']
GEV_Pvals["MS3","RS1"]       = summary(Fullmod)$coefficients['ArmMFS:GEV','Pr(>|t|)']
GFP_Pvals["MS1","RS1"]       = summary(Fullmod)$coefficients['GFP:ArmMFS','Pr(>|t|)']
TP1_Pvals["MS3","RS1"]       = summary(Fullmod)$coefficients['ArmMFS:TP1','Pr(>|t|)']
TP2_Pvals["MS3","RS1"]       = summary(Fullmod)$coefficients['ArmMFS:TP2','Pr(>|t|)']
TP4_Pvals["MS3","RS1"]       = summary(Fullmod)$coefficients['ArmMFS:TP4','Pr(>|t|)']
TP5_Pvals["MS3","RS1"]       = summary(Fullmod)$coefficients['ArmMFS:TP5','Pr(>|t|)']

Fullmod <- lmer(UnpleasantnessRatings ~ UnpleasantnessT0 + TP1*Arm + (1|Subject), data = MS3RS1,REML = FALSE)
summary(Fullmod)
Intplot <- ggpredict(Fullmod,terms = c('Arm','TP1'))
plot(Intplot,add.data = TRUE)

# # RS2
# Fullmod <- lmer(UnpleasantnessRatings ~ UnpleasantnessT0 + GFP:Arm + Duration:Arm + Occurance:Arm + Coverage:Arm + GEV:Arm + TP2:Arm + TP1:Arm + TP4:Arm + TP5:Arm + (1|Subject), data = MS3RS2,REML = FALSE)
# anova(Fullmod)
# Duration_Pvals["MS3","RS2"]  = summary(Fullmod)$coefficients['Arm:Duration','Pr(>|t|)']
# Occurance_Pvals["MS3","RS2"] = summary(Fullmod)$coefficients['Arm:Occurance','Pr(>|t|)']
# Coverage_Pvals["MS3","RS2"]  = summary(Fullmod)$coefficients['Arm:Coverage','Pr(>|t|)']
# GEV_Pvals["MS3","RS2"]       = summary(Fullmod)$coefficients['Arm:GEV','Pr(>|t|)']
# TP1_Pvals["MS3","RS2"]       = summary(Fullmod)$coefficients['Arm:TP1','Pr(>|t|)']
# TP2_Pvals["MS3","RS2"]       = summary(Fullmod)$coefficients['Arm:TP2','Pr(>|t|)']
# TP3_Pvals["MS3","RS2"]       = summary(Fullmod)$coefficients['Arm:TP3','Pr(>|t|)']
# TP4_Pvals["MS3","RS2"]       = summary(Fullmod)$coefficients['Arm:TP4','Pr(>|t|)']
# TP5_Pvals["MS3","RS2"]       = summary(Fullmod)$coefficients['Arm:TP5','Pr(>|t|)']
# 
# # Diff
# Fullmod <- lmer(UnpleasantnessRatings ~ UnpleasantnessT0 + GFP:Arm + Duration:Arm + Occurance:Arm + Coverage:Arm + GEV:Arm + TP2:Arm + TP1:Arm + TP4:Arm + TP5:Arm + (1|Subject), data = MS3diff,REML = FALSE)
# anova(Fullmod)
# Duration_Pvals["MS3","Diff"]  = summary(Fullmod)$coefficients['Arm:Duration','Pr(>|t|)']
# Occurance_Pvals["MS3","Diff"] = summary(Fullmod)$coefficients['Arm:Occurance','Pr(>|t|)']
# Coverage_Pvals["MS3","Diff"]  = summary(Fullmod)$coefficients['Arm:Coverage','Pr(>|t|)']
# GEV_Pvals["MS3","Diff"]       = summary(Fullmod)$coefficients['Arm:GEV','Pr(>|t|)']
# TP1_Pvals["MS3","Diff"]       = summary(Fullmod)$coefficients['Arm:TP1','Pr(>|t|)']
# TP2_Pvals["MS3","Diff"]       = summary(Fullmod)$coefficients['Arm:TP2','Pr(>|t|)']
# TP3_Pvals["MS3","Diff"]       = summary(Fullmod)$coefficients['Arm:TP3','Pr(>|t|)']
# TP4_Pvals["MS3","Diff"]       = summary(Fullmod)$coefficients['Arm:TP4','Pr(>|t|)']
# TP5_Pvals["MS3","Diff"]       = summary(Fullmod)$coefficients['Arm:TP5','Pr(>|t|)']



## MS4
# RS1
Fullmod <- lmer(UnpleasantnessRatings ~ UnpleasantnessT0 + GFP:Arm + Duration:Arm + Occurance:Arm + Coverage:Arm + GEV:Arm + TP2:Arm + TP3:Arm + TP1:Arm + TP5:Arm + (1|Subject), data = MS4RS1,REML = FALSE)
anova(Fullmod)
Duration_Pvals["MS4","RS1"]  = summary(Fullmod)$coefficients['ArmMFS:Duration','Pr(>|t|)']
Occurance_Pvals["MS4","RS1"] = summary(Fullmod)$coefficients['ArmMFS:Occurance','Pr(>|t|)']
Coverage_Pvals["MS4","RS1"]  = summary(Fullmod)$coefficients['ArmMFS:Coverage','Pr(>|t|)']
GEV_Pvals["MS4","RS1"]       = summary(Fullmod)$coefficients['ArmMFS:GEV','Pr(>|t|)']
GFP_Pvals["MS1","RS1"]       = summary(Fullmod)$coefficients['GFP:ArmMFS','Pr(>|t|)']
TP1_Pvals["MS4","RS1"]       = summary(Fullmod)$coefficients['ArmMFS:TP1','Pr(>|t|)']
TP2_Pvals["MS4","RS1"]       = summary(Fullmod)$coefficients['ArmMFS:TP2','Pr(>|t|)']
TP3_Pvals["MS4","RS1"]       = summary(Fullmod)$coefficients['ArmMFS:TP3','Pr(>|t|)']
TP5_Pvals["MS4","RS1"]       = summary(Fullmod)$coefficients['ArmMFS:TP5','Pr(>|t|)']

# # RS2
# Fullmod <- lmer(UnpleasantnessRatings ~ UnpleasantnessT0 + GFP:Arm + Duration:Arm + Occurance:Arm + Coverage:Arm + GEV:Arm + TP2:Arm + TP3:Arm + TP1:Arm + TP5:Arm + (1|Subject), data = MS4RS2,REML = FALSE)
# anova(Fullmod)
# Duration_Pvals["MS4","RS2"]  = summary(Fullmod)$coefficients['Arm:Duration','Pr(>|t|)']
# Occurance_Pvals["MS4","RS2"] = summary(Fullmod)$coefficients['Arm:Occurance','Pr(>|t|)']
# Coverage_Pvals["MS4","RS2"]  = summary(Fullmod)$coefficients['Arm:Coverage','Pr(>|t|)']
# GEV_Pvals["MS4","RS2"]       = summary(Fullmod)$coefficients['Arm:GEV','Pr(>|t|)']
# TP1_Pvals["MS4","RS2"]       = summary(Fullmod)$coefficients['Arm:TP1','Pr(>|t|)']
# TP2_Pvals["MS4","RS2"]       = summary(Fullmod)$coefficients['Arm:TP2','Pr(>|t|)']
# TP3_Pvals["MS4","RS2"]       = summary(Fullmod)$coefficients['Arm:TP3','Pr(>|t|)']
# TP4_Pvals["MS4","RS2"]       = summary(Fullmod)$coefficients['Arm:TP4','Pr(>|t|)']
# TP5_Pvals["MS4","RS2"]       = summary(Fullmod)$coefficients['Arm:TP5','Pr(>|t|)']
# 
# # Diff
# Fullmod <- lmer(UnpleasantnessRatings ~ UnpleasantnessT0 + GFP:Arm + Duration:Arm + Occurance:Arm + Coverage:Arm + GEV:Arm + TP2:Arm + TP3:Arm + TP1:Arm + TP5:Arm + (1|Subject), data = MS4diff,REML = FALSE)
# anova(Fullmod)
# Duration_Pvals["MS4","Diff"]  = summary(Fullmod)$coefficients['Arm:Duration','Pr(>|t|)']
# Occurance_Pvals["MS4","Diff"] = summary(Fullmod)$coefficients['Arm:Occurance','Pr(>|t|)']
# Coverage_Pvals["MS4","Diff"]  = summary(Fullmod)$coefficients['Arm:Coverage','Pr(>|t|)']
# GEV_Pvals["MS4","Diff"]       = summary(Fullmod)$coefficients['Arm:GEV','Pr(>|t|)']
# TP1_Pvals["MS4","Diff"]       = summary(Fullmod)$coefficients['Arm:TP1','Pr(>|t|)']
# TP2_Pvals["MS4","Diff"]       = summary(Fullmod)$coefficients['Arm:TP2','Pr(>|t|)']
# TP3_Pvals["MS4","Diff"]       = summary(Fullmod)$coefficients['Arm:TP3','Pr(>|t|)']
# TP4_Pvals["MS4","Diff"]       = summary(Fullmod)$coefficients['Arm:TP4','Pr(>|t|)']
# TP5_Pvals["MS4","Diff"]       = summary(Fullmod)$coefficients['Arm:TP5','Pr(>|t|)']

## MS5
# RS1
Fullmod <- lmer(UnpleasantnessRatings ~ UnpleasantnessT0 + GFP*Arm + Duration*Arm + Occurance*Arm + Coverage*Arm + GEV*Arm + TP2*Arm + TP3*Arm + TP4*Arm + TP1*Arm + (1|Subject), data = MS5RS1,REML = FALSE)
anova(Fullmod)
Duration_Pvals["MS5","RS1"]  = summary(Fullmod)$coefficients['ArmMFS:Duration','Pr(>|t|)']
Occurance_Pvals["MS5","RS1"] = summary(Fullmod)$coefficients['ArmMFS:Occurance','Pr(>|t|)']
Coverage_Pvals["MS5","RS1"]  = summary(Fullmod)$coefficients['ArmMFS:Coverage','Pr(>|t|)']
GEV_Pvals["MS5","RS1"]       = summary(Fullmod)$coefficients['ArmMFS:GEV','Pr(>|t|)']
GFP_Pvals["MS1","RS1"]       = summary(Fullmod)$coefficients['GFP:ArmMFS','Pr(>|t|)']
TP1_Pvals["MS5","RS1"]       = summary(Fullmod)$coefficients['ArmMFS:TP1','Pr(>|t|)']
TP2_Pvals["MS5","RS1"]       = summary(Fullmod)$coefficients['ArmMFS:TP2','Pr(>|t|)']
TP3_Pvals["MS5","RS1"]       = summary(Fullmod)$coefficients['ArmMFS:TP3','Pr(>|t|)']
TP4_Pvals["MS5","RS1"]       = summary(Fullmod)$coefficients['ArmMFS:TP4','Pr(>|t|)']

# # RS2
# Fullmod <- lmer(UnpleasantnessRatings ~ UnpleasantnessT0 + GFP:Arm + Duration:Arm + Occurance:Arm + Coverage:Arm + GEV:Arm + TP2:Arm + TP3:Arm + TP4:Arm + TP1:Arm + (1|Subject), data = MS5RS2,REML = FALSE)
# anova(Fullmod)
# Duration_Pvals["MS5","RS2"]  = summary(Fullmod)$coefficients['Arm:Duration','Pr(>|t|)']
# Occurance_Pvals["MS5","RS2"] = summary(Fullmod)$coefficients['Arm:Occurance','Pr(>|t|)']
# Coverage_Pvals["MS5","RS2"]  = summary(Fullmod)$coefficients['Arm:Coverage','Pr(>|t|)']
# GEV_Pvals["MS5","RS2"]       = summary(Fullmod)$coefficients['Arm:GEV','Pr(>|t|)']
# TP1_Pvals["MS5","RS2"]       = summary(Fullmod)$coefficients['Arm:TP1','Pr(>|t|)']
# TP2_Pvals["MS5","RS2"]       = summary(Fullmod)$coefficients['Arm:TP2','Pr(>|t|)']
# TP3_Pvals["MS5","RS2"]       = summary(Fullmod)$coefficients['Arm:TP3','Pr(>|t|)']
# TP4_Pvals["MS5","RS2"]       = summary(Fullmod)$coefficients['Arm:TP4','Pr(>|t|)']
# TP5_Pvals["MS5","RS2"]       = summary(Fullmod)$coefficients['Arm:TP5','Pr(>|t|)']
# 
# # Diff
# Fullmod <- lmer(UnpleasantnessRatings ~ UnpleasantnessT0 + GFP:Arm + Duration:Arm + Occurance:Arm + Coverage:Arm + GEV:Arm + TP2:Arm + TP3:Arm + TP4:Arm + TP1:Arm + (1|Subject), data = MS5diff,REML = FALSE)
# anova(Fullmod)
# Duration_Pvals["MS5","Diff"]  = summary(Fullmod)$coefficients['Arm:Duration','Pr(>|t|)']
# Occurance_Pvals["MS5","Diff"] = summary(Fullmod)$coefficients['Arm:Occurance','Pr(>|t|)']
# Coverage_Pvals["MS5","Diff"]  = summary(Fullmod)$coefficients['Arm:Coverage','Pr(>|t|)']
# GEV_Pvals["MS5","Diff"]       = summary(Fullmod)$coefficients['Arm:GEV','Pr(>|t|)']
# TP1_Pvals["MS5","Diff"]       = summary(Fullmod)$coefficients['Arm:TP1','Pr(>|t|)']
# TP2_Pvals["MS5","Diff"]       = summary(Fullmod)$coefficients['Arm:TP2','Pr(>|t|)']
# TP3_Pvals["MS5","Diff"]       = summary(Fullmod)$coefficients['Arm:TP3','Pr(>|t|)']
# TP4_Pvals["MS5","Diff"]       = summary(Fullmod)$coefficients['Arm:TP4','Pr(>|t|)']
# TP5_Pvals["MS5","Diff"]       = summary(Fullmod)$coefficients['Arm:TP5','Pr(>|t|)']

# Check for significance after FDR correction
Pvalues <- unlist(TP5_Pvals)
non_na_indices <- which(!is.na(Pvalues))
All_adj_p_values <- p.adjust(Pvalues[non_na_indices], method = "fdr")
Pvalues[non_na_indices] <- All_adj_p_values
Pvalues[Pvalues<.05]


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
non_na_indices <- which(!is.na(Pvalues))
All_adj_p_values <- p.adjust(Pvalues[non_na_indices], method = "fdr")
Pvalues[non_na_indices] <- All_adj_p_values
Pvalues[Pvalues<.05]



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



## Double Subtracted ----
## Microstates x PPRatings ----
MicrostatsxPP <- read_excel("Path_to_Microstate_PPRatings_data/MicrostatsxPPratings.xlsx")
View(MicrostatsxPP)

HFS  = MicrostatsxPP[MicrostatsxPP$Arm == "MFS",]
Cont = MicrostatsxPP[MicrostatsxPP$Arm == "Cont",]

library(pracma)
ind = linspace(1,1991,200)
temp = HFS[ind,]
temp$DSInt = NA
temp$DSUnp = NA

index = 1
# Initialize a for loop
for (i in ind) {
  
  HFSIntT1 = mean(HFS$IntensityRatings[i:(i+9)])
  HFSUnpT1 = mean(HFS$UnpleasantnessRatings[i:(i+9)])
  ContIntT1 = mean(Cont$IntensityRatings[i:(i+9)])
  ContUnpT1 = mean(Cont$UnpleasantnessRatings[i:(i+9)])
  
  HFSIntT0 = mean(HFS$IntensityT0[i:(i+9)])
  HFSUnpT0 = mean(HFS$UnpleasantnessT0[i:(i+9)])
  ContIntT0 = mean(Cont$IntensityT0[i:(i+9)])
  ContUnpT0 = mean(Cont$UnpleasantnessT0[i:(i+9)])
  
  DSInt = (HFSIntT1 - HFSIntT0) - (ContIntT1 - ContIntT0)
  DSUnp = (HFSUnpT1 - HFSUnpT0) - (ContUnpT1 - ContUnpT0)
  
  temp$DSInt[index] = DSInt
  temp$DSUnp[index] = DSUnp
  
  index = index + 1
  
}

MicrostatsxPP = temp

# Standardize predictors
MicrostatsxPP$GFP <- as.numeric(scale(MicrostatsxPP$GFP))
MicrostatsxPP$Duration <- as.numeric(scale(MicrostatsxPP$Duration))
MicrostatsxPP$Occurance <- as.numeric(scale(MicrostatsxPP$Occurance))

# Change the sign for participant 21 since there was a mistake with control vs HFS
MicrostatsxPP$DSInt[MicrostatsxPP$Subject==21] = MicrostatsxPP$DSInt[MicrostatsxPP$Subject==21]*-1
MicrostatsxPP$DSUnp[MicrostatsxPP$Subject==21] = MicrostatsxPP$DSUnp[MicrostatsxPP$Subject==21]*-1

# 
# # Create the double subtracted data
# MicrostatsxPP$Intsubtr = MicrostatsxPP$IntensityRatings - MicrostatsxPP$IntensityT0
# MicrostatsxPP$Unpsubtr = MicrostatsxPP$UnpleasantnessRatings - MicrostatsxPP$UnpleasantnessT0
# 
# HFS$DSint = HFS$Intsubtr - Cont$Intsubtr
# HFS$DSUnp = HFS$Unpsubtr - Cont$Unpsubtr
# 
# MicrostatsxPP = HFS

# Function to compute p-values for individual predictors
Fit_Models_DS <- function(predictor_names, dependent_variable, data,MS) {
  # Initialize a named vector to store p-values
  p_values <- numeric(length(predictor_names))
  new_predictor_names <- paste0(MS, predictor_names)
  p_values <- setNames(c(p_values), c(new_predictor_names))

  # Loop over each predictor
  for (i in 1:length(predictor_names)) {
    # Create the formula for the model
    formula <- formula(paste(dependent_variable, "~", predictor_names[i]))
    
    # Fit the linear model
    model <- lm(formula, data = data)
    print(summary(model))
    # Extract the p-value for the predictor and store it with the predictor name
    print(summary(model)$coefficients[predictor_names[i], "Pr(>|t|)"])
    p_values[new_predictor_names[i]] <- summary(model)$coefficients[2, "Pr(>|t|)"]
  }
  # Return the named vector of p-values
  print(p_values)
  return(p_values)
}


# MicrostatsxPP[MicrostatsxPP$DSInt < 0,] = NA
# MicrostatsxPP$DSInt = log(MicrostatsxPP$DSInt)
# hist(MicrostatsxPP$DSInt)

MicrostatsxPP[MicrostatsxPP$DSUnp < 0,] = NA
MicrostatsxPP$DSUnp = log(MicrostatsxPP$DSUnp)
hist(MicrostatsxPP$DSUnp)

MS1RS1 = MicrostatsxPP[MicrostatsxPP$MS == 'MS1' & MicrostatsxPP$RS == 'RS1',]
MS2RS1 = MicrostatsxPP[MicrostatsxPP$MS == 'MS2' & MicrostatsxPP$RS == 'RS1',]
MS3RS1 = MicrostatsxPP[MicrostatsxPP$MS == 'MS3' & MicrostatsxPP$RS == 'RS1',]
MS4RS1 = MicrostatsxPP[MicrostatsxPP$MS == 'MS4' & MicrostatsxPP$RS == 'RS1',]
MS5RS1 = MicrostatsxPP[MicrostatsxPP$MS == 'MS5' & MicrostatsxPP$RS == 'RS1',]

# MS1RS2 = MicrostatsxPP[MicrostatsxPP$MS == 'MS1' & MicrostatsxPP$RS == 'RS2',]
# MS2RS2 = MicrostatsxPP[MicrostatsxPP$MS == 'MS2' & MicrostatsxPP$RS == 'RS2',]
# MS3RS2 = MicrostatsxPP[MicrostatsxPP$MS == 'MS3' & MicrostatsxPP$RS == 'RS2',]
# MS4RS2 = MicrostatsxPP[MicrostatsxPP$MS == 'MS4' & MicrostatsxPP$RS == 'RS2',]
# MS5RS2 = MicrostatsxPP[MicrostatsxPP$MS == 'MS5' & MicrostatsxPP$RS == 'RS2',]
# 
# MS1diff = MS1RS1
# MS1diff[,c(6:16)] = MS1RS1[,c(6:16)] - MS1RS2[,c(6:16)]
# MS2diff = MS2RS1
# MS2diff[,c(6:16)] = MS2RS1[,c(6:16)] - MS2RS2[,c(6:16)]
# MS3diff = MS3RS1
# MS3diff[,c(6:16)] = MS3RS1[,c(6:16)] - MS3RS2[,c(6:16)]
# MS4diff = MS4RS1
# MS4diff[,c(6:16)] = MS4RS1[,c(6:16)] - MS4RS2[,c(6:16)]
# MS5diff = MS5RS1
# MS5diff[,c(6:16)] = MS5RS1[,c(6:16)] - MS5RS2[,c(6:16)]

### Intensity ----
# Initialize variables

# Create column and row names
col_names <- c("RS1","RS2","Diff")
row_names <- c("MS1", "MS2", "MS3", "MS4", "MS5")
# Create empty data frames
Duration_Pvals <- data.frame(matrix(ncol = length(col_names), nrow = length(row_names)))
colnames(Duration_Pvals) <- col_names
rownames(Duration_Pvals) <- row_names
Occurance_Pvals = Duration_Pvals
Coverage_Pvals = Duration_Pvals
GEV_Pvals = Duration_Pvals
GFP_Pvals = Duration_Pvals
TP1_Pvals = Duration_Pvals
TP2_Pvals = Duration_Pvals
TP3_Pvals = Duration_Pvals
TP4_Pvals = Duration_Pvals
TP5_Pvals = Duration_Pvals


## MS1
# RS1
# Fullmod <- lmer(DSInt ~ GFP  + Duration  + Occurance  + Coverage  + GEV  + TP2  + TP3  + TP4  + TP5  + (1|Subject), data = MS1RS1,REML = FALSE)
Fullmod <- lm(DSInt ~ GFP  + Duration  + Occurance  + Coverage  + GEV  + TP2  + TP3  + TP4  + TP5, data = MS1RS1)
summary(Fullmod)
Duration_Pvals["MS1","RS1"]  = summary(Fullmod)$coefficients['Duration','Pr(>|t|)']
Occurance_Pvals["MS1","RS1"] = summary(Fullmod)$coefficients['Occurance','Pr(>|t|)']
Coverage_Pvals["MS1","RS1"]  = summary(Fullmod)$coefficients['Coverage','Pr(>|t|)']
GEV_Pvals["MS1","RS1"]       = summary(Fullmod)$coefficients['GEV','Pr(>|t|)']
GFP_Pvals["MS1","RS1"]       = summary(Fullmod)$coefficients['GFP','Pr(>|t|)']
TP2_Pvals["MS1","RS1"]       = summary(Fullmod)$coefficients['TP2','Pr(>|t|)']
TP3_Pvals["MS1","RS1"]       = summary(Fullmod)$coefficients['TP3','Pr(>|t|)']
TP4_Pvals["MS1","RS1"]       = summary(Fullmod)$coefficients['TP4','Pr(>|t|)']
TP5_Pvals["MS1","RS1"]       = summary(Fullmod)$coefficients['TP5','Pr(>|t|)']


# # RS2
# Fullmod <- lmer(DSInt ~ GFP  + Duration  + Occurance  + Coverage  + GEV  + TP2  + TP3  + TP4  + TP5  + (1|Subject), data = MS1RS2,REML = FALSE)
# anova(Fullmod)
# Duration_Pvals["MS1","RS2"]  = summary(Fullmod)$coefficients['Duration','Pr(>|t|)']
# Occurance_Pvals["MS1","RS2"] = summary(Fullmod)$coefficients['Occurance','Pr(>|t|)']
# Coverage_Pvals["MS1","RS2"]  = summary(Fullmod)$coefficients['Coverage','Pr(>|t|)']
# GEV_Pvals["MS1","RS2"]       = summary(Fullmod)$coefficients['GEV','Pr(>|t|)']
# TP1_Pvals["MS1","RS2"]       = summary(Fullmod)$coefficients['TP1','Pr(>|t|)']
# TP2_Pvals["MS1","RS2"]       = summary(Fullmod)$coefficients['TP2','Pr(>|t|)']
# TP3_Pvals["MS1","RS2"]       = summary(Fullmod)$coefficients['TP3','Pr(>|t|)']
# TP4_Pvals["MS1","RS2"]       = summary(Fullmod)$coefficients['TP4','Pr(>|t|)']
# TP5_Pvals["MS1","RS2"]       = summary(Fullmod)$coefficients['TP5','Pr(>|t|)']
# # Diff
# Fullmod <- lmer(DSInt ~ GFP  + Duration  + Occurance  + Coverage  + GEV  + TP2  + TP3  + TP4  + TP5  + (1|Subject), data = MS1diff,REML = FALSE)
# anova(Fullmod)
# Duration_Pvals["MS1","Diff"]  = summary(Fullmod)$coefficients['Duration','Pr(>|t|)']
# Occurance_Pvals["MS1","Diff"] = summary(Fullmod)$coefficients['Occurance','Pr(>|t|)']
# Coverage_Pvals["MS1","Diff"]  = summary(Fullmod)$coefficients['Coverage','Pr(>|t|)']
# GEV_Pvals["MS1","Diff"]       = summary(Fullmod)$coefficients['GEV','Pr(>|t|)']
# TP1_Pvals["MS1","Diff"]       = summary(Fullmod)$coefficients['TP1','Pr(>|t|)']
# TP2_Pvals["MS1","Diff"]       = summary(Fullmod)$coefficients['TP2','Pr(>|t|)']
# TP3_Pvals["MS1","Diff"]       = summary(Fullmod)$coefficients['TP3','Pr(>|t|)']
# TP4_Pvals["MS1","Diff"]       = summary(Fullmod)$coefficients['TP4','Pr(>|t|)']
# TP5_Pvals["MS1","Diff"]       = summary(Fullmod)$coefficients['TP5','Pr(>|t|)']


## MS2
# RS1
# Fullmod <- lmer(DSInt ~ GFP  + Duration  + Occurance  + Coverage  + GEV  + TP1  + TP3  + TP4  + TP5  + (1|Subject), data = MS2RS1,REML = FALSE)
Fullmod <- lm(DSInt ~ GFP + Duration  + Occurance  + Coverage  + GEV  + TP1  + TP3  + TP4  + TP5, data = MS2RS1)
summary(Fullmod)
anova(Fullmod)
Duration_Pvals["MS2","RS1"]  = summary(Fullmod)$coefficients['Duration','Pr(>|t|)']
Occurance_Pvals["MS2","RS1"] = summary(Fullmod)$coefficients['Occurance','Pr(>|t|)']
Coverage_Pvals["MS2","RS1"]  = summary(Fullmod)$coefficients['Coverage','Pr(>|t|)']
GEV_Pvals["MS2","RS1"]       = summary(Fullmod)$coefficients['GEV','Pr(>|t|)']
GFP_Pvals["MS2","RS1"]       = summary(Fullmod)$coefficients['GFP','Pr(>|t|)']
TP1_Pvals["MS2","RS1"]       = summary(Fullmod)$coefficients['TP1','Pr(>|t|)']
TP3_Pvals["MS2","RS1"]       = summary(Fullmod)$coefficients['TP3','Pr(>|t|)']
TP4_Pvals["MS2","RS1"]       = summary(Fullmod)$coefficients['TP4','Pr(>|t|)']
TP5_Pvals["MS2","RS1"]       = summary(Fullmod)$coefficients['TP5','Pr(>|t|)']

Fullmod <- lm(DSInt ~  TP3, data = MS2RS1)
anova(Fullmod)
Intplot <- ggpredict(Fullmod,terms = c('TP3'))
plot(Intplot,add.data = TRUE)

# # RS2
# Fullmod <- lmer(DSInt ~ GFP  + Duration  + Occurance  + Coverage  + GEV  + TP1  + TP3  + TP4  + TP5  + (1|Subject), data = MS2RS2,REML = FALSE)
# anova(Fullmod)
# Duration_Pvals["MS2","RS2"]  = summary(Fullmod)$coefficients['Duration','Pr(>|t|)']
# Occurance_Pvals["MS2","RS2"] = summary(Fullmod)$coefficients['Occurance','Pr(>|t|)']
# Coverage_Pvals["MS2","RS2"]  = summary(Fullmod)$coefficients['Coverage','Pr(>|t|)']
# GEV_Pvals["MS2","RS2"]       = summary(Fullmod)$coefficients['GEV','Pr(>|t|)']
# TP1_Pvals["MS2","RS2"]       = summary(Fullmod)$coefficients['TP1','Pr(>|t|)']
# TP2_Pvals["MS2","RS2"]       = summary(Fullmod)$coefficients['TP2','Pr(>|t|)']
# TP3_Pvals["MS2","RS2"]       = summary(Fullmod)$coefficients['TP3','Pr(>|t|)']
# TP4_Pvals["MS2","RS2"]       = summary(Fullmod)$coefficients['TP4','Pr(>|t|)']
# TP5_Pvals["MS2","RS2"]       = summary(Fullmod)$coefficients['TP5','Pr(>|t|)']
# 
# # Diff
# Fullmod <- lmer(DSInt ~ GFP  + Duration  + Occurance  + Coverage  + GEV  + TP1  + TP3  + TP4  + TP5  + (1|Subject), data = MS2diff,REML = FALSE)
# anova(Fullmod)
# Duration_Pvals["MS2","Diff"]  = summary(Fullmod)$coefficients['Duration','Pr(>|t|)']
# Occurance_Pvals["MS2","Diff"] = summary(Fullmod)$coefficients['Occurance','Pr(>|t|)']
# Coverage_Pvals["MS2","Diff"]  = summary(Fullmod)$coefficients['Coverage','Pr(>|t|)']
# GEV_Pvals["MS2","Diff"]       = summary(Fullmod)$coefficients['GEV','Pr(>|t|)']
# TP1_Pvals["MS2","Diff"]       = summary(Fullmod)$coefficients['TP1','Pr(>|t|)']
# TP2_Pvals["MS2","Diff"]       = summary(Fullmod)$coefficients['TP2','Pr(>|t|)']
# TP3_Pvals["MS2","Diff"]       = summary(Fullmod)$coefficients['TP3','Pr(>|t|)']
# TP4_Pvals["MS2","Diff"]       = summary(Fullmod)$coefficients['TP4','Pr(>|t|)']
# TP5_Pvals["MS2","Diff"]       = summary(Fullmod)$coefficients['TP5','Pr(>|t|)']


## MS3
# RS1
# Fullmod <- lmer(DSInt ~ GFP  + Duration  + Occurance  + Coverage  + GEV  + TP2  + TP1  + TP4  + TP5  + (1|Subject), data = MS3RS1,REML = FALSE)
Fullmod <- lm(DSInt ~ GFP  + Duration  + Occurance  + Coverage  + GEV  + TP2  + TP1  + TP4  + TP5, data = MS3RS1)
summary(Fullmod)
Duration_Pvals["MS3","RS1"]  = summary(Fullmod)$coefficients['Duration','Pr(>|t|)']
Occurance_Pvals["MS3","RS1"] = summary(Fullmod)$coefficients['Occurance','Pr(>|t|)']
Coverage_Pvals["MS3","RS1"]  = summary(Fullmod)$coefficients['Coverage','Pr(>|t|)']
GEV_Pvals["MS3","RS1"]       = summary(Fullmod)$coefficients['GEV','Pr(>|t|)']
GFP_Pvals["MS3","RS1"]       = summary(Fullmod)$coefficients['GFP','Pr(>|t|)']
TP1_Pvals["MS3","RS1"]       = summary(Fullmod)$coefficients['TP1','Pr(>|t|)']
TP2_Pvals["MS3","RS1"]       = summary(Fullmod)$coefficients['TP2','Pr(>|t|)']
TP4_Pvals["MS3","RS1"]       = summary(Fullmod)$coefficients['TP4','Pr(>|t|)']
TP5_Pvals["MS3","RS1"]       = summary(Fullmod)$coefficients['TP5','Pr(>|t|)']

# # RS2
# Fullmod <- lmer(DSInt ~ GFP  + Duration  + Occurance  + Coverage  + GEV  + TP2  + TP1  + TP4  + TP5  + (1|Subject), data = MS3RS2,REML = FALSE)
# anova(Fullmod)
# Duration_Pvals["MS3","RS2"]  = summary(Fullmod)$coefficients['Duration','Pr(>|t|)']
# Occurance_Pvals["MS3","RS2"] = summary(Fullmod)$coefficients['Occurance','Pr(>|t|)']
# Coverage_Pvals["MS3","RS2"]  = summary(Fullmod)$coefficients['Coverage','Pr(>|t|)']
# GEV_Pvals["MS3","RS2"]       = summary(Fullmod)$coefficients['GEV','Pr(>|t|)']
# TP1_Pvals["MS3","RS2"]       = summary(Fullmod)$coefficients['TP1','Pr(>|t|)']
# TP2_Pvals["MS3","RS2"]       = summary(Fullmod)$coefficients['TP2','Pr(>|t|)']
# TP3_Pvals["MS3","RS2"]       = summary(Fullmod)$coefficients['TP3','Pr(>|t|)']
# TP4_Pvals["MS3","RS2"]       = summary(Fullmod)$coefficients['TP4','Pr(>|t|)']
# TP5_Pvals["MS3","RS2"]       = summary(Fullmod)$coefficients['TP5','Pr(>|t|)']
# 
# # Diff
# Fullmod <- lmer(DSInt ~ GFP  + Duration  + Occurance  + Coverage  + GEV  + TP2  + TP1  + TP4  + TP5  + (1|Subject), data = MS3diff,REML = FALSE)
# anova(Fullmod)
# Duration_Pvals["MS3","Diff"]  = summary(Fullmod)$coefficients['Duration','Pr(>|t|)']
# Occurance_Pvals["MS3","Diff"] = summary(Fullmod)$coefficients['Occurance','Pr(>|t|)']
# Coverage_Pvals["MS3","Diff"]  = summary(Fullmod)$coefficients['Coverage','Pr(>|t|)']
# GEV_Pvals["MS3","Diff"]       = summary(Fullmod)$coefficients['GEV','Pr(>|t|)']
# TP1_Pvals["MS3","Diff"]       = summary(Fullmod)$coefficients['TP1','Pr(>|t|)']
# TP2_Pvals["MS3","Diff"]       = summary(Fullmod)$coefficients['TP2','Pr(>|t|)']
# TP3_Pvals["MS3","Diff"]       = summary(Fullmod)$coefficients['TP3','Pr(>|t|)']
# TP4_Pvals["MS3","Diff"]       = summary(Fullmod)$coefficients['TP4','Pr(>|t|)']
# TP5_Pvals["MS3","Diff"]       = summary(Fullmod)$coefficients['TP5','Pr(>|t|)']


## MS4
# RS1
# Fullmod <- lmer(DSInt ~ GFP  + Duration  + Occurance  + Coverage  + GEV  + TP2  + TP3  + TP1  + TP5  + (1|Subject), data = MS4RS1,REML = FALSE)
Fullmod <- lm(DSInt ~ GFP  + Duration  + Occurance  + Coverage  + GEV  + TP2  + TP3  + TP1  + TP5, data = MS4RS1)
summary(Fullmod)
Duration_Pvals["MS4","RS1"]  = summary(Fullmod)$coefficients['Duration','Pr(>|t|)']
Occurance_Pvals["MS4","RS1"] = summary(Fullmod)$coefficients['Occurance','Pr(>|t|)']
Coverage_Pvals["MS4","RS1"]  = summary(Fullmod)$coefficients['Coverage','Pr(>|t|)']
GEV_Pvals["MS4","RS1"]       = summary(Fullmod)$coefficients['GEV','Pr(>|t|)']
GFP_Pvals["MS4","RS1"]       = summary(Fullmod)$coefficients['GFP','Pr(>|t|)']
TP1_Pvals["MS4","RS1"]       = summary(Fullmod)$coefficients['TP1','Pr(>|t|)']
TP2_Pvals["MS4","RS1"]       = summary(Fullmod)$coefficients['TP2','Pr(>|t|)']
TP3_Pvals["MS4","RS1"]       = summary(Fullmod)$coefficients['TP3','Pr(>|t|)']
TP5_Pvals["MS4","RS1"]       = summary(Fullmod)$coefficients['TP5','Pr(>|t|)']

# Fullmod <- lm(DSInt ~  Coverage, data = MS4RS1)
# Intplot <- ggpredict(Fullmod,terms = c('Coverage'))
# plot(Intplot,add.data = TRUE)

# # RS2
# Fullmod <- lmer(DSInt ~ GFP  + Duration  + Occurance  + Coverage  + GEV  + TP2  + TP3  + TP1  + TP5  + (1|Subject), data = MS4RS2,REML = FALSE)
# anova(Fullmod)
# Duration_Pvals["MS4","RS2"]  = summary(Fullmod)$coefficients['Duration','Pr(>|t|)']
# Occurance_Pvals["MS4","RS2"] = summary(Fullmod)$coefficients['Occurance','Pr(>|t|)']
# Coverage_Pvals["MS4","RS2"]  = summary(Fullmod)$coefficients['Coverage','Pr(>|t|)']
# GEV_Pvals["MS4","RS2"]       = summary(Fullmod)$coefficients['GEV','Pr(>|t|)']
# TP1_Pvals["MS4","RS2"]       = summary(Fullmod)$coefficients['TP1','Pr(>|t|)']
# TP2_Pvals["MS4","RS2"]       = summary(Fullmod)$coefficients['TP2','Pr(>|t|)']
# TP3_Pvals["MS4","RS2"]       = summary(Fullmod)$coefficients['TP3','Pr(>|t|)']
# TP4_Pvals["MS4","RS2"]       = summary(Fullmod)$coefficients['TP4','Pr(>|t|)']
# TP5_Pvals["MS4","RS2"]       = summary(Fullmod)$coefficients['TP5','Pr(>|t|)']
# 
# # Diff
# Fullmod <- lmer(DSInt ~ GFP  + Duration  + Occurance  + Coverage  + GEV  + TP2  + TP3  + TP1  + TP5  + (1|Subject), data = MS4diff,REML = FALSE)
# anova(Fullmod)
# Duration_Pvals["MS4","Diff"]  = summary(Fullmod)$coefficients['Duration','Pr(>|t|)']
# Occurance_Pvals["MS4","Diff"] = summary(Fullmod)$coefficients['Occurance','Pr(>|t|)']
# Coverage_Pvals["MS4","Diff"]  = summary(Fullmod)$coefficients['Coverage','Pr(>|t|)']
# GEV_Pvals["MS4","Diff"]       = summary(Fullmod)$coefficients['GEV','Pr(>|t|)']
# TP1_Pvals["MS4","Diff"]       = summary(Fullmod)$coefficients['TP1','Pr(>|t|)']
# TP2_Pvals["MS4","Diff"]       = summary(Fullmod)$coefficients['TP2','Pr(>|t|)']
# TP3_Pvals["MS4","Diff"]       = summary(Fullmod)$coefficients['TP3','Pr(>|t|)']
# TP4_Pvals["MS4","Diff"]       = summary(Fullmod)$coefficients['TP4','Pr(>|t|)']
# TP5_Pvals["MS4","Diff"]       = summary(Fullmod)$coefficients['TP5','Pr(>|t|)']


## MS5
# RS1
# Fullmod <- lmer(DSInt ~ GFP  + Duration  + Occurance  + Coverage  + GEV  + TP2  + TP3  + TP4  + TP1  + (1|Subject), data = MS5RS1,REML = FALSE)
Fullmod <- lm(DSInt ~ GFP  + Duration  + Occurance  + Coverage  + GEV  + TP2  + TP3  + TP4  + TP1, data = MS5RS1)
summary(Fullmod)
Duration_Pvals["MS5","RS1"]  = summary(Fullmod)$coefficients['Duration','Pr(>|t|)']
Occurance_Pvals["MS5","RS1"] = summary(Fullmod)$coefficients['Occurance','Pr(>|t|)']
Coverage_Pvals["MS5","RS1"]  = summary(Fullmod)$coefficients['Coverage','Pr(>|t|)']
GEV_Pvals["MS5","RS1"]       = summary(Fullmod)$coefficients['GEV','Pr(>|t|)']
GFP_Pvals["MS5","RS1"]       = summary(Fullmod)$coefficients['GFP','Pr(>|t|)']
TP1_Pvals["MS5","RS1"]       = summary(Fullmod)$coefficients['TP1','Pr(>|t|)']
TP2_Pvals["MS5","RS1"]       = summary(Fullmod)$coefficients['TP2','Pr(>|t|)']
TP3_Pvals["MS5","RS1"]       = summary(Fullmod)$coefficients['TP3','Pr(>|t|)']
TP4_Pvals["MS5","RS1"]       = summary(Fullmod)$coefficients['TP4','Pr(>|t|)']

# # RS2
# Fullmod <- lmer(DSInt ~ GFP  + Duration  + Occurance  + Coverage  + GEV  + TP2  + TP3  + TP4  + TP1  + (1|Subject), data = MS5RS2,REML = FALSE)
# anova(Fullmod)
# Duration_Pvals["MS5","RS2"]  = summary(Fullmod)$coefficients['Duration','Pr(>|t|)']
# Occurance_Pvals["MS5","RS2"] = summary(Fullmod)$coefficients['Occurance','Pr(>|t|)']
# Coverage_Pvals["MS5","RS2"]  = summary(Fullmod)$coefficients['Coverage','Pr(>|t|)']
# GEV_Pvals["MS5","RS2"]       = summary(Fullmod)$coefficients['GEV','Pr(>|t|)']
# TP1_Pvals["MS5","RS2"]       = summary(Fullmod)$coefficients['TP1','Pr(>|t|)']
# TP2_Pvals["MS5","RS2"]       = summary(Fullmod)$coefficients['TP2','Pr(>|t|)']
# TP3_Pvals["MS5","RS2"]       = summary(Fullmod)$coefficients['TP3','Pr(>|t|)']
# TP4_Pvals["MS5","RS2"]       = summary(Fullmod)$coefficients['TP4','Pr(>|t|)']
# TP5_Pvals["MS5","RS2"]       = summary(Fullmod)$coefficients['TP5','Pr(>|t|)']
# 
# # Diff
# Fullmod <- lmer(DSInt ~ GFP  + Duration  + Occurance  + Coverage  + GEV  + TP2  + TP3  + TP4  + TP1  + (1|Subject), data = MS5diff,REML = FALSE)
# anova(Fullmod)
# Duration_Pvals["MS5","Diff"]  = summary(Fullmod)$coefficients['Duration','Pr(>|t|)']
# Occurance_Pvals["MS5","Diff"] = summary(Fullmod)$coefficients['Occurance','Pr(>|t|)']
# Coverage_Pvals["MS5","Diff"]  = summary(Fullmod)$coefficients['Coverage','Pr(>|t|)']
# GEV_Pvals["MS5","Diff"]       = summary(Fullmod)$coefficients['GEV','Pr(>|t|)']
# TP1_Pvals["MS5","Diff"]       = summary(Fullmod)$coefficients['TP1','Pr(>|t|)']
# TP2_Pvals["MS5","Diff"]       = summary(Fullmod)$coefficients['TP2','Pr(>|t|)']
# TP3_Pvals["MS5","Diff"]       = summary(Fullmod)$coefficients['TP3','Pr(>|t|)']
# TP4_Pvals["MS5","Diff"]       = summary(Fullmod)$coefficients['TP4','Pr(>|t|)']
# TP5_Pvals["MS5","Diff"]       = summary(Fullmod)$coefficients['TP5','Pr(>|t|)']

# Check for significance after FDR correction
Pvalues <- unlist(TP4_Pvals)
non_na_indices <- which(!is.na(Pvalues))
All_adj_p_values <- p.adjust(Pvalues[non_na_indices], method = "fdr")
Pvalues[non_na_indices] <- All_adj_p_values
Pvalues[Pvalues<.05]



# Fit predictors separately to avoid the confusion of the other predictors
predictors <- c("GFP", "Duration", "Occurance", "Coverage", "GEV", "TP2", "TP3", "TP4","TP5")
MS1_DSInt <- Fit_Models_DS(predictors, "DSInt", MS1RS1,'MS1')

predictors <- c("GFP", "Duration", "Occurance", "Coverage", "GEV", "TP1", "TP3", "TP4","TP5")
MS2_DSInt <- Fit_Models_DS(predictors, "DSInt", MS2RS1,'MS2')

predictors <- c("GFP", "Duration", "Occurance", "Coverage", "GEV", "TP2", "TP1", "TP4","TP5")
MS3_DSInt <- Fit_Models_DS(predictors, "DSInt", MS3RS1,'MS3')

predictors <- c("GFP", "Duration", "Occurance", "Coverage", "GEV", "TP2", "TP3", "TP1","TP5")
MS4_DSInt <- Fit_Models_DS(predictors, "DSInt", MS4RS1,'MS4')

predictors <- c("GFP", "Duration", "Occurance", "Coverage", "GEV", "TP1", "TP2", "TP3", "TP4")
MS5_DSInt <- Fit_Models_DS(predictors, "DSInt", MS5RS1,'MS5')

## Multiple comparisons correction for Area Length
Pvalues <- unlist(c(MS1_DSInt,MS2_DSInt,MS3_DSInt,MS4_DSInt,MS5_DSInt))
non_na_indices <- which(!is.na(Pvalues))
All_adj_p_values <- p.adjust(Pvalues[non_na_indices], method = "fdr")
Pvalues[non_na_indices] <- All_adj_p_values
Pvalues[Pvalues<.05]


### Unpleasantness ----
# Initialize variables

# Create column and row names
col_names <- c("RS1","RS2","Diff")
row_names <- c("MS1", "MS2", "MS3", "MS4", "MS5")
# Create empty data frames
Duration_Pvals <- data.frame(matrix(ncol = length(col_names), nrow = length(row_names)))
colnames(Duration_Pvals) <- col_names
rownames(Duration_Pvals) <- row_names
Occurance_Pvals = Duration_Pvals
Coverage_Pvals = Duration_Pvals
GEV_Pvals = Duration_Pvals
GFP_Pvals = Duration_Pvals
TP1_Pvals = Duration_Pvals
TP2_Pvals = Duration_Pvals
TP3_Pvals = Duration_Pvals
TP4_Pvals = Duration_Pvals
TP5_Pvals = Duration_Pvals

## MS1
# RS1
Fullmod <- lm(DSUnp ~ GFP  + Duration  + Occurance  + Coverage  + GEV  + TP2  + TP3  + TP4  + TP5, data = MS1RS1)
summary(Fullmod)
Duration_Pvals["MS1","RS1"]  = summary(Fullmod)$coefficients['Duration','Pr(>|t|)']
Occurance_Pvals["MS1","RS1"] = summary(Fullmod)$coefficients['Occurance','Pr(>|t|)']
Coverage_Pvals["MS1","RS1"]  = summary(Fullmod)$coefficients['Coverage','Pr(>|t|)']
GEV_Pvals["MS1","RS1"]       = summary(Fullmod)$coefficients['GEV','Pr(>|t|)']
GFP_Pvals["MS1","RS1"]       = summary(Fullmod)$coefficients['GFP','Pr(>|t|)']
TP2_Pvals["MS1","RS1"]       = summary(Fullmod)$coefficients['TP2','Pr(>|t|)']
TP3_Pvals["MS1","RS1"]       = summary(Fullmod)$coefficients['TP3','Pr(>|t|)']
TP4_Pvals["MS1","RS1"]       = summary(Fullmod)$coefficients['TP4','Pr(>|t|)']
TP5_Pvals["MS1","RS1"]       = summary(Fullmod)$coefficients['TP5','Pr(>|t|)']

# # RS2
# Fullmod <- lmer(DSUnp ~ GFP  + Duration  + Occurance  + Coverage  + GEV  + TP2  + TP3  + TP4  + TP5  + (1|Subject), data = MS1RS2,REML = FALSE)
# anova(Fullmod)
# Duration_Pvals["MS1","RS2"]  = summary(Fullmod)$coefficients['Duration','Pr(>|t|)']
# Occurance_Pvals["MS1","RS2"] = summary(Fullmod)$coefficients['Occurance','Pr(>|t|)']
# Coverage_Pvals["MS1","RS2"]  = summary(Fullmod)$coefficients['Coverage','Pr(>|t|)']
# GEV_Pvals["MS1","RS2"]       = summary(Fullmod)$coefficients['GEV','Pr(>|t|)']
# TP1_Pvals["MS1","RS2"]       = summary(Fullmod)$coefficients['TP1','Pr(>|t|)']
# TP2_Pvals["MS1","RS2"]       = summary(Fullmod)$coefficients['TP2','Pr(>|t|)']
# TP3_Pvals["MS1","RS2"]       = summary(Fullmod)$coefficients['TP3','Pr(>|t|)']
# TP4_Pvals["MS1","RS2"]       = summary(Fullmod)$coefficients['TP4','Pr(>|t|)']
# TP5_Pvals["MS1","RS2"]       = summary(Fullmod)$coefficients['TP5','Pr(>|t|)']
# 
# # Diff
# Fullmod <- lmer(DSUnp ~ GFP  + Duration  + Occurance  + Coverage  + GEV  + TP2  + TP3  + TP4  + TP5  + (1|Subject), data = MS1diff,REML = FALSE)
# anova(Fullmod)
# Duration_Pvals["MS1","Diff"]  = summary(Fullmod)$coefficients['Duration','Pr(>|t|)']
# Occurance_Pvals["MS1","Diff"] = summary(Fullmod)$coefficients['Occurance','Pr(>|t|)']
# Coverage_Pvals["MS1","Diff"]  = summary(Fullmod)$coefficients['Coverage','Pr(>|t|)']
# GEV_Pvals["MS1","Diff"]       = summary(Fullmod)$coefficients['GEV','Pr(>|t|)']
# TP1_Pvals["MS1","Diff"]       = summary(Fullmod)$coefficients['TP1','Pr(>|t|)']
# TP2_Pvals["MS1","Diff"]       = summary(Fullmod)$coefficients['TP2','Pr(>|t|)']
# TP3_Pvals["MS1","Diff"]       = summary(Fullmod)$coefficients['TP3','Pr(>|t|)']
# TP4_Pvals["MS1","Diff"]       = summary(Fullmod)$coefficients['TP4','Pr(>|t|)']
# TP5_Pvals["MS1","Diff"]       = summary(Fullmod)$coefficients['TP5','Pr(>|t|)']


## MS2
# RS1
Fullmod <- lm(DSUnp ~ GFP  + Duration  + Occurance  + Coverage  + GEV  + TP1  + TP3  + TP4  + TP5, data = MS2RS1)
summary(Fullmod)
Duration_Pvals["MS2","RS1"]  = summary(Fullmod)$coefficients['Duration','Pr(>|t|)']
Occurance_Pvals["MS2","RS1"] = summary(Fullmod)$coefficients['Occurance','Pr(>|t|)']
Coverage_Pvals["MS2","RS1"]  = summary(Fullmod)$coefficients['Coverage','Pr(>|t|)']
GEV_Pvals["MS2","RS1"]       = summary(Fullmod)$coefficients['GEV','Pr(>|t|)']
GFP_Pvals["MS2","RS1"]       = summary(Fullmod)$coefficients['GFP','Pr(>|t|)']
TP1_Pvals["MS2","RS1"]       = summary(Fullmod)$coefficients['TP1','Pr(>|t|)']
TP3_Pvals["MS2","RS1"]       = summary(Fullmod)$coefficients['TP3','Pr(>|t|)']
TP4_Pvals["MS2","RS1"]       = summary(Fullmod)$coefficients['TP4','Pr(>|t|)']
TP5_Pvals["MS2","RS1"]       = summary(Fullmod)$coefficients['TP5','Pr(>|t|)']

# # RS2
# Fullmod <- lmer(DSUnp ~ GFP  + Duration  + Occurance  + Coverage  + GEV  + TP1  + TP3  + TP4  + TP5  + (1|Subject), data = MS2RS2,REML = FALSE)
# anova(Fullmod)
# Duration_Pvals["MS2","RS2"]  = summary(Fullmod)$coefficients['Duration','Pr(>|t|)']
# Occurance_Pvals["MS2","RS2"] = summary(Fullmod)$coefficients['Occurance','Pr(>|t|)']
# Coverage_Pvals["MS2","RS2"]  = summary(Fullmod)$coefficients['Coverage','Pr(>|t|)']
# GEV_Pvals["MS2","RS2"]       = summary(Fullmod)$coefficients['GEV','Pr(>|t|)']
# TP1_Pvals["MS2","RS2"]       = summary(Fullmod)$coefficients['TP1','Pr(>|t|)']
# TP2_Pvals["MS2","RS2"]       = summary(Fullmod)$coefficients['TP2','Pr(>|t|)']
# TP3_Pvals["MS2","RS2"]       = summary(Fullmod)$coefficients['TP3','Pr(>|t|)']
# TP4_Pvals["MS2","RS2"]       = summary(Fullmod)$coefficients['TP4','Pr(>|t|)']
# TP5_Pvals["MS2","RS2"]       = summary(Fullmod)$coefficients['TP5','Pr(>|t|)']
# 
# # Diff
# Fullmod <- lmer(DSUnp ~ GFP  + Duration  + Occurance  + Coverage  + GEV  + TP1  + TP3  + TP4  + TP5  + (1|Subject), data = MS2diff,REML = FALSE)
# anova(Fullmod)
# Duration_Pvals["MS2","Diff"]  = summary(Fullmod)$coefficients['Duration','Pr(>|t|)']
# Occurance_Pvals["MS2","Diff"] = summary(Fullmod)$coefficients['Occurance','Pr(>|t|)']
# Coverage_Pvals["MS2","Diff"]  = summary(Fullmod)$coefficients['Coverage','Pr(>|t|)']
# GEV_Pvals["MS2","Diff"]       = summary(Fullmod)$coefficients['GEV','Pr(>|t|)']
# TP1_Pvals["MS2","Diff"]       = summary(Fullmod)$coefficients['TP1','Pr(>|t|)']
# TP2_Pvals["MS2","Diff"]       = summary(Fullmod)$coefficients['TP2','Pr(>|t|)']
# TP3_Pvals["MS2","Diff"]       = summary(Fullmod)$coefficients['TP3','Pr(>|t|)']
# TP4_Pvals["MS2","Diff"]       = summary(Fullmod)$coefficients['TP4','Pr(>|t|)']
# TP5_Pvals["MS2","Diff"]       = summary(Fullmod)$coefficients['TP5','Pr(>|t|)']


## MS3
# RS1
Fullmod <- lm(DSUnp ~ GFP  + Duration  + Occurance  + Coverage  + GEV  + TP2  + TP1  + TP4  + TP5, data = MS3RS1)
summary(Fullmod)
Duration_Pvals["MS3","RS1"]  = summary(Fullmod)$coefficients['Duration','Pr(>|t|)']
Occurance_Pvals["MS3","RS1"] = summary(Fullmod)$coefficients['Occurance','Pr(>|t|)']
Coverage_Pvals["MS3","RS1"]  = summary(Fullmod)$coefficients['Coverage','Pr(>|t|)']
GEV_Pvals["MS3","RS1"]       = summary(Fullmod)$coefficients['GEV','Pr(>|t|)']
GFP_Pvals["MS3","RS1"]       = summary(Fullmod)$coefficients['GFP','Pr(>|t|)']
TP1_Pvals["MS3","RS1"]       = summary(Fullmod)$coefficients['TP1','Pr(>|t|)']
TP2_Pvals["MS3","RS1"]       = summary(Fullmod)$coefficients['TP2','Pr(>|t|)']
TP4_Pvals["MS3","RS1"]       = summary(Fullmod)$coefficients['TP4','Pr(>|t|)']
TP5_Pvals["MS3","RS1"]       = summary(Fullmod)$coefficients['TP5','Pr(>|t|)']

Fullmod <- lm(DSUnp ~  TP1, data = MS3RS1)
summary(Fullmod)
Intplot <- ggpredict(Fullmod,terms = c('TP1'))
plot(Intplot,add.data = TRUE)

# # RS2
# Fullmod <- lmer(DSUnp ~ GFP  + Duration  + Occurance  + Coverage  + GEV  + TP2  + TP1  + TP4  + TP5  + (1|Subject), data = MS3RS2,REML = FALSE)
# anova(Fullmod)
# Duration_Pvals["MS3","RS2"]  = summary(Fullmod)$coefficients['Duration','Pr(>|t|)']
# Occurance_Pvals["MS3","RS2"] = summary(Fullmod)$coefficients['Occurance','Pr(>|t|)']
# Coverage_Pvals["MS3","RS2"]  = summary(Fullmod)$coefficients['Coverage','Pr(>|t|)']
# GEV_Pvals["MS3","RS2"]       = summary(Fullmod)$coefficients['GEV','Pr(>|t|)']
# TP1_Pvals["MS3","RS2"]       = summary(Fullmod)$coefficients['TP1','Pr(>|t|)']
# TP2_Pvals["MS3","RS2"]       = summary(Fullmod)$coefficients['TP2','Pr(>|t|)']
# TP3_Pvals["MS3","RS2"]       = summary(Fullmod)$coefficients['TP3','Pr(>|t|)']
# TP4_Pvals["MS3","RS2"]       = summary(Fullmod)$coefficients['TP4','Pr(>|t|)']
# TP5_Pvals["MS3","RS2"]       = summary(Fullmod)$coefficients['TP5','Pr(>|t|)']
# 
# # Diff
# Fullmod <- lmer(DSUnp ~ GFP  + Duration  + Occurance  + Coverage  + GEV  + TP2  + TP1  + TP4  + TP5  + (1|Subject), data = MS3diff,REML = FALSE)
# anova(Fullmod)
# Duration_Pvals["MS3","Diff"]  = summary(Fullmod)$coefficients['Duration','Pr(>|t|)']
# Occurance_Pvals["MS3","Diff"] = summary(Fullmod)$coefficients['Occurance','Pr(>|t|)']
# Coverage_Pvals["MS3","Diff"]  = summary(Fullmod)$coefficients['Coverage','Pr(>|t|)']
# GEV_Pvals["MS3","Diff"]       = summary(Fullmod)$coefficients['GEV','Pr(>|t|)']
# TP1_Pvals["MS3","Diff"]       = summary(Fullmod)$coefficients['TP1','Pr(>|t|)']
# TP2_Pvals["MS3","Diff"]       = summary(Fullmod)$coefficients['TP2','Pr(>|t|)']
# TP3_Pvals["MS3","Diff"]       = summary(Fullmod)$coefficients['TP3','Pr(>|t|)']
# TP4_Pvals["MS3","Diff"]       = summary(Fullmod)$coefficients['TP4','Pr(>|t|)']
# TP5_Pvals["MS3","Diff"]       = summary(Fullmod)$coefficients['TP5','Pr(>|t|)']



## MS4
# RS1
Fullmod <- lm(DSUnp ~ GFP  + Duration  + Occurance  + Coverage  + GEV  + TP2  + TP3  + TP1  + TP5, data = MS4RS1)
summary(Fullmod)
Duration_Pvals["MS4","RS1"]  = summary(Fullmod)$coefficients['Duration','Pr(>|t|)']
Occurance_Pvals["MS4","RS1"] = summary(Fullmod)$coefficients['Occurance','Pr(>|t|)']
Coverage_Pvals["MS4","RS1"]  = summary(Fullmod)$coefficients['Coverage','Pr(>|t|)']
GEV_Pvals["MS4","RS1"]       = summary(Fullmod)$coefficients['GEV','Pr(>|t|)']
GFP_Pvals["MS4","RS1"]       = summary(Fullmod)$coefficients['GFP','Pr(>|t|)']
TP1_Pvals["MS4","RS1"]       = summary(Fullmod)$coefficients['TP1','Pr(>|t|)']
TP2_Pvals["MS4","RS1"]       = summary(Fullmod)$coefficients['TP2','Pr(>|t|)']
TP3_Pvals["MS4","RS1"]       = summary(Fullmod)$coefficients['TP3','Pr(>|t|)']
TP5_Pvals["MS4","RS1"]       = summary(Fullmod)$coefficients['TP5','Pr(>|t|)']

# # RS2
# Fullmod <- lmer(DSUnp ~ GFP  + Duration  + Occurance  + Coverage  + GEV  + TP2  + TP3  + TP1  + TP5  + (1|Subject), data = MS4RS2,REML = FALSE)
# anova(Fullmod)
# Duration_Pvals["MS4","RS2"]  = summary(Fullmod)$coefficients['Duration','Pr(>|t|)']
# Occurance_Pvals["MS4","RS2"] = summary(Fullmod)$coefficients['Occurance','Pr(>|t|)']
# Coverage_Pvals["MS4","RS2"]  = summary(Fullmod)$coefficients['Coverage','Pr(>|t|)']
# GEV_Pvals["MS4","RS2"]       = summary(Fullmod)$coefficients['GEV','Pr(>|t|)']
# TP1_Pvals["MS4","RS2"]       = summary(Fullmod)$coefficients['TP1','Pr(>|t|)']
# TP2_Pvals["MS4","RS2"]       = summary(Fullmod)$coefficients['TP2','Pr(>|t|)']
# TP3_Pvals["MS4","RS2"]       = summary(Fullmod)$coefficients['TP3','Pr(>|t|)']
# TP4_Pvals["MS4","RS2"]       = summary(Fullmod)$coefficients['TP4','Pr(>|t|)']
# TP5_Pvals["MS4","RS2"]       = summary(Fullmod)$coefficients['TP5','Pr(>|t|)']
# 
# # Diff
# Fullmod <- lmer(DSUnp ~ GFP  + Duration  + Occurance  + Coverage  + GEV  + TP2  + TP3  + TP1  + TP5  + (1|Subject), data = MS4diff,REML = FALSE)
# anova(Fullmod)
# Duration_Pvals["MS4","Diff"]  = summary(Fullmod)$coefficients['Duration','Pr(>|t|)']
# Occurance_Pvals["MS4","Diff"] = summary(Fullmod)$coefficients['Occurance','Pr(>|t|)']
# Coverage_Pvals["MS4","Diff"]  = summary(Fullmod)$coefficients['Coverage','Pr(>|t|)']
# GEV_Pvals["MS4","Diff"]       = summary(Fullmod)$coefficients['GEV','Pr(>|t|)']
# TP1_Pvals["MS4","Diff"]       = summary(Fullmod)$coefficients['TP1','Pr(>|t|)']
# TP2_Pvals["MS4","Diff"]       = summary(Fullmod)$coefficients['TP2','Pr(>|t|)']
# TP3_Pvals["MS4","Diff"]       = summary(Fullmod)$coefficients['TP3','Pr(>|t|)']
# TP4_Pvals["MS4","Diff"]       = summary(Fullmod)$coefficients['TP4','Pr(>|t|)']
# TP5_Pvals["MS4","Diff"]       = summary(Fullmod)$coefficients['TP5','Pr(>|t|)']

## MS5
# RS1
Fullmod <- lm(DSUnp ~ GFP  + Duration  + Occurance  + Coverage  + GEV  + TP2  + TP3  + TP4  + TP1, data = MS5RS1)
summary(Fullmod)
Duration_Pvals["MS5","RS1"]  = summary(Fullmod)$coefficients['Duration','Pr(>|t|)']
Occurance_Pvals["MS5","RS1"] = summary(Fullmod)$coefficients['Occurance','Pr(>|t|)']
Coverage_Pvals["MS5","RS1"]  = summary(Fullmod)$coefficients['Coverage','Pr(>|t|)']
GEV_Pvals["MS5","RS1"]       = summary(Fullmod)$coefficients['GEV','Pr(>|t|)']
GFP_Pvals["MS5","RS1"]       = summary(Fullmod)$coefficients['GFP','Pr(>|t|)']
TP1_Pvals["MS5","RS1"]       = summary(Fullmod)$coefficients['TP1','Pr(>|t|)']
TP2_Pvals["MS5","RS1"]       = summary(Fullmod)$coefficients['TP2','Pr(>|t|)']
TP3_Pvals["MS5","RS1"]       = summary(Fullmod)$coefficients['TP3','Pr(>|t|)']
TP4_Pvals["MS5","RS1"]       = summary(Fullmod)$coefficients['TP4','Pr(>|t|)']

# # RS2
# Fullmod <- lmer(DSUnp ~ GFP  + Duration  + Occurance  + Coverage  + GEV  + TP2  + TP3  + TP4  + TP1  + (1|Subject), data = MS5RS2,REML = FALSE)
# anova(Fullmod)
# Duration_Pvals["MS5","RS2"]  = summary(Fullmod)$coefficients['Duration','Pr(>|t|)']
# Occurance_Pvals["MS5","RS2"] = summary(Fullmod)$coefficients['Occurance','Pr(>|t|)']
# Coverage_Pvals["MS5","RS2"]  = summary(Fullmod)$coefficients['Coverage','Pr(>|t|)']
# GEV_Pvals["MS5","RS2"]       = summary(Fullmod)$coefficients['GEV','Pr(>|t|)']
# TP1_Pvals["MS5","RS2"]       = summary(Fullmod)$coefficients['TP1','Pr(>|t|)']
# TP2_Pvals["MS5","RS2"]       = summary(Fullmod)$coefficients['TP2','Pr(>|t|)']
# TP3_Pvals["MS5","RS2"]       = summary(Fullmod)$coefficients['TP3','Pr(>|t|)']
# TP4_Pvals["MS5","RS2"]       = summary(Fullmod)$coefficients['TP4','Pr(>|t|)']
# TP5_Pvals["MS5","RS2"]       = summary(Fullmod)$coefficients['TP5','Pr(>|t|)']
# 
# # Diff
# Fullmod <- lmer(DSUnp ~ GFP  + Duration  + Occurance  + Coverage  + GEV  + TP2  + TP3  + TP4  + TP1  + (1|Subject), data = MS5diff,REML = FALSE)
# anova(Fullmod)
# Duration_Pvals["MS5","Diff"]  = summary(Fullmod)$coefficients['Duration','Pr(>|t|)']
# Occurance_Pvals["MS5","Diff"] = summary(Fullmod)$coefficients['Occurance','Pr(>|t|)']
# Coverage_Pvals["MS5","Diff"]  = summary(Fullmod)$coefficients['Coverage','Pr(>|t|)']
# GEV_Pvals["MS5","Diff"]       = summary(Fullmod)$coefficients['GEV','Pr(>|t|)']
# TP1_Pvals["MS5","Diff"]       = summary(Fullmod)$coefficients['TP1','Pr(>|t|)']
# TP2_Pvals["MS5","Diff"]       = summary(Fullmod)$coefficients['TP2','Pr(>|t|)']
# TP3_Pvals["MS5","Diff"]       = summary(Fullmod)$coefficients['TP3','Pr(>|t|)']
# TP4_Pvals["MS5","Diff"]       = summary(Fullmod)$coefficients['TP4','Pr(>|t|)']
# TP5_Pvals["MS5","Diff"]       = summary(Fullmod)$coefficients['TP5','Pr(>|t|)']

# Check for significance after FDR correction
Pvalues <- unlist(TP5_Pvals)
non_na_indices <- which(!is.na(Pvalues))
All_adj_p_values <- p.adjust(Pvalues[non_na_indices], method = "fdr")
Pvalues[non_na_indices] <- All_adj_p_values
Pvalues[Pvalues<.05]

# Fit predictors separately to avoid the confusion of the other predictors
predictors <- c("GFP", "Duration", "Occurance", "Coverage", "GEV", "TP2", "TP3", "TP4","TP5")
MS1_DSUnp <- Fit_Models_DS(predictors, "DSUnp", MS1RS1,'MS1')

predictors <- c("GFP", "Duration", "Occurance", "Coverage", "GEV", "TP1", "TP3", "TP4","TP5")
MS2_DSUnp <- Fit_Models_DS(predictors, "DSUnp", MS2RS1,'MS2')

predictors <- c("GFP", "Duration", "Occurance", "Coverage", "GEV", "TP2", "TP1", "TP4","TP5")
MS3_DSUnp <- Fit_Models_DS(predictors, "DSUnp", MS3RS1,'MS3')

predictors <- c("GFP", "Duration", "Occurance", "Coverage", "GEV", "TP2", "TP3", "TP1","TP5")
MS4_DSUnp <- Fit_Models_DS(predictors, "DSUnp", MS4RS1,'MS4')

predictors <- c("GFP", "Duration", "Occurance", "Coverage", "GEV", "TP1", "TP2", "TP3", "TP4")
MS5_DSUnp <- Fit_Models_DS(predictors, "DSUnp", MS5RS1,'MS5')

## Multiple comparisons correction for Area Length
Pvalues <- unlist(c(MS1_DSUnp,MS2_DSUnp,MS3_DSUnp,MS4_DSUnp,MS5_DSUnp))
non_na_indices <- which(!is.na(Pvalues))
All_adj_p_values <- p.adjust(Pvalues[non_na_indices], method = "fdr")
Pvalues[non_na_indices] <- All_adj_p_values
Pvalues[Pvalues<.05]



## Averaged but not subtracted----
## Microstates x PPRatings ----
MicrostatsxPP <- read_excel("Path_to_Microstate_PPRatings_data/MicrostatsxPPratings.xlsx")
View(MicrostatsxPP)

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


HFS  = MicrostatsxPP[MicrostatsxPP$Arm == "MFS",]
Cont = MicrostatsxPP[MicrostatsxPP$Arm == "Cont",]

library(pracma)
ind = linspace(1,1991,200)
temp = HFS[ind,]

temp = MicrostatsxPP[MicrostatsxPP$Trial == 1,]

tempHFS = temp[temp$Arm == 'MFS',]
tempCont = temp[temp$Arm == 'Cont',]

index = 1
# Initialize a for loop
for (i in ind) {
  
  HFSIntT1 = mean(HFS$IntensityRatings[i:(i+9)])
  HFSUnpT1 = mean(HFS$UnpleasantnessRatings[i:(i+9)])
  ContIntT1 = mean(Cont$IntensityRatings[i:(i+9)])
  ContUnpT1 = mean(Cont$UnpleasantnessRatings[i:(i+9)])
  
  HFSIntT0 = mean(HFS$IntensityT0[i:(i+9)])
  HFSUnpT0 = mean(HFS$UnpleasantnessT0[i:(i+9)])
  ContIntT0 = mean(Cont$IntensityT0[i:(i+9)])
  ContUnpT0 = mean(Cont$UnpleasantnessT0[i:(i+9)])
  
  tempHFS$IntensityRatings[index] = HFSIntT1
  tempHFS$IntensityT0[index] = HFSIntT0
  tempCont$IntensityRatings[index] = ContIntT1
  tempCont$IntensityT0[index] = ContIntT0
  
  tempHFS$UnpleasantnessRatings[index] = HFSUnpT1
  tempHFS$UnpleasantnessT0[index] = HFSUnpT0
  tempCont$UnpleasantnessRatings[index] = ContUnpT1
  tempCont$UnpleasantnessT0[index] = ContUnpT0
  
  index = index + 1
  
}

temp <- rbind(tempHFS, tempCont)

MicrostatsxPP = temp

# Standardize predictors
MicrostatsxPP$GFP <- as.numeric(scale(MicrostatsxPP$GFP))
MicrostatsxPP$Duration <- as.numeric(scale(MicrostatsxPP$Duration))
MicrostatsxPP$Occurance <- as.numeric(scale(MicrostatsxPP$Occurance))

# 
# # Create the double subtracted data
# MicrostatsxPP$Intsubtr = MicrostatsxPP$IntensityRatings - MicrostatsxPP$IntensityT0
# MicrostatsxPP$Unpsubtr = MicrostatsxPP$UnpleasantnessRatings - MicrostatsxPP$UnpleasantnessT0
# 
# HFS$DSint = HFS$Intsubtr - Cont$Intsubtr
# HFS$DSUnp = HFS$Unpsubtr - Cont$Unpsubtr
# 
# MicrostatsxPP = HFS

MS1RS1 = MicrostatsxPP[MicrostatsxPP$MS == 'MS1' & MicrostatsxPP$RS == 'RS1',]
MS2RS1 = MicrostatsxPP[MicrostatsxPP$MS == 'MS2' & MicrostatsxPP$RS == 'RS1',]
MS3RS1 = MicrostatsxPP[MicrostatsxPP$MS == 'MS3' & MicrostatsxPP$RS == 'RS1',]
MS4RS1 = MicrostatsxPP[MicrostatsxPP$MS == 'MS4' & MicrostatsxPP$RS == 'RS1',]
MS5RS1 = MicrostatsxPP[MicrostatsxPP$MS == 'MS5' & MicrostatsxPP$RS == 'RS1',]

# MS1RS2 = MicrostatsxPP[MicrostatsxPP$MS == 'MS1' & MicrostatsxPP$RS == 'RS2',]
# MS2RS2 = MicrostatsxPP[MicrostatsxPP$MS == 'MS2' & MicrostatsxPP$RS == 'RS2',]
# MS3RS2 = MicrostatsxPP[MicrostatsxPP$MS == 'MS3' & MicrostatsxPP$RS == 'RS2',]
# MS4RS2 = MicrostatsxPP[MicrostatsxPP$MS == 'MS4' & MicrostatsxPP$RS == 'RS2',]
# MS5RS2 = MicrostatsxPP[MicrostatsxPP$MS == 'MS5' & MicrostatsxPP$RS == 'RS2',]
# 
# MS1diff = MS1RS1
# MS1diff[,c(6:16)] = MS1RS1[,c(6:16)] - MS1RS2[,c(6:16)]
# MS2diff = MS2RS1
# MS2diff[,c(6:16)] = MS2RS1[,c(6:16)] - MS2RS2[,c(6:16)]
# MS3diff = MS3RS1
# MS3diff[,c(6:16)] = MS3RS1[,c(6:16)] - MS3RS2[,c(6:16)]
# MS4diff = MS4RS1
# MS4diff[,c(6:16)] = MS4RS1[,c(6:16)] - MS4RS2[,c(6:16)]
# MS5diff = MS5RS1
# MS5diff[,c(6:16)] = MS5RS1[,c(6:16)] - MS5RS2[,c(6:16)]

### Intensity ----
# Initialize variables

# Create column and row names
col_names <- c("RS1","RS2","Diff")
row_names <- c("MS1", "MS2", "MS3", "MS4", "MS5")
# Create empty data frames
Duration_Pvals <- data.frame(matrix(ncol = length(col_names), nrow = length(row_names)))
colnames(Duration_Pvals) <- col_names
rownames(Duration_Pvals) <- row_names
Occurance_Pvals = Duration_Pvals
Coverage_Pvals = Duration_Pvals
GEV_Pvals = Duration_Pvals
GFP_Pvals = Duration_Pvals
TP1_Pvals = Duration_Pvals
TP2_Pvals = Duration_Pvals
TP3_Pvals = Duration_Pvals
TP4_Pvals = Duration_Pvals
TP5_Pvals = Duration_Pvals


## MS1
# RS1
# Fullmod <- lmer(DSInt ~ GFP  + Duration  + Occurance  + Coverage  + GEV  + TP2  + TP3  + TP4  + TP5  + (1|Subject), data = MS1RS1,REML = FALSE)
Fullmod <- lm(DSInt ~ GFP  + Duration  + Occurance  + Coverage  + GEV  + TP2  + TP3  + TP4  + TP5, data = MS1RS1)
summary(Fullmod)
Duration_Pvals["MS1","RS1"]  = summary(Fullmod)$coefficients['Duration','Pr(>|t|)']
Occurance_Pvals["MS1","RS1"] = summary(Fullmod)$coefficients['Occurance','Pr(>|t|)']
Coverage_Pvals["MS1","RS1"]  = summary(Fullmod)$coefficients['Coverage','Pr(>|t|)']
GEV_Pvals["MS1","RS1"]       = summary(Fullmod)$coefficients['GEV','Pr(>|t|)']
GFP_Pvals["MS1","RS1"]       = summary(Fullmod)$coefficients['GFP','Pr(>|t|)']
TP2_Pvals["MS1","RS1"]       = summary(Fullmod)$coefficients['TP2','Pr(>|t|)']
TP3_Pvals["MS1","RS1"]       = summary(Fullmod)$coefficients['TP3','Pr(>|t|)']
TP4_Pvals["MS1","RS1"]       = summary(Fullmod)$coefficients['TP4','Pr(>|t|)']
TP5_Pvals["MS1","RS1"]       = summary(Fullmod)$coefficients['TP5','Pr(>|t|)']


# # RS2
# Fullmod <- lmer(DSInt ~ GFP  + Duration  + Occurance  + Coverage  + GEV  + TP2  + TP3  + TP4  + TP5  + (1|Subject), data = MS1RS2,REML = FALSE)
# anova(Fullmod)
# Duration_Pvals["MS1","RS2"]  = summary(Fullmod)$coefficients['Duration','Pr(>|t|)']
# Occurance_Pvals["MS1","RS2"] = summary(Fullmod)$coefficients['Occurance','Pr(>|t|)']
# Coverage_Pvals["MS1","RS2"]  = summary(Fullmod)$coefficients['Coverage','Pr(>|t|)']
# GEV_Pvals["MS1","RS2"]       = summary(Fullmod)$coefficients['GEV','Pr(>|t|)']
# TP1_Pvals["MS1","RS2"]       = summary(Fullmod)$coefficients['TP1','Pr(>|t|)']
# TP2_Pvals["MS1","RS2"]       = summary(Fullmod)$coefficients['TP2','Pr(>|t|)']
# TP3_Pvals["MS1","RS2"]       = summary(Fullmod)$coefficients['TP3','Pr(>|t|)']
# TP4_Pvals["MS1","RS2"]       = summary(Fullmod)$coefficients['TP4','Pr(>|t|)']
# TP5_Pvals["MS1","RS2"]       = summary(Fullmod)$coefficients['TP5','Pr(>|t|)']
# # Diff
# Fullmod <- lmer(DSInt ~ GFP  + Duration  + Occurance  + Coverage  + GEV  + TP2  + TP3  + TP4  + TP5  + (1|Subject), data = MS1diff,REML = FALSE)
# anova(Fullmod)
# Duration_Pvals["MS1","Diff"]  = summary(Fullmod)$coefficients['Duration','Pr(>|t|)']
# Occurance_Pvals["MS1","Diff"] = summary(Fullmod)$coefficients['Occurance','Pr(>|t|)']
# Coverage_Pvals["MS1","Diff"]  = summary(Fullmod)$coefficients['Coverage','Pr(>|t|)']
# GEV_Pvals["MS1","Diff"]       = summary(Fullmod)$coefficients['GEV','Pr(>|t|)']
# TP1_Pvals["MS1","Diff"]       = summary(Fullmod)$coefficients['TP1','Pr(>|t|)']
# TP2_Pvals["MS1","Diff"]       = summary(Fullmod)$coefficients['TP2','Pr(>|t|)']
# TP3_Pvals["MS1","Diff"]       = summary(Fullmod)$coefficients['TP3','Pr(>|t|)']
# TP4_Pvals["MS1","Diff"]       = summary(Fullmod)$coefficients['TP4','Pr(>|t|)']
# TP5_Pvals["MS1","Diff"]       = summary(Fullmod)$coefficients['TP5','Pr(>|t|)']


## MS2
# RS1
# Fullmod <- lmer(DSInt ~ GFP  + Duration  + Occurance  + Coverage  + GEV  + TP1  + TP3  + TP4  + TP5  + (1|Subject), data = MS2RS1,REML = FALSE)
Fullmod <- lm(DSInt ~ GFP + Duration  + Occurance  + Coverage  + GEV  + TP1  + TP3  + TP4  + TP5, data = MS2RS1)
summary(Fullmod)
anova(Fullmod)
Duration_Pvals["MS2","RS1"]  = summary(Fullmod)$coefficients['Duration','Pr(>|t|)']
Occurance_Pvals["MS2","RS1"] = summary(Fullmod)$coefficients['Occurance','Pr(>|t|)']
Coverage_Pvals["MS2","RS1"]  = summary(Fullmod)$coefficients['Coverage','Pr(>|t|)']
GEV_Pvals["MS2","RS1"]       = summary(Fullmod)$coefficients['GEV','Pr(>|t|)']
GFP_Pvals["MS2","RS1"]       = summary(Fullmod)$coefficients['GFP','Pr(>|t|)']
TP1_Pvals["MS2","RS1"]       = summary(Fullmod)$coefficients['TP1','Pr(>|t|)']
TP3_Pvals["MS2","RS1"]       = summary(Fullmod)$coefficients['TP3','Pr(>|t|)']
TP4_Pvals["MS2","RS1"]       = summary(Fullmod)$coefficients['TP4','Pr(>|t|)']
TP5_Pvals["MS2","RS1"]       = summary(Fullmod)$coefficients['TP5','Pr(>|t|)']

Fullmod <- lm(DSInt ~  TP3, data = MS2RS1)
anova(Fullmod)
Intplot <- ggpredict(Fullmod,terms = c('TP3'))
plot(Intplot,add.data = TRUE)

# # RS2
# Fullmod <- lmer(DSInt ~ GFP  + Duration  + Occurance  + Coverage  + GEV  + TP1  + TP3  + TP4  + TP5  + (1|Subject), data = MS2RS2,REML = FALSE)
# anova(Fullmod)
# Duration_Pvals["MS2","RS2"]  = summary(Fullmod)$coefficients['Duration','Pr(>|t|)']
# Occurance_Pvals["MS2","RS2"] = summary(Fullmod)$coefficients['Occurance','Pr(>|t|)']
# Coverage_Pvals["MS2","RS2"]  = summary(Fullmod)$coefficients['Coverage','Pr(>|t|)']
# GEV_Pvals["MS2","RS2"]       = summary(Fullmod)$coefficients['GEV','Pr(>|t|)']
# TP1_Pvals["MS2","RS2"]       = summary(Fullmod)$coefficients['TP1','Pr(>|t|)']
# TP2_Pvals["MS2","RS2"]       = summary(Fullmod)$coefficients['TP2','Pr(>|t|)']
# TP3_Pvals["MS2","RS2"]       = summary(Fullmod)$coefficients['TP3','Pr(>|t|)']
# TP4_Pvals["MS2","RS2"]       = summary(Fullmod)$coefficients['TP4','Pr(>|t|)']
# TP5_Pvals["MS2","RS2"]       = summary(Fullmod)$coefficients['TP5','Pr(>|t|)']
# 
# # Diff
# Fullmod <- lmer(DSInt ~ GFP  + Duration  + Occurance  + Coverage  + GEV  + TP1  + TP3  + TP4  + TP5  + (1|Subject), data = MS2diff,REML = FALSE)
# anova(Fullmod)
# Duration_Pvals["MS2","Diff"]  = summary(Fullmod)$coefficients['Duration','Pr(>|t|)']
# Occurance_Pvals["MS2","Diff"] = summary(Fullmod)$coefficients['Occurance','Pr(>|t|)']
# Coverage_Pvals["MS2","Diff"]  = summary(Fullmod)$coefficients['Coverage','Pr(>|t|)']
# GEV_Pvals["MS2","Diff"]       = summary(Fullmod)$coefficients['GEV','Pr(>|t|)']
# TP1_Pvals["MS2","Diff"]       = summary(Fullmod)$coefficients['TP1','Pr(>|t|)']
# TP2_Pvals["MS2","Diff"]       = summary(Fullmod)$coefficients['TP2','Pr(>|t|)']
# TP3_Pvals["MS2","Diff"]       = summary(Fullmod)$coefficients['TP3','Pr(>|t|)']
# TP4_Pvals["MS2","Diff"]       = summary(Fullmod)$coefficients['TP4','Pr(>|t|)']
# TP5_Pvals["MS2","Diff"]       = summary(Fullmod)$coefficients['TP5','Pr(>|t|)']


## MS3
# RS1
# Fullmod <- lmer(DSInt ~ GFP  + Duration  + Occurance  + Coverage  + GEV  + TP2  + TP1  + TP4  + TP5  + (1|Subject), data = MS3RS1,REML = FALSE)
Fullmod <- lm(DSInt ~ GFP  + Duration  + Occurance  + Coverage  + GEV  + TP2  + TP1  + TP4  + TP5, data = MS3RS1)
summary(Fullmod)
Duration_Pvals["MS3","RS1"]  = summary(Fullmod)$coefficients['Duration','Pr(>|t|)']
Occurance_Pvals["MS3","RS1"] = summary(Fullmod)$coefficients['Occurance','Pr(>|t|)']
Coverage_Pvals["MS3","RS1"]  = summary(Fullmod)$coefficients['Coverage','Pr(>|t|)']
GEV_Pvals["MS3","RS1"]       = summary(Fullmod)$coefficients['GEV','Pr(>|t|)']
GFP_Pvals["MS3","RS1"]       = summary(Fullmod)$coefficients['GFP','Pr(>|t|)']
TP1_Pvals["MS3","RS1"]       = summary(Fullmod)$coefficients['TP1','Pr(>|t|)']
TP2_Pvals["MS3","RS1"]       = summary(Fullmod)$coefficients['TP2','Pr(>|t|)']
TP4_Pvals["MS3","RS1"]       = summary(Fullmod)$coefficients['TP4','Pr(>|t|)']
TP5_Pvals["MS3","RS1"]       = summary(Fullmod)$coefficients['TP5','Pr(>|t|)']

# # RS2
# Fullmod <- lmer(DSInt ~ GFP  + Duration  + Occurance  + Coverage  + GEV  + TP2  + TP1  + TP4  + TP5  + (1|Subject), data = MS3RS2,REML = FALSE)
# anova(Fullmod)
# Duration_Pvals["MS3","RS2"]  = summary(Fullmod)$coefficients['Duration','Pr(>|t|)']
# Occurance_Pvals["MS3","RS2"] = summary(Fullmod)$coefficients['Occurance','Pr(>|t|)']
# Coverage_Pvals["MS3","RS2"]  = summary(Fullmod)$coefficients['Coverage','Pr(>|t|)']
# GEV_Pvals["MS3","RS2"]       = summary(Fullmod)$coefficients['GEV','Pr(>|t|)']
# TP1_Pvals["MS3","RS2"]       = summary(Fullmod)$coefficients['TP1','Pr(>|t|)']
# TP2_Pvals["MS3","RS2"]       = summary(Fullmod)$coefficients['TP2','Pr(>|t|)']
# TP3_Pvals["MS3","RS2"]       = summary(Fullmod)$coefficients['TP3','Pr(>|t|)']
# TP4_Pvals["MS3","RS2"]       = summary(Fullmod)$coefficients['TP4','Pr(>|t|)']
# TP5_Pvals["MS3","RS2"]       = summary(Fullmod)$coefficients['TP5','Pr(>|t|)']
# 
# # Diff
# Fullmod <- lmer(DSInt ~ GFP  + Duration  + Occurance  + Coverage  + GEV  + TP2  + TP1  + TP4  + TP5  + (1|Subject), data = MS3diff,REML = FALSE)
# anova(Fullmod)
# Duration_Pvals["MS3","Diff"]  = summary(Fullmod)$coefficients['Duration','Pr(>|t|)']
# Occurance_Pvals["MS3","Diff"] = summary(Fullmod)$coefficients['Occurance','Pr(>|t|)']
# Coverage_Pvals["MS3","Diff"]  = summary(Fullmod)$coefficients['Coverage','Pr(>|t|)']
# GEV_Pvals["MS3","Diff"]       = summary(Fullmod)$coefficients['GEV','Pr(>|t|)']
# TP1_Pvals["MS3","Diff"]       = summary(Fullmod)$coefficients['TP1','Pr(>|t|)']
# TP2_Pvals["MS3","Diff"]       = summary(Fullmod)$coefficients['TP2','Pr(>|t|)']
# TP3_Pvals["MS3","Diff"]       = summary(Fullmod)$coefficients['TP3','Pr(>|t|)']
# TP4_Pvals["MS3","Diff"]       = summary(Fullmod)$coefficients['TP4','Pr(>|t|)']
# TP5_Pvals["MS3","Diff"]       = summary(Fullmod)$coefficients['TP5','Pr(>|t|)']


## MS4
# RS1
# Fullmod <- lmer(DSInt ~ GFP  + Duration  + Occurance  + Coverage  + GEV  + TP2  + TP3  + TP1  + TP5  + (1|Subject), data = MS4RS1,REML = FALSE)
Fullmod <- lm(DSInt ~ GFP  + Duration  + Occurance  + Coverage  + GEV  + TP2  + TP3  + TP1  + TP5, data = MS4RS1)
summary(Fullmod)
Duration_Pvals["MS4","RS1"]  = summary(Fullmod)$coefficients['Duration','Pr(>|t|)']
Occurance_Pvals["MS4","RS1"] = summary(Fullmod)$coefficients['Occurance','Pr(>|t|)']
Coverage_Pvals["MS4","RS1"]  = summary(Fullmod)$coefficients['Coverage','Pr(>|t|)']
GEV_Pvals["MS4","RS1"]       = summary(Fullmod)$coefficients['GEV','Pr(>|t|)']
GFP_Pvals["MS4","RS1"]       = summary(Fullmod)$coefficients['GFP','Pr(>|t|)']
TP1_Pvals["MS4","RS1"]       = summary(Fullmod)$coefficients['TP1','Pr(>|t|)']
TP2_Pvals["MS4","RS1"]       = summary(Fullmod)$coefficients['TP2','Pr(>|t|)']
TP3_Pvals["MS4","RS1"]       = summary(Fullmod)$coefficients['TP3','Pr(>|t|)']
TP5_Pvals["MS4","RS1"]       = summary(Fullmod)$coefficients['TP5','Pr(>|t|)']

# Fullmod <- lm(DSInt ~  Coverage, data = MS4RS1)
# Intplot <- ggpredict(Fullmod,terms = c('Coverage'))
# plot(Intplot,add.data = TRUE)

# # RS2
# Fullmod <- lmer(DSInt ~ GFP  + Duration  + Occurance  + Coverage  + GEV  + TP2  + TP3  + TP1  + TP5  + (1|Subject), data = MS4RS2,REML = FALSE)
# anova(Fullmod)
# Duration_Pvals["MS4","RS2"]  = summary(Fullmod)$coefficients['Duration','Pr(>|t|)']
# Occurance_Pvals["MS4","RS2"] = summary(Fullmod)$coefficients['Occurance','Pr(>|t|)']
# Coverage_Pvals["MS4","RS2"]  = summary(Fullmod)$coefficients['Coverage','Pr(>|t|)']
# GEV_Pvals["MS4","RS2"]       = summary(Fullmod)$coefficients['GEV','Pr(>|t|)']
# TP1_Pvals["MS4","RS2"]       = summary(Fullmod)$coefficients['TP1','Pr(>|t|)']
# TP2_Pvals["MS4","RS2"]       = summary(Fullmod)$coefficients['TP2','Pr(>|t|)']
# TP3_Pvals["MS4","RS2"]       = summary(Fullmod)$coefficients['TP3','Pr(>|t|)']
# TP4_Pvals["MS4","RS2"]       = summary(Fullmod)$coefficients['TP4','Pr(>|t|)']
# TP5_Pvals["MS4","RS2"]       = summary(Fullmod)$coefficients['TP5','Pr(>|t|)']
# 
# # Diff
# Fullmod <- lmer(DSInt ~ GFP  + Duration  + Occurance  + Coverage  + GEV  + TP2  + TP3  + TP1  + TP5  + (1|Subject), data = MS4diff,REML = FALSE)
# anova(Fullmod)
# Duration_Pvals["MS4","Diff"]  = summary(Fullmod)$coefficients['Duration','Pr(>|t|)']
# Occurance_Pvals["MS4","Diff"] = summary(Fullmod)$coefficients['Occurance','Pr(>|t|)']
# Coverage_Pvals["MS4","Diff"]  = summary(Fullmod)$coefficients['Coverage','Pr(>|t|)']
# GEV_Pvals["MS4","Diff"]       = summary(Fullmod)$coefficients['GEV','Pr(>|t|)']
# TP1_Pvals["MS4","Diff"]       = summary(Fullmod)$coefficients['TP1','Pr(>|t|)']
# TP2_Pvals["MS4","Diff"]       = summary(Fullmod)$coefficients['TP2','Pr(>|t|)']
# TP3_Pvals["MS4","Diff"]       = summary(Fullmod)$coefficients['TP3','Pr(>|t|)']
# TP4_Pvals["MS4","Diff"]       = summary(Fullmod)$coefficients['TP4','Pr(>|t|)']
# TP5_Pvals["MS4","Diff"]       = summary(Fullmod)$coefficients['TP5','Pr(>|t|)']


## MS5
# RS1
# Fullmod <- lmer(DSInt ~ GFP  + Duration  + Occurance  + Coverage  + GEV  + TP2  + TP3  + TP4  + TP1  + (1|Subject), data = MS5RS1,REML = FALSE)
Fullmod <- lm(DSInt ~ GFP  + Duration  + Occurance  + Coverage  + GEV  + TP2  + TP3  + TP4  + TP1, data = MS5RS1)
summary(Fullmod)
Duration_Pvals["MS5","RS1"]  = summary(Fullmod)$coefficients['Duration','Pr(>|t|)']
Occurance_Pvals["MS5","RS1"] = summary(Fullmod)$coefficients['Occurance','Pr(>|t|)']
Coverage_Pvals["MS5","RS1"]  = summary(Fullmod)$coefficients['Coverage','Pr(>|t|)']
GEV_Pvals["MS5","RS1"]       = summary(Fullmod)$coefficients['GEV','Pr(>|t|)']
GFP_Pvals["MS5","RS1"]       = summary(Fullmod)$coefficients['GFP','Pr(>|t|)']
TP1_Pvals["MS5","RS1"]       = summary(Fullmod)$coefficients['TP1','Pr(>|t|)']
TP2_Pvals["MS5","RS1"]       = summary(Fullmod)$coefficients['TP2','Pr(>|t|)']
TP3_Pvals["MS5","RS1"]       = summary(Fullmod)$coefficients['TP3','Pr(>|t|)']
TP4_Pvals["MS5","RS1"]       = summary(Fullmod)$coefficients['TP4','Pr(>|t|)']

# # RS2
# Fullmod <- lmer(DSInt ~ GFP  + Duration  + Occurance  + Coverage  + GEV  + TP2  + TP3  + TP4  + TP1  + (1|Subject), data = MS5RS2,REML = FALSE)
# anova(Fullmod)
# Duration_Pvals["MS5","RS2"]  = summary(Fullmod)$coefficients['Duration','Pr(>|t|)']
# Occurance_Pvals["MS5","RS2"] = summary(Fullmod)$coefficients['Occurance','Pr(>|t|)']
# Coverage_Pvals["MS5","RS2"]  = summary(Fullmod)$coefficients['Coverage','Pr(>|t|)']
# GEV_Pvals["MS5","RS2"]       = summary(Fullmod)$coefficients['GEV','Pr(>|t|)']
# TP1_Pvals["MS5","RS2"]       = summary(Fullmod)$coefficients['TP1','Pr(>|t|)']
# TP2_Pvals["MS5","RS2"]       = summary(Fullmod)$coefficients['TP2','Pr(>|t|)']
# TP3_Pvals["MS5","RS2"]       = summary(Fullmod)$coefficients['TP3','Pr(>|t|)']
# TP4_Pvals["MS5","RS2"]       = summary(Fullmod)$coefficients['TP4','Pr(>|t|)']
# TP5_Pvals["MS5","RS2"]       = summary(Fullmod)$coefficients['TP5','Pr(>|t|)']
# 
# # Diff
# Fullmod <- lmer(DSInt ~ GFP  + Duration  + Occurance  + Coverage  + GEV  + TP2  + TP3  + TP4  + TP1  + (1|Subject), data = MS5diff,REML = FALSE)
# anova(Fullmod)
# Duration_Pvals["MS5","Diff"]  = summary(Fullmod)$coefficients['Duration','Pr(>|t|)']
# Occurance_Pvals["MS5","Diff"] = summary(Fullmod)$coefficients['Occurance','Pr(>|t|)']
# Coverage_Pvals["MS5","Diff"]  = summary(Fullmod)$coefficients['Coverage','Pr(>|t|)']
# GEV_Pvals["MS5","Diff"]       = summary(Fullmod)$coefficients['GEV','Pr(>|t|)']
# TP1_Pvals["MS5","Diff"]       = summary(Fullmod)$coefficients['TP1','Pr(>|t|)']
# TP2_Pvals["MS5","Diff"]       = summary(Fullmod)$coefficients['TP2','Pr(>|t|)']
# TP3_Pvals["MS5","Diff"]       = summary(Fullmod)$coefficients['TP3','Pr(>|t|)']
# TP4_Pvals["MS5","Diff"]       = summary(Fullmod)$coefficients['TP4','Pr(>|t|)']
# TP5_Pvals["MS5","Diff"]       = summary(Fullmod)$coefficients['TP5','Pr(>|t|)']

# Check for significance after FDR correction
Pvalues <- unlist(TP4_Pvals)
non_na_indices <- which(!is.na(Pvalues))
All_adj_p_values <- p.adjust(Pvalues[non_na_indices], method = "fdr")
Pvalues[non_na_indices] <- All_adj_p_values
Pvalues[Pvalues<.05]

# Function to compute p-values for individual predictors
Fit_Models_Averaged <- function(predictor_names, dependent_variable, data,MS) {
  # Initialize a named vector to store p-values
  p_values <- numeric(length(predictor_names))
  new_predictor_names <- paste0(MS, predictor_names)
  p_values <- setNames(c(p_values), c(new_predictor_names))
  
  # Loop over each predictor
  for (i in 1:length(predictor_names)) {
    # Create the formula for the model
    formula <- formula(paste(dependent_variable, "~ IntensityT0 + Arm*", predictor_names[i]))
    
    # Fit the linear model
    model <- lm(formula, data = data)
    print(summary(model))
    # Extract the p-value for the predictor and store it with the predictor name
    print(summary(model)$coefficients[5, "Pr(>|t|)"])
    p_values[new_predictor_names[i]] <- summary(model)$coefficients[5, "Pr(>|t|)"]
  }
  # Return the named vector of p-values
  print(p_values)
  return(p_values)
}


# Fit predictors separately to avoid the confusion of the other predictors
predictors <- c("GFP", "Duration", "Occurance", "Coverage", "GEV", "TP2", "TP3", "TP4","TP5")
MS1_Avgint <- Fit_Models_Averaged(predictors, "IntensityRatings", MS1RS1,'MS1')

predictors <- c("GFP", "Duration", "Occurance", "Coverage", "GEV", "TP1", "TP3", "TP4","TP5")
MS2_Avgint <- Fit_Models_Averaged(predictors, "IntensityRatings", MS2RS1,'MS2')

predictors <- c("GFP", "Duration", "Occurance", "Coverage", "GEV", "TP2", "TP1", "TP4","TP5")
MS3_Avgint <- Fit_Models_Averaged(predictors, "IntensityRatings", MS3RS1,'MS3')

predictors <- c("GFP", "Duration", "Occurance", "Coverage", "GEV", "TP2", "TP3", "TP1","TP5")
MS4_Avgint <- Fit_Models_Averaged(predictors, "IntensityRatings", MS4RS1,'MS4')

predictors <- c("GFP", "Duration", "Occurance", "Coverage", "GEV", "TP1", "TP2", "TP3", "TP4")
MS5_Avgint <- Fit_Models_Averaged(predictors, "IntensityRatings", MS5RS1,'MS5')

## Multiple comparisons correction for Area Length
Pvalues <- unlist(c(MS1_Avgint,MS2_Avgint,MS3_Avgint,MS4_Avgint,MS5_Avgint))
non_na_indices <- which(!is.na(Pvalues))
All_adj_p_values <- p.adjust(Pvalues[non_na_indices], method = "fdr")
Pvalues[non_na_indices] <- All_adj_p_values
Pvalues[Pvalues<.05]


### Unpleasantness ----
# Initialize variables

# Create column and row names
col_names <- c("RS1","RS2","Diff")
row_names <- c("MS1", "MS2", "MS3", "MS4", "MS5")
# Create empty data frames
Duration_Pvals <- data.frame(matrix(ncol = length(col_names), nrow = length(row_names)))
colnames(Duration_Pvals) <- col_names
rownames(Duration_Pvals) <- row_names
Occurance_Pvals = Duration_Pvals
Coverage_Pvals = Duration_Pvals
GEV_Pvals = Duration_Pvals
GFP_Pvals = Duration_Pvals
TP1_Pvals = Duration_Pvals
TP2_Pvals = Duration_Pvals
TP3_Pvals = Duration_Pvals
TP4_Pvals = Duration_Pvals
TP5_Pvals = Duration_Pvals

## MS1
# RS1
Fullmod <- lm(DSUnp ~ GFP  + Duration  + Occurance  + Coverage  + GEV  + TP2  + TP3  + TP4  + TP5, data = MS1RS1)
summary(Fullmod)
Duration_Pvals["MS1","RS1"]  = summary(Fullmod)$coefficients['Duration','Pr(>|t|)']
Occurance_Pvals["MS1","RS1"] = summary(Fullmod)$coefficients['Occurance','Pr(>|t|)']
Coverage_Pvals["MS1","RS1"]  = summary(Fullmod)$coefficients['Coverage','Pr(>|t|)']
GEV_Pvals["MS1","RS1"]       = summary(Fullmod)$coefficients['GEV','Pr(>|t|)']
GFP_Pvals["MS1","RS1"]       = summary(Fullmod)$coefficients['GFP','Pr(>|t|)']
TP2_Pvals["MS1","RS1"]       = summary(Fullmod)$coefficients['TP2','Pr(>|t|)']
TP3_Pvals["MS1","RS1"]       = summary(Fullmod)$coefficients['TP3','Pr(>|t|)']
TP4_Pvals["MS1","RS1"]       = summary(Fullmod)$coefficients['TP4','Pr(>|t|)']
TP5_Pvals["MS1","RS1"]       = summary(Fullmod)$coefficients['TP5','Pr(>|t|)']

# # RS2
# Fullmod <- lmer(DSUnp ~ GFP  + Duration  + Occurance  + Coverage  + GEV  + TP2  + TP3  + TP4  + TP5  + (1|Subject), data = MS1RS2,REML = FALSE)
# anova(Fullmod)
# Duration_Pvals["MS1","RS2"]  = summary(Fullmod)$coefficients['Duration','Pr(>|t|)']
# Occurance_Pvals["MS1","RS2"] = summary(Fullmod)$coefficients['Occurance','Pr(>|t|)']
# Coverage_Pvals["MS1","RS2"]  = summary(Fullmod)$coefficients['Coverage','Pr(>|t|)']
# GEV_Pvals["MS1","RS2"]       = summary(Fullmod)$coefficients['GEV','Pr(>|t|)']
# TP1_Pvals["MS1","RS2"]       = summary(Fullmod)$coefficients['TP1','Pr(>|t|)']
# TP2_Pvals["MS1","RS2"]       = summary(Fullmod)$coefficients['TP2','Pr(>|t|)']
# TP3_Pvals["MS1","RS2"]       = summary(Fullmod)$coefficients['TP3','Pr(>|t|)']
# TP4_Pvals["MS1","RS2"]       = summary(Fullmod)$coefficients['TP4','Pr(>|t|)']
# TP5_Pvals["MS1","RS2"]       = summary(Fullmod)$coefficients['TP5','Pr(>|t|)']
# 
# # Diff
# Fullmod <- lmer(DSUnp ~ GFP  + Duration  + Occurance  + Coverage  + GEV  + TP2  + TP3  + TP4  + TP5  + (1|Subject), data = MS1diff,REML = FALSE)
# anova(Fullmod)
# Duration_Pvals["MS1","Diff"]  = summary(Fullmod)$coefficients['Duration','Pr(>|t|)']
# Occurance_Pvals["MS1","Diff"] = summary(Fullmod)$coefficients['Occurance','Pr(>|t|)']
# Coverage_Pvals["MS1","Diff"]  = summary(Fullmod)$coefficients['Coverage','Pr(>|t|)']
# GEV_Pvals["MS1","Diff"]       = summary(Fullmod)$coefficients['GEV','Pr(>|t|)']
# TP1_Pvals["MS1","Diff"]       = summary(Fullmod)$coefficients['TP1','Pr(>|t|)']
# TP2_Pvals["MS1","Diff"]       = summary(Fullmod)$coefficients['TP2','Pr(>|t|)']
# TP3_Pvals["MS1","Diff"]       = summary(Fullmod)$coefficients['TP3','Pr(>|t|)']
# TP4_Pvals["MS1","Diff"]       = summary(Fullmod)$coefficients['TP4','Pr(>|t|)']
# TP5_Pvals["MS1","Diff"]       = summary(Fullmod)$coefficients['TP5','Pr(>|t|)']


## MS2
# RS1
Fullmod <- lm(DSUnp ~ GFP  + Duration  + Occurance  + Coverage  + GEV  + TP1  + TP3  + TP4  + TP5, data = MS2RS1)
summary(Fullmod)
Duration_Pvals["MS2","RS1"]  = summary(Fullmod)$coefficients['Duration','Pr(>|t|)']
Occurance_Pvals["MS2","RS1"] = summary(Fullmod)$coefficients['Occurance','Pr(>|t|)']
Coverage_Pvals["MS2","RS1"]  = summary(Fullmod)$coefficients['Coverage','Pr(>|t|)']
GEV_Pvals["MS2","RS1"]       = summary(Fullmod)$coefficients['GEV','Pr(>|t|)']
GFP_Pvals["MS2","RS1"]       = summary(Fullmod)$coefficients['GFP','Pr(>|t|)']
TP1_Pvals["MS2","RS1"]       = summary(Fullmod)$coefficients['TP1','Pr(>|t|)']
TP3_Pvals["MS2","RS1"]       = summary(Fullmod)$coefficients['TP3','Pr(>|t|)']
TP4_Pvals["MS2","RS1"]       = summary(Fullmod)$coefficients['TP4','Pr(>|t|)']
TP5_Pvals["MS2","RS1"]       = summary(Fullmod)$coefficients['TP5','Pr(>|t|)']

# # RS2
# Fullmod <- lmer(DSUnp ~ GFP  + Duration  + Occurance  + Coverage  + GEV  + TP1  + TP3  + TP4  + TP5  + (1|Subject), data = MS2RS2,REML = FALSE)
# anova(Fullmod)
# Duration_Pvals["MS2","RS2"]  = summary(Fullmod)$coefficients['Duration','Pr(>|t|)']
# Occurance_Pvals["MS2","RS2"] = summary(Fullmod)$coefficients['Occurance','Pr(>|t|)']
# Coverage_Pvals["MS2","RS2"]  = summary(Fullmod)$coefficients['Coverage','Pr(>|t|)']
# GEV_Pvals["MS2","RS2"]       = summary(Fullmod)$coefficients['GEV','Pr(>|t|)']
# TP1_Pvals["MS2","RS2"]       = summary(Fullmod)$coefficients['TP1','Pr(>|t|)']
# TP2_Pvals["MS2","RS2"]       = summary(Fullmod)$coefficients['TP2','Pr(>|t|)']
# TP3_Pvals["MS2","RS2"]       = summary(Fullmod)$coefficients['TP3','Pr(>|t|)']
# TP4_Pvals["MS2","RS2"]       = summary(Fullmod)$coefficients['TP4','Pr(>|t|)']
# TP5_Pvals["MS2","RS2"]       = summary(Fullmod)$coefficients['TP5','Pr(>|t|)']
# 
# # Diff
# Fullmod <- lmer(DSUnp ~ GFP  + Duration  + Occurance  + Coverage  + GEV  + TP1  + TP3  + TP4  + TP5  + (1|Subject), data = MS2diff,REML = FALSE)
# anova(Fullmod)
# Duration_Pvals["MS2","Diff"]  = summary(Fullmod)$coefficients['Duration','Pr(>|t|)']
# Occurance_Pvals["MS2","Diff"] = summary(Fullmod)$coefficients['Occurance','Pr(>|t|)']
# Coverage_Pvals["MS2","Diff"]  = summary(Fullmod)$coefficients['Coverage','Pr(>|t|)']
# GEV_Pvals["MS2","Diff"]       = summary(Fullmod)$coefficients['GEV','Pr(>|t|)']
# TP1_Pvals["MS2","Diff"]       = summary(Fullmod)$coefficients['TP1','Pr(>|t|)']
# TP2_Pvals["MS2","Diff"]       = summary(Fullmod)$coefficients['TP2','Pr(>|t|)']
# TP3_Pvals["MS2","Diff"]       = summary(Fullmod)$coefficients['TP3','Pr(>|t|)']
# TP4_Pvals["MS2","Diff"]       = summary(Fullmod)$coefficients['TP4','Pr(>|t|)']
# TP5_Pvals["MS2","Diff"]       = summary(Fullmod)$coefficients['TP5','Pr(>|t|)']


## MS3
# RS1
Fullmod <- lm(DSUnp ~ GFP  + Duration  + Occurance  + Coverage  + GEV  + TP2  + TP1  + TP4  + TP5, data = MS3RS1)
summary(Fullmod)
Duration_Pvals["MS3","RS1"]  = summary(Fullmod)$coefficients['Duration','Pr(>|t|)']
Occurance_Pvals["MS3","RS1"] = summary(Fullmod)$coefficients['Occurance','Pr(>|t|)']
Coverage_Pvals["MS3","RS1"]  = summary(Fullmod)$coefficients['Coverage','Pr(>|t|)']
GEV_Pvals["MS3","RS1"]       = summary(Fullmod)$coefficients['GEV','Pr(>|t|)']
GFP_Pvals["MS3","RS1"]       = summary(Fullmod)$coefficients['GFP','Pr(>|t|)']
TP1_Pvals["MS3","RS1"]       = summary(Fullmod)$coefficients['TP1','Pr(>|t|)']
TP2_Pvals["MS3","RS1"]       = summary(Fullmod)$coefficients['TP2','Pr(>|t|)']
TP4_Pvals["MS3","RS1"]       = summary(Fullmod)$coefficients['TP4','Pr(>|t|)']
TP5_Pvals["MS3","RS1"]       = summary(Fullmod)$coefficients['TP5','Pr(>|t|)']

Fullmod <- lm(DSUnp ~  TP1, data = MS3RS1)
summary(Fullmod)
Intplot <- ggpredict(Fullmod,terms = c('TP1'))
plot(Intplot,add.data = TRUE)

# # RS2
# Fullmod <- lmer(DSUnp ~ GFP  + Duration  + Occurance  + Coverage  + GEV  + TP2  + TP1  + TP4  + TP5  + (1|Subject), data = MS3RS2,REML = FALSE)
# anova(Fullmod)
# Duration_Pvals["MS3","RS2"]  = summary(Fullmod)$coefficients['Duration','Pr(>|t|)']
# Occurance_Pvals["MS3","RS2"] = summary(Fullmod)$coefficients['Occurance','Pr(>|t|)']
# Coverage_Pvals["MS3","RS2"]  = summary(Fullmod)$coefficients['Coverage','Pr(>|t|)']
# GEV_Pvals["MS3","RS2"]       = summary(Fullmod)$coefficients['GEV','Pr(>|t|)']
# TP1_Pvals["MS3","RS2"]       = summary(Fullmod)$coefficients['TP1','Pr(>|t|)']
# TP2_Pvals["MS3","RS2"]       = summary(Fullmod)$coefficients['TP2','Pr(>|t|)']
# TP3_Pvals["MS3","RS2"]       = summary(Fullmod)$coefficients['TP3','Pr(>|t|)']
# TP4_Pvals["MS3","RS2"]       = summary(Fullmod)$coefficients['TP4','Pr(>|t|)']
# TP5_Pvals["MS3","RS2"]       = summary(Fullmod)$coefficients['TP5','Pr(>|t|)']
# 
# # Diff
# Fullmod <- lmer(DSUnp ~ GFP  + Duration  + Occurance  + Coverage  + GEV  + TP2  + TP1  + TP4  + TP5  + (1|Subject), data = MS3diff,REML = FALSE)
# anova(Fullmod)
# Duration_Pvals["MS3","Diff"]  = summary(Fullmod)$coefficients['Duration','Pr(>|t|)']
# Occurance_Pvals["MS3","Diff"] = summary(Fullmod)$coefficients['Occurance','Pr(>|t|)']
# Coverage_Pvals["MS3","Diff"]  = summary(Fullmod)$coefficients['Coverage','Pr(>|t|)']
# GEV_Pvals["MS3","Diff"]       = summary(Fullmod)$coefficients['GEV','Pr(>|t|)']
# TP1_Pvals["MS3","Diff"]       = summary(Fullmod)$coefficients['TP1','Pr(>|t|)']
# TP2_Pvals["MS3","Diff"]       = summary(Fullmod)$coefficients['TP2','Pr(>|t|)']
# TP3_Pvals["MS3","Diff"]       = summary(Fullmod)$coefficients['TP3','Pr(>|t|)']
# TP4_Pvals["MS3","Diff"]       = summary(Fullmod)$coefficients['TP4','Pr(>|t|)']
# TP5_Pvals["MS3","Diff"]       = summary(Fullmod)$coefficients['TP5','Pr(>|t|)']



## MS4
# RS1
Fullmod <- lm(DSUnp ~ GFP  + Duration  + Occurance  + Coverage  + GEV  + TP2  + TP3  + TP1  + TP5, data = MS4RS1)
summary(Fullmod)
Duration_Pvals["MS4","RS1"]  = summary(Fullmod)$coefficients['Duration','Pr(>|t|)']
Occurance_Pvals["MS4","RS1"] = summary(Fullmod)$coefficients['Occurance','Pr(>|t|)']
Coverage_Pvals["MS4","RS1"]  = summary(Fullmod)$coefficients['Coverage','Pr(>|t|)']
GEV_Pvals["MS4","RS1"]       = summary(Fullmod)$coefficients['GEV','Pr(>|t|)']
GFP_Pvals["MS4","RS1"]       = summary(Fullmod)$coefficients['GFP','Pr(>|t|)']
TP1_Pvals["MS4","RS1"]       = summary(Fullmod)$coefficients['TP1','Pr(>|t|)']
TP2_Pvals["MS4","RS1"]       = summary(Fullmod)$coefficients['TP2','Pr(>|t|)']
TP3_Pvals["MS4","RS1"]       = summary(Fullmod)$coefficients['TP3','Pr(>|t|)']
TP5_Pvals["MS4","RS1"]       = summary(Fullmod)$coefficients['TP5','Pr(>|t|)']

# # RS2
# Fullmod <- lmer(DSUnp ~ GFP  + Duration  + Occurance  + Coverage  + GEV  + TP2  + TP3  + TP1  + TP5  + (1|Subject), data = MS4RS2,REML = FALSE)
# anova(Fullmod)
# Duration_Pvals["MS4","RS2"]  = summary(Fullmod)$coefficients['Duration','Pr(>|t|)']
# Occurance_Pvals["MS4","RS2"] = summary(Fullmod)$coefficients['Occurance','Pr(>|t|)']
# Coverage_Pvals["MS4","RS2"]  = summary(Fullmod)$coefficients['Coverage','Pr(>|t|)']
# GEV_Pvals["MS4","RS2"]       = summary(Fullmod)$coefficients['GEV','Pr(>|t|)']
# TP1_Pvals["MS4","RS2"]       = summary(Fullmod)$coefficients['TP1','Pr(>|t|)']
# TP2_Pvals["MS4","RS2"]       = summary(Fullmod)$coefficients['TP2','Pr(>|t|)']
# TP3_Pvals["MS4","RS2"]       = summary(Fullmod)$coefficients['TP3','Pr(>|t|)']
# TP4_Pvals["MS4","RS2"]       = summary(Fullmod)$coefficients['TP4','Pr(>|t|)']
# TP5_Pvals["MS4","RS2"]       = summary(Fullmod)$coefficients['TP5','Pr(>|t|)']
# 
# # Diff
# Fullmod <- lmer(DSUnp ~ GFP  + Duration  + Occurance  + Coverage  + GEV  + TP2  + TP3  + TP1  + TP5  + (1|Subject), data = MS4diff,REML = FALSE)
# anova(Fullmod)
# Duration_Pvals["MS4","Diff"]  = summary(Fullmod)$coefficients['Duration','Pr(>|t|)']
# Occurance_Pvals["MS4","Diff"] = summary(Fullmod)$coefficients['Occurance','Pr(>|t|)']
# Coverage_Pvals["MS4","Diff"]  = summary(Fullmod)$coefficients['Coverage','Pr(>|t|)']
# GEV_Pvals["MS4","Diff"]       = summary(Fullmod)$coefficients['GEV','Pr(>|t|)']
# TP1_Pvals["MS4","Diff"]       = summary(Fullmod)$coefficients['TP1','Pr(>|t|)']
# TP2_Pvals["MS4","Diff"]       = summary(Fullmod)$coefficients['TP2','Pr(>|t|)']
# TP3_Pvals["MS4","Diff"]       = summary(Fullmod)$coefficients['TP3','Pr(>|t|)']
# TP4_Pvals["MS4","Diff"]       = summary(Fullmod)$coefficients['TP4','Pr(>|t|)']
# TP5_Pvals["MS4","Diff"]       = summary(Fullmod)$coefficients['TP5','Pr(>|t|)']

## MS5
# RS1
Fullmod <- lm(DSUnp ~ GFP  + Duration  + Occurance  + Coverage  + GEV  + TP2  + TP3  + TP4  + TP1, data = MS5RS1)
summary(Fullmod)
Duration_Pvals["MS5","RS1"]  = summary(Fullmod)$coefficients['Duration','Pr(>|t|)']
Occurance_Pvals["MS5","RS1"] = summary(Fullmod)$coefficients['Occurance','Pr(>|t|)']
Coverage_Pvals["MS5","RS1"]  = summary(Fullmod)$coefficients['Coverage','Pr(>|t|)']
GEV_Pvals["MS5","RS1"]       = summary(Fullmod)$coefficients['GEV','Pr(>|t|)']
GFP_Pvals["MS5","RS1"]       = summary(Fullmod)$coefficients['GFP','Pr(>|t|)']
TP1_Pvals["MS5","RS1"]       = summary(Fullmod)$coefficients['TP1','Pr(>|t|)']
TP2_Pvals["MS5","RS1"]       = summary(Fullmod)$coefficients['TP2','Pr(>|t|)']
TP3_Pvals["MS5","RS1"]       = summary(Fullmod)$coefficients['TP3','Pr(>|t|)']
TP4_Pvals["MS5","RS1"]       = summary(Fullmod)$coefficients['TP4','Pr(>|t|)']

# # RS2
# Fullmod <- lmer(DSUnp ~ GFP  + Duration  + Occurance  + Coverage  + GEV  + TP2  + TP3  + TP4  + TP1  + (1|Subject), data = MS5RS2,REML = FALSE)
# anova(Fullmod)
# Duration_Pvals["MS5","RS2"]  = summary(Fullmod)$coefficients['Duration','Pr(>|t|)']
# Occurance_Pvals["MS5","RS2"] = summary(Fullmod)$coefficients['Occurance','Pr(>|t|)']
# Coverage_Pvals["MS5","RS2"]  = summary(Fullmod)$coefficients['Coverage','Pr(>|t|)']
# GEV_Pvals["MS5","RS2"]       = summary(Fullmod)$coefficients['GEV','Pr(>|t|)']
# TP1_Pvals["MS5","RS2"]       = summary(Fullmod)$coefficients['TP1','Pr(>|t|)']
# TP2_Pvals["MS5","RS2"]       = summary(Fullmod)$coefficients['TP2','Pr(>|t|)']
# TP3_Pvals["MS5","RS2"]       = summary(Fullmod)$coefficients['TP3','Pr(>|t|)']
# TP4_Pvals["MS5","RS2"]       = summary(Fullmod)$coefficients['TP4','Pr(>|t|)']
# TP5_Pvals["MS5","RS2"]       = summary(Fullmod)$coefficients['TP5','Pr(>|t|)']
# 
# # Diff
# Fullmod <- lmer(DSUnp ~ GFP  + Duration  + Occurance  + Coverage  + GEV  + TP2  + TP3  + TP4  + TP1  + (1|Subject), data = MS5diff,REML = FALSE)
# anova(Fullmod)
# Duration_Pvals["MS5","Diff"]  = summary(Fullmod)$coefficients['Duration','Pr(>|t|)']
# Occurance_Pvals["MS5","Diff"] = summary(Fullmod)$coefficients['Occurance','Pr(>|t|)']
# Coverage_Pvals["MS5","Diff"]  = summary(Fullmod)$coefficients['Coverage','Pr(>|t|)']
# GEV_Pvals["MS5","Diff"]       = summary(Fullmod)$coefficients['GEV','Pr(>|t|)']
# TP1_Pvals["MS5","Diff"]       = summary(Fullmod)$coefficients['TP1','Pr(>|t|)']
# TP2_Pvals["MS5","Diff"]       = summary(Fullmod)$coefficients['TP2','Pr(>|t|)']
# TP3_Pvals["MS5","Diff"]       = summary(Fullmod)$coefficients['TP3','Pr(>|t|)']
# TP4_Pvals["MS5","Diff"]       = summary(Fullmod)$coefficients['TP4','Pr(>|t|)']
# TP5_Pvals["MS5","Diff"]       = summary(Fullmod)$coefficients['TP5','Pr(>|t|)']

# Check for significance after FDR correction
Pvalues <- unlist(TP5_Pvals)
non_na_indices <- which(!is.na(Pvalues))
All_adj_p_values <- p.adjust(Pvalues[non_na_indices], method = "fdr")
Pvalues[non_na_indices] <- All_adj_p_values
Pvalues[Pvalues<.05]


Fit_Models_Averaged <- function(predictor_names, dependent_variable, data,MS) {
  # Initialize a named vector to store p-values
  p_values <- numeric(length(predictor_names))
  new_predictor_names <- paste0(MS, predictor_names)
  p_values <- setNames(c(p_values), c(new_predictor_names))
  
  # Loop over each predictor
  for (i in 1:length(predictor_names)) {
    # Create the formula for the model
    formula <- formula(paste(dependent_variable, "~ UnpleasantnessT0 + Arm*", predictor_names[i]))
    
    # Fit the linear model
    model <- lm(formula, data = data)
    print(summary(model))
    # Extract the p-value for the predictor and store it with the predictor name
    print(summary(model)$coefficients[5, "Pr(>|t|)"])
    p_values[new_predictor_names[i]] <- summary(model)$coefficients[5, "Pr(>|t|)"]
  }
  # Return the named vector of p-values
  print(p_values)
  return(p_values)
}

# Fit predictors separately to avoid the confusion of the other predictors
predictors <- c("GFP", "Duration", "Occurance", "Coverage", "GEV", "TP2", "TP3", "TP4","TP5")
MS1_AvgUnp <- Fit_Models_Averaged(predictors, "UnpleasantnessRatings", MS1RS1,'MS1')

predictors <- c("GFP", "Duration", "Occurance", "Coverage", "GEV", "TP1", "TP3", "TP4","TP5")
MS2_AvgUnp <- Fit_Models_Averaged(predictors, "UnpleasantnessRatings", MS2RS1,'MS2')

predictors <- c("GFP", "Duration", "Occurance", "Coverage", "GEV", "TP2", "TP1", "TP4","TP5")
MS3_AvgUnp <- Fit_Models_Averaged(predictors, "UnpleasantnessRatings", MS3RS1,'MS3')

predictors <- c("GFP", "Duration", "Occurance", "Coverage", "GEV", "TP2", "TP3", "TP1","TP5")
MS4_AvgUnp <- Fit_Models_Averaged(predictors, "UnpleasantnessRatings", MS4RS1,'MS4')

predictors <- c("GFP", "Duration", "Occurance", "Coverage", "GEV", "TP1", "TP2", "TP3", "TP4")
MS5_AvgUnp <- Fit_Models_Averaged(predictors, "UnpleasantnessRatings", MS5RS1,'MS5')

## Multiple comparisons correction for Area Length
Pvalues <- unlist(c(MS1_AvgUnp,MS2_AvgUnp,MS3_AvgUnp,MS4_AvgUnp,MS5_AvgUnp))
non_na_indices <- which(!is.na(Pvalues))
All_adj_p_values <- p.adjust(Pvalues[non_na_indices], method = "fdr")
Pvalues[non_na_indices] <- All_adj_p_values
Pvalues[Pvalues<.05]





## Plotting
Microstats <- read_excel("Path_to_Microstate_results/Microstats.xlsx")


Microstats$MS[Microstats$MS == 'MS4'] <- "A"
Microstats$MS[Microstats$MS == 'MS3'] <- "B"
Microstats$MS[Microstats$MS == 'MS2'] <- "C"
Microstats$MS[Microstats$MS == 'MS1'] <- "D"
Microstats$MS[Microstats$MS == 'MS5'] <- "E"

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


