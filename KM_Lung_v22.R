##################################################################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##################################################################################################
# Start
#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
# Dr Reza Rafiee, August 2017
# Research Fellow, Queen's University Belfast
# Survival analysis (including OS and DFS) on lung cancer data (LUAD and LUSC)
# version 2.2, 1st of September 2017
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Loading all survival and dependent libraries
library(survival)
library(survivalROC)
library(pec)
library(pROC)
library(rms)
#============================================
setwd("/home/reza/Documents/Live/")
#============================================
##############################   LUAD  ################################
# Loading a csv file
# Pathfolder <- "C:/Users/3052092/Documents/Live/"   # initialise the path of the csv file 
Pathfolder <- "/home/reza/Documents/Live/"           # initialise the path of the csv file (from my Linux) 
csvfilename <- "LUAD_Samples_DDRD_ordered_23June17.csv" # LUAD
csvfilename <- "LUSC_Samples_DDRD_ordered_23June17.csv" # LUSC
pathcsvfile <- paste(Pathfolder,csvfilename,sep = "")
datafileObject <- read.csv(pathcsvfile, header=T)  # LUAD samples, TCGA
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Use only specific stage (II or datafileObject_stageII_III <- datafileObject[-c(which(datafileObject$AJCC_PATHOLOGIC_TUMOR_STAGE == "Stage IV"), which(datafileObject$AJCC_PATHOLOGIC_TUMOR_STAGE == "Stage I"), which(datafileObject$AJCC_PATHOLOGIC_TUMOR_STAGE == "Stage IA"),which(datafileObject$AJCC_PATHOLOGIC_TUMOR_STAGE == "Stage IB") ),]III, etc.) (possible chemotherapy) samples
length(which(datafileObject$AJCC_PATHOLOGIC_TUMOR_STAGE == "Stage I"))  # n=5 (LUAD), n= 3 (LUSC) 
length(which(datafileObject$AJCC_PATHOLOGIC_TUMOR_STAGE == "Stage II"))  # n=1 (LUAD), n=3 (LUSC)
length(which(datafileObject$AJCC_PATHOLOGIC_TUMOR_STAGE == "Stage IIA"))  # n=50 (LUAD), n=65 (LUSC)
length(which(datafileObject$AJCC_PATHOLOGIC_TUMOR_STAGE == "Stage IIB"))  # n=71 (LUAD), n=94 (LUSC)
length(which(datafileObject$AJCC_PATHOLOGIC_TUMOR_STAGE == "Stage IIIA"))  # n=73 (LUAD), n=63 (LUSC)
length(which(datafileObject$AJCC_PATHOLOGIC_TUMOR_STAGE == "Stage IIIB"))  # n=11 (LUAD), n= 18 (LUSC)
length(which(datafileObject$AJCC_PATHOLOGIC_TUMOR_STAGE == "Stage IV"))  # n=25 (LUAD), n=7 (LUSC)
length(which(datafileObject$AJCC_PATHOLOGIC_TUMOR_STAGE == "Stage IA"))  # n=132 (LUAD), n=89 (LUSC)
length(which(datafileObject$AJCC_PATHOLOGIC_TUMOR_STAGE == "Stage IB"))  # n=140 (LUAD), n=151 (LUSC)
length(which(datafileObject$AJCC_PATHOLOGIC_TUMOR_STAGE == ""))  # n= (LUAD), n=5 (LUSC)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Use only stage II, stage III and IV (possible chemotherapy) samples
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
datafileObject_stageII_III <- datafileObject[-c(which(datafileObject$AJCC_PATHOLOGIC_TUMOR_STAGE == "Stage I"), 
                                                which(datafileObject$AJCC_PATHOLOGIC_TUMOR_STAGE == "Stage IA"),
                                                which(datafileObject$AJCC_PATHOLOGIC_TUMOR_STAGE == "Stage IB")),]
length(which(datafileObject$AJCC_PATHOLOGIC_TUMOR_STAGE == ""))  # n=9, number of NA for the stage
datafileObject_stageII_III <- datafileObject_stageII_III[-c(which(datafileObject_stageII_III$AJCC_PATHOLOGIC_TUMOR_STAGE == "")),]
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Use only stage II and stage III (possible chemotherapy) samples
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
datafileObject_stageII_III <- datafileObject[-c(which(datafileObject$AJCC_PATHOLOGIC_TUMOR_STAGE == "Stage IV"), which(datafileObject$AJCC_PATHOLOGIC_TUMOR_STAGE == "Stage I"), which(datafileObject$AJCC_PATHOLOGIC_TUMOR_STAGE == "Stage IA"),which(datafileObject$AJCC_PATHOLOGIC_TUMOR_STAGE == "Stage IB") ),]
length(which(datafileObject$AJCC_PATHOLOGIC_TUMOR_STAGE == ""))  # n=9 (LUAD), number of NA for the stage, n=5 (LUSC)
datafileObject_stageII_III <- datafileObject_stageII_III[-c(which(datafileObject_stageII_III$AJCC_PATHOLOGIC_TUMOR_STAGE == "")),]
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Exclude patients who received Targeted Molecular Therapy
length(which(datafileObject_stageII_III$TARGETED_MOLECULAR_THERAPY == "NO")) #n=27 (LUAD),  n=45 (LUSC)
length(which(datafileObject_stageII_III$TARGETED_MOLECULAR_THERAPY == "YES")) # n=31 (LUAD), n=41 (LUSC)
length(which(datafileObject_stageII_III$TARGETED_MOLECULAR_THERAPY == "")) # n=148 (LUAD), n=160 (LUSC)
length(c(which(datafileObject_stageII_III$TARGETED_MOLECULAR_THERAPY == "YES"),which(datafileObject_stageII_III$TARGETED_MOLECULAR_THERAPY == ""))) # n= 179 (LUAD) n=201 (LUSC)
datafileObject_stageII_III <- datafileObject_stageII_III[-c(which(datafileObject_stageII_III$TARGETED_MOLECULAR_THERAPY == "YES"),which(datafileObject_stageII_III$TARGETED_MOLECULAR_THERAPY == "")),]
# n=27 (LUAD), n=45 (LUSC)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Use only stage II (possible chemotherapy) samples
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
datafileObject_stageII_III <- datafileObject[-c(which(datafileObject$AJCC_PATHOLOGIC_TUMOR_STAGE == "Stage IV"),
                                                which(datafileObject$AJCC_PATHOLOGIC_TUMOR_STAGE == "Stage I"),
                                                which(datafileObject$AJCC_PATHOLOGIC_TUMOR_STAGE == "Stage IA"),
                                                which(datafileObject$AJCC_PATHOLOGIC_TUMOR_STAGE == "Stage IB"),
                                                which(datafileObject$AJCC_PATHOLOGIC_TUMOR_STAGE == "Stage IIIA"),
                                                which(datafileObject$AJCC_PATHOLOGIC_TUMOR_STAGE == "Stage IIIB")),]

length(which(datafileObject$AJCC_PATHOLOGIC_TUMOR_STAGE == ""))  # n=9 (LUAD), number of NA for the stage, n= 5 (LUSC)
datafileObject_stageII_III <- datafileObject_stageII_III[-c(which(datafileObject_stageII_III$AJCC_PATHOLOGIC_TUMOR_STAGE == "")),] # n=122 (LUAd), n= 165 (LUSC)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
length(which(datafileObject$DFS_STATUS == ""))  # n=81 (LUAD), number of NA for the stage, n=126 (LUSC)
length(which(datafileObject_stageII_III$DFS_STATUS == ""))  # n=17, 40, number of NA for the stage, n=53, 57 (LUSC)
# LUAD
tmp_all <- datafileObject_stageII_III[,c(19,20,51,52,77,82,83)] # "DFS_MONTHS","DFS_STATUS","OS_MONTHS","OS_STATUS","TREATMENT_OUTCOME_FIRST_COURSE","DDRD.score","DDRD.call"       
# LUSC
tmp_all <- datafileObject_stageII_III[,c(19,20,50,51,75,80,81)] # "DFS_MONTHS","DFS_STATUS","OS_MONTHS","OS_STATUS","TREATMENT_OUTCOME_FIRST_COURSE","DDRD.score","DDRD.call"       
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# only those which have treatment field values
length(which(tmp_all$TREATMENT_OUTCOME_FIRST_COURSE == "")) # n=85, 171, 151 (LUAD), n= 174, 179 (LUSC)
tmp_all_NATratment <- tmp_all[-c(which(tmp_all$TREATMENT_OUTCOME_FIRST_COURSE == "")),]
#tmp_all_NATratment <- tmp_all
colnames(tmp_all_NATratment)[7] <- "DDRD.call" # changing the colname DDRD.quartile to DDRD.call

write.csv(tmp_all_NATratment, file="OS_DFS_Treatment_4DDRDSubgroups_LUSC_StageII_III_070917.csv")
write.csv(tmp_all_NATratment, file="OS_DFS_Treatment_4DDRDSubgroups_LUSC_StageII_III_IV_080917.csv")

write.csv(datafileObject_stageII_III, file="ObjectFile_4DDRDSubgroups_LUSC_StageII_III_080917.csv")
write.csv(datafileObject_stageII_III, file="ObjectFile_4DDRDSubgroups_LUSC_StageII_III_IV_080917.csv")

write.csv(datafileObject_stageII_III, file="ObjectFile_4DDRDSubgroups_LUAD_StageII_III_IV_080917.csv")
write.csv(datafileObject_stageII_III, file="ObjectFile_4DDRDSubgroups_LUAD_StageII_III_080917.csv")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# length(which(tmp_all$TREATMENT_OUTCOME_FIRST_COURSE == "Complete Remission/Response")) # n= (LUAD), n=62  (LUSC)
# length(which(tmp_all$TREATMENT_OUTCOME_FIRST_COURSE == "Progressive Disease")) # n= (LUAD), n=4  (LUSC)
# length(which(tmp_all$TREATMENT_OUTCOME_FIRST_COURSE == "Stable Disease")) # n= (LUAD), n=6  (LUSC)
# tmp_all_NATratment <- tmp_all[-c(which(tmp_all$TREATMENT_OUTCOME_FIRST_COURSE == "Progressive Disease"),which(tmp_all$TREATMENT_OUTCOME_FIRST_COURSE == "Stable Disease"),which(tmp_all$TREATMENT_OUTCOME_FIRST_COURSE == "")),]
# colnames(tmp_all_NATratment)[7] <- "DDRD.call" # changing the colname DDRD.quartile to DDRD.call
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# OS 
tmp <- tmp_all_NATratment[,c(3,4,7)]
#tmp <- tmp_all[,c(3,4,7)] # only "OS_MONTHS","OS_STATUS","DDRD.call"
length(which(is.na(tmp$OS_MONTHS))) #n=5 (LUAD), number of NA, n=1 (LUSC)
# removing NA 
tmp <- tmp[-c(which(is.na(tmp$OS_MONTHS))),]
# -----------
tmp$OS_STATUS <- ifelse(tmp$OS_STATUS == "LIVING", 0, ifelse(tmp$OS_STATUS == "DECEASED", 1, NA))
tmp$OS_MONTHS <- tmp$OS_MONTHS/12

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# DFS 
tmp <- tmp_all_NATratment[,c(1,2,7)]
length(which(is.na(tmp$DFS_MONTHS))) #n=8 (LUAD), n=7,8 (LUSC)
# removing NA 
tmp <- tmp[-c(which(is.na(tmp$DFS_MONTHS))),]
# -----------
tmp$DFS_STATUS <- ifelse(tmp$DFS_STATUS == "DiseaseFree", 0, ifelse(tmp$DFS_STATUS == "Recurred/Progressed",1, NA))
tmp$DFS_MONTHS <- tmp$DFS_MONTHS/12
#######################################################################
# Only two groups: "12" and "34"
idx_call_12 <- c(which(tmp$DDRD.call == "1"),which(tmp$DDRD.call == "2"))
idx_call_34 <- c(which(tmp$DDRD.call == "3"),which(tmp$DDRD.call == "4"))
tmp$DDRD.call[idx_call_12] <- "12"
tmp$DDRD.call[idx_call_34] <- "34"
#######################################################################
# Only two groups: "2" and "134"
idx_call_134 <- c(which(tmp$DDRD.call == "1"),which(tmp$DDRD.call == "3"), which(tmp$DDRD.call == "4"))
tmp$DDRD.call[idx_call_134] <- "134"
#######################################################################
# Only two groups: "123" and "4"
idx_call_123 <- c(which(tmp$DDRD.call == "1"),which(tmp$DDRD.call == "2"), which(tmp$DDRD.call == "3"))
tmp$DDRD.call[idx_call_123] <- "123"
#######################################################################
# Only two groups: "124" and "3"
idx_call_123 <- c(which(tmp$DDRD.call == "1"),which(tmp$DDRD.call == "2"), which(tmp$DDRD.call == "4"))
tmp$DDRD.call[idx_call_123] <- "124"
#######################################################################
# Only two groups: "14" and "23"
idx_call_14 <- c(which(tmp$DDRD.call == "1"),which(tmp$DDRD.call == "4"))
idx_call_23 <- c(which(tmp$DDRD.call == "2"),which(tmp$DDRD.call == "3"))
tmp$DDRD.call[idx_call_14] <- "14"
tmp$DDRD.call[idx_call_23] <- "23"
#######################################################################
# Only three groups: "12", "3" and "4"
idx_call_12 <- c(which(tmp$DDRD.call == "1"),which(tmp$DDRD.call == "2"))
tmp$DDRD.call[idx_call_12] <- "12"
#######################################################################

### OS
f <- npsurv(Surv(OS_MONTHS, OS_STATUS) ~ DDRD.call , data=tmp)  # OS
survdiff(Surv(OS_MONTHS, OS_STATUS) ~ DDRD.call , data=tmp, rho = 0)  # OS
cox1 <- coxph(Surv(OS_MONTHS, OS_STATUS) ~ as.character(DDRD.call), method = "exact", data = tmp) # OS
print(summary(cox1))
### DFS
f <- npsurv(Surv(DFS_MONTHS, DFS_STATUS) ~ DDRD.call, data=tmp) # DFS
survdiff(Surv(DFS_MONTHS, DFS_STATUS) ~ DDRD.call, data=tmp, rho = 0) # DFS
cox1 <- coxph(Surv(DFS_MONTHS, DFS_STATUS) ~ as.character(DDRD.call), method = "exact", data = tmp) # DFS use as.character
print(summary(cox1))
#######################################################################
dev.off()
No_of_samples <- nrow(tmp)
# Colour transparent and then add lines, with censor points.
survplot(f, conf="none",  n.risk=F, lty=1, lwd=2, col="#FFFFFF00", 
         xlab="Time (years)", ylab="Overall survival", label.curves=T,  time.inc=2)

survplot(f, conf="none",  n.risk=F, lty=1, lwd=2, col="#FFFFFF00", 
         xlab="Time (years)", ylab="Disease-free survival", label.curves=T,  time.inc=2)

lines(f, mark.time=T, lwd=2, col=c("blue", "red","orange","darkgreen")) # DDRD.call=1 (blue), DDRD.call=2 (red), DDRD.call=3 (orange), DDRD.call=4 (darkgreen)
leg.txt <- c("DDRD.call=1", "DDRD.call=2","DDRD.call=3","DDRD.call=4") 
legend(13, 1.03, leg.txt, lty = 1, lwd = 2, col = c("blue", "red","orange","darkgreen")) # LUAD, CIMP-high
title(paste("LUAD, stage II & III, n=",No_of_samples,sep = "")) # LUAD
title(paste("LUAD, stage II, n=",No_of_samples,sep = "")) # LUAD
title(paste("LUAD, stage II & III, no targeted molecular therapy, n=",No_of_samples,sep = "")) # LUSC
title(paste("LUSC, stage II & III, no targeted molecular therapy, n=",No_of_samples,sep = "")) # LUSC
title(paste("LUSC, stage II, n=",No_of_samples,sep = "")) # LUSC
title(paste("LUSC, stage II, III & IV, n=",No_of_samples,sep = "")) # LUSC

# when using 2 subgroups (12 vs. 34)
lines(f, mark.time=T, lwd=2, col=c("purple","green")) # DDRD.call=12 (purple), DDRD.call=34 (green)
leg.txt <- c("DDRD.call=12", "DDRD.call=34") 
legend(12, 0.9, leg.txt, lty = 1, lwd = 2, col = c("purple","green"))
title(paste("LUAD, stage II & III, n=",No_of_samples,sep = "")) # LUAD
#*****************************
# when using 2 subgroups (2 vs. 134)
lines(f, mark.time=T, lwd=2, col=c("darkolivegreen2","red")) # DDRD.call=12 (purple), DDRD.call=34 (green)
leg.txt <- c("DDRD.call=134", "DDRD.call=2") 
legend(12, 0.4, leg.txt, lty = 1, lwd = 2, col = c("darkolivegreen2","red")) # LUAD
title(paste("LUAD, stage II & III, n=",No_of_samples,sep = "")) # LUAD
title(paste("LUAD, stage II, n=",No_of_samples,sep = "")) # LUAD
#*****************************
# when using 2 subgroups (123 vs. 4)
lines(f, mark.time=T, lwd=2, col=c("black","darkgreen")) # DDRD.call=123 (black), DDRD.call=4 (darkgreen)
leg.txt <- c("DDRD.call=123", "DDRD.call=4") 
legend(8.5, 1.03, leg.txt, lty = 1, lwd = 2, col = c("black","darkgreen")) # LUSC
title(paste("LUSC, stage II & III, n=",No_of_samples,sep = "")) # LUSC
title(paste("LUSC, stage II, n=",No_of_samples,sep = "")) # LUSC
title(paste("LUSC, stage II, III & IV, n=",No_of_samples,sep = "")) # LUSC
#*****************************
# when using 2 subgroups (124 vs. 3)
lines(f, mark.time=T, lwd=2, col=c("grey","orange")) # DDRD.call=124 (grey)
leg.txt <- c("DDRD.call=124", "DDRD.call=3") 
legend(8.5, 1.03, leg.txt, lty = 1, lwd = 2, col = c("grey","orange")) # LUSC
title(paste("LUSC, stage II & III, n=",No_of_samples,sep = "")) # LUSC
title(paste("LUSC, stage II, n=",No_of_samples,sep = "")) # LUSC
title(paste("LUSC, stage II, III & IV, n=",No_of_samples,sep = "")) # LUSC
#*****************************
# when using 2 subgroups (14 vs. 23)
lines(f, mark.time=T, lwd=2, col=c("cyan","brown1")) # DDRD.call=12 (purple), DDRD.call=34 (green)
leg.txt <- c("DDRD.call=14", "DDRD.call=23") 
legend(9, 1.03, leg.txt, lty = 1, lwd = 2, col = c("cyan","brown1"))
title(paste("LUSC, stage II & III, n=",No_of_samples,sep = "")) # LUSC
title(paste("LUSC, stage II, III & IV, n=",No_of_samples,sep = "")) # LUSC
#*****************************
# when using 3 subgroups (12, 3 and 4)
lines(f, mark.time=T, lwd=2, col=c("purple","orange","darkgreen")) # DDRD.call=12 (purple), DDRD.call=34 (green)
leg.txt <- c("DDRD.call=12", "DDRD.call=3", "DDRD.call=4") 
legend(9, 1.03, leg.txt, lty = 1, lwd = 2, col = c("purple","orange","darkgreen"))
title(paste("LUSC, stage II & III, no targeted molecular therapy, n=",No_of_samples,sep = "")) # LUSC
title(paste("LUSC, stage II & III, n=",No_of_samples,sep = "")) # LUSC
title(paste("LUSC, stage II, III & IV, n=",No_of_samples,sep = "")) # LUSC
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#length(which(datafileObject_stageII_III$TREATMENT_OUTCOME_FIRST_COURSE == "")) # n=151 out of 206
# LUSC_StageII_III_50_OS_4DDRDCall_p00884
# coxph(formula = Surv(OS_MONTHS, OS_STATUS) ~ DDRD.call, data = tmp, 
#       method = "exact")
# 
# n= 50, number of events= 12 
# 
# coef exp(coef) se(coef)     z Pr(>|z|)  
# DDRD.call 0.5411    1.7179   0.3287 1.646   0.0997 .
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# exp(coef) exp(-coef) lower .95 upper .95
# DDRD.call     1.718     0.5821     0.902     3.272
# 
# Rsquare= 0.06   (max possible= 0.806 )
# Likelihood ratio test= 3.1  on 1 df,   p=0.07808
# Wald test            = 2.71  on 1 df,   p=0.09973
# Score (logrank) test = 2.9  on 1 df,   p=0.08844
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# LUSC_StageII_III_50_OS_123vs4DDRDCall_p00228
# coxph(formula = Surv(OS_MONTHS, OS_STATUS) ~ DDRD.call, data = tmp, 
#       method = "exact")
# 
# n= 50, number of events= 12 
# 
# coef exp(coef) se(coef)     z Pr(>|z|)  
# DDRD.call4 1.3060    3.6913   0.6144 2.125   0.0335 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# exp(coef) exp(-coef) lower .95 upper .95
# DDRD.call4     3.691     0.2709     1.107     12.31
# 
# Rsquare= 0.093   (max possible= 0.806 )
# Likelihood ratio test= 4.88  on 1 df,   p=0.02709
# Wald test            = 4.52  on 1 df,   p=0.03355
# Score (logrank) test = 5.18  on 1 df,   p=0.02282
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# LUSC_StageII_III_47_DFS_4DDRDCall_p00534
# Call:survdiff(Surv(DFS_MONTHS, DFS_STATUS) ~ DDRD.call, data=tmp, rho = 0) # DFS
#   survdiff(formula = Surv(DFS_MONTHS, DFS_STATUS) ~ DDRD.call, 
#            data = tmp, rho = 0)
# 
# N Observed Expected (O-E)^2/E (O-E)^2/V
# DDRD.call=1  7        1     1.76     0.328     0.373
# DDRD.call=2 14        4     5.52     0.419     0.645
# DDRD.call=3  7        1     3.76     2.028     2.727
# DDRD.call=4 19       10     4.96     5.133     7.504
# 
# Chisq= 8  on 3 degrees of freedom, p= 0.0457  
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# coxph(formula = Surv(DFS_MONTHS, DFS_STATUS) ~ DDRD.call, data = tmp, 
#       method = "exact")
# 
# n= 47, number of events= 16 
# 
# coef exp(coef) se(coef)     z Pr(>|z|)  
# DDRD.call 0.5362    1.7095   0.2776 1.932   0.0534 .
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# exp(coef) exp(-coef) lower .95 upper .95
# DDRD.call     1.709      0.585    0.9922     2.945
# 
# Rsquare= 0.085   (max possible= 0.901 )
# Likelihood ratio test= 4.16  on 1 df,   p=0.04127
# Wald test            = 3.73  on 1 df,   p=0.05341
# Score (logrank) test = 3.98  on 1 df,   p=0.04615
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# LUSC_StageII_III_47_DFS_123vs4DDRDCall_p00107
# Call:
#   survdiff(formula = Surv(DFS_MONTHS, DFS_STATUS) ~ DDRD.call, 
#            data = tmp, rho = 0)
# 
# N Observed Expected (O-E)^2/E (O-E)^2/V
# DDRD.call=123 28        6    11.04      2.30       7.5
# DDRD.call=4   19       10     4.96      5.13       7.5
# 
# Chisq= 7.5  on 1 degrees of freedom, p= 0.00616 
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# Call:
#   coxph(formula = Surv(DFS_MONTHS, DFS_STATUS) ~ DDRD.call, data = tmp, 
#         method = "exact")
# 
# n= 47, number of events= 16 
# 
# coef exp(coef) se(coef)     z Pr(>|z|)  
# DDRD.call4 1.3261    3.7662   0.5192 2.554   0.0107 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# exp(coef) exp(-coef) lower .95 upper .95
# DDRD.call4     3.766     0.2655     1.361     10.42
# 
# Rsquare= 0.134   (max possible= 0.901 )
# Likelihood ratio test= 6.79  on 1 df,   p=0.009182
# Wald test            = 6.52  on 1 df,   p=0.01065
# Score (logrank) test = 7.5  on 1 df,   p=0.006157
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# LUSC_StageII_50_OS_4DDRDCall_p00884
# Call:
#   survdiff(formula = Surv(OS_MONTHS, OS_STATUS) ~ DDRD.call, data = tmp, 
#            rho = 0)
# 
# N Observed Expected (O-E)^2/E (O-E)^2/V
# DDRD.call=1  8        1     1.22    0.0396    0.0449
# DDRD.call=2 13        2     3.85    0.8878    1.3186
# DDRD.call=3  8        1     2.69    1.0615    1.4020
# DDRD.call=4 21        8     4.24    3.3291    5.1825
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
# Chisq= 5.4  on 3 degrees of freedom, p= 0.146 
# Call:
#   coxph(formula = Surv(OS_MONTHS, OS_STATUS) ~ DDRD.call, data = tmp, 
#         method = "exact")
# 
# n= 50, number of events= 12 
# 
# coef exp(coef) se(coef)     z Pr(>|z|)  
# DDRD.call 0.5411    1.7179   0.3287 1.646   0.0997 .
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# exp(coef) exp(-coef) lower .95 upper .95
# DDRD.call     1.718     0.5821     0.902     3.272
# 
# Rsquare= 0.06   (max possible= 0.806 )
# Likelihood ratio test= 3.1  on 1 df,   p=0.07808
# Wald test            = 2.71  on 1 df,   p=0.09973
# Score (logrank) test = 2.9  on 1 df,   p=0.08844
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# LUSC_StageII_50_OS_123vs4DDRDCall_p00335
# survdiff(formula = Surv(OS_MONTHS, OS_STATUS) ~ DDRD.call, data = tmp, 
#          rho = 0)
# 
# N Observed Expected (O-E)^2/E (O-E)^2/V
# DDRD.call=123 29        4     7.76      1.82      5.18
# DDRD.call=4   21        8     4.24      3.33      5.18
# 
# Chisq= 5.2  on 1 degrees of freedom, p= 0.0228 
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# Call:
#   coxph(formula = Surv(OS_MONTHS, OS_STATUS) ~ DDRD.call, data = tmp, 
#         method = "exact")
# 
# n= 50, number of events= 12 
# 
# coef exp(coef) se(coef)     z Pr(>|z|)  
# DDRD.call4 1.3060    3.6913   0.6144 2.125   0.0335 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# exp(coef) exp(-coef) lower .95 upper .95
# DDRD.call4     3.691     0.2709     1.107     12.31
# 
# Rsquare= 0.093   (max possible= 0.806 )
# Likelihood ratio test= 4.88  on 1 df,   p=0.02709
# Wald test            = 4.52  on 1 df,   p=0.03355
# Score (logrank) test = 5.18  on 1 df,   p=0.02282
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# LUSC_StageII_47_DFS_4DDRDCall_p00534
# Call:
#   survdiff(formula = Surv(DFS_MONTHS, DFS_STATUS) ~ DDRD.call, 
#            data = tmp, rho = 0)
# 
# N Observed Expected (O-E)^2/E (O-E)^2/V
# DDRD.call=1  7        1     1.76     0.328     0.373
# DDRD.call=2 14        4     5.52     0.419     0.645
# DDRD.call=3  7        1     3.76     2.028     2.727
# DDRD.call=4 19       10     4.96     5.133     7.504
# 
# Chisq= 8  on 3 degrees of freedom, p= 0.0457 
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# Call:
#   coxph(formula = Surv(DFS_MONTHS, DFS_STATUS) ~ DDRD.call, data = tmp, 
#         method = "exact")
# 
# n= 47, number of events= 16 
# 
# coef exp(coef) se(coef)     z Pr(>|z|)  
# DDRD.call 0.5362    1.7095   0.2776 1.932   0.0534 .
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# exp(coef) exp(-coef) lower .95 upper .95
# DDRD.call     1.709      0.585    0.9922     2.945
# 
# Rsquare= 0.085   (max possible= 0.901 )
# Likelihood ratio test= 4.16  on 1 df,   p=0.04127
# Wald test            = 3.73  on 1 df,   p=0.05341
# Score (logrank) test = 3.98  on 1 df,   p=0.04615
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# LUSC_StageII_47_DFS_123vs4DDRDCall_p00107
# Call:
#   survdiff(formula = Surv(DFS_MONTHS, DFS_STATUS) ~ DDRD.call, 
#            data = tmp, rho = 0)
# 
# N Observed Expected (O-E)^2/E (O-E)^2/V
# DDRD.call=123 28        6    11.04      2.30       7.5
# DDRD.call=4   19       10     4.96      5.13       7.5
# 
# Chisq= 7.5  on 1 degrees of freedom, p= 0.00616 
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# Call:
#   coxph(formula = Surv(DFS_MONTHS, DFS_STATUS) ~ DDRD.call, data = tmp, 
#         method = "exact")
# 
# n= 47, number of events= 16 
# 
# coef exp(coef) se(coef)     z Pr(>|z|)  
# DDRD.call4 1.3261    3.7662   0.5192 2.554   0.0107 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# exp(coef) exp(-coef) lower .95 upper .95
# DDRD.call4     3.766     0.2655     1.361     10.42
# 
# Rsquare= 0.134   (max possible= 0.901 )
# Likelihood ratio test= 6.79  on 1 df,   p=0.009182
# Wald test            = 6.52  on 1 df,   p=0.01065
# Score (logrank) test = 7.5  on 1 df,   p=0.006157
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
