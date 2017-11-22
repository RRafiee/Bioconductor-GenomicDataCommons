##################################################################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##################################################################################################
# Start
#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
# Dr Reza Rafiee, August 2017
# Research Fellow, Queen's University Belfast
# Survival analysis (including OS and DFS) on lung cancer data (LUAD and LUSC)
# version 2.3, 17th November 2017
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
##############################   LUAD/LUSC  ################################
# Loading a csv file
# Pathfolder <- "C:/Users/3052092/Documents/Live/"   # initialise the path of the csv file 
Pathfolder <- "/home/reza/Documents/Live/"           # initialise the path of the csv file (from my Linux) 

#16/11/2017

csvfilename <- "LUSC_Stage_II_III_Samples_Ordered_192_Chemo_161117.csv" # LUSC, 4DDRDs, stage II& III, Chemo, with/without alterations
pathcsvfile <- paste(Pathfolder,csvfilename,sep = "")
datafileObject <- as.data.frame(read.csv(pathcsvfile, header=T, row.names = 1,stringsAsFactors=FALSE))  # LUAD samples, TCGA

# csvfilename <- "Statistics_Data_Drugs_LUAD_Chemo_221117.csv"
# pathcsvfile <- paste(Pathfolder,csvfilename,sep = "")
# datafileObject <- as.data.frame(read.csv(pathcsvfile, header=T, stringsAsFactors=FALSE))  # LUAD samples, TCGA with Chemo data


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 16/11/2017
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Use only received chemo samples
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
datafileObject_stageII_III <- datafileObject[which(datafileObject$Chemotherapy != "No Chemo"),]  # petients who received chemotherapy
datafileObject_stageII_III <- datafileObject_stageII_III[-c(grep("Radiation",as.character(datafileObject_stageII_III$Chemotherapy))),]
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#tmp_all_NATratment <- tmp_all_NATratment[which(tmp_all_NATratment$DDRD.call == "4"),] # only group 4

tmp_all_NATratment_withAlterations <- datafileObject_stageII_III[which(datafileObject_stageII_III$BRCA_FA.Alterations == "1"),] # with alterations, n=35
tmp_all_NATratment <- tmp_all_NATratment_withAlterations
tmp_all_NATratment_withoutAlterations <- datafileObject_stageII_III[which(datafileObject_stageII_III$BRCA_FA.Alterations == "0"),] # without alterations, n=37
tmp_all_NATratment <- tmp_all_NATratment_withoutAlterations
# 
# # only BRCA1 alterations, 13 October 2017
# tmp_all_NATratment_withAlterations <- datafileObject[which(datafileObject$BRCA1.Alterations == "1"),] # with alterations, n=35
# tmp_all_NATratment <- tmp_all_NATratment_withAlterations
# tmp_all_NATratment_withoutAlterations <- datafileObject[which(datafileObject$BRCA1.Alterations == "0"),] # without alterations, n=37
# tmp_all_NATratment <- tmp_all_NATratment_withoutAlterations

datafileObject_stageII_III$OS_MONTHS <- as.numeric(as.character(datafileObject_stageII_III$OS_MONTHS))
tmp_all_NATratment <- datafileObject_stageII_III # 16/11/2017

write.csv(tmp_all_NATratment, file="OS_DFS_Treatment_4DDRDSubgroups_LUSC_StageII_III_ChemoTherapy_Prescribed_171117.csv")

#-------------------------------------------------------------
tmp <- tmp_all_NATratment[,c(4,5,8)] #"OS_MONTHS" "OS_STATUS" "DDRD.call":8/"BRCA_FA.Alterations":10
tmp <- tmp_all_NATratment[,c(2,3,8)] #"DFS_MONTHS" "DFS_STATUS" "DDRD.call":8/"BRCA_FA.Alterations":10
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# OS 
tmp <- tmp_all_NATratment[,c(4,5,8)] # 16/11/2017
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
length(which(is.na(tmp$DFS_MONTHS))) #n=8 (LUAD), n=7,8 (LUSC), 21 (LUSC_Untreated_NatPaper)
# removing NA 
tmp <- tmp[-c(which(is.na(tmp$DFS_MONTHS))),]
# -----------
tmp$DFS_STATUS <- ifelse(tmp$DFS_STATUS == "DiseaseFree", 0, ifelse(tmp$DFS_STATUS == "Recurred/Progressed",1, NA))
tmp$DFS_MONTHS <- tmp$DFS_MONTHS/12
# #######################################################################
# Hazard Ratio
# library(survival)
# data.survdiff <- survdiff(Surv(time, status) ~ group)
# p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
# HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
# up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
# low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
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
legend(3.2, 0.35, leg.txt, lty = 1, lwd = 2, col = c("blue", "red","orange","darkgreen")) # LUAD, CIMP-high
title(paste("LUSC, stage II & III, treated, n=",No_of_samples,sep = "")) # LUSC, 

##################################################################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##################################################################################################
Pathfolder <- "/home/reza/Documents/Live/"           # initialise the path of the csv file (from my Linux) 
csvfilename <- "Statistics_Data_Drugs_LUAD_Chemo_221117.csv"
pathcsvfile <- paste(Pathfolder,csvfilename,sep = "")
datafileObject <- as.data.frame(read.csv(pathcsvfile, header=T, stringsAsFactors=FALSE))  # LUAD samples, TCGA with Chemo data
#........................................................................................
LUAD_Chemo_Received_Without_Duplicate <- datafileObject[!duplicated(datafileObject$Bcr.Patient.Barcode),] # excluding duplicated samples 
# Adding "-01" to sample ID
LUAD_Chemo_Received_Without_Duplicate$Bcr.Patient.Barcode <- paste(LUAD_Chemo_Received_Without_Duplicate$Bcr.Patient.Barcode,"-01",sep = "")
# removing only Radiation as well as "No drug prescribed"
LUAD_Chemo_Received_Without_Duplicate <- LUAD_Chemo_Received_Without_Duplicate[-c(which(LUAD_Chemo_Received_Without_Duplicate$Therapy.Drug.Names == "No drug prescribed"), which(LUAD_Chemo_Received_Without_Duplicate$Therapy.Drug.Names == "Radiation")),] 
# Replace "-" with "."
LUAD_Chemo_Received_Without_Duplicate$Bcr.Patient.Barcode <- gsub("-", "\\.",LUAD_Chemo_Received_Without_Duplicate$Bcr.Patient.Barcode)
# Match the samples with main dataset for survival data
csvfilename <- "LUAD_Samples_DDRD_ordered_23June17.csv" # LUAD
pathcsvfile <- paste(Pathfolder,csvfilename,sep = "")
datafileObject_LUAD <- read.csv(pathcsvfile, header=T, row.names = 1)  # for LUAD or LUSC
#........................................................................................
LUAD_Chemo_Received_Without_Duplicate <- LUAD_Chemo_Received_Without_Duplicate[order(LUAD_Chemo_Received_Without_Duplicate$Bcr.Patient.Barcode),] 
datafileObject_LUAD <- datafileObject_LUAD[order(datafileObject_LUAD$Sample_ID),]

#LUAD_Chemo_Received_Without_Duplicate$Bcr.Patient.Barcode
rownames(LUAD_Chemo_Received_Without_Duplicate) <- LUAD_Chemo_Received_Without_Duplicate$Bcr.Patient.Barcode
rownames(datafileObject_LUAD) <- as.character(datafileObject_LUAD$Sample_ID)

# datafileObject_LUAD[rownames(LUAD_Chemo_Received_Without_Duplicate) %in% rownames(datafileObject_LUAD),]

datafileObject_LUAD_matched <- datafileObject_LUAD[which(rownames(datafileObject_LUAD) %in% rownames(LUAD_Chemo_Received_Without_Duplicate)),]

which(rownames(LUAD_Chemo_Received_Without_Duplicate) %in% rownames(datafileObject_LUAD_matched) == FALSE) # 26 in LUAD_Chemo_Received_Without_Duplicate need to be removed

LUAD_Chemo_Received_Without_Duplicate <- LUAD_Chemo_Received_Without_Duplicate[-26,]

LUAD_Chemo_Received_Without_Duplicate_WithClinicalInfo <- cbind(datafileObject_LUAD_matched,LUAD_Chemo_Received_Without_Duplicate)

write.csv(LUAD_Chemo_Received_Without_Duplicate_WithClinicalInfo,"LUAD_Chemo_Received_Without_Duplicate_WithClinicalInfo_152_AllStages_221117.csv")

as.character(LUAD_Chemo_Received_Without_Duplicate_WithClinicalInfo$Sample_ID)
as.character(LUAD_Chemo_Received_Without_Duplicate_WithClinicalInfo$Bcr.Patient.Barcode)
