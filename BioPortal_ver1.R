##################################################################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##################################################################################################
# Start
#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
# Dr Reza Rafiee, August 2017
# Research Fellow, Queen's University Belfast
# Survival analysis (including OS and DFS) on lung cancer data (LUAD and LUSC)
# version 1.0, September 2017
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setwd("/home/reza/Documents/Live/")
##############################   LUAD/LUSC  ################################
# Loading a csv file
# Pathfolder <- "C:/Users/3052092/Documents/Live/"   # initialise the path of the csv file 
Pathfolder <- "/home/reza/Documents/Live/"           # initialise the path of the csv file (from my Linux) 

# csvfilename <- "LUAD_Samples_DDRD_ordered_23June17.csv" # LUAD
# csvfilename <- "LUSC_Samples_DDRD_ordered_23June17.csv" # LUSC

# Using the new dataset including EMT scores, 29/09/2017
csvfilename <- "LUAD_Samples_DDRD_ordered_EMT_29Sep17.csv" # LUAD
csvfilename <- "LUSC_Samples_DDRD_ordered_EMT_29Sep17.csv" # LUSC
#
csvfilename <- "TCGA_LUSC_Clinical_Info_178_NatureMethod_withDDRDCall_13September17.csv" # LUSC, n=178 Nature Method Paper
#-----------------------------------------------------------------------------------------------------------------------
csvfilename <- "TCGA_LUAD_Clinical_Info_230_NatureMethod_withDDRDCalls_16August17.csv" # LUAD, n=230 Nature Method Paper
# #-----------------------------------------------------------------------------------------------------------------------
# csvfilename2 <- "TCGA_Lung_EMTScores_LUSC.csv" # All LUSC pateients with EMT scores, 29/09/2017
# csvfilename2 <- "TCGA_Lung_EMTScores_LUAD.csv" # All LUAD pateients with EMT scores, 29/09/2017
# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pathcsvfile <- paste(Pathfolder,csvfilename,sep = "")
datafileObject <- read.csv(pathcsvfile, header=T, row.names = 1)  # for LUAD or LUSC

# # for EMT dataset, 29/09/2017
# pathcsvfile2 <- paste(Pathfolder,csvfilename2,sep = "")
# datafileObject2 <- read.csv(pathcsvfile2, header=T, row.names = 1)  # LUSC & LUAD samples, TCGA
# rownames(datafileObject2) <- gsub("\\-", ".",rownames(datafileObject2)) # same sample ID matched with LUSC ones
# datafileObject2 <- cbind(rownames(datafileObject2),datafileObject2)
# # Order Sample IDs 
# colnames(datafileObject2)[1] <- "Sample_ID"
# datafileObject <- datafileObject[order(datafileObject$Sample_ID),]
# datafileObject2 <- datafileObject2[order(datafileObject2$Sample_ID),]
# # Add a new column for EMT
# datafileObject <- cbind(datafileObject,datafileObject2$EMT.15.Gene.Score)
# colnames(datafileObject)[81] <- "EMT.score" # LUSC
# colnames(datafileObject)[83] <- "EMT.score" # LUAD
# write.csv(datafileObject, file="LUSC_Samples_DDRD_ordered_EMT_29Sep17.csv")
# write.csv(datafileObject, file="LUAD_Samples_DDRD_ordered_EMT_29Sep17.csv")


#------------------------------------------- Only for n=230 ------------------------------------------------------------
colnames(datafileObject)[6] <- "AJCC_PATHOLOGIC_TUMOR_STAGE"
colnames(datafileObject)[1] <- "Sample_ID"
which(datafileObject$AJCC_PATHOLOGIC_TUMOR_STAGE == "[Not Available]")
datafileObject <- datafileObject[-c(which(datafileObject$AJCC_PATHOLOGIC_TUMOR_STAGE == "[Not Available]")),]
#-----------------------------------------------------------------------------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Use only stage II and stage III (possible chemotherapy) samples
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
datafileObject_stageII_III <- datafileObject[-c(which(datafileObject$AJCC_PATHOLOGIC_TUMOR_STAGE == "Stage IV"), which(datafileObject$AJCC_PATHOLOGIC_TUMOR_STAGE == "Stage I"), which(datafileObject$AJCC_PATHOLOGIC_TUMOR_STAGE == "Stage IA"),which(datafileObject$AJCC_PATHOLOGIC_TUMOR_STAGE == "Stage IB") ),]
length(which(datafileObject$AJCC_PATHOLOGIC_TUMOR_STAGE == ""))  # n=9 (LUAD), number of NA for the stage, n=5 (LUSC)

datafileObject_stageII_III <- datafileObject_stageII_III[-c(which(datafileObject_stageII_III$AJCC_PATHOLOGIC_TUMOR_STAGE == "")),]
#-----------------------------------------------------------------------------------------------------------------------
# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Use only patients who have not any kind of metastasis (it means only N0 but N1, N2, N3)
length(which(datafileObject$AJCC_NODES_PATHOLOGIC_PN == "N0")) #n=318
length(which(datafileObject$AJCC_NODES_PATHOLOGIC_PN == "N1")) #n=131
length(which(datafileObject$AJCC_NODES_PATHOLOGIC_PN == "N2")) #n=40
length(which(datafileObject$AJCC_NODES_PATHOLOGIC_PN == "N3")) #n=5


# 09/10/2017
###########################################################################
# Choosing patients who received chemotherapy
#tmp_all <- datafileObject_stageII_III

length(which(tmp_all$TREATMENT_OUTCOME_FIRST_COURSE == "")) # n=85, 171, 151 (LUAD), n= 174, 179 (LUSC)

if (length(which(tmp_all$TREATMENT_OUTCOME_FIRST_COURSE == "")) == 0)
{
  tmp_all_NATratment <- tmp_all
} else 
{
  tmp_all_NATratment <- tmp_all[-c(which(tmp_all$TREATMENT_OUTCOME_FIRST_COURSE == "")),]  
}
# only Group 4 DDRD
tmp_all_NATratment <- tmp_all_NATratment[which(tmp_all_NATratment$DDRD.quartile == "4"),]

###########################################################################

#tmp_all <- datafileObject_stageII_III[,c(18,19,50,51,76,81,82)] # For LUAD,"DFS_MONTHS","DFS_STATUS","OS_MONTHS","OS_STATUS","TREATMENT_OUTCOME_FIRST_COURSE","DDRD.score","DDRD.call"       

tmp_all <- datafileObject_stageII_III[,c(1,18,19,50,51,76,81,82,83)] # For LUAD with EMT scores # [,c(1,18,19,50,51,76,81,82,83)]
tmp_all <- datafileObject_stageII_III[,c(1,18,19,49,50,74,79,80,81)] # For LUSC with EMT scores # [,c(1,18,19,49,50,74,79,80)]

# LUAD, n=230
tmp_all <- datafileObject_stageII_III

# if treated then next line
 tmp_all_NATratment <- tmp_all[-c(which(tmp_all$TREATMENT_OUTCOME_FIRST_COURSE == "")),] 
# else (untreated)
 tmp_all_NATratment <- tmp_all # when using untreated (n=178)
# only for LUSC
colnames(tmp_all_NATratment)[8] <- "DDRD.call" # changing the colname DDRD.quartile to DDRD.call

# order based on DDRD scores
tmp_all_NATratment <- tmp_all_NATratment[order(tmp_all_NATratment$DDRD.score),]
# order based on DDRD call
tmp_all_NATratment <- tmp_all_NATratment[order(tmp_all_NATratment$DDRD.call),]
# order based on EMT scores
tmp_all_NATratment <- tmp_all_NATratment[order(tmp_all_NATratment$EMT.score),]
# order based on Sample_ID
tmp_all_NATratment <- tmp_all_NATratment[order(tmp_all_NATratment$Sample_ID),]


# Using "-" instead of "." in sample ID
tmp_all_NATratment$Sample_ID <- gsub("\\.", "-",tmp_all_NATratment$Sample_ID)
#tmp_all_NATratment$Sample_ID <- gsub("-01", "",tmp_all_NATratment$Sample_ID)

write.csv(tmp_all_NATratment, file="LUAD_4DDRDSubgroups_StageII_III_OrderbyDDRDScore_withEMT_206_101017.csv")
write.csv(tmp_all_NATratment, file="LUAD_4DDRDSubgroups_StageII_III_OrderbySample_ID_withEMT_206_101017.csv")



write.csv(tmp_all_NATratment, file="LUSC_4DDRDSubgroups_StageII_III_Treated_150917.csv")
write.csv(tmp_all_NATratment, file="LUAD_4DDRDSubgroups_StageII_III_Treated_200917.csv")
write.csv(tmp_all_NATratment, file="LUSC_4DDRDSubgroups_StageII_III_Untreated_250917.csv")
write.csv(tmp_all_NATratment, file="LUAD_4DDRDSubgroups_StageII_III_Untreated_260917.csv")
write.csv(tmp_all_NATratment, file="LUSC_4DDRDSubgroups_StageII_III_270917.csv") # 27/09/2017

# ordered based on EMT scores
write.csv(tmp_all_NATratment, file="LUSC_4DDRDSubgroups_StageII_III_Treated_EMTOrdered_290917.csv")
# order based on DDRD scores
write.csv(tmp_all_NATratment, file="LUSC_4DDRDSubgroups_StageII_III_Treated_DDRDScoreOrdered_290917.csv")

length(which(tmp_all_NATratment$DDRD.call == "4"))
length(which(tmp_all_NATratment$DDRD.call == "3"))
length(which(tmp_all_NATratment$DDRD.call == "2"))
length(which(tmp_all_NATratment$DDRD.call == "1"))

tmp_all_NATratment_DDRD34 <- tmp_all_NATratment[c(which(tmp_all_NATratment$DDRD.call == "3"),which(tmp_all_NATratment$DDRD.call == "4")),]
tmp_all_NATratment_DDRD12 <- tmp_all_NATratment[c(which(tmp_all_NATratment$DDRD.call == "1"),which(tmp_all_NATratment$DDRD.call == "2")),]

tmp_all_NATratment_DDRD34 <- tmp_all_NATratment_DDRD34[order(tmp_all_NATratment_DDRD34$EMT.score),]
tmp_all_NATratment_DDRD12 <- tmp_all_NATratment_DDRD12[order(tmp_all_NATratment_DDRD12$EMT.score),]

# ordered based on EMT scores
write.csv(tmp_all_NATratment_DDRD34, file="LUSC_34DDRDSubgroups_StageII_III_Treated_EMTOrdered_290917.csv")
write.csv(tmp_all_NATratment_DDRD12, file="LUSC_12DDRDSubgroups_StageII_III_Treated_EMTOrdered_290917.csv")

#boxplot(tmp_all_NATratment$DDRD.score,tmp_all_NATratment$EMT.score, main= "", ylab="EMT score",col=c("red","pink"),las=1,notch = FALSE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plotting the DDRD score of 4 subgroups
#library(ggplot2)
x <- 1:length(tmp_all_NATratment$DDRD.call)
f_x <- tmp_all_NATratment$DDRD.score
plot(x, f_x , xlim=range(x), ylim=range(f_x), xlab="Samples", ylab="DDRD score", 
     main = "DDRD score",pch=19, bg= "green" ,lwd=2, cex=0.5)

DDRD_1 <- which(tmp_all_NATratment$DDRD.call == "1")
DDRD_2 <- which(tmp_all_NATratment$DDRD.call == "2")
DDRD_3 <- which(tmp_all_NATratment$DDRD.call == "3")
DDRD_4 <- which(tmp_all_NATratment$DDRD.call == "4")

#DDRD_4[1]:DDRD_4[length(DDRD_4)]

lines(DDRD_1[1]:DDRD_1[length(DDRD_1)], f_x[DDRD_1], xlim=range(x), ylim=range(f_x), pch=19, col="blue", bg="blue", lwd=1)
lines(DDRD_2[1]:DDRD_2[length(DDRD_2)], f_x[DDRD_2], xlim=range(x), ylim=range(f_x), pch=19, col="red", bg="red", lwd=1)
lines(DDRD_3[1]:DDRD_3[length(DDRD_3)], f_x[DDRD_3], xlim=range(x), ylim=range(f_x), pch=19, col="orange", bg= "orange", lwd=1)
lines(DDRD_4[1]:DDRD_4[length(DDRD_4)], f_x[DDRD_4], xlim=range(x), ylim=range(f_x), pch=19, col="darkgreen", bg="darkgreen", lwd=1)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
DDRD_Assay_Genelist_44 <- c("CXCL10","MX1","IDO1","IFI44L","CD2","GBP5","PRAME","ITGAL","LRP4","APOL3","CDR1","FYB","TSPAN7",
                         "RAC2","KLHDC7B","GRB14","AC138128.1","KIF26A","CD274","CD109","ETV7","MFAP5","OLFM4","PI15",
                         "FOSB","FAM19A5","NLRC5","PRICKLE1","EGR1","CLDN10","ADAMTS4","SP140L","ANXA1","RSAD2","ESR1",
                         "IKZF3","OR2I1P","EGFR","NAT1","LATS2","CYP2B6","PTPRC","PPP1R1A","AL137218.1")
DDRD_Assay_Genelist_41 <- c("CXCL10","MX1","IDO1","IFI44L","CD2","GBP5","PRAME","ITGAL","LRP4","APOL3","CDR1","FYB","TSPAN7","RAC2",    
                          "KLHDC7B","GRB14","KIF26A","CD274","CD109","ETV7","MFAP5","OLFM4","PI15","FOSB","FAM19A5","NLRC5","PRICKLE1","EGR1",    
                          "CLDN10","ADAMTS4","SP140L","ANXA1","RSAD2","ESR1","IKZF3","EGFR","NAT1","LATS2","CYP2B6","PTPRC","PPP1R1A") 
DDRD_Assay_Genelist_All <- cbind(DDRD_Assay_Genelist_44,DDRD_Assay_Genelist_41)
write.csv(DDRD_Assay_Genelist_All, file="LUSC_DDRD_Assay_Genelist_All_150917.csv")
# CXCL10 MX1 IDO1 IFI44L CD2 GBP5 PRAME ITGAL LRP4 APOL3 CDR1 FYB TSPAN7 RAC2 KLHDC7B GRB14 KIF26A CD274 CD109 ETV7 MFAP5 OLFM4 PI15 FOSB FAM19A5 NLRC5 PRICKLE1 EGR1 CLDN10 ADAMTS4 SP140L ANXA1 RSAD2 ESR1 IKZF3 EGFR NAT1 LATS2 CYP2B6 PTPRC PPP1R1A
# Group2 genes which not up or down regulated (z.score: 1.8): CXCL10 MX1 IDO1 IFI44L CD2 GBP5 ITGAL APOL3 CDR1 FYB KLHDC7B GRB14  ETV7 MFAP5 OLFM4  NLRC5 EGR1 CLDN10 SP140L ANXA1 RSAD2 IKZF3 NAT1 PTPRC
# Group3 genes which are up or down regulated: EGFR MX1 GBP5 MFAP5 SP140L RSAD2 ANXA1 CXCL10 IFI44L PRAME FYB KLHDC7B GRB14 CD274 CD109  PI15 LATS2 NLRC5 PRICKLE1  ESR1 NAT1 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# For all 4 DDRD subgroups, LUSC
# http://www.cbioportal.org/index.do?session_id=59bc2361498e5df2e2955b1a&show_samples=false&heatmap_track_groups=lusc_tcga_rna_seq_v2_mrna_median_Zscores%2CCXCL10%2CMX1%2CIDO1%2CIFI44L%2CCD2%2CGBP5%2CPRAME%2CITGAL%2CLRP4%2CAPOL3%2CCDR1%2CFYB%2CRAC2%2CKLHDC7B%2CGRB14%2CKIF26A%2CCD274%2CCD109%2CETV7%2CMFAP5%2COLFM4%2CPI15%2CFOSB%2CFAM19A5%2CNLRC5%2CPRICKLE1%2CEGR1%2CCLDN10%2CADAMTS4%2CSP140L%2CANXA1%2CRSAD2%2CESR1%2CIKZF3%2CEGFR%2CNAT1%2CLATS2%2CCYP2B6%2CPTPRC%2CPPP1R1A&clinicallist=AJCC_METASTASIS_PATHOLOGIC_PM,AJCC_TUMOR_PATHOLOGIC_PT,DFS_STATUS,DFS_MONTHS,KRAS_MUTATION,LOCATION_LUNG_PARENCHYMA,MUTATION_STATUS,OS_MONTHS,OS_STATUS,TREATMENT_OUTCOME_FIRST_COURSE,SMOKING_PACK_YEARS,FRACTION_GENOME_ALTERED,PERFORMANCE_STATUS_TIMING,GENDER,PRIMARY_SITE,RESIDUAL_TUMOR,AJCC_PATHOLOGIC_TUMOR_STAGE,TARGETED_MOLECULAR_THERAPY
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# For DDRD 4, LUSC
#http://www.cbioportal.org/index.do?session_id=59bc3329498e5df2e2955c12&show_samples=false&heatmap_track_groups=lusc_tcga_rna_seq_v2_mrna_median_Zscores%2CCXCL10%2CMX1%2CIDO1%2CIFI44L%2CCD2%2CGBP5%2CPRAME%2CITGAL%2CLRP4%2CAPOL3%2CCDR1%2CFYB%2CRAC2%2CKLHDC7B%2CGRB14%2CKIF26A%2CCD274%2CCD109%2CETV7%2CMFAP5%2COLFM4%2CPI15%2CFOSB%2CFAM19A5%2CNLRC5%2CPRICKLE1%2CEGR1%2CCLDN10%2CADAMTS4%2CSP140L%2CANXA1%2CRSAD2%2CESR1%2CIKZF3%2CEGFR%2CNAT1%2CLATS2%2CCYP2B6%2CPTPRC%2CPPP1R1A&
# Genes I used in Bioportal: CXCL10 MX1 IDO1 IFI44L CD2 GBP5 PRAME ITGAL LRP4 APOL3 CDR1 FYB TSPAN7 RAC2 KLHDC7B GRB14 KIF26A CD274 CD109 ETV7 MFAP5 OLFM4 PI15 FOSB FAM19A5 NLRC5 PRICKLE1 EGR1 CLDN10 ADAMTS4 SP140L ANXA1 RSAD2 ESR1 IKZF3 EGFR NAT1 LATS2 CYP2B6 PTPRC PPP1R1A
# FA/BRCA pathway: BRCA1 BRCA2 BRIP1 FANCA FANCB FANCC FANCD2 FANCE FANCF FANCG FANCL FANCM
# EMT genes: CDH11 RAB31 COL5A1 COL10A1 VCAN FAP FN1 ANGPTL2 GJB2 INHBA MMP14 PLAU THBS1 THBS2 GFPT2 
#--------------------------------------
#--------------------------------------
# LUSC treated: DDRD 3 & 4 (n=39)
# http://www.cbioportal.org/index.do?session_id=59c5453a498e5df2e295a6a3&show_samples=false&&heatmap_track_groups=lusc_tcga_rna_seq_v2_mrna_median_Zscores%2CBRCA1%2CBRCA2%2CBRIP1%2CFANCA%2CFANCB%2CFANCC%2CFANCD2%2CFANCE%2CFANCF%2CFANCG%2CFANCL%2CFANCM&
# LUSC Untreated: DDRD 3 & 4
# http://www.cbioportal.org/index.do?session_id=59c8d615498e5df2e295b758&show_samples=false&heatmap_track_groups=lusc_tcga_rna_seq_v2_mrna_median_Zscores%2CBRCA1%2CBRCA2%2CBRIP1%2CFANCA%2CFANCB%2CFANCC%2CFANCD2%2CFANCE%2CFANCF%2CFANCG%2CFANCL%2CFANCM&
# LUSC Untreated: DDRD 1 & 2
# http://www.cbioportal.org/index.do?session_id=59c8dd66498e5df2e295b79a&show_samples=false&heatmap_track_groups=lusc_tcga_rna_seq_v2_mrna_median_Zscores%2CBRCA1%2CBRCA2%2CBRIP1%2CFANCA%2CFANCB%2CFANCC%2CFANCD2%2CFANCE%2CFANCF%2CFANCG%2CFANCL%2CFANCM&
# LUSC Treated: DDRD 1 & 2
# http://www.cbioportal.org/index.do?session_id=59c8e160498e5df2e295b7b5&show_samples=false&heatmap_track_groups=lusc_tcga_rna_seq_v2_mrna_median_Zscores%2CBRCA1%2CBRCA2%2CBRIP1%2CFANCA%2CFANCB%2CFANCC%2CFANCD2%2CFANCE%2CFANCF%2CFANCG%2CFANCL%2CFANCM&
# -------------------------------------
# -------------------------------------
# for all 4 DDRD subgroup, treated,LUAD, n=55
#http://www.cbioportal.org/index.do?session_id=59c296b9498e5df2e2958a76&show_samples=false&&heatmap_track_groups=luad_tcga_rna_seq_v2_mrna_median_Zscores%2CCXCL10%2CMX1%2CIDO1%2CIFI44L%2CCD2%2CGBP5%2CPRAME%2CITGAL%2CLRP4%2CAPOL3%2CCDR1%2CFYB%2CRAC2%2CKLHDC7B%2CGRB14%2CKIF26A%2CCD274%2CCD109%2CETV7%2CMFAP5%2COLFM4%2CPI15%2CFOSB%2CFAM19A5%2CNLRC5%2CPRICKLE1%2CEGR1%2CCLDN10%2CADAMTS4%2CSP140L%2CANXA1%2CRSAD2%2CESR1%2CIKZF3%2CEGFR%2CNAT1%2CLATS2%2CCYP2B6%2CPTPRC%2CPPP1R1A&
#http://www.cbioportal.org/index.do?session_id=59c296b9498e5df2e2958a76&show_samples=false&&&clinicallist=AJCC_METASTASIS_PATHOLOGIC_PM,AJCC_TUMOR_PATHOLOGIC_PT,DFS_STATUS,LOCATION_LUNG_PARENCHYMA,DFS_MONTHS,KRAS_MUTATION,MUTATION_STATUS,OS_MONTHS,OS_STATUS,TREATMENT_OUTCOME_FIRST_COURSE,SMOKING_PACK_YEARS,PERFORMANCE_STATUS_TIMING,GENDER,FRACTION_GENOME_ALTERED,PRIMARY_SITE,RESIDUAL_TUMOR,AJCC_PATHOLOGIC_TUMOR_STAGE,TARGETED_MOLECULAR_THERAPY,ALK_TRANSLOCATION_STATUS&heatmap_track_groups=luad_tcga_rna_seq_v2_mrna_median_Zscores%2CCXCL10%2CMX1%2CIDO1%2CIFI44L%2CCD2%2CGBP5%2CPRAME%2CITGAL%2CLRP4%2CAPOL3%2CCDR1%2CFYB%2CRAC2%2CKLHDC7B%2CGRB14%2CKIF26A%2CCD274%2CCD109%2CETV7%2CMFAP5%2COLFM4%2CPI15%2CFOSB%2CFAM19A5%2CNLRC5%2CPRICKLE1%2CEGR1%2CCLDN10%2CADAMTS4%2CSP140L%2CANXA1%2CRSAD2%2CESR1%2CIKZF3%2CEGFR%2CNAT1%2CLATS2%2CCYP2B6%2CPTPRC%2CPPP1R1A
# LUAD, Treated, DDRD 3 & 4, n=28
# http://www.cbioportal.org/index.do?session_id=59c9109a498e5df2e295b99b&show_samples=false&&&heatmap_track_groups=luad_tcga_rna_seq_v2_mrna_median_Zscores%2CBRCA1%2CBRCA2%2CBRIP1%2CFANCA%2CFANCB%2CFANCC%2CFANCD2%2CFANCE%2CFANCF%2CFANCG%2CFANCL%2CFANCM&clinicallist=TREATMENT_OUTCOME_FIRST_COURSE
# LUAD, Treated, DDRD 1 & 2, n=27
# http://www.cbioportal.org/index.do?session_id=59c9154f498e5df2e295b9e5&show_samples=false&heatmap_track_groups=luad_tcga_rna_seq_v2_mrna_median_Zscores%2CBRCA1%2CBRCA2%2CBRIP1%2CFANCA%2CFANCB%2CFANCC%2CFANCD2%2CFANCE%2CFANCF%2CFANCG%2CFANCL%2CFANCM&clinicallist=TREATMENT_OUTCOME_FIRST_COURSE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LUSC, DDRD 1 & 2, All stage II & III, (Treated as well as those patients who have missing values in treatement outcome columns)
# http://www.cbioportal.org/index.do?session_id=59cbbbfb498e5df2e295daf4&show_samples=false&heatmap_track_groups=lusc_tcga_rna_seq_v2_mrna_median_Zscores%2CBRCA1%2CBRCA2%2CBRIP1%2CFANCA%2CFANCB%2CFANCC%2CFANCD2%2CFANCE%2CFANCF%2CFANCG%2CFANCL%2CFANCM&
# LUSC, DDRD 3 & 4, All stage II & III, (Treated as well as those patients who have missing values in treatement outcome columns)
# http://www.cbioportal.org/index.do?session_id=59cbce35498e5df2e295dc89&show_samples=false&&heatmap_track_groups=lusc_tcga_rna_seq_v2_mrna_median_Zscores%2CBRCA1%2CBRCA2%2CBRIP1%2CFANCA%2CFANCB%2CFANCC%2CFANCD2%2CFANCE%2CFANCF%2CFANCG%2CFANCL%2CFANCM&clinicallist=TREATMENT_OUTCOME_FIRST_COURSE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LUSC
# DDRD group 1, n=16
# Ordered by DDRD score (low to high)
TCGA-43-A474-01
TCGA-NC-A5HF-01
TCGA-NC-A5HT-01
TCGA-56-A4BX-01
TCGA-77-A5G6-01
TCGA-94-8035-01
TCGA-33-A5GW-01
TCGA-63-A5MI-01
TCGA-58-8391-01
TCGA-56-8308-01
TCGA-85-8352-01
TCGA-94-A5I4-01
TCGA-85-8052-01
TCGA-58-8390-01
TCGA-43-A56V-01
TCGA-92-8063-01
# Group 2, n=17
# Ordered by DDRD score (low to high)
TCGA-NK-A5D1-01
TCGA-68-A59I-01
TCGA-85-A4CN-01
TCGA-98-A53D-01
TCGA-85-8351-01
TCGA-85-8353-01
TCGA-77-8139-01
TCGA-96-A4JK-01
TCGA-77-A5G3-01
TCGA-58-A46J-01
TCGA-63-A5MU-01
TCGA-J1-A4AH-01
TCGA-NC-A5HO-01
TCGA-85-A4JC-01
TCGA-63-A5MV-01
TCGA-22-A5C4-01
TCGA-56-8629-01
# Group 3, n=14
# Ordered by DDRD score (low to high)
TCGA-77-A5G1-01
TCGA-34-A5IX-01
TCGA-37-A5EN-01
TCGA-NC-A5HQ-01
TCGA-NC-A5HR-01
TCGA-85-8288-01
TCGA-58-A46K-01
TCGA-85-A4JB-01
TCGA-63-A5MJ-01
TCGA-37-A5EM-01
TCGA-85-8276-01
TCGA-NC-A5HN-01
TCGA-NC-A5HG-01
TCGA-58-A46L-01
# Group 4, n=25
# Ordered by DDRD score (low to high)
TCGA-NK-A5CX-01
TCGA-85-A511-01
TCGA-94-8491-01
TCGA-63-A5MM-01
TCGA-77-A5GF-01
TCGA-NK-A7XE-01
TCGA-77-A5FZ-01
TCGA-85-A510-01
TCGA-63-A5MP-01
TCGA-92-8065-01
TCGA-63-A5MT-01
TCGA-56-8625-01
TCGA-34-8454-01
TCGA-NC-A5HK-01
TCGA-33-AASI-01
TCGA-58-8387-01
TCGA-90-A4EE-01
TCGA-58-A46M-01
TCGA-56-8307-01
TCGA-63-A5MN-01
TCGA-NC-A5HJ-01
TCGA-96-A4JL-01
TCGA-NC-A5HE-01
TCGA-90-A59Q-01
TCGA-43-A475-01
# Two groups in group 4: based on upregulated genes
#**************************************************************************************
#**************************************************************************************
# LUSC only 25 DDRD Group 4
#http://www.cbioportal.org/index.do?session_id=59c5316b498e5df2e295a56e&show_samples=false&heatmap_track_groups=lusc_tcga_rna_seq_v2_mrna_median_Zscores%2CBRCA1%2CBRCA2%2CBRIP1%2CFANCA%2CFANCB%2CFANCC%2CFANCD2%2CFANCE%2CFANCF%2CFANCG%2CFANCL%2CFANCM&
# Normal (Cases without Alteration(s) in Query Gene(s))
# 1)TCGA-NK-A5CX-01
# 3)TCGA-94-8491-01
# 4)TCGA-63-A5MM-01
# 5)TCGA-77-A5GF-01
# 6)TCGA-NK-A7XE-01
# 7)TCGA-77-A5FZ-01
# 10)TCGA-92-8065-01
# 12)TCGA-56-8625-01
# 15)TCGA-33-AASI-01
# 17)TCGA-90-A4EE-01
# 20)TCGA-63-A5MN-01
# 21)TCGA-NC-A5HJ-01
# 23)TCGA-NC-A5HE-01
# 24)TCGA-90-A59Q-01

###############
# Upregulated ones (Cases with Alteration(s) in genes)
# 2)TCGA-85-A511-01
# 8)TCGA-85-A510-01
# 9)TCGA-63-A5MP-01
# 11)TCGA-63-A5MT-01
# 13)TCGA-34-8454-01
# 14)TCGA-NC-A5HK-01
# 16)TCGA-58-8387-01
# 18)TCGA-58-A46M-01
# 19)TCGA-56-8307-01
# 22)TCGA-96-A4JL-01
# 25)TCGA-43-A475-01

# Overall Survival Kaplan-Meier Estimate
# 
# Cases with Alteration(s) in Query Gene(s)
# Case ID	Number at Risk	Status	Survival Rate	Time (months)
# TCGA-NC-A5HK	11	censored	1	4.2
# TCGA-43-A475	10	censored	1	9.72
# TCGA-58-8387	9	deceased	0.8888888888888888	13.24
# TCGA-85-A511	8	deceased	0.7777777777777777	14.95
# TCGA-85-A510	7	censored	0.7777777777777777	15.83
# TCGA-63-A5MT	6	censored	0.7777777777777777	16.36
# TCGA-63-A5MP	5	censored	0.7777777777777777	25.26
# TCGA-56-8307	4	censored	0.7777777777777777	26.87
# TCGA-96-A4JL	3	censored	0.7777777777777777	27.66
# TCGA-58-A46M	2	censored	0.7777777777777777	35.22
# TCGA-34-8454	1	censored	0.7777777777777777	38.76
# 
# Cases without Alteration(s) in Query Gene(s)
# Case ID	Number at Risk	Status	Survival Rate	Time (months)
# TCGA-NK-A7XE	14	censored	1	0.43
# TCGA-92-8065	13	censored	1	2.3
# TCGA-NK-A5CX	12	censored	1	3.65
# TCGA-56-8625	11	deceased	0.9090909090909091	10.35
# TCGA-90-A59Q	10	deceased	0.8181818181818181	10.58
# TCGA-63-A5MN	9	deceased	0.7272727272727272	11.33
# TCGA-NC-A5HJ	8	deceased	0.6363636363636362	13.73
# TCGA-63-A5MM	7	deceased	0.5454545454545453	14.98
# TCGA-90-A4EE	6	censored	0.5454545454545453	22.6
# TCGA-94-8491	5	censored	0.5454545454545453	26.61
# TCGA-77-A5GF	4	deceased	0.40909090909090895	27.6
# TCGA-33-AASI	3	deceased	0.2727272727272726	44.15
# TCGA-NC-A5HE	2	censored	0.2727272727272726	76.74
# TCGA-77-A5FZ	1	deceased	0	126.08

#----------------------------------------------------------
# Disease Free Survival Kaplan-Meier Estimate
# 
# Cases with Alteration(s) in Query Gene(s)
# Case ID	Number at Risk	Status	Survival Rate	Time (months)
# TCGA-NC-A5HK	10	censored	1	4.2
# TCGA-43-A475	9	censored	1	9.72
# TCGA-85-A510	8	relapsed	0.875	11.93
# TCGA-85-A511	7	relapsed	0.75	12.81
# TCGA-63-A5MT	6	relapsed	0.625	14.03
# TCGA-63-A5MP	5	relapsed	0.5	16.79
# TCGA-56-8307	4	censored	0.5	26.87
# TCGA-96-A4JL	3	censored	0.5	27.66
# TCGA-58-A46M	2	censored	0.5	35.22
# TCGA-34-8454	1	censored	0.5	38.76
# 
# Cases without Alteration(s) in Query Gene(s)
# Case ID	Number at Risk	Status	Survival Rate	Time (months)
# TCGA-NK-A7XE	12	censored	1	0.43
# TCGA-92-8065	11	censored	1	2.3
# TCGA-NK-A5CX	10	censored	1	3.65
# TCGA-63-A5MM	9	relapsed	0.8888888888888888	7.16
# TCGA-56-8625	8	relapsed	0.7777777777777777	8.94
# TCGA-90-A59Q	7	relapsed	0.6666666666666665	9
# TCGA-63-A5MN	6	relapsed	0.5555555555555555	10.51
# TCGA-NC-A5HJ	5	relapsed	0.4444444444444444	10.78
# TCGA-94-8491	4	relapsed	0.3333333333333333	21.06
# TCGA-90-A4EE	3	censored	0.3333333333333333	22.6
# TCGA-77-A5GF	2	relapsed	0.16666666666666666	23.72
# TCGA-NC-A5HE	1	censored	0.16666666666666666	76.74
#**************************************************************************************
#**************************************************************************************
# http://www.cbioportal.org/index.do?session_id=59c5453a498e5df2e295a6a3&show_samples=false&
# When using DDRD group 3&4
# Overall Survival Kaplan-Meier Estimate
# 
# Cases with Alteration(s) in Query Gene(s)
# Case ID	Number at Risk	Status	Survival Rate	Time (months)
# TCGA-NC-A5HK	20	censored	1	4.2
# TCGA-43-A475	19	censored	1	9.72
# TCGA-58-8387	18	deceased	0.9444444444444444	13.24
# TCGA-NC-A5HQ	17	deceased	0.8888888888888888	14.72
# TCGA-85-A511	16	deceased	0.8333333333333333	14.95
# TCGA-85-A510	15	censored	0.8333333333333333	15.83
# TCGA-63-A5MT	14	censored	0.8333333333333333	16.36
# TCGA-37-A5EN	13	censored	0.8333333333333333	21.68
# TCGA-63-A5MP	12	censored	0.8333333333333333	25.26
# TCGA-56-8307	11	censored	0.8333333333333333	26.87
# TCGA-96-A4JL	10	censored	0.8333333333333333	27.66
# TCGA-37-A5EM	9	censored	0.8333333333333333	28.48
# TCGA-85-A4JB	8	censored	0.8333333333333333	30.95
# TCGA-58-A46M	7	censored	0.8333333333333333	35.22
# TCGA-34-8454	6	censored	0.8333333333333333	38.76
# TCGA-NC-A5HR	5	censored	0.8333333333333333	40.87
# TCGA-NC-A5HN	4	censored	0.8333333333333333	49.24
# TCGA-58-A46L	3	censored	0.8333333333333333	56.6
# TCGA-63-A5MJ	2	censored	0.8333333333333333	59.92
# TCGA-NC-A5HG	1	censored	0.8333333333333333	64.49
# 
# Cases without Alteration(s) in Query Gene(s)
# Case ID	Number at Risk	Status	Survival Rate	Time (months)
# TCGA-NK-A7XE	19	censored	1	0.43
# TCGA-92-8065	18	censored	1	2.3
# TCGA-NK-A5CX	17	censored	1	3.65
# TCGA-56-8625	16	deceased	0.9375	10.35
# TCGA-90-A59Q	15	deceased	0.875	10.58
# TCGA-63-A5MN	14	deceased	0.8125	11.33
# TCGA-85-8288	13	deceased	0.75	13.21
# TCGA-NC-A5HJ	12	deceased	0.6875	13.73
# TCGA-63-A5MM	11	deceased	0.625	14.98
# TCGA-90-A4EE	10	censored	0.625	22.6
# TCGA-94-8491	9	censored	0.625	26.61
# TCGA-77-A5GF	8	deceased	0.546875	27.6
# TCGA-34-A5IX	7	censored	0.546875	33.87
# TCGA-58-A46K	6	deceased	0.4557291666666667	34.33
# TCGA-85-8276	5	censored	0.4557291666666667	34.49
# TCGA-33-AASI	4	deceased	0.341796875	44.15
# TCGA-NC-A5HE	3	censored	0.341796875	76.74
# TCGA-77-A5FZ	2	deceased	0.1708984375	126.08
# TCGA-77-A5G1	1	censored	0.1708984375	132.26
#-------------------------------------------------------------------------
# Disease Free Survival Kaplan-Meier Estimate
# 
# Cases with Alteration(s) in Query Gene(s)
# Case ID	Number at Risk	Status	Survival Rate	Time (months)
# TCGA-NC-A5HK	18	censored	1	4.2
# TCGA-43-A475	17	censored	1	9.72
# TCGA-85-A510	16	relapsed	0.9375	11.93
# TCGA-85-A511	15	relapsed	0.875	12.81
# TCGA-63-A5MT	14	relapsed	0.8125	14.03
# TCGA-63-A5MP	13	relapsed	0.75	16.79
# TCGA-37-A5EN	12	censored	0.75	21.68
# TCGA-56-8307	11	censored	0.75	26.87
# TCGA-96-A4JL	10	censored	0.75	27.66
# TCGA-37-A5EM	9	censored	0.75	28.48
# TCGA-85-A4JB	8	censored	0.75	30.95
# TCGA-58-A46M	7	censored	0.75	35.22
# TCGA-34-8454	6	censored	0.75	38.76
# TCGA-NC-A5HR	5	censored	0.75	40.87
# TCGA-NC-A5HN	4	censored	0.75	49.24
# TCGA-58-A46L	3	censored	0.75	56.6
# TCGA-63-A5MJ	2	censored	0.75	59.92
# TCGA-NC-A5HG	1	censored	0.75	64.49
# 
# Cases without Alteration(s) in Query Gene(s)
# Case ID	Number at Risk	Status	Survival Rate	Time (months)
# TCGA-NK-A7XE	16	censored	1	0.43
# TCGA-92-8065	15	censored	1	2.3
# TCGA-NK-A5CX	14	censored	1	3.65
# TCGA-63-A5MM	13	relapsed	0.9230769230769231	7.16
# TCGA-56-8625	12	relapsed	0.8461538461538461	8.94
# TCGA-90-A59Q	11	relapsed	0.7692307692307692	9
# TCGA-63-A5MN	10	relapsed	0.6923076923076923	10.51
# TCGA-NC-A5HJ	9	relapsed	0.6153846153846153	10.78
# TCGA-94-8491	8	relapsed	0.5384615384615384	21.06
# TCGA-90-A4EE	7	censored	0.5384615384615384	22.6
# TCGA-77-A5GF	6	relapsed	0.44871794871794873	23.72
# TCGA-58-A46K	5	relapsed	0.35897435897435903	24.15
# TCGA-34-A5IX	4	censored	0.35897435897435903	33.87
# TCGA-85-8276	3	relapsed	0.23931623931623935	34.43
# TCGA-77-A5G1	2	relapsed	0.11965811965811968	72.17
# TCGA-NC-A5HE	1	censored	0.11965811965811968	76.74

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
######################################################
# Sample ID, LUAD, Untreated, n=95 (out of 230, Nature Method)
# DDRD 3 & 4
# http://www.cbioportal.org/index.do?session_id=59ca42ef498e5df2e295c5bc&show_samples=false&heatmap_track_groups=luad_tcga_pub_mrna_median_Zscores%2CBRCA1%2CBRCA2%2CBRIP1%2CFANCA%2CFANCB%2CFANCC%2CFANCD2%2CFANCE%2CFANCF%2CFANCG%2CFANCL%2CFANCM&
# DDRD 1 & 2
# http://www.cbioportal.org/index.do?session_id=59ca447d498e5df2e295c5ce&show_samples=false&heatmap_track_groups=luad_tcga_pub_mrna_median_Zscores%2CBRCA1%2CBRCA2%2CBRIP1%2CFANCA%2CFANCB%2CFANCC%2CFANCD2%2CFANCE%2CFANCF%2CFANCG%2CFANCL%2CFANCM&
# ----------------------------------------------------
# DDRD Group 1
TCGA-95-7039-01
TCGA-91-6849-01
TCGA-86-6562-01
TCGA-78-7161-01
TCGA-78-7154-01
TCGA-78-7150-01
TCGA-78-7149-01
TCGA-78-7148-01
TCGA-75-7030-01
TCGA-75-6214-01
TCGA-75-6203-01
TCGA-73-7498-01
TCGA-73-4676-01
TCGA-73-4675-01
TCGA-73-4659-01
TCGA-69-7760-01
TCGA-64-5775-01
TCGA-64-1679-01
TCGA-55-6983-01
TCGA-55-6981-01
TCGA-50-6593-01
TCGA-50-5936-01
TCGA-49-4510-01
TCGA-44-7670-01
TCGA-44-6774-01
TCGA-44-2665-01
TCGA-38-6178-01
TCGA-38-4627-01
TCGA-38-4626-01
TCGA-05-4415-01
TCGA-05-4396-01
TCGA-05-4384-01
# DDRD Group 2
TCGA-97-7554-01
TCGA-86-7714-01
TCGA-80-5607-01
TCGA-78-7166-01
TCGA-78-7158-01
TCGA-75-6212-01
TCGA-64-1678-01
TCGA-55-7576-01
TCGA-55-1594-01
TCGA-53-7626-01
TCGA-50-5933-01
TCGA-50-5932-01
TCGA-50-5072-01
TCGA-50-5051-01
TCGA-49-6742-01
TCGA-49-4512-01
TCGA-49-4490-01
TCGA-44-6146-01
TCGA-38-4628-01
# DDRD Group 3
TCGA-91-7771-01
TCGA-78-7536-01
TCGA-75-6207-01
TCGA-67-6217-01
TCGA-64-1677-01
TCGA-55-7907-01
TCGA-55-6982-01
TCGA-55-6970-01
TCGA-50-5044-01
TCGA-49-6761-01
TCGA-49-6745-01
TCGA-49-4506-01
TCGA-49-4505-01
TCGA-44-3396-01
TCGA-44-2659-01
TCGA-05-5429-01
TCGA-05-5428-01
TCGA-05-5423-01
TCGA-05-4432-01
TCGA-05-4424-01
TCGA-05-4398-01
TCGA-05-4395-01
# DDRD Group 4
TCGA-99-7458-01
TCGA-86-6851-01
TCGA-78-7539-01
TCGA-78-7147-01
TCGA-75-5126-01
TCGA-55-7914-01
TCGA-55-7727-01
TCGA-55-6979-01
TCGA-55-6978-01
TCGA-55-1596-01
TCGA-53-7813-01
TCGA-50-6595-01
TCGA-50-5941-01
TCGA-50-5068-01
TCGA-50-5055-01
TCGA-49-6767-01
TCGA-49-6744-01
TCGA-49-4507-01
TCGA-49-4494-01
TCGA-44-6779-01
TCGA-05-5420-01
TCGA-05-4418-01
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
######################################################
# All LUSC treated stage II & III, 4 DDRD paients, n=71
LUSC_4DDRD_Treated_withAlterations <- read.table("OncoPrintPatients_withAlterations_LUSC_4DDRD_OrderedByDDRDSCore.txt")
colnames(LUSC_4DDRD_Treated_withAlterations)[1] <- "Sample_ID"
colnames(LUSC_4DDRD_Treated_withAlterations)[2] <- "BRCA_FA.Alterations"
colnames(LUSC_4DDRD_Treated_withAlterations)[3] <- "DDRD.call"
LUSC_4DDRD_Treated_withAlterations$Sample_ID <- paste(LUSC_4DDRD_Treated_withAlterations$Sample_ID,"-01", sep = "") 
write.csv(LUSC_4DDRD_Treated_withAlterations, file="OncoPrintPatients_withAlterations_LUSC_4DDRD_OrderedByDDRDSCore_031017.csv")
######################################################
# All LUSC stage II & III, 4 DDRD paients, n=246
#http://www.cbioportal.org/index.do?session_id=59d3bc29498e5df2e2961fdd&show_samples=false&
# ordered by Sample ID
# http://www.cbioportal.org/index.do?session_id=59d4b013498e5df2e2962877&show_samples=false&
LUSC_4DDRD_stageII_III_withAlterations <- read.table("LUSC_Sample_ID_withAlterations.txt")

Pathfolder <- "/home/reza/Documents/Live/"           # initialise the path of the csv file (from my Linux) 
csvfilename <- "LUSC_4DDRDSubgroups_StageII_III_OrderedBySampleID_040917.csv"
pathcsvfile <- paste(Pathfolder,csvfilename,sep = "")
datafileObject <- read.csv(pathcsvfile, header=T, row.names = 1)  # for LUAD or LUSC
datafileObject <- cbind(datafileObject,LUSC_4DDRD_stageII_III_withAlterations)
datafileObject <- datafileObject[,-9]
colnames(datafileObject)[9] <- "BRCA_FA.Alterations"
datafileObject <- datafileObject[order(datafileObject$DDRD.score),] # order by DDRD scores from low to high
write.csv(datafileObject, file="OncoPrintPatients_withAlterations_LUSC_stageII_III_4DDRD_OrderedByDDRDSCore_041017.csv")


######################################################
# All LUSC 4 DDRD subgroup patients, Treated, n=72, ordered by DDRD score (from low to high)
TCGA-43-A474-01
TCGA-NC-A5HF-01
TCGA-NC-A5HT-01
TCGA-56-A4BX-01
TCGA-77-A5G6-01
TCGA-94-8035-01
TCGA-33-A5GW-01
TCGA-63-A5MI-01
TCGA-58-8391-01
TCGA-56-8308-01
TCGA-85-8352-01
TCGA-94-A5I4-01
TCGA-85-8052-01
TCGA-58-8390-01
TCGA-43-A56V-01
TCGA-92-8063-01
TCGA-NK-A5D1-01
TCGA-68-A59I-01
TCGA-85-A4CN-01
TCGA-98-A53D-01
TCGA-85-8351-01
TCGA-85-8353-01
TCGA-77-8139-01
TCGA-96-A4JK-01
TCGA-77-A5G3-01
TCGA-58-A46J-01
TCGA-63-A5MU-01
TCGA-J1-A4AH-01
TCGA-NC-A5HO-01
TCGA-85-A4JC-01
TCGA-63-A5MV-01
TCGA-22-A5C4-01
TCGA-56-8629-01
TCGA-77-A5G1-01
TCGA-34-A5IX-01
TCGA-37-A5EN-01
TCGA-NC-A5HQ-01
TCGA-NC-A5HR-01
TCGA-85-8288-01
TCGA-58-A46K-01
TCGA-85-A4JB-01
TCGA-63-A5MJ-01
TCGA-37-A5EM-01
TCGA-85-8276-01
TCGA-NC-A5HN-01
TCGA-NC-A5HG-01
TCGA-58-A46L-01
TCGA-NK-A5CX-01
TCGA-85-A511-01
TCGA-94-8491-01
TCGA-63-A5MM-01
TCGA-77-A5GF-01
TCGA-NK-A7XE-01
TCGA-77-A5FZ-01
TCGA-85-A510-01
TCGA-63-A5MP-01
TCGA-92-8065-01
TCGA-63-A5MT-01
TCGA-56-8625-01
TCGA-34-8454-01
TCGA-NC-A5HK-01
TCGA-33-AASI-01
TCGA-58-8387-01
TCGA-90-A4EE-01
TCGA-58-A46M-01
TCGA-56-8307-01
TCGA-63-A5MN-01
TCGA-NC-A5HJ-01
TCGA-96-A4JL-01
TCGA-NC-A5HE-01
TCGA-90-A59Q-01
TCGA-43-A475-01

# from above list:
# Altered cases (DDRD 1 to DDRD 4)

# ordered by subgroups from 1 to 4

# DDRD 1

########
# DDRD 2


########
# DDRD 3


########
# DDRD 4


########
# Case matrix: (1= Case harbors alteration in one of the input genes)
# TCGA-22-A5C4-01	1
# TCGA-33-A5GW-01	0
# TCGA-33-AASI-01	0
# TCGA-34-8454-01	1
# TCGA-34-A5IX-01	0
# TCGA-37-A5EM-01	1
# TCGA-37-A5EN-01	1
# TCGA-43-A474-01	0
# TCGA-43-A475-01	1
# TCGA-43-A56V-01	1
# TCGA-56-8307-01	1
# TCGA-56-8308-01	0
# TCGA-56-8625-01	0
# TCGA-56-8629-01	0
# TCGA-56-A4BX-01	1
# TCGA-58-8387-01	1
# TCGA-58-8390-01	0
# TCGA-58-8391-01	1
# TCGA-58-A46J-01	1
# TCGA-58-A46K-01	0
# TCGA-58-A46L-01	1
# TCGA-58-A46M-01	1
# TCGA-63-A5MI-01	1
# TCGA-63-A5MJ-01	1
# TCGA-63-A5MM-01	0
# TCGA-63-A5MN-01	0
# TCGA-63-A5MP-01	1
# TCGA-63-A5MT-01	1
# TCGA-63-A5MU-01	1
# TCGA-63-A5MV-01	0
# TCGA-68-A59I-01	0
# TCGA-77-8139-01	1
# TCGA-77-A5FZ-01	0
# TCGA-77-A5G1-01	0
# TCGA-77-A5G3-01	1
# TCGA-77-A5G6-01	1
# TCGA-77-A5GF-01	0
# TCGA-85-8052-01	0
# TCGA-85-8276-01	0
# TCGA-85-8288-01	0
# TCGA-85-8351-01	0
# TCGA-85-8352-01	0
# TCGA-85-8353-01	0
# TCGA-85-A4CN-01	0
# TCGA-85-A4JB-01	1
# TCGA-85-A4JC-01	0
# TCGA-85-A510-01	1
# TCGA-85-A511-01	1
# TCGA-90-A4EE-01	0
# TCGA-90-A59Q-01	0
# TCGA-92-8063-01	0
# TCGA-92-8065-01	0
# TCGA-94-8035-01	0
# TCGA-94-8491-01	0
# TCGA-94-A5I4-01	1
# TCGA-96-A4JK-01	1
# TCGA-96-A4JL-01	1
# TCGA-98-A53D-01	0
# TCGA-J1-A4AH-01	1
# TCGA-NC-A5HE-01	0
# TCGA-NC-A5HF-01	0
# TCGA-NC-A5HG-01	1
# TCGA-NC-A5HJ-01	0
# TCGA-NC-A5HK-01	1
# TCGA-NC-A5HN-01	1
# TCGA-NC-A5HO-01	1
# TCGA-NC-A5HQ-01	1
# TCGA-NC-A5HR-01	1
# TCGA-NC-A5HT-01	0
# TCGA-NK-A5CX-01	0
# TCGA-NK-A5D1-01	1
# TCGA-NK-A7XE-01	0

#############################################################
# Group , LUSC, Treated, Dataset for boxplot in: http://shiny.chemgrid.org/boxplotr/

LUSC_4DDRD_stageII_III <- tmp_all_NATratment

G1 <- matrix(nrow=length(which(tmp_all_NATratment$DDRD.call == "1")),ncol = 1)
G2 <- matrix(nrow=length(which(tmp_all_NATratment$DDRD.call == "2")),ncol = 1)
G3 <- matrix(nrow=length(which(tmp_all_NATratment$DDRD.call == "3")),ncol = 1)
G4 <- matrix(nrow=length(which(tmp_all_NATratment$DDRD.call == "4")),ncol = 1)
G1234 <- data.frame(matrix(nrow=length(which(tmp_all_NATratment$DDRD.call == "4")),ncol = 4))

G1234[1:16,1] <- tmp_all_NATratment[which(tmp_all_NATratment$DDRD.call == "1"),7]
G1234[1:17,2] <- tmp_all_NATratment[which(tmp_all_NATratment$DDRD.call == "2"),7]
G1234[1:14,3] <- tmp_all_NATratment[which(tmp_all_NATratment$DDRD.call == "3"),7]
G1234[1:25,4] <- tmp_all_NATratment[which(tmp_all_NATratment$DDRD.call == "4"),7]

colnames(G1234)[1] <- "DDRD.1"
colnames(G1234)[2] <- "DDRD.2"
colnames(G1234)[3] <- "DDRD.3"
colnames(G1234)[4] <- "DDRD.4"

write.csv(G1234, file="G1234_DDRDScore_LUSC_Treated.csv",row.names=FALSE)
write.csv(G1234, file="G1234_EMT_LUSC_Treated.csv",row.names=FALSE)
######################################################
######################################################
######################################################
# LUAD, n=206, stage II & III, with BRCA/FA genes alterations from Bioportal, EMT scores and DDRD scores
LUAD_4DDRD_stageII_III_withAlterations <- read.table("LUAD_BRCAFA_Alterations_206.txt")

Pathfolder <- "/home/reza/Documents/Live/"           # initialise the path of the csv file (from my Linux) 
csvfilename <- "LUAD_4DDRDSubgroups_StageII_III_OrderbySample_ID_withEMT_206_101017.csv"
pathcsvfile <- paste(Pathfolder,csvfilename,sep = "")
datafileObject <- read.csv(pathcsvfile, header=T, row.names = 1)  # for LUAD or LUSC
datafileObject <- cbind(datafileObject,LUAD_4DDRD_stageII_III_withAlterations)
datafileObject <- datafileObject[,-10]
colnames(datafileObject)[10] <- "BRCA_FA.Alterations"
datafileObject <- datafileObject[order(datafileObject$DDRD.score),] # order by DDRD scores from low to high
write.csv(datafileObject, file="OncoPrintPatients_withAlterations_LUAD_stageII_III_4DDRD_OrderedByDDRDSCore_101017.csv")

# LUAD, stage II & III, n=55 treated, ordered by DDRD score from low to high 
# OncoPrintPatients_withAlterations_LUAD_stageII_III_4DDRD_OrderedByDDRDSCore_101017.csv
# http://www.cbioportal.org/index.do?session_id=59dcc7d8498e5df2e2966b22&show_samples=false&heatmap_track_groups=luad_tcga_rna_seq_v2_mrna_median_Zscores%2CBRCA1%2CBRCA2%2CBRIP1%2CFANCA%2CFANCB%2CFANCC%2CFANCD2%2CFANCE%2CFANCF%2CFANCG%2CFANCL%2CFANCM&clinicallist=TREATMENT_OUTCOME_FIRST_COURSE,GENDER,AJCC_METASTASIS_PATHOLOGIC_PM

######################################################
######################################################
# LUSC, stage II & III, Treated
# 11 DDRD group4 LUSC with Recurred/Progressed
# Grp4_1 (Recurred)
TCGA-85-A511-01
TCGA-85-A510-01
TCGA-63-A5MP-01
TCGA-63-A5MT-01
TCGA-94-8491-01
TCGA-63-A5MM-01
TCGA-77-A5GF-01
TCGA-56-8625-01
TCGA-63-A5MN-01
TCGA-NC-A5HJ-01
TCGA-90-A59Q-01
# Grp4_2 (DiseaseFree)
TCGA-NK-A5CX-01
TCGA-NK-A7XE-01
TCGA-92-8065-01
TCGA-34-8454-01
TCGA-NC-A5HK-01
TCGA-90-A4EE-01
TCGA-58-A46M-01
TCGA-56-8307-01
TCGA-96-A4JL-01
TCGA-NC-A5HE-01
TCGA-43-A475-01
