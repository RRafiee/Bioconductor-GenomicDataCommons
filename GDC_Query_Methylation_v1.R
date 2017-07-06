#################################################################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##################################################################################################
# Start
#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
# Dr Reza Rafiee, June 2017
# Research Fellow, Queen's University Belfast
# Genomic Data common (GDC) query script for lung cancer - Methylation data (Illumina 450K and 27K platforms)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# GenomicDataCommons R-Package
# The National Cancer Institute (NCI) Genomic Data Commons provides the cancer research community
# with an open and unified repository for sharing and accessing data across numerous cancer studies
# and projects via a high-performance data transfer and query infrastructure. The Bioconductor project
# is an open source and open development software project built on the R statistical programming 
# environment. A major goal of the Bioconductor project is to facilitate the use, analysis, 
# and comprehension of genomic data. The GenomicDataCommons Bioconductor package provides basic 
# infrastructure for querying, accessing, and mining genomic datasets available from the GDC.
# We expect that Bioconductor developer and bioinformatics community will build on the GenomicDataCommons 
# package to add higher-level functionality and expose cancer genomics data to many state-of-the-art
# bioinformatics methods available in Bioconductor.

#setwd("~/ICGC/TargetExomeInputFiles/Fastq")

#Folder_path <- "C:/Users/3052092/Documents/GDCdata/TCGA-LUAD/legacy/DNA_methylation/Methylation_beta_value"
Folder_path <- "C:/Users/3052092/Documents/GDCdata/TCGA-LUSC/legacy/DNA_methylation/Methylation_beta_value"

# reading the list of methylation folders and corresponding files 
temp_folder <- list.dirs(path = Folder_path, full.names = TRUE, recursive = TRUE)

table1 <- as.data.frame.array(matrix(nrow = length(temp_folder)-1, ncol = 3,0))
colnames(table1)[1] <- "Path_Filename" 
colnames(table1)[2] <- "Sample_ID"  # technically this field is a barcode
colnames(table1)[3] <- "Flag"  # Flag showing the match between sample names of the clinical data and methylation

for (i in 2:length(temp_folder))
{
  #i <- 3
  temp_folder[i]
  txt_file_i <- list.files(path = temp_folder[i], pattern = NULL, all.files = FALSE,
             full.names = FALSE, recursive = FALSE,
             ignore.case = FALSE, include.dirs = FALSE)
  LUXX_Samples_DDRD <- read.table(paste(temp_folder[i],"/",txt_file_i,sep = ""),header=T, sep="\t")  # loading the single file 
  # reading the sample name/ID which is in the header of second column
  #colnames(LUXX_Samples_DDRD)[2]
  # 
  table1[i-1,1] <- paste(temp_folder[i],"/",txt_file_i,sep = "")
  table1[i-1,2] <- colnames(LUXX_Samples_DDRD)[2]
  table1[i-1,3] <- 0 # not yet matched with the clinical data (my reference)
}
#write.table(table1,"Index_table_LUAD_TCGA_Methylation_5July17.txt",sep="\t",row.names=TRUE)
write.table(table1,"Index_table_LUSC_TCGA_Methylation_5July17.txt",sep="\t",row.names=TRUE)

# match the Sample_ID with the reference (those which we have DDRD score for them)

#table1 <- read.table("Index_table_LUAD_TCGA_Methylation_5July17.txt",header=T, sep="\t")  # loading the single file 
table1 <- read.table("Index_table_LUSC_TCGA_Methylation_5July17.txt",header=T, sep="\t")  # loading the single file 

table1_methylation <- cbind(table1,0,0) # adding two extra columns in table1 for keeping the short version of Sample_IDs and platform
colnames(table1_methylation)[4] <- "Short_Sample_ID"
colnames(table1_methylation)[5] <- "Platform"

# 
for (i in 1:nrow(table1_methylation))
{
  #i <- 1
  pos11 <- regexpr('TCGA', table1_methylation[i,2])
  N1 <- substr(table1_methylation[i,2], pos11[1],pos11[1]+14)
  #N1 <- gsub("-", ".", N1)
  table1_methylation[i,4] <- N1
}

#load("df_TCGA_517_LUAD_41ordered_matched_Gene_CorrectGeneName_DDRDScore_3rdJuly17.rda")
#rownames(combinmatrix1_GroupsOrdered14)[1]

load("df_TCGA_501_LUSC_41ordered_matched_Gene_CorrectGeneName_DDRDScore_3rdJuly17.rda")
#rownames(combinmatrix1_GroupsOrdered14)[1]


length(which(table1_methylation$Short_Sample_ID %in% rownames(combinmatrix1_GroupsOrdered14) == FALSE))
length(which(table1_methylation$Short_Sample_ID %in% rownames(combinmatrix1_GroupsOrdered14) == TRUE)) # 535

table1_methylation[which(table1_methylation$Short_Sample_ID %in% rownames(combinmatrix1_GroupsOrdered14) == TRUE),3] <- "1"

# 450K or 27K?
for (i in 1:nrow(table1_methylation))
{
  #i <- 7
  #if (as.character(table1_methylation$Flag[i]) == "1")
  #{
    pos11 <- regexpr('HumanMethylation',as.character(table1_methylation$Path_Filename[i]))
    N1 <- substr(as.character(table1_methylation$Path_Filename[i]), pos11[1],pos11[1]+18)
    table1_methylation$Platform[i] <- N1
  #}
}

#write.table(table1_methylation,"Index_table_LUAD_TCGA_Methylation450K27K_FlagMatched_05July17.txt",sep="\t",row.names=TRUE)
write.table(table1_methylation,"Index_table_LUSC_TCGA_Methylation450K27K_FlagMatched_05July17.txt",sep="\t",row.names=TRUE)

#length(which(as.character(table1_methylation$Flag)=="1")) # for LUSC, n=535 while my reference is n=517
length(which(as.character(table1_methylation$Flag)=="1")) # for LUSC, n=501 and covered exactly all my reference samples

#load("LUAD_DNAMethylation450k_dataframe_28June17.rda") #Error in View : cannot allocate vector of size 3.7 Mb
#colnames(data)
#head(rownames(data))
#as.character(LUAD_Samples_DDRD$Sample_ID[1])

#-----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------
# We have multiple methylation data (LUAD and not LUSC) for a specific sample name (Sample_ID)
# This part handle this issue
table1_methylation <- read.table("Index_table_LUAD_TCGA_Methylation450K27K_FlagMatched_05July17.txt",header=T, sep="\t")  # loading the single file 
length(which(as.character(table1_methylation$Flag)=="1"))
table1_methylation_matched_withDuplicates <- table1_methylation[which(as.character(table1_methylation$Flag)=="1"),] 
table1_methylation_matched_withDuplicates <- cbind(table1_methylation_matched_withDuplicates,"0")
colnames(table1_methylation_matched_withDuplicates)[6] <- "Duplicate"

# order sample names column
table1_methylation_matched_withDuplicates_ordered <- as.data.frame(table1_methylation_matched_withDuplicates[order(table1_methylation_matched_withDuplicates$Short_Sample_ID),])

# DO NOT use R functions like duplicate, unique, etc.
UPi <- nrow(table1_methylation_matched_withDuplicates_ordered) - 1 
i <- 1
while (i <= UPi)
 {
  j <- i + 1 
   while (j <= nrow(table1_methylation_matched_withDuplicates_ordered))
   {
     ix1 <- as.character(table1_methylation_matched_withDuplicates_ordered$Short_Sample_ID[i])
     ix2 <- as.character(table1_methylation_matched_withDuplicates_ordered$Short_Sample_ID[j])
     if (ix1 == ix2)
     {
       #print("Matched")
       table1_methylation_matched_withDuplicates_ordered$Duplicate[i] <- "NA"
       table1_methylation_matched_withDuplicates_ordered$Duplicate[j] <- "NA"
     }
     j <- j + 1
   }
  i <- i + 1
 }
  
length(which(is.na(table1_methylation_matched_withDuplicates_ordered$Duplicate))) # n=30 (each sample has two duplicates)
table1_methylation_matched_withDuplicates_ordered[which(is.na(table1_methylation_matched_withDuplicates_ordered$Duplicate)),c(2,4,5)]
write.csv(table1_methylation_matched_withDuplicates_ordered,"Index_table_LUAD_TCGA_Methylation450K27K_FlagMatched_duplicates_06July17.csv")
write.table(table1_methylation_matched_withDuplicates_ordered,"Index_table_LUAD_TCGA_Methylation450K27K_FlagMatched_duplicates_06July17.txt",sep="\t",row.names=TRUE)

length(which(table1_methylation_matched_withDuplicates_ordered$Duplicate == 0)) # n=505 unique samples
# I will choose 10 samples from duplicates so totally n=515 unique methylaion data for LUAD

# Remove the two follwoing samples from the reference (n=515 in stead of 517)
rownames(combinmatrix1_GroupsOrdered14)[314] # TCGA.67.3770.01
rownames(combinmatrix1_GroupsOrdered14)[369] # TCGA.44.2659.01
#length(which(rownames(combinmatrix1_GroupsOrdered14) %in% unique(as.character(table1_methylation$Short_Sample_ID))))



#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
# End
##################################################################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##################################################################################################
