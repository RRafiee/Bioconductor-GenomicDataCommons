##################################################################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##################################################################################################
# Start
#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
# Dr Reza Rafiee, June 2017
# Research Fellow, Queen's University Belfast
# Genomic Data common (GDC) query script for lung cancer
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

## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("GenomicDataCommons")
install_all_packages_automatic("magrittr")

library(magrittr)
library(GenomicDataCommons)
GenomicDataCommons::status()

#############################################################
#Check basic functionality
# GenomicDataCommons::status()
# 
# library(magrittr)
# ge_manifest = files() %>% 
#   filter( ~ cases.project.project_id == 'TCGA-OV' &
#             type == 'gene_expression' &
#             analysis.workflow_type == 'HTSeq - Counts') %>%
#   manifest()
# 
# ge_manifest = files() %>% 
#   filter( ~ cases.project.project_id == '	TCGA-LUAD' &
#             type == 'Gene expression quantification' &
#             analysis.workflow_type == 'HiSeq - Counts') %>%
#   manifest()
#############################################################

#############################################################
# install_all_packages_automatic("mnormt")
# install_all_packages_automatic("TCGAbiolinks")

# check for the default fields
# so that we can use one of them to build a filter
# pQuery <- projects()
# default_fields(pQuery)
# pQuery <- filter(pQuery,~ project_id == 'TCGA-LUAD')
# get_filter(pQuery)

##############################################################################
library(mnormt)
library(TCGAbiolinks)
# Gene expression - RNA_Seq parameter setting ################################
TCGAprojectname <- "TCGA-LUAD" #"TCGA-LUSC"
TCGAdatacategory <- "Gene expression"
TCGAdatatype <- "Gene expression quantification"
TCGAplatform <- "Illumina HiSeq"
TCGAfiletype <- "normalized_results"
TCGAexperimentalstratgy <- "RNA-Seq"
##############################################################################
# DNA methylation - 450k parameter setting ###################################
TCGAprojectname <- "TCGA-LUSC" #"TCGA-LUAD" 
TCGAdatacategory <- "DNA Methylation"
TCGAdatatype <- "Methylation Beta Value"
TCGAplatform <- "Illumina Human Methylation 450"
TCGAfiletype <- ""
TCGAexperimentalstratgy <- "Methylation Array"
##############################################################################
# DNA methylation - 27k parameter setting ####################################
TCGAprojectname <- "TCGA-LUSC" #"TCGA-LUAD" #
TCGAdatacategory <- "DNA Methylation"
TCGAdatatype <- "Methylation beta value"
TCGAplatform <- "Illumina Human Methylation 27"
TCGAfiletype <- ""
TCGAexperimentalstratgy <- "Methylation Array"
##############################################################################
# Make a query using GDCquery function
TCGA_query <- GDCquery(project = TCGAprojectname, 
                       legacy = TRUE,
                       data.category = TCGAdatacategory, 
                       data.type = TCGAdatatype,
                       platform = TCGAplatform,
                       #file.type = TCGAfiletype, # should be commented for the Methylation array
                       experimental.strategy = TCGAexperimentalstratgy)

GDCdownload(TCGA_query)
df <- GDCprepare(TCGA_query, 
                 save=TRUE,
                 save.filename = paste(TCGAprojectname,"_",TCGAplatform,".rda",sep = ""),
                 summarizedExperiment = FALSE)
#######################################################################
# query_LUAD <- GDCquery(project = "TCGA-LUAD", 
#                   legacy = TRUE,
#                   data.category = "Gene expression", 
#                   data.type = "Gene expression quantification",
#                   platform = "Illumina HiSeq",
#                   file.type = "normalized_results", 
#                   experimental.strategy = "RNA-Seq")
# 
# 
# query_LUAD <- GDCquery(project = "TCGA-LUAD",
#                        legacy = TRUE,
#                        data.category = "DNA methylation",
#                        data.type = "Methylation beta value",
#                        experimental.strategy = "Methylation Array",
#                        platform = "Illumina Human Methylation 450")
# 
# GDCdownload(query_LUAD)
# df <- GDCprepare(query_LUAD, 
#                  save=TRUE,
#                  save.filename = "LUAD_DNAMethylation450k_dataframe_26June17.rda",
#                  summarizedExperiment = FALSE)
# 
# 
# 
# query_LUSC <- GDCquery(project = "TCGA-LUSC", 
#                   legacy = TRUE,
#                   data.category = "Gene expression", 
#                   data.type = "Gene expression quantification",
#                   platform = "Illumina HiSeq",
#                   file.type = "normalized_results", 
#                   experimental.strategy = "RNA-Seq")
# 
# 
# GDCdownload(query_LUAD)
# GDCdownload(query_LUSC)
# 
# query <- query_LUSC
# query <- query_LUAD
# 
# # Case 1: Dataframe as output
# df <- GDCprepare(query, 
#                  save=TRUE,
#                  save.filename = "LUAD_geneExp_dataframe.rda",
#                  summarizedExperiment = FALSE)
# 
# # Case 1: Dataframe as output
# df <- GDCprepare(query, 
#                  save=TRUE,
#                  save.filename = "LUSC_geneExp_dataframe.rda",
#                  summarizedExperiment = FALSE)
# 
# # # Case 2: SummarizedExperiment object
# # se <- GDCprepare(query, 
# #                  save=TRUE,
# #                  save.filename = "LUAD_geneExp_dataframe_TRUE.rda",
# #                  summarizedExperiment = TRUE)
# # Error in the above command when summarizedExperiment is TRUE: lexical error: invalid char in json text.
# # <?xml version="1.0" ?> <respons
# # (right here) ------^
# ## get gene Expression values
# #geneExp <- SummarizedExperiment::assay(se)
# #############################################################

# # Loading data
# load("LUAD_geneExp_dataframe.rda")
# df_LUAD_TCGA <- data
# df <- df_LUAD_TCGA
# 
# 
# load("LUSC_geneExp_dataframe.rda")
# df_LUSC_TCGA <- data
# df <- df_LUSC_TCGA


#############################################################
# This part aims to match colnames of df with the DDRD table

LUAD_Samples_DDRD <- read.table("tcga_luad_ClinicalData_cbioportal_incDDRD.txt",header=T, sep="\t")  # 
LUSC_Samples_DDRD <- read.table("tcga_lusc_ClinicalData_cbioportal_incDDRD.txt",header=T, sep="\t")  # 

head(LUAD_Samples_DDRD$AJCC_METASTASIS_PATHOLOGIC_PM)
head(LUAD_Samples_DDRD$AJCC_NODES_PATHOLOGIC_PN)

#LUXX_Samples_DDRD <- LUAD_Samples_DDRD
LUXX_Samples_DDRD <- LUSC_Samples_DDRD

# extract the name of samples matched with the DDRD table
for (i in 1:ncol(df))
{
  pos11 <- regexpr('TCGA', colnames(df)[i])
  N1 <- substr(colnames(df)[i], pos11[1],pos11[1]+14)
  N1 <- gsub("-", ".", N1)
  colnames(df)[i] <- N1
}

# which(colnames(df) == as.character(LUXX_Samples_DDRD[,1]))

# order the df table based on column names and LUXX_Samples_DDRS table based on the first row
LUXX_Samples_DDRD_ordered <- LUXX_Samples_DDRD[order(LUXX_Samples_DDRD[,1], decreasing = TRUE),]
df_ordered <- df[,order(colnames(df), decreasing = TRUE)]

matchtable <- match(colnames(df_ordered), LUXX_Samples_DDRD_ordered[,1])
which(matchtable[]!="NA")

df_ordered_matched <- df_ordered[,which(matchtable[]!="NA")]

# DDRD assay gene symbols
DDRD_Assay_Genelist <- c("CXCL10","MX1","IDO1","IFI44L","CD2","GBP5","PRAME","ITGAL","LRP4","APOL3","CDR1","FYB","TSPAN7",
"RAC2","KLHDC7B","GRB14","AC138128.1","KIF26A","CD274","CD109","ETV7","MFAP5","OLFM4","PI15",
"FOSB","FAM19A5","NLRC5","PRICKLE1","EGR1","CLDN10","ADAMTS4","SP140L","ANXA1","RSAD2","ESR1",
"IKZF3","OR2I1P","EGFR","NAT1","LATS2","CYP2B6","PTPRC","PPP1R1A","AL137218.1")

matchtable2 <- vector()
matchtable21<- matchtable2
for (j in 1:length(DDRD_Assay_Genelist))
{
  
  #j <- 2
  matchtable2 <- grep(DDRD_Assay_Genelist[j],rownames(df_ordered_matched))
  if (length(matchtable2) > 1)
  {
    for (k in 1:length(matchtable2))
    {
      #k <- 4
      gene_name_str <- rownames(df_ordered_matched)[matchtable2[k]]
      part_before_pipe <- unlist(strsplit(gene_name_str, "|", fixed=TRUE))[1] # first part of gene name EMX1 of EMX1|2016
      idx_match <- grep(paste("^",DDRD_Assay_Genelist[j],"$",sep = ""),part_before_pipe)
      if (length(idx_match) != 0)
      {
        matchtable2 <- matchtable2[k]
      }
    }
  }
  matchtable21 <- c(matchtable21,matchtable2)
}

# order the matrix
df_ordered_matched_Genes <- as.data.frame(df_ordered_matched[matchtable21,])

# There is no any matches for the following genes
matchtable3 <- grep("AC138128.1",rownames(df_ordered_matched)) # no gene name in downloaded table
matchtable3 <- grep("OR2I1P",rownames(df_ordered_matched))     # no gene name ...
matchtable3 <- grep("AL137218.1",rownames(df_ordered_matched)) # no gene name ...

# Transpose df_ordered_matched_Genes table
rn_Samples <- colnames(df_ordered_matched_Genes)
df_ordered_matched_Genes <- t(df_ordered_matched_Genes)
rownames(df_ordered_matched_Genes) <- rn_Samples

# save tables
save(df_ordered_matched_Genes,file="df_TCGA_LUAD_ordered_matched_Gene_23June17.rda")
#save(df_ordered_matched_Genes,file="df_TCGA_LUSC_ordered_matched_Gene_23June17.rda")
save(LUXX_Samples_DDRD_ordered,file="LUAD_Samples_DDRD_ordered_23June17.rda")
#save(LUXX_Samples_DDRD_ordered,file="LUSC_Samples_DDRD_ordered_23June17.rda")


# ## try http:// if https:// URLs are not supported
 source("https://bioconductor.org/biocLite.R")
 biocLite("TCGAbiolinksGUI")
 library(TCGAbiolinksGUI)
# #install_all_packages_automatic("shinydashboard")
# install.packages("shiny")
# install.packages("dplyr")
 library(shiny)
 library(dplyr)
 library(TCGAbiolinksGUI)  
 TCGAbiolinksGUI(run = TRUE)
#---------------------------------------------------
 
 
 #LUSC_Methylation_DDRD <- read.table("jhu-usc.edu_LUSC.HumanMethylation450.3.lvl-3.TCGA-22-5489-11A-01D-1633-05.txt",header=T, sep="\t")  # 
 #load("LUAD_DNAMethylation450k_dataframe_28June17.rda")
 #head(data$`TCGA-95-7043-01A-11D-1947-05`)
 
#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
# End
##################################################################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##################################################################################################