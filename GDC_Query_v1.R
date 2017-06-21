##################################################################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##################################################################################################
# Start
#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
# Dr Reza Rafiee, June 2017
# Research Fellow, Queen's University Belfast
# Genomic Data common (GDC) query script for lung cancer
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

install_all_packages_automatic <- function(x) {
  x <- as.character(substitute(x))
  if(isTRUE(x %in% .packages(all.available=TRUE))) {
    eval(parse(text = sprintf("require(\"%s\")", x)))
  } else {
    #update.packages(ask= FALSE) #update installed packages.
    eval(parse(text = sprintf("install.packages(\"%s\", dependencies = TRUE)", x)))
  }
  if(isTRUE(x %in% .packages(all.available=TRUE))) {
    eval(parse(text = sprintf("require(\"%s\")", x)))
  } else {
    source("http://bioconductor.org/biocLite.R")
    #biocLite(character(), ask=FALSE) #update installed packages.
    eval(parse(text = sprintf("biocLite(\"%s\")", x)))
    eval(parse(text = sprintf("require(\"%s\")", x)))
  }
}

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
library(GenomicDataCommons)

#############################################################
#Check basic functionality
GenomicDataCommons::status()

library(magrittr)
ge_manifest = files() %>% 
  filter( ~ cases.project.project_id == 'TCGA-OV' &
            type == 'gene_expression' &
            analysis.workflow_type == 'HTSeq - Counts') %>%
  manifest()

ge_manifest = files() %>% 
  filter( ~ cases.project.project_id == 'TCGA-LUAD' &
            type == 'gene_expression' &
            analysis.workflow_type == 'HTSeq - Counts') %>%
  manifest()
#############################################################

#############################################################
install_all_packages_automatic("mnormt")
install_all_packages_automatic("TCGAbiolinks")
library(TCGAbiolinks)

query <- GDCquery(project = "TCGA-LUAD", 
                  legacy = TRUE,
                  data.category = "Gene expression", 
                  data.type = "Gene expression quantification",
                  platform = "Illumina HiSeq",
                  file.type = "normalized_results", 
                  experimental.strategy = "RNA-Seq")
GDCdownload(query)
# Case 1: Dataframe as output
df <- GDCprepare(query, 
                 save=TRUE,
                 save.filename = "LUAD_geneExp_dataframe.rda",
                 summarizedExperiment = FALSE)

# Case 2: SummarizedExperiment object
se <- GDCprepare(query, 
                 save=TRUE,
                 save.filename = "LUAD_geneExp_dataframe.rda",
                 summarizedExperiment = TRUE)

# Error in the above command when summarizedExperiment is TRUE: lexical error: invalid char in json text.
# <?xml version="1.0" ?> <respons
# (right here) ------^

## get gene Expression values
geneExp <- SummarizedExperiment::assay(se)
#############################################################

# Loading data
load("LUAD_geneExp_dataframe.rda")
df_LUAD_TCGA <- data







#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
# End
##################################################################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##################################################################################################
