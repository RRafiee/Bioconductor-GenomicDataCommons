#################################################################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##################################################################################################
# Start
#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
# Dr Reza Rafiee, July 2017
# Research Fellow, Queen's University Belfast
# Applying Similarity Network Fusion (SNF) to different data types:
# 1) mRNA expression
# 2) Methylation 
# 3) CNVs
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


#################################################
#install_all_packages_automatic(SNFtool)
library(SNFtool)
#################################################
load("df_TCGA_517_LUAD_41ordered_matched_Gene_CorrectGeneName_DDRDScore_3rdJuly17.rda")
load("df_TCGA_501_LUSC_41ordered_matched_Gene_CorrectGeneName_DDRDScore_3rdJuly17.rda")

Data1 <- as.data.frame(combinmatrix1_GroupsOrdered14[,c(1:6,8,10,12,14,18,20,26,33,35,40)]) # Only 16 genes
#Data1 <- as.data.frame(combinmatrix1_GroupsOrdered14[,1:41])
truelabel = as.vector(combinmatrix1_GroupsOrdered14[,42]) ## the ground truth of the simulated data
# truelabel[which(truelabel=="4")] <- "3" grouping 3 and 4 into 3

## though it is optional depending on the data the users want to use.
Data1 = standardNormalization(Data1)
## First, set all the parameters:
K = 30; # number of neighbors, usually (10~30)
alpha = 0.8; # hyperparameter, usually (0.3~0.8)
T1 = 20; # Number of Iterations, usually (10~20)
## Data1 is of size n x d_1,
## where n is the number of patients, d_1 is the number of genes,
## Data2 is of size n x d_2,
## where n is the number of patients, d_2 is the number of methylation
# data(Data1)
# data(Data2)
## Here, the simulation data (SNFdata) has two data types. They are complementary to each other.
## And two data types have the same number of points.
## The first half data belongs to the first cluster; the rest belongs to the second cluster.
# truelabel = c(matrix(1,100,1),matrix(2,100,1)); ## the ground truth of the simulated data
## Calculate distance matrices
## (here we calculate Euclidean Distance, you can use other distance, e.g,correlation)
## If the data are all continuous values, we recommend the users to perform
## standard normalization before using SNF,
## Calculate the pair-wise distance;
## If the data is continuous, we recommend to use the function "dist2" as follows
Dist1 = dist2(as.matrix(Data1),as.matrix(Data1))

#Dist1 = chiDist2(as.matrix(Data1), as.matrix(Data1)) # is not proper
#Dist1 = dist(as.matrix(Data1), method = "euclidean",diag = TRUE, upper = TRUE, p = 2) #minkowski, euclidean

# Dist2 = dist2(as.matrix(Data2),as.matrix(Data2))

## next, construct similarity graphs: Computes affinity matrix from a generic distance matrix
# alpha is the Variance for local model, K is the Number of nearest neighbors and Dist1 is the Distance matrix
W1 = affinityMatrix(Dist1, K, alpha)
# W2 = affinityMatrix(Dist2, K, alpha)

## These similarity graphs have complementary information about clusters.
displayClusters(W1, truelabel)
#displayClusters(W2, truelabel)

W = SNF(list(W1,W1), K, T1)

## You can display clusters in the data by the following function
## where C is the number of clusters.
estimationResult = estimateNumberOfClustersGivenGraph(W, 2:5);
C = 4 # number of clusters which is 4
group = spectralClustering(W,C); # the final subtypes information

## Get a matrix containing the group information
## for the samples such as the SpectralClustering result and the True label
M_label=cbind(group,truelabel)
colnames(M_label)=c("spectralClustering","TrueLabel")

## Use the getColorsForGroups function to assign a color to each group
## NB is more than 8 groups, you will have to input a vector
## of colors into the getColorsForGroups function
M_label_colors=t(apply(M_label,1,getColorsForGroups))

## or choose you own colors for each label, for example:
M_label_colors=cbind("spectralClustering"=getColorsForGroups(M_label[,"spectralClustering"],
                                                             colors=c("blue","red","orange","darkgreen")),"TrueLabel"=getColorsForGroups(M_label[,"TrueLabel"],
                                                                                                                      colors=c("blue","red","Orange","darkgreen")))

## Visualize the clusters present in the given similarity matrix
## as well as some sample information
## In this presentation no clustering method is ran the samples
## are ordered in function of their group label present in the group arguments
displayClustersWithHeatmap(W, group, M_label_colors[,"spectralClustering"]) #Kernel NMF, spectralClustering
displayClustersWithHeatmap(W, group, M_label_colors) #ColSideColors=NULL

## Here we provide two ways to estimate the number of clusters. Note that,
## these two methods cannot guarantee the accuracy of esstimated number of
## clusters, but just to offer two insights about the datasets.

G1 <- "1"
LG_C <- which(M_label[,1] == G1) # spectral clustering
LG_T <- which(M_label[,2] == G1) # true labels
length(LG_T) #G4:129, G3:129
length(LG_C) #G4:97,  G3:88
#
# Combining data types.
# SNF can be used to incorporate arbitrary types of discrete (binary or categorical) 
# and continuous data. For integration of discrete data, we recommend the use of chi-squared 
# distance as the similarity measure. Compatibility of data sources can be checked via 
# normalized mutual information (NMI). If the patient similarity obtained from different data
# sources is completely discordant; NMI can help to clarify which data should and
# which should not be combined.

# Before applying our SNF, we performed three steps of preprocessing: 
# outlier removal, missing-data imputation and normalization.
# If a patient had more than 20% missing data in a certain data type,
# we did not consider this patient. Similarly, if a certain biological feature (for example, mRNA expression)
# had more than 20% of missing values across patients, we filtered out this feature. 
# Also, for missing data, we used K nearest neighbor (KNN) imputation21,
# where the number of neighbors is the same with K value used in our method 

#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
# End
##################################################################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##################################################################################################




