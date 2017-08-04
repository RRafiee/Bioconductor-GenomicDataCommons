##################################################################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##################################################################################################
# Start
#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
# Dr Reza Rafiee, August 2017
# Research Fellow, Queen's University Belfast
# Survival analysis (including OS and DFS) on lung cancer data (LUAD and LUSC)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##############################   LUAD  ################################
# Loading a csv file
# Pathfolder <- "C:/Users/3052092/Documents/Live/"   # initialise the path of the csv file 
Pathfolder <- "/home/reza/Documents/Live/"           # initialise the path of the csv file (from my Linux) 
csvfilename <- "LUAD_Samples_DDRD_ordered_23June17.csv"
pathcsvfile <- paste(Pathfolder,csvfilename,sep = "")
datafileObject <- read.csv(pathcsvfile, header=T)  # LUSC samples, TCGA

# Pre-processing of clinical data: binarization and so on
length(which(datafileObject$TREATMENT_OUTCOME_FIRST_COURSE == ""))  #368
length(which(datafileObject$TREATMENT_OUTCOME_FIRST_COURSE != ""))  #149
which(datafileObject$TREATMENT_OUTCOME_FIRST_COURSE != "")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ LUAD ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Select variables which you need to find association with subgroups
# Gender, ALK mutation,
tmp <- datafileObject[,c(10,20,29,42,46,47,77,77,77,77,79,82,83)] 
tmp1 <- tmp
tmp1$GENDER <- ifelse(tmp1$GENDER == "FEMALE", 0, ifelse(tmp1$GENDER == "MALE", 1, NA))
tmp1$ALK_TRANSLOCATION_STATUS <- ifelse(tmp1$ALK_TRANSLOCATION_STATUS == "NO", 0, ifelse(tmp1$ALK_TRANSLOCATION_STATUS == "YES", 1, NA))
tmp1$MUTATION_STATUS <- ifelse(tmp1$MUTATION_STATUS == "NO", 0, ifelse(tmp1$MUTATION_STATUS == "YES", 1, NA))
tmp1$TUMOR_STATUS  <- ifelse(tmp1$TUMOR_STATUS == "TUMOR FREE", 0, ifelse(tmp1$TUMOR_STATUS == "WITH TUMOR", 1, NA))
tmp1$DFS_STATUS <- ifelse(tmp1$DFS_STATUS == "DiseaseFree", 0, ifelse(tmp1$TUMOR_STATUS == "Recurred/Progressed", 1, NA))


# stable disease to describe a tumor that is neither growing nor shrinking. 
# Stable disease also means that no new tumors have developed and that 
# the cancer has not spread to any new regions of the body
# (the cancer is not getting better or worse and has not metastasized further).

# The term used for the absence of all detectable cancer after your treatment is complete response (CR). 
# Complete response doesn't necessarily mean that you are cured, but it is the best result that can be reported. 
# It means the cancerous tumor is now gone and there is no evidence of disease.

# Complete Remission/Response, Partial Remission/Response, Progressive Disease, Stable Disease
tmp1$TREATMENT_OUTCOME_FIRST_COURSE <- ifelse(tmp1$TREATMENT_OUTCOME_FIRST_COURSE == "Stable Disease", 1, 
                                              ifelse(tmp1$TREATMENT_OUTCOME_FIRST_COURSE == "Progressive Disease", 0,
                                                     ifelse(tmp1$TREATMENT_OUTCOME_FIRST_COURSE == "Partial Remission/Response", 0,
                                                            ifelse(tmp1$TREATMENT_OUTCOME_FIRST_COURSE == "Complete Remission/Response", 0, NA ))))

tmp1$TREATMENT_OUTCOME_FIRST_COURSE.1 <- ifelse(tmp1$TREATMENT_OUTCOME_FIRST_COURSE.1 == "Stable Disease", 0, 
                                              ifelse(tmp1$TREATMENT_OUTCOME_FIRST_COURSE.1 == "Progressive Disease", 1,
                                                     ifelse(tmp1$TREATMENT_OUTCOME_FIRST_COURSE.1 == "Partial Remission/Response", 0,
                                                            ifelse(tmp1$TREATMENT_OUTCOME_FIRST_COURSE.1 == "Complete Remission/Response", 0, NA ))))

tmp1$TREATMENT_OUTCOME_FIRST_COURSE.2 <- ifelse(tmp1$TREATMENT_OUTCOME_FIRST_COURSE.2 == "Stable Disease", 0, 
                                              ifelse(tmp1$TREATMENT_OUTCOME_FIRST_COURSE.2 == "Progressive Disease", 0,
                                                     ifelse(tmp1$TREATMENT_OUTCOME_FIRST_COURSE.2 == "Partial Remission/Response", 1,
                                                            ifelse(tmp1$TREATMENT_OUTCOME_FIRST_COURSE.2 == "Complete Remission/Response", 0, NA ))))

tmp1$TREATMENT_OUTCOME_FIRST_COURSE.3 <- ifelse(tmp1$TREATMENT_OUTCOME_FIRST_COURSE.3 == "Stable Disease", 0, 
                                                ifelse(tmp1$TREATMENT_OUTCOME_FIRST_COURSE.3 == "Progressive Disease", 0,
                                                       ifelse(tmp1$TREATMENT_OUTCOME_FIRST_COURSE.3 == "Partial Remission/Response", 0,
                                                              ifelse(tmp1$TREATMENT_OUTCOME_FIRST_COURSE.3 == "Complete Remission/Response", 1, NA ))))

tmp1$KRAS_MUTATION <- ifelse(tmp1$KRAS_MUTATION == "NO", 0, ifelse(tmp1$KRAS_MUTATION == "YES", 1, NA))  # KRAS mutation

# Only for EGFR mutation (Exon 19 Deletion, L858R, L861Q, Other, T790M) 
tmp1$MUTATION_TYPE <- ifelse(tmp1$MUTATION_TYPE == "Other", 0, ifelse(tmp1$MUTATION_TYPE == "Exon 19 Deletion", 1,   
                                                                      ifelse(tmp1$MUTATION_TYPE == "T790M", 1,
                                                                             ifelse(tmp1$MUTATION_TYPE == "L858R", 1,
                                                                                    ifelse(tmp1$MUTATION_TYPE == "L861Q", 1, NA))))) # Double check this one: L861Q 

##################################################################################################
#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
# End
##################################################################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##################################################################################################
