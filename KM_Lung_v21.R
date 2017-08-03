
##################################################################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##################################################################################################
# Start
#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
# Dr Reza Rafiee, August 2017
# Research Fellow, Queen's University Belfast
# Survival analysis on lung cancer data (LUAD and LUSC)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# install_all_packages_automatic <- function(x) {
#   x <- as.character(substitute(x))
#   if(isTRUE(x %in% .packages(all.available=TRUE))) {
#     eval(parse(text = sprintf("require(\"%s\")", x)))
#   } else {
#     #update.packages(ask= FALSE) #update installed packages.
#     eval(parse(text = sprintf("install.packages(\"%s\", dependencies = TRUE)", x)))
#   }
#   if(isTRUE(x %in% .packages(all.available=TRUE))) {
#     eval(parse(text = sprintf("require(\"%s\")", x)))
#   } else {
#     source("http://bioconductor.org/biocLite.R")
#     #biocLite(character(), ask=FALSE) #update installed packages.
#     eval(parse(text = sprintf("biocLite(\"%s\")", x)))
#     eval(parse(text = sprintf("require(\"%s\")", x)))
#   }
# }
# 
# Install all required packages as well as dependencies
# install_all_packages_automatic(pec)
# install_all_packages_automatic(pROC)
# install_all_packages_automatic(survival)
# install_all_packages_automatic(survivalROC)
# install_all_packages_automatic("rms")

# Loading all survival and dependent libraries
library(survival)
library(survivalROC)
library(pec)
library(pROC)
library(rms)

##############################   LUAD  ################################
# Loading a csv file
Pathfolder <- "C:/Users/3052092/Documents/Live/"  # initialise the path of the csv file
csvfilename <- "LUAD_Samples_DDRD_ordered_23June17.csv"
pathcsvfile <- paste(Pathfolder,csvfilename,sep = "")
datafileObject <- read.csv(pathcsvfile, header=T)  # LUSC samples, TCGA

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##############################   LUSC  ################################
# Loading a csv file
Pathfolder <- "C:/Users/3052092/Documents/Live/"  # initialise the path of the csv file
csvfilename <- "LUSC_Samples_DDRD_ordered_23June17.csv"
pathcsvfile <- paste(Pathfolder,csvfilename,sep = "")
datafileObject <- read.csv(pathcsvfile, header=T)  # LUSC samples, TCGA
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Overall Survival using KM to get an estimate of the survival curve
tmp <- datafileObject[,c(51,52,83)] 
tmp$OS_STATUS <- ifelse(tmp$OS_STATUS == "LIVING", 0, ifelse(tmp$OS_STATUS == "DECEASED", 1, NA))
tmp$OS_MONTHS <- tmp$OS_MONTHS/12
#######################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Overall Survival using KM to get an estimate of the survival curve
tmp <- datafileObject[,c(50,51,81)] 
tmp$OS_STATUS <- ifelse(tmp$OS_STATUS == "LIVING", 0, ifelse(tmp$OS_STATUS == "DECEASED", 1, NA))
tmp$OS_MONTHS <- tmp$OS_MONTHS/12
colnames(tmp)[3] <- "DDRD.call"  # change "DDRD.quartile" to "DDRD.call" 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#######################################################################

# Only two groups: "123" and "4"
idx_call_123 <- c(which(tmp$DDRD.call == "1"),which(tmp$DDRD.call == "2"), which(tmp$DDRD.call == "3"))
tmp$DDRD.call[idx_call_123] <- "123"

# Only two groups: "12" and "34"
idx_call_12 <- c(which(tmp$DDRD.call == "1"),which(tmp$DDRD.call == "2"))
idx_call_34 <- c(which(tmp$DDRD.call == "3"),which(tmp$DDRD.call == "4"))
tmp$DDRD.call[idx_call_12] <- "12"
tmp$DDRD.call[idx_call_34] <- "34"

# Only two groups: "1" and "234"
idx_call_234 <- c(which(tmp$DDRD.call == "2"),which(tmp$DDRD.call == "3"), which(tmp$DDRD.call == "4"))
tmp$DDRD.call[idx_call_234] <- "234"


# tmp <- tmp[-c(which(tmp$DDRD.call == "4")),] # excluding group 4
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#######################################################################

f <- npsurv(Surv(OS_MONTHS, OS_STATUS) ~ DDRD.call, data=tmp)
# f$type: "right"
# datafileObject[as.integer(f$na.action),51]

# Get pvalue and see if there is a survival difference between subgorups
# H0 (in rho=0): no difference in survival functionbetween groups

survdiff(Surv(OS_MONTHS, OS_STATUS) ~ DDRD.call, data=tmp, rho = 0) # p= 0.765 (we don't reject the null hypothesis) so there is no difference in survival among all 4 subgroups   

# 4 subgroups (LUAD)
# n=507, 10 observations deleted due to missingness.
# N Observed Expected (O-E)^2/E (O-E)^2/V
# DDRD.call=1 126       46     48.2    0.1000    0.1369
# DDRD.call=2 129       45     44.0    0.0222    0.0294
# DDRD.call=3 127       46     40.7    0.6835    0.8863
# DDRD.call=4 125       46     50.1    0.3308    0.4578
# Chisq= 1.1  on 3 degrees of freedom, p= 0.765 
# index of missing samples- LUAD,OS (datafileObject[as.integer(f$na.action),51]): 150 184 185 190 191 193 194 197 199 517
# as.character(datafileObject[as.integer(f$na.action),2])
# [1] "TCGA.80.5607.01" "TCGA.75.7031.01" "TCGA.75.7030.01" "TCGA.75.6211.01" "TCGA.75.6207.01" "TCGA.75.6205.01"
# [7] "TCGA.75.6203.01" "TCGA.75.5126.01" "TCGA.75.5122.01" "TCGA.05.4244.01"


# 4 subgroups (LUSC)
# n=494, 7 observations deleted due to missingness.
# N Observed Expected (O-E)^2/E (O-E)^2/V
# DDRD.quartile=1 134       58     58.2  0.000811   0.00112
# DDRD.quartile=2 116       56     42.8  4.079454   5.15259
# DDRD.quartile=3 112       45     53.4  1.314110   1.77067
# DDRD.quartile=4 132       52     56.6  0.376898   0.52703
# Chisq= 5.8  on 3 degrees of freedom, p= 0.121 

# 2 subgroups (LUAD: 123 vs. 4)
# n=507, 10 observations deleted due to missingness.
# N Observed Expected (O-E)^2/E (O-E)^2/V
# DDRD.call=123 382      137    132.9     0.125     0.458
# DDRD.call=4   125       46     50.1     0.331     0.458
# Chisq= 0.5  on 1 degrees of freedom, p= 0.499 

# 2 subgroups (LUSC: 123 vs. 4)
# n=494, 7 observations deleted due to missingness.
# N Observed Expected (O-E)^2/E (O-E)^2/V
# DDRD.call=123 362      159    154.4     0.138     0.527
# DDRD.call=4   132       52     56.6     0.377     0.527
# Chisq= 0.5  on 1 degrees of freedom, p= 0.468 

# 2 subgroups (LUAD: 12 vs. 34)
# n=507, 10 observations deleted due to missingness.
# N Observed Expected (O-E)^2/E (O-E)^2/V
# DDRD.call=12 255       91     92.2    0.0158    0.0319
# DDRD.call=34 252       92     90.8    0.0160    0.0319
# Chisq= 0  on 1 degrees of freedom, p= 0.858 

# 2 subgroups (LUSC: 12 vs. 34)
# n=494, 7 observations deleted due to missingness.
# N Observed Expected (O-E)^2/E (O-E)^2/V
# DDRD.call=12 250      114      101      1.67      3.24
# DDRD.call=34 244       97      110      1.54      3.24
# Chisq= 3.2  on 1 degrees of freedom, p= 0.072 

# 2 subgorups (LUAD: 1 vs. 234)
# n=494, 7 observations deleted due to missingness.
# N Observed Expected (O-E)^2/E (O-E)^2/V
# DDRD.call=1   134       58     58.2  0.000811   0.00112
# DDRD.call=234 360      153    152.8  0.000309   0.00112
# Chisq= 0  on 1 degrees of freedom, p= 0.973 
 
# 2 subgroups (LUSC: 1 vs. 234)
# n=507, 10 observations deleted due to missingness.
# N Observed Expected (O-E)^2/E (O-E)^2/V
# DDRD.call=1   126       46     48.2    0.1000     0.137
# DDRD.call=234 381      137    134.8    0.0358     0.137
# Chisq= 0.1  on 1 degrees of freedom, p= 0.711 

cox1 <- coxph(Surv(OS_MONTHS, OS_STATUS) ~ DDRD.call, method = "exact", data = tmp)
print(summary(cox1))

# 4 subgroups (LUAD)
# n= 507, number of events= 183 
# (10 observations deleted due to missingness)
# coef exp(coef)  se(coef)      z Pr(>|z|)
# DDRD.call -0.002774  0.997230  0.064404 -0.043    0.966
# exp(coef) exp(-coef) lower .95 upper .95
# DDRD.call    0.9972      1.003     0.879     1.131
# Rsquare= 0   (max possible= 0.979 )
# Likelihood ratio test= 0  on 1 df,   p=0.9656
# Wald test            = 0  on 1 df,   p=0.9656
# Score (logrank) test = 0  on 1 df,   p=0.9656

# 4 subgroups (LUSC)
# n= 494, number of events= 211 
# (7 observations deleted due to missingness)
# coef exp(coef) se(coef)      z Pr(>|z|)
# DDRD.quartile -0.06236   0.93955  0.05989 -1.041    0.298
# exp(coef) exp(-coef) lower .95 upper .95
# DDRD.quartile    0.9395      1.064    0.8355     1.057
# Rsquare= 0.002   (max possible= 0.989 )
# Likelihood ratio test= 1.08  on 1 df,   p=0.2976
# Wald test            = 1.08  on 1 df,   p=0.2978
# Score (logrank) test = 1.09  on 1 df,   p=0.2975


# 2 subgorups (LUAD: 123 vs. 4)
# n= 507, number of events= 183 
# (10 observations deleted due to missingness)
# coef exp(coef) se(coef)      z Pr(>|z|)
# DDRD.call4 -0.1155    0.8909   0.1708 -0.676    0.499
# exp(coef) exp(-coef) lower .95 upper .95
# DDRD.call4    0.8909      1.122    0.6374     1.245
# Rsquare= 0.001   (max possible= 0.979 )
# Likelihood ratio test= 0.47  on 1 df,   p=0.4949
# Wald test            = 0.46  on 1 df,   p=0.4989
# Score (logrank) test = 0.46  on 1 df,   p=0.4987

# 2 subgorups (LUSC: 123 vs. 4)
# n= 494, number of events= 211 
# (7 observations deleted due to missingness)
# coef exp(coef) se(coef)      z Pr(>|z|)
# DDRD.call4 -0.1173    0.8893   0.1617 -0.726    0.468
# exp(coef) exp(-coef) lower .95 upper .95
# DDRD.call4    0.8893      1.125    0.6477     1.221
# Rsquare= 0.001   (max possible= 0.989 )
# Likelihood ratio test= 0.54  on 1 df,   p=0.4637
# Wald test            = 0.53  on 1 df,   p=0.4681
# Score (logrank) test = 0.53  on 1 df,   p=0.4679

# 2 subgorups (LUAD: 12 vs. 34)
# n= 507, number of events= 183 
# (10 observations deleted due to missingness)
# coef exp(coef) se(coef)     z Pr(>|z|)
# DDRD.call34 0.02647   1.02683  0.14816 0.179    0.858
# exp(coef) exp(-coef) lower .95 upper .95
# DDRD.call34     1.027     0.9739     0.768     1.373
# Rsquare= 0   (max possible= 0.979 )
# Likelihood ratio test= 0.03  on 1 df,   p=0.8582
# Wald test            = 0.03  on 1 df,   p=0.8582
# Score (logrank) test = 0.03  on 1 df,   p=0.8582

# 2 subgorups (LUSC: 12 vs. 34)
# coef exp(coef) se(coef)      z Pr(>|z|)  
# DDRD.call34 -0.2491    0.7795   0.1389 -1.794   0.0728 .
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# exp(coef) exp(-coef) lower .95 upper .95
# DDRD.call34    0.7795      1.283    0.5938     1.023
# Rsquare= 0.007   (max possible= 0.989 )
# Likelihood ratio test= 3.23  on 1 df,   p=0.07212
# Wald test            = 3.22  on 1 df,   p=0.07277
# Score (logrank) test = 3.24  on 1 df,   p=0.07204

# 2 subgorups (LUAD: 1 vs. 234)
# n= 507, number of events= 183 
# (10 observations deleted due to missingness)
# coef exp(coef) se(coef)    z Pr(>|z|)
# DDRD.call234 0.06331   1.06535  0.17112 0.37    0.711
# exp(coef) exp(-coef) lower .95 upper .95
# DDRD.call234     1.065     0.9387    0.7618      1.49
# Rsquare= 0   (max possible= 0.979 )
# Likelihood ratio test= 0.14  on 1 df,   p=0.71
# Wald test            = 0.14  on 1 df,   p=0.7114
# Score (logrank) test = 0.14  on 1 df,   p=0.7114

# 2 subgorups (LUSC: 1 vs. 234)
# n= 494, number of events= 211 
# (7 observations deleted due to missingness)
# coef exp(coef) se(coef)     z Pr(>|z|)
# DDRD.call234 0.00518   1.00519  0.15447 0.034    0.973
# exp(coef) exp(-coef) lower .95 upper .95
# DDRD.call234     1.005     0.9948    0.7426     1.361
# Rsquare= 0   (max possible= 0.989 )
# Likelihood ratio test= 0  on 1 df,   p=0.9732
# Wald test            = 0  on 1 df,   p=0.9732
# Score (logrank) test = 0  on 1 df,   p=0.9732


dev.off()
# Colour transparent and then add lines, with censor points.
survplot(f, conf="none",  n.risk=TRUE, lty=1, lwd=2, col="#FFFFFF00", 
         xlab="Time (years)", ylab="Overall survival", label.curves=T,  time.inc=2)

lines(f, mark.time=T, lwd=2, col=c("blue", "red","orange","darkgreen")) # DDRD.call=1 (blue), DDRD.call=2 (red), DDRD.call=3 (orange), DDRD.call=4 (darkgreen)
leg.txt <- c("DDRD.call=1", "DDRD.call=2","DDRD.call=3","DDRD.call=4") 

legend(17, 0.9, leg.txt, lty = 1, lwd = 2, col = c("blue", "red","orange","darkgreen")) # LUAD
legend(11, 0.9, leg.txt, lty = 1, lwd = 2, col = c("blue", "red","orange","darkgreen")) # LUSC

title("LUAD, n=507, 10 observations deleted due to missingness") # LUAD
title("LUSC, n=494, 7 observations deleted due to missingness")  # LUSC

# when using 2 subgroups (123 vs. 4)
lines(f, mark.time=T, lwd=2, col=c("black","darkgreen")) # DDRD.call=123 (black), DDRD.call=4 (darkgreen)
leg.txt <- c("DDRD.call=123", "DDRD.call=4") 
legend(16, 0.9, leg.txt, lty = 1, lwd = 2, col = c("black","darkgreen")) # LUAD
legend(10, 0.9, leg.txt, lty = 1, lwd = 2, col = c("black","darkgreen")) # LUSC
title("LUAD, n=507, 10 observations deleted due to missingness") # LUAD
title("LUSC, n=494, 7 observations deleted due to missingness") # LUAD

# when using 2 subgroups (12 vs. 34)
lines(f, mark.time=T, lwd=2, col=c("purple","green")) # DDRD.call=12 (purple), DDRD.call=34 (green)
leg.txt <- c("DDRD.call=12", "DDRD.call=34") 

legend(16, 0.9, leg.txt, lty = 1, lwd = 2, col = c("purple","green"))
legend(10, 0.9, leg.txt, lty = 1, lwd = 2, col = c("purple","green"))

title("LUAD, n=507, 10 observations deleted due to missingness")
title("LUSC, n=494, 7 observations deleted due to missingness")

# when using 2 subgroups (1 vs. 234)
lines(f, mark.time=T, lwd=2, col=c("blue","gold")) # DDRD.call=123 (black), DDRD.call=4 (darkgreen)
leg.txt <- c("DDRD.call=1", "DDRD.call=234") 
legend(16, 0.9, leg.txt, lty = 1, lwd = 2, col = c("blue","gold")) # LUAD
legend(10, 0.9, leg.txt, lty = 1, lwd = 2, col = c("blue","gold")) # LUSC
title("LUAD, n=507, 10 observations deleted due to missingness")
title("LUSC, n=494, 7 observations deleted due to missingness")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Disease-free Survival - LUAD
tmp <- datafileObject[,c(19,20,83)] # Overall survival
tmp$DFS_STATUS <- ifelse(tmp$DFS_STATUS == "DiseaseFree", 0, ifelse(tmp$DFS_STATUS == "Recurred/Progressed",1, NA))
tmp$DFS_MONTHS <- tmp$DFS_MONTHS/12
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Disease-free Survival - LUSC
tmp <- datafileObject[,c(19,20,81)] # Overall survival
tmp$DFS_STATUS <- ifelse(tmp$DFS_STATUS == "DiseaseFree", 0, ifelse(tmp$DFS_STATUS == "Recurred/Progressed",1, NA))
tmp$DFS_MONTHS <- tmp$DFS_MONTHS/12
colnames(tmp)[3] <- "DDRD.call"  # change "DDRD.quartile" to "DDRD.call" 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Only two groups: "123" and "4"
idx_call_123 <- c(which(tmp$DDRD.call == "1"),which(tmp$DDRD.call == "2"), which(tmp$DDRD.call == "3"))
tmp$DDRD.call[idx_call_123] <- "123"

# Only two groups: "12" and "34"
idx_call_12 <- c(which(tmp$DDRD.call == "1"),which(tmp$DDRD.call == "2"))
idx_call_34 <- c(which(tmp$DDRD.call == "3"),which(tmp$DDRD.call == "4"))
tmp$DDRD.call[idx_call_12] <- "12"
tmp$DDRD.call[idx_call_34] <- "34"

# Only two groups: "1" and "234"
idx_call_234 <- c(which(tmp$DDRD.call == "2"),which(tmp$DDRD.call == "3"), which(tmp$DDRD.call == "4"))
tmp$DDRD.call[idx_call_234] <- "234"


#tmp <- tmp[-c(which(tmp$DDRD.call == "234")),] # excluding group 4

f <- npsurv(Surv(DFS_MONTHS, DFS_STATUS) ~ DDRD.call, data=tmp)

# Get pvalue and see if there is a survival difference between subgorups
# H0 (in rho=0): no difference in survival functionbetween groups

survdiff(Surv(DFS_MONTHS, DFS_STATUS) ~ DDRD.call, data=tmp, rho = 0) # p= 0.674 (we don't reject the null hypothesis) so there is no difference in survival among all 4 subgroups   

# 4 subgroups (LUAD)
# n=431, 86 observations deleted due to missingness.
# 
# N Observed Expected (O-E)^2/E (O-E)^2/V
# DDRD.call=1 106       43     48.0   0.51643   0.70608
# DDRD.call=2 108       45     45.3   0.00264   0.00352
# DDRD.call=3 112       51     44.4   0.96666   1.28275
# DDRD.call=4 105       46     47.2   0.03207   0.04319
# 
# Chisq= 1.5  on 3 degrees of freedom, p= 0.674 

# 4 subgroups (LUSC)
# n=374, 127 observations deleted due to missingness.
# N Observed Expected (O-E)^2/E (O-E)^2/V
# DDRD.call=1 103       43     33.3  2.831130   3.81413
# DDRD.call=2  83       29     29.2  0.000836   0.00108
# DDRD.call=3  83       25     34.8  2.760128   3.79987
# DDRD.call=4 105       33     32.8  0.001885   0.00253
# Chisq= 5.6  on 3 degrees of freedom, p= 0.131 


# 2 subgroups (LUAD: 123 vs. 4) 
# n=431, 86 observations deleted due to missingness.
# N Observed Expected (O-E)^2/E (O-E)^2/V
# DDRD.call=123 326      139    137.8    0.0110    0.0432
# DDRD.call=4   105       46     47.2    0.0321    0.0432
# Chisq= 0  on 1 degrees of freedom, p= 0.835 

# 2 subgroups (LUSC: 123 vs. 4) 
# n=374, 127 observations deleted due to missingness.
# N Observed Expected (O-E)^2/E (O-E)^2/V
# DDRD.call=123 269       97     97.2  0.000635   0.00253
# DDRD.call=4   105       33     32.8  0.001885   0.00253
# Chisq= 0  on 1 degrees of freedom, p= 0.96 

# 2 subgroups (LUAD: 12 vs. 34) 
# n=431, 86 observations deleted due to missingness.
# N Observed Expected (O-E)^2/E (O-E)^2/V
# DDRD.call=12 214       88     93.3     0.304     0.615
# DDRD.call=34 217       97     91.7     0.309     0.615
# Chisq= 0.6  on 1 degrees of freedom, p= 0.433 

# 2 subgroups (LUSC: 12 vs. 34) 
# n=374, 127 observations deleted due to missingness.
# N Observed Expected (O-E)^2/E (O-E)^2/V
# DDRD.call=12 186       72     62.4      1.46      2.82
# DDRD.call=34 188       58     67.6      1.35      2.82
# Chisq= 2.8  on 1 degrees of freedom, p= 0.0929 

# 2 subgroups (LUAD: 1 vs. 234) 
# n=431, 86 observations deleted due to missingness.
# N Observed Expected (O-E)^2/E (O-E)^2/V
# DDRD.call=1   106       43       48     0.516     0.706
# DDRD.call=234 325      142      137     0.181     0.706
# Chisq= 0.7  on 1 degrees of freedom, p= 0.401 

# 2 subgroups (LUSC: 1 vs. 234) 
# n=374, 127 observations deleted due to missingness.
# N Observed Expected (O-E)^2/E (O-E)^2/V
# DDRD.call=1   103       43     33.3     2.831      3.81
# DDRD.call=234 271       87     96.7     0.975      3.81
# Chisq= 3.8  on 1 degrees of freedom, p= 0.0508 

cox1 <- coxph(Surv(DFS_MONTHS, DFS_STATUS) ~ DDRD.call, method = "exact", data = tmp)
print(summary(cox1))

# 4 subgroup (LUAD)
# n= 431, number of events= 185 
# (86 observations deleted due to missingness)
# coef exp(coef) se(coef)     z Pr(>|z|)
# DDRD.call 0.03851   1.03926  0.06517 0.591    0.555
# exp(coef) exp(-coef) lower .95 upper .95
# DDRD.call     1.039     0.9622    0.9146     1.181
# Rsquare= 0.001   (max possible= 0.988 )
# Likelihood ratio test= 0.35  on 1 df,   p=0.5545
# Wald test            = 0.35  on 1 df,   p=0.5546
# Score (logrank) test = 0.35  on 1 df,   p=0.5545

# 4 subgroup (LUSC)
# n= 374, number of events= 130 
# (127 observations deleted due to missingness)
# coef exp(coef) se(coef)      z Pr(>|z|)
# DDRD.call -0.11589   0.89057  0.07826 -1.481    0.139
# exp(coef) exp(-coef) lower .95 upper .95
# DDRD.call    0.8906      1.123    0.7639     1.038
# Rsquare= 0.006   (max possible= 0.971 )
# Likelihood ratio test= 2.2  on 1 df,   p=0.1379
# Wald test            = 2.19  on 1 df,   p=0.1387
# Score (logrank) test = 2.2  on 1 df,   p=0.1379

# 2 subgroups (LUAD: 123 vs. 4)
# n= 431, number of events= 185 
# (86 observations deleted due to missingness)
# coef exp(coef) se(coef)      z Pr(>|z|)
# DDRD.call4 -0.0354    0.9652   0.1703 -0.208    0.835
# exp(coef) exp(-coef) lower .95 upper .95
# DDRD.call4    0.9652      1.036    0.6912     1.348
# Rsquare= 0   (max possible= 0.988 )
# Likelihood ratio test= 0.04  on 1 df,   p=0.8349
# Wald test            = 0.04  on 1 df,   p=0.8354
# Score (logrank) test = 0.04  on 1 df,   p=0.8354

# 2 subgroups (LUSC: 123 vs. 4)
# n= 374, number of events= 130 
# (127 observations deleted due to missingness)
# coef exp(coef) se(coef)    z Pr(>|z|)
# DDRD.call4 0.01015   1.01020  0.20182 0.05     0.96
# exp(coef) exp(-coef) lower .95 upper .95
# DDRD.call4      1.01     0.9899    0.6802       1.5
# Rsquare= 0   (max possible= 0.971 )
# Likelihood ratio test= 0  on 1 df,   p=0.9599
# Wald test            = 0  on 1 df,   p=0.9599
# Score (logrank) test = 0  on 1 df,   p=0.9599

# 2 subgroups (LUAD: 12 vs. 34)
# n= 431, number of events= 185 
# (86 observations deleted due to missingness)
# coef exp(coef) se(coef)     z Pr(>|z|)
# DDRD.call34 0.1157    1.1226   0.1475 0.784    0.433
# exp(coef) exp(-coef) lower .95 upper .95
# DDRD.call34     1.123     0.8908    0.8407     1.499
# Rsquare= 0.001   (max possible= 0.988 )
# Likelihood ratio test= 0.62  on 1 df,   p=0.4327
# Wald test            = 0.61  on 1 df,   p=0.433
# Score (logrank) test = 0.62  on 1 df,   p=0.4327

# 2 subgroups (LUSC: 12 vs. 34)
# n= 374, number of events= 130 
# (127 observations deleted due to missingness)
# coef exp(coef) se(coef)      z Pr(>|z|)  
# DDRD.call34 -0.2959    0.7439   0.1768 -1.674   0.0941 .
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# exp(coef) exp(-coef) lower .95 upper .95
# DDRD.call34    0.7439      1.344    0.5261     1.052
# Rsquare= 0.008   (max possible= 0.971 )
# Likelihood ratio test= 2.82  on 1 df,   p=0.09299
# Wald test            = 2.8  on 1 df,   p=0.09413
# Score (logrank) test = 2.82  on 1 df,   p=0.09294

# 2 subgroups (LUAD: 1 vs. 234)
# coef exp(coef) se(coef)    z Pr(>|z|)
# DDRD.call234 0.1470    1.1584   0.1751 0.84    0.401
# exp(coef) exp(-coef) lower .95 upper .95
# DDRD.call234     1.158     0.8633    0.8219     1.633
# Rsquare= 0.002   (max possible= 0.988 )
# Likelihood ratio test= 0.72  on 1 df,   p=0.3952
# Wald test            = 0.7  on 1 df,   p=0.4012
# Score (logrank) test = 0.71  on 1 df,   p=0.4007

# 2 subgroups (LUSC: 1 vs. 234)
# coef exp(coef) se(coef)      z Pr(>|z|)  
# DDRD.call234 -0.3625    0.6959   0.1867 -1.942   0.0521 .
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# exp(coef) exp(-coef) lower .95 upper .95
# DDRD.call234    0.6959      1.437    0.4827     1.003
# Rsquare= 0.01   (max possible= 0.971 )
# Likelihood ratio test= 3.61  on 1 df,   p=0.05752
# Wald test            = 3.77  on 1 df,   p=0.05209
# Score (logrank) test = 3.81  on 1 df,   p=0.05082

dev.off()
# Colour transparent and then add lines, with censor points.
survplot(f, conf="none",  n.risk=TRUE, lty=1, lwd=2, col="#FFFFFF00", 
         xlab="Time (years)", ylab="Disease-free survival", label.curves=T,  time.inc=2)

# 4 subgroups
lines(f, mark.time=T, lwd=2, col=c("blue","red", "orange" ,"darkgreen"))
leg.txt <- c("DDRD.call=1", "DDRD.call=2","DDRD.call=3","DDRD.call=4") 

legend(17, 0.9, leg.txt, lty = 1, lwd = 2, col = c("blue", "red","orange","darkgreen")) # LUAD
legend(10, 0.9, leg.txt, lty = 1, lwd = 2, col = c("blue", "red","orange","darkgreen")) # LUSC
title("LUAD, n=431, 86 observations deleted due to missingness.")
title("LUSC, n=374, 127 observations deleted due to missingness.")

# when using 2 subgroups (123 vs. 4)
lines(f, mark.time=T, lwd=2, col=c("black","darkgreen"))
leg.txt <- c("DDRD.call=123", "DDRD.call=4") 
legend(16, 0.9, leg.txt, lty = 1, lwd = 2, col = c("black","darkgreen")) #LUAD
legend(10, 0.9, leg.txt, lty = 1, lwd = 2, col = c("black","darkgreen")) #LUSC

title("LUAD, n=431, 86 observations deleted due to missingness.")
title("LUSC, n=374, 127 observations deleted due to missingness.")

# when using 2 subgroups (12 vs. 34)
lines(f, mark.time=T, lwd=2, col=c("purple","green"))
leg.txt <- c("DDRD.call=12", "DDRD.call=34") 
legend(16, 0.9, leg.txt, lty = 1, lwd = 2, col = c("purple","green")) # LUAD
legend(10, 0.9, leg.txt, lty = 1, lwd = 2, col = c("purple","green")) # LUSC
title("LUAD, n=431, 86 observations deleted due to missingness.")
title("LUSC, n=374, 127 observations deleted due to missingness.")

# when using 2 subgroups (1 vs. 234)
lines(f, mark.time=T, lwd=2, col=c("blue","gold"))
leg.txt <- c("DDRD.call=1", "DDRD.call=234") 
legend(16, 0.9, leg.txt, lty = 1, lwd = 2, col = c("blue","gold")) # LUAD
legend(10, 0.9, leg.txt, lty = 1, lwd = 2, col = c("blue","gold")) # LUSC
title("LUAD, n=431, 86 observations deleted due to missingness.")
title("LUSC, n=374, 127 observations deleted due to missingness.")


# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # Survival modelling
# Pathfolder <- "C:/Users/3052092/Documents/Live/"  # initialise the path of the csv file
# csvfilename <- "km_luad_OS_Reza.csv"
# pathcsvfile <- paste(Pathfolder,csvfilename,sep = "")
# datafileObject <- read.csv(pathcsvfile, header=T)  # LUAD samples, TCGA
# 
# tmp <- datafileObject[,c(20,21,27)] # Overall survival
# tmp$OS_STATUS <- ifelse(tmp$OS_STATUS == "LIVING", 0, 1)
# tmp$OS_MONTHS <- tmp$OS_MONTHS/12
# 
# # Checking the prorpotionality of hazard
# #var1_ph <- cox.zph(coxph(Surv(OS_Time, OS_YN) ~ Variable1, data = datafileObject))
# var1_ph <- cox.zph(coxph(Surv(OS_MONTHS, OS_STATUS) ~ DDRD.call, data = tmp))
# var1_ph_pvalue <- var1_ph$table[3] # p-value: [1] 0.02712499 is significant ===> reject the null hypothesis (Doesn't pass the test for the proportionality of hazard)
# 
# #var2_ph <- cox.zph(coxph(Surv(OS_MONTHS, OS_STATUS) ~ DDRD.call.1, data = tmp))
# #var2_ph_pvalue <- var2_ph$table[3] # p-value: [1] 0.7164459 is not significant ==> doesn't reject the null hypothesis (Pass the test for the proportionality of hazard)
# #par(mfrow=c(2,2))
# 
# plot(main= "proportional hazard test", cox.zph(coxph(Surv(OS_MONTHS, OS_STATUS) ~ DDRD.call, data = tmp))) # variable1
# legend(bty="n","topleft",lty=0.1, col="blue", paste("p-value: ",round(var1_ph_pvalue,digits = 3)))
# #plot(main= "proportional hazard test", cox.zph(coxph(Surv(OS_Time, OS_YN) ~ DDRD.call.1, data = datafileObject))) # variable2
# #legend(bty="n","topleft",lty=0.1, col="blue", paste("p-value: ",round(var2_ph_pvalue,digits = 3)))
# 
# # Building a survival model (Variable1 is DN_MBEN and Variable2 is TRvs.NonTR)
# # Survival_model <- coxph(Surv(OS_Time, OS_YN) ~ Variable1 +  Variable2, data = datafileObject)
# Survival_model <- coxph(Surv(OS_MONTHS, OS_STATUS) ~ DDRD.call, data = tmp)
# summary(Survival_model)

#######################################################################

#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
# End
##################################################################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##################################################################################################