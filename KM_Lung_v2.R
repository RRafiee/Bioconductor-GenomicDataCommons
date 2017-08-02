
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

#######################################################################
# Loading a csv file
Pathfolder <- "C:/Users/3052092/Documents/Live/"  # initialise the path of the csv file
csvfilename <- "km_luad_OS_Reza.csv"
pathcsvfile <- paste(Pathfolder,csvfilename,sep = "")
datafileObject <- read.csv(pathcsvfile, header=T)  # LUAD samples, TCGA

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Overall Survival using KM to get an estimate of the survival curve

tmp <- datafileObject[,c(20,21,27)] 
tmp$OS_STATUS <- ifelse(tmp$OS_STATUS == "LIVING", 0, 1)
tmp$OS_MONTHS <- tmp$OS_MONTHS/12

# Only two groups: "123" and "4"
idx_call_123 <- c(which(tmp$DDRD.call == "1"),which(tmp$DDRD.call == "2"), which(tmp$DDRD.call == "3"))
tmp$DDRD.call[idx_call_123] <- "123"

# Only two groups: "12" and "34"
idx_call_12 <- c(which(tmp$DDRD.call == "1"),which(tmp$DDRD.call == "2"))
idx_call_34 <- c(which(tmp$DDRD.call == "3"),which(tmp$DDRD.call == "4"))
tmp$DDRD.call[idx_call_12] <- "12"
tmp$DDRD.call[idx_call_34] <- "34"

# tmp <- tmp[-c(which(tmp$DDRD.call == "4")),] # excluding group 4
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
f <- npsurv(Surv(OS_MONTHS, OS_STATUS) ~ DDRD.call, data=tmp)

# Get pvalue and see if there is a survival difference between subgorups
# H0 (in rho=0): no difference in survival functionbetween groups

survdiff(Surv(OS_MONTHS, OS_STATUS) ~ DDRD.call, data=tmp, rho = 0) # p= 0.765 (we don't reject the null hypothesis) so there is no difference in survival among all 4 subgroups   

# 4 subgroups
# n=507, 10 observations deleted due to missingness.
# N Observed Expected (O-E)^2/E (O-E)^2/V
# DDRD.call=1 126       46     48.2    0.1000    0.1369
# DDRD.call=2 129       45     44.0    0.0222    0.0294
# DDRD.call=3 127       46     40.7    0.6835    0.8863
# DDRD.call=4 125       46     50.1    0.3308    0.4578
# Chisq= 1.1  on 3 degrees of freedom, p= 0.765 

# 2 subgroups (123 vs. 4)
# n=507, 10 observations deleted due to missingness.
# N Observed Expected (O-E)^2/E (O-E)^2/V
# DDRD.call=123 382      137    132.9     0.125     0.458
# DDRD.call=4   125       46     50.1     0.331     0.458
# Chisq= 0.5  on 1 degrees of freedom, p= 0.499 

# 2 subgroups (12 vs. 34)
# n=507, 10 observations deleted due to missingness.
# N Observed Expected (O-E)^2/E (O-E)^2/V
# DDRD.call=12 255       91     92.2    0.0158    0.0319
# DDRD.call=34 252       92     90.8    0.0160    0.0319
# Chisq= 0  on 1 degrees of freedom, p= 0.858 

cox1 <- coxph(Surv(OS_MONTHS, OS_STATUS) ~ DDRD.call, method = "exact", data = tmp)
print(summary(cox1))

# 4 subgroups
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

# 2 subgorups (123 vs. 4)
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

# 2 subgorups (12 vs. 34)
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

dev.off()
# Colour transparent and then add lines, with censor points.
survplot(f, conf="none",  n.risk=TRUE, lty=1, lwd=2, col="#FFFFFF00", 
         xlab="Time (years)", ylab="Overall survival", label.curves=T,  time.inc=2)

lines(f, mark.time=T, lwd=2, col=c("blue", "red","orange","darkgreen")) # DDRD.call=1 (blue), DDRD.call=2 (red), DDRD.call=3 (orange), DDRD.call=4 (darkgreen)
leg.txt <- c("DDRD.call=1", "DDRD.call=2","DDRD.call=3","DDRD.call=4") 
legend(17, 0.9, leg.txt, lty = 1, lwd = 2, col = c("blue", "red","orange","darkgreen"))
title("LUAD, n=507, 10 observations deleted due to missingness")

# when using 2 subgroups (123 vs. 4)
lines(f, mark.time=T, lwd=2, col=c("black","darkgreen")) # DDRD.call=123 (black), DDRD.call=4 (darkgreen)
leg.txt <- c("DDRD.call=123", "DDRD.call=4") 
legend(16, 0.9, leg.txt, lty = 1, lwd = 2, col = c("black","darkgreen"))
title("LUAD, n=507, 10 observations deleted due to missingness")

# when using 2 subgroups (12 vs. 34)
lines(f, mark.time=T, lwd=2, col=c("purple","green")) # DDRD.call=12 (purple), DDRD.call=34 (green)
leg.txt <- c("DDRD.call=12", "DDRD.call=34") 
legend(16, 0.9, leg.txt, lty = 1, lwd = 2, col = c("purple","green"))
title("LUAD, n=507, 10 observations deleted due to missingness")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Disease-free Survival

tmp <- datafileObject[,c(12,13,27)] # Overall survival
tmp$DFS_STATUS <- ifelse(tmp$DFS_STATUS == "DiseaseFree", 0, 1)
tmp$DFS_MONTHS <- tmp$DFS_MONTHS/12

# Only two groups: "123" and "4"
idx_call_123 <- c(which(tmp$DDRD.call == "1"),which(tmp$DDRD.call == "2"), which(tmp$DDRD.call == "3"))
tmp$DDRD.call[idx_call_123] <- "123"

# Only two groups: "12" and "34"
idx_call_12 <- c(which(tmp$DDRD.call == "1"),which(tmp$DDRD.call == "2"))
idx_call_34 <- c(which(tmp$DDRD.call == "3"),which(tmp$DDRD.call == "4"))
tmp$DDRD.call[idx_call_12] <- "12"
tmp$DDRD.call[idx_call_34] <- "34"

# tmp <- tmp[-c(which(tmp$DDRD.call == "4")),] # excluding group 4

f <- npsurv(Surv(DFS_MONTHS, DFS_STATUS) ~ DDRD.call, data=tmp)

# Get pvalue and see if there is a survival difference between subgorups
# H0 (in rho=0): no difference in survival functionbetween groups

survdiff(Surv(DFS_MONTHS, DFS_STATUS) ~ DDRD.call, data=tmp, rho = 0) # p= 0.674 (we don't reject the null hypothesis) so there is no difference in survival among all 4 subgroups   

# 4 subgroups
# n=431, 86 observations deleted due to missingness.
# 
# N Observed Expected (O-E)^2/E (O-E)^2/V
# DDRD.call=1 106       43     48.0   0.51643   0.70608
# DDRD.call=2 108       45     45.3   0.00264   0.00352
# DDRD.call=3 112       51     44.4   0.96666   1.28275
# DDRD.call=4 105       46     47.2   0.03207   0.04319
# 
# Chisq= 1.5  on 3 degrees of freedom, p= 0.674 

# 2 subgroups (123 vs. 4) 
# n=431, 86 observations deleted due to missingness.
# N Observed Expected (O-E)^2/E (O-E)^2/V
# DDRD.call=123 326      139    137.8    0.0110    0.0432
# DDRD.call=4   105       46     47.2    0.0321    0.0432
# Chisq= 0  on 1 degrees of freedom, p= 0.835 

# 2 subgroups (12 vs. 34) 
# n=431, 86 observations deleted due to missingness.
# N Observed Expected (O-E)^2/E (O-E)^2/V
# DDRD.call=12 214       88     93.3     0.304     0.615
# DDRD.call=34 217       97     91.7     0.309     0.615
# Chisq= 0.6  on 1 degrees of freedom, p= 0.433 

cox1 <- coxph(Surv(DFS_MONTHS, DFS_STATUS) ~ DDRD.call, method = "exact", data = tmp)
print(summary(cox1))

# 4 subgroup
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

# 2 subgroups (123 vs. 4)
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

# 2 subgroups (12 vs. 34)
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

dev.off()
# Colour transparent and then add lines, with censor points.
survplot(f, conf="none",  n.risk=TRUE, lty=1, lwd=2, col="#FFFFFF00", 
         xlab="Time (years)", ylab="Disease-free survival", label.curves=T,  time.inc=2)

lines(f, mark.time=T, lwd=2, col=c("orange","darkgreen","red","blue"))
leg.txt <- c("DDRD.call=1", "DDRD.call=2","DDRD.call=3","DDRD.call=4") 
legend(17, 0.9, leg.txt, lty = 1, lwd = 2, col = c("blue", "red","orange","darkgreen"))
title("LUAD, n=431, 86 observations deleted due to missingness.")

# when using 2 subgroups (123 vs. 4)
lines(f, mark.time=T, lwd=2, col=c("black","darkgreen"))
leg.txt <- c("DDRD.call=123", "DDRD.call=4") 
legend(16, 0.9, leg.txt, lty = 1, lwd = 2, col = c("black","darkgreen"))
title("LUAD, n=431, 86 observations deleted due to missingness.")

# when using 2 subgroups (12 vs. 34)
lines(f, mark.time=T, lwd=2, col=c("purple","green"))
leg.txt <- c("DDRD.call=12", "DDRD.call=34") 
legend(16, 0.9, leg.txt, lty = 1, lwd = 2, col = c("purple","green"))
title("LUAD, n=431, 86 observations deleted due to missingness.")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Survival modelling

Pathfolder <- "C:/Users/3052092/Documents/Live/"  # initialise the path of the csv file
csvfilename <- "km_luad_OS_Reza.csv"
pathcsvfile <- paste(Pathfolder,csvfilename,sep = "")
datafileObject <- read.csv(pathcsvfile, header=T)  # LUAD samples, TCGA

tmp <- datafileObject[,c(20,21,27)] # Overall survival
tmp$OS_STATUS <- ifelse(tmp$OS_STATUS == "LIVING", 0, 1)
tmp$OS_MONTHS <- tmp$OS_MONTHS/12

# Checking the prorpotionality of hazard
#var1_ph <- cox.zph(coxph(Surv(OS_Time, OS_YN) ~ Variable1, data = datafileObject))
var1_ph <- cox.zph(coxph(Surv(OS_MONTHS, OS_STATUS) ~ DDRD.call, data = tmp))
var1_ph_pvalue <- var1_ph$table[3] # p-value: [1] 0.02712499 is significant ===> reject the null hypothesis (Doesn't pass the test for the proportionality of hazard)

#var2_ph <- cox.zph(coxph(Surv(OS_MONTHS, OS_STATUS) ~ DDRD.call.1, data = tmp))
#var2_ph_pvalue <- var2_ph$table[3] # p-value: [1] 0.7164459 is not significant ==> doesn't reject the null hypothesis (Pass the test for the proportionality of hazard)
#par(mfrow=c(2,2))

plot(main= "proportional hazard test", cox.zph(coxph(Surv(OS_MONTHS, OS_STATUS) ~ DDRD.call, data = tmp))) # variable1
legend(bty="n","topleft",lty=0.1, col="blue", paste("p-value: ",round(var1_ph_pvalue,digits = 3)))
#plot(main= "proportional hazard test", cox.zph(coxph(Surv(OS_Time, OS_YN) ~ DDRD.call.1, data = datafileObject))) # variable2
#legend(bty="n","topleft",lty=0.1, col="blue", paste("p-value: ",round(var2_ph_pvalue,digits = 3)))


# Building a survival model (Variable1 is DN_MBEN and Variable2 is TRvs.NonTR)
# Survival_model <- coxph(Surv(OS_Time, OS_YN) ~ Variable1 +  Variable2, data = datafileObject)
Survival_model <- coxph(Surv(OS_MONTHS, OS_STATUS) ~ DDRD.call, data = tmp)
summary(Survival_model)


#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
# End
##################################################################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##################################################################################################