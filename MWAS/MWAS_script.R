#MWAS script
##The method we're gonna use is Bayesian Ridge Regression
###The covariates are the class effect of management group and covariates of expected heterozygosity and breed proportions
##Follow: http://morotalab.org/apsc5984-2020/day41/day41.html
setwd("/work/wtp/alakamp/Subsetting/MWAS")
library(data.table)
library(dplyr)
library(stringr)
library(BGLR)
library(coda)

set.seed(3746)
#Read in phenotype data and ORF relative abundance table
data = fread("Rumen_data_genotyped.txt")
#data$AnimalID = as.character(data$AnimalID)

#Feed to gain as our efficiency trait
data$FtG = data$ADDMI_keep / data$ADG_keep

ORF = fread("ORF_RA_table_1.0.csv", header = T) %>%
  filter(., Sample_ID %in% data$AnimalID)

animal_names = ORF[,1]
ORF[,1] = NULL
ORF = as.matrix(ORF)
rownames(ORF) = animal_names$Sample_ID
ORF = ORF[match(data$AnimalID, rownames(ORF)), ]

total_m2 = matrix(nrow = 3) %>%
  as.data.frame()
rownames(total_m2) = c("lower_limit", "mean", "upper_limit")

#All available data####
#Center and scale ORF matrix
##Centered by Gaussian standard deviation bc 
ORFcs = scale(ORF, center = FALSE, scale = apply(ORF, 2, sd, na.rm = TRUE))

##ADDMI####
full_ADDMI_ETA<-list( list(~factor(Management_Group) + heterosis +
                perANS + perHHS + perARS +
                perSHS + perDSS + perBMS + 
                perBRS + perBNS + perSGS + 
                perBVS + perCHS + perCAS + 
                perGVS + perLMS + perMAS +        
                perSAS + perSMS + perTAS + 
                perBVSo + perCHSo + perLMSo +
                perHH +  perAN +  perSM +  
                perCH +  perM2 +  perM3 + 
                perRS +  perBV +  perBR +  
                perXB +  perRO ,
                data=data,model="FIXED"),
           list(X = ORFcs, model="BRR")
)

full_ADDMI_model<-BGLR(y=data$ADDMI_keep,ETA=full_ADDMI_ETA, 
                       nIter=100000, burnIn=5000, thin = 20, saveAt = "full_ADDMI_")

#Diagnostics
list.files()
# Residual variance
varE_full_ADDMI <-scan("full_ADDMI_varE.dat")
plot(varE_full_ADDMI,type="o",col=2,cex=.5,
     ylab=expression(sigma[epsilon]^2),
     xlab="Sample",main="Residual Variance");
abline(h=full_ADDMI_model$varE,col=4,lwd=2);
abline(v=full_ADDMI_model$burnIn/full_ADDMI_model$thin,col=4)
effectiveSize(mcmc(varE_full_ADDMI))

# Metagenomic Variance across chain
lambda_full_ADDMI <-scan("full_ADDMI_ETA_2_varB.dat")
plot(lambda_full_ADDMI, type="o", col=2, cex=.5,
     xlab="Sample",ylab=expression(lambda),
     main="Metagenome Variance");
abline(h=full_ADDMI_model$ETA[[2]]$varB,col=4,lwd=2);
abline(v=full_ADDMI_model$burnIn/full_ADDMI_model$thin,col=4)
effectiveSize(mcmc(lambda_full_ADDMI))

#microbiome variance
sigma2m_full_ADDMI <- ncol(ORFcs) * full_ADDMI_model$ETA[[2]]$varB
# m2
sigma2m_full_ADDMI / (sigma2m_full_ADDMI + full_ADDMI_model$varE) #0.17

#Mean and sd m^2
m2_full_ADDMI =  (ncol(ORFcs) * lambda_full_ADDMI) / ( (ncol(ORFcs) * lambda_full_ADDMI) + varE_full_ADDMI) 
total_m2$Full_ADDMI = c(round(mean(m2_full_ADDMI) - 2*sd(m2_full_ADDMI), 2),
                      round(mean(m2_full_ADDMI), 2),
                      round(mean(m2_full_ADDMI) + 2*sd(m2_full_ADDMI), 2))
total_m2$V1 = NULL

#ORF effects
full_ADDMI_b <- full_ADDMI_model$ETA[[2]]$b
full_ADDMI_SDb = full_ADDMI_model$ETA[[2]]$SD.b
plot(full_ADDMI_b ^2, ylab="Estimated Squared-ORF Effect",
     type="o",cex=.5,col="red",main="Full ADDMI ORF Effects",
     xlab="ORF")
points(full_ADDMI_b^2,cex=0.5,col="blue")

full_ADDMI_ORF_effects = full_ADDMI_model$ETA[[2]]$b %>%
  as.data.frame()
colnames(full_ADDMI_ORF_effects) = "Effect"
full_ADDMI_ORF_effects$Effect_sq = (full_ADDMI_ORF_effects$Effect)^2
full_ADDMI_ORF_effects = full_ADDMI_ORF_effects %>%
  arrange(desc(Effect_sq))
full_ADDMI_ORF_effects$ORF_ID = rownames(full_ADDMI_ORF_effects)
full_ADDMI_ORF_effects = full_ADDMI_ORF_effects %>%
  relocate(ORF_ID, .before = Effect)

##ADG####
full_ADG_ETA<-list( list(~factor(Management_Group) + heterosis +
                             perANS + perHHS + perARS +
                             perSHS + perDSS + perBMS + 
                             perBRS + perBNS + perSGS + 
                             perBVS + perCHS + perCAS + 
                             perGVS + perLMS + perMAS +        
                             perSAS + perSMS + perTAS + 
                             perBVSo + perCHSo + perLMSo +
                             perHH +  perAN +  perSM +  
                             perCH +  perM2 +  perM3 + 
                             perRS +  perBV +  perBR +  
                             perXB +  perRO ,
                           data=data,model="FIXED"),
                      list(X = ORFcs, model="BRR")
)

full_ADG_model<-BGLR(y=data$ADG_keep,ETA=full_ADG_ETA, 
                     nIter=100000, burnIn=5000, thin = 20, saveAt = "full_ADG_")

#Diagnostics
list.files()
# Residual variance
varE_full_ADG <-scan("full_ADG_varE.dat")
plot(varE_full_ADG, type="o", col=2, cex=.5,
     ylab=expression(sigma[epsilon]^2),
     xlab="Sample",main="Residual Variance");
abline(h=full_ADG_model$varE,col=4,lwd=2);
abline(v=full_ADG_model$burnIn/full_ADG_model$thin,col=4)
effectiveSize(mcmc(varE_full_ADG))

# Metagenomic Variance across chain
lambda_full_ADG <-scan("full_ADG_ETA_2_varB.dat")
plot(lambda_full_ADG, type="o", col=2, cex=.5,
     xlab="Sample",ylab=expression(lambda),
     main="Metagenome Variance");
abline(h=full_ADG_model$ETA[[2]]$varB,col=4,lwd=2);
abline(v=full_ADG_model$burnIn/full_ADG_model$thin,col=4)
effectiveSize(mcmc(lambda_full_ADG))

#microbiome variance
sigma2m_full_ADG <- ncol(ORFcs) * full_ADG_model$ETA[[2]]$varB
# m2
sigma2m_full_ADG / (sigma2m_full_ADG + full_ADG_model$varE) #0.34

#Mean and sd m^2
m2_full_ADG =  (ncol(ORFcs) * lambda_full_ADG) / ( (ncol(ORFcs) * lambda_full_ADG) + varE_full_ADG) 
total_m2$Full_ADG = c(round(mean(m2_full_ADG) - 2*sd(m2_full_ADG), 2),
                      round(mean(m2_full_ADG), 2),
                      round(mean(m2_full_ADG) + 2*sd(m2_full_ADG), 2))

#ORF Effects
full_ADG_b <- full_ADG_model$ETA[[2]]$b
full_ADG_SDb = full_ADG_model$ETA[[2]]$SD.b
plot(full_ADG_b ^2, ylab="Estimated Squared-ORF Effect",
     type="o",cex=.5,col="red",main="Full ADG ORF Effects",
     xlab="ORF")
points(full_ADG_b^2,cex=0.5,col="blue")

full_ADG_ORF_effects = full_ADG_model$ETA[[2]]$b %>%
  as.data.frame()
colnames(full_ADG_ORF_effects) = "Effect"
full_ADG_ORF_effects$Effect_sq = (full_ADG_ORF_effects$Effect)^2
full_ADG_ORF_effects = full_ADG_ORF_effects %>%
  arrange(desc(Effect_sq))

##FtG####
full_FtG_ETA<-list( list(~factor(Management_Group) + heterosis +
                           perANS + perHHS + perARS +
                           perSHS + perDSS + perBMS + 
                           perBRS + perBNS + perSGS + 
                           perBVS + perCHS + perCAS + 
                           perGVS + perLMS + perMAS +        
                           perSAS + perSMS + perTAS + 
                           perBVSo + perCHSo + perLMSo +
                           perHH +  perAN +  perSM +  
                           perCH +  perM2 +  perM3 + 
                           perRS +  perBV +  perBR +  
                           perXB +  perRO ,
                         data=data,model="FIXED"),
                    list(X = ORFcs, model="BRR")
)

full_FtG_model<-BGLR(y=data$FtG,ETA=full_FtG_ETA, 
                     nIter=100000, burnIn=5000, thin = 20, saveAt = "full_FtG_")

#Diagnostics
list.files()
# Residual variance
varE_full_FtG <-scan("full_FtG_varE.dat")
plot(varE_full_FtG, type="o", col=2, cex=.5,
     ylab=expression(sigma[epsilon]^2),
     xlab="Sample",main="Residual Variance");
abline(h=full_FtG_model$varE,col=4,lwd=2);
abline(v=full_FtG_model$burnIn/full_FtG_model$thin,col=4)
effectiveSize(mcmc(varE_full_FtG))

# Metagenomic Variance across chain
lambda_full_FtG <-scan("full_FtG_ETA_2_varB.dat")
plot(lambda_full_FtG, type="o", col=2, cex=.5,
     xlab="Sample",ylab=expression(lambda),
     main="Metagenome Variance");
abline(h=full_FtG_model$ETA[[2]]$varB,col=4,lwd=2);
abline(v=full_FtG_model$burnIn/full_FtG_model$thin,col=4)
effectiveSize(mcmc(lambda_full_FtG))

#microbiome variance
sigma2m_full_FtG <- ncol(ORFcs) * full_FtG_model$ETA[[2]]$varB
# m2
sigma2m_full_FtG / (sigma2m_full_FtG + full_FtG_model$varE) #0.17

#Mean and sd m^2
m2_full_FtG =  (ncol(ORFcs) * lambda_full_FtG) / ( (ncol(ORFcs) * lambda_full_FtG) + varE_full_FtG) 
total_m2$Full_FtG = c(round(mean(m2_full_FtG) - 2*sd(m2_full_FtG), 2),
                      round(mean(m2_full_FtG), 2),
                      round(mean(m2_full_FtG) + 2*sd(m2_full_FtG), 2))

#ORF Effects
full_FtG_b <- full_FtG_model$ETA[[2]]$b
full_FtG_SDb = full_FtG_model$ETA[[2]]$SD.b
plot(full_FtG_b ^2, ylab="Estimated Squared-ORF Effect",
     type="o",cex=.5,col="red",main="Full FtG ORF Effects",
     xlab="ORF")
points(full_FtG_b^2,cex=0.5,col="blue")

full_FtG_ORF_effects = full_FtG_model$ETA[[2]]$b %>%
  as.data.frame()
colnames(full_FtG_ORF_effects) = "Effect"
full_FtG_ORF_effects$Effect_sq = (full_FtG_ORF_effects$Effect)^2
full_FtG_ORF_effects = full_FtG_ORF_effects %>%
  arrange(desc(Effect_sq))
full_FtG_ORF_effects$ORF_ID = rownames(full_FtG_ORF_effects)
full_FtG_ORF_effects = full_FtG_ORF_effects %>%
  relocate(ORF_ID, .before = Effect)

#Cocentrate data####
#Subset matrix to concentrate animals (steers)
concentrate_data = data %>%
  filter(., Gender == "Steer")

#Get ORF of only animals in this diet group
ORFcs = ORF[rownames(ORF) %in% concentrate_data$AnimalID,]

#Remove any ORF which is NOT present in this diet group
ORFcs = ORFcs[,colSums(ORFcs) != 0]

#Center and scale remaining ORF
ORFcs = scale(ORFcs, center = FALSE, scale = apply(ORFcs, 2, sd, na.rm = TRUE))

#Logic checks
if(!(identical(rownames(ORFcs), as.character(concentrate_data$AnimalID)))){
  stop("Data order does NOT match!")
}

if(!(identical(nrow(ORFcs), length(concentrate_data$AnimalID)))){
  stop("Data order does NOT match!")
}

##ADDMI####
concentrate_ADDMI_ETA<-list( list(~factor(Management_Group) + heterosis +
                             perANS + perHHS + perARS +
                             perSHS + perDSS + perBMS + 
                             perBRS + perBNS + perSGS + 
                             perBVS + perCHS + perCAS + 
                             perGVS + perLMS + perMAS +        
                             perSAS + perSMS + perTAS + 
                             perBVSo + perCHSo + perLMSo +
                             perHH +  perAN +  perSM +  
                             perCH +  perM2 +  perM3 + 
                             perRS +  perBV +  perBR +  
                             perXB +  perRO ,
                           data=concentrate_data, model="FIXED"),
                      list(X = ORFcs, model="BRR")
)

concentrate_ADDMI_model<-BGLR(y=concentrate_data$ADDMI_keep,ETA=concentrate_ADDMI_ETA, 
                       nIter=100000, burnIn=5000, thin = 20, saveAt = "concentrate_ADDMI_")

#Diagnostics
list.files()
# Residual variance
varE_concentrate_ADDMI <-scan("concentrate_ADDMI_varE.dat")
plot(varE_concentrate_ADDMI,type="o",col=2,cex=.5,
     ylab=expression(sigma[epsilon]^2),
     xlab="Sample",main="Residual Variance");
abline(h=concentrate_ADDMI_model$varE,col=4,lwd=2);
abline(v=concentrate_ADDMI_model$burnIn/concentrate_ADDMI_model$thin,col=4)
effectiveSize(mcmc(varE_concentrate_ADDMI))

# Metagenomic Variance across chain
lambda_concentrate_ADDMI <-scan("concentrate_ADDMI_ETA_2_varB.dat")
plot(lambda_concentrate_ADDMI, type="o", col=2, cex=.5,
     xlab="Sample",ylab=expression(lambda),
     main="Metagenome Variance");
abline(h=concentrate_ADDMI_model$ETA[[2]]$varB,col=4,lwd=2);
abline(v=concentrate_ADDMI_model$burnIn/concentrate_ADDMI_model$thin,col=4)
effectiveSize(mcmc(lambda_concentrate_ADDMI))

#microbiome variance
sigma2m_concentrate_ADDMI <- ncol(ORFcs) * concentrate_ADDMI_model$ETA[[2]]$varB
# m2
sigma2m_concentrate_ADDMI / (sigma2m_concentrate_ADDMI + concentrate_ADDMI_model$varE) #0.29

#Mean and sd m^2
m2_concentrate_ADDMI =  (ncol(ORFcs) * lambda_concentrate_ADDMI) / ( (ncol(ORFcs) * lambda_concentrate_ADDMI) + varE_concentrate_ADDMI) 
total_m2$concentrate_ADDMI = c(round(mean(m2_concentrate_ADDMI) - 2*sd(m2_concentrate_ADDMI), 2),
                        round(mean(m2_concentrate_ADDMI), 2),
                        round(mean(m2_concentrate_ADDMI) + 2*sd(m2_concentrate_ADDMI), 2))

#ORF effects
concentrate_ADDMI_b <- concentrate_ADDMI_model$ETA[[2]]$b
concentrate_ADDMI_SDb = concentrate_ADDMI_model$ETA[[2]]$SD.b
plot(concentrate_ADDMI_b ^2, ylab="Estimated Squared-ORF Effect",
     type="o",cex=.5,col="red",main="concentrate ADDMI ORF Effects",
     xlab="ORF")
points(concentrate_ADDMI_b^2,cex=0.5,col="blue")

concentrate_ADDMI_ORF_effects = concentrate_ADDMI_model$ETA[[2]]$b %>%
  as.data.frame()
colnames(concentrate_ADDMI_ORF_effects) = "Effect"
concentrate_ADDMI_ORF_effects$Effect_sq = (concentrate_ADDMI_ORF_effects$Effect)^2
concentrate_ADDMI_ORF_effects = concentrate_ADDMI_ORF_effects %>%
  arrange(desc(Effect_sq))
concentrate_ADDMI_ORF_effects$ORF_ID = rownames(concentrate_ADDMI_ORF_effects)
concentrate_ADDMI_ORF_effects = concentrate_ADDMI_ORF_effects %>%
  relocate(ORF_ID, .before = Effect)

##ADG####
concentrate_ADG_ETA<-list( list(~factor(Management_Group) + heterosis +
                           perANS + perHHS + perARS +
                           perSHS + perDSS + perBMS + 
                           perBRS + perBNS + perSGS + 
                           perBVS + perCHS + perCAS + 
                           perGVS + perLMS + perMAS +        
                           perSAS + perSMS + perTAS + 
                           perBVSo + perCHSo + perLMSo +
                           perHH +  perAN +  perSM +  
                           perCH +  perM2 +  perM3 + 
                           perRS +  perBV +  perBR +  
                           perXB +  perRO ,
                         data=concentrate_data, model="FIXED"),
                    list(X = ORFcs, model="BRR")
)

concentrate_ADG_model<-BGLR(y=concentrate_data$ADG_keep,ETA=concentrate_ADG_ETA, 
                            nIter=100000, burnIn=5000, thin = 20, saveAt = "concentrate_ADG_")

#Diagnostics
list.files()
# Residual variance
varE_concentrate_ADG <-scan("concentrate_ADG_varE.dat")
plot(varE_concentrate_ADG, type="o", col=2, cex=.5,
     ylab=expression(sigma[epsilon]^2),
     xlab="Sample",main="Residual Variance");
abline(h=concentrate_ADG_model$varE,col=4,lwd=2);
abline(v=concentrate_ADG_model$burnIn/concentrate_ADG_model$thin,col=4)
effectiveSize(mcmc(varE_concentrate_ADG))

# Metagenomic Variance across chain
lambda_concentrate_ADG <-scan("concentrate_ADG_ETA_2_varB.dat")
plot(lambda_concentrate_ADG, type="o", col=2, cex=.5,
     xlab="Sample",ylab=expression(lambda),
     main="Metagenome Variance");
abline(h=concentrate_ADG_model$ETA[[2]]$varB,col=4,lwd=2);
abline(v=concentrate_ADG_model$burnIn/concentrate_ADG_model$thin,col=4)
effectiveSize(mcmc(lambda_concentrate_ADG))

#microbiome variance
sigma2m_concentrate_ADG <- ncol(ORFcs) * concentrate_ADG_model$ETA[[2]]$varB
# m2
sigma2m_concentrate_ADG / (sigma2m_concentrate_ADG + concentrate_ADG_model$varE) #0.34

#Mean and sd m^2
m2_concentrate_ADG =  (ncol(ORFcs) * lambda_concentrate_ADG) / ( (ncol(ORFcs) * lambda_concentrate_ADG) + varE_concentrate_ADG) 
total_m2$concentrate_ADG = c(round(mean(m2_concentrate_ADG) - 2*sd(m2_concentrate_ADG), 2),
                      round(mean(m2_concentrate_ADG), 2),
                      round(mean(m2_concentrate_ADG) + 2*sd(m2_concentrate_ADG), 2))

#ORF Effects
concentrate_ADG_b <- concentrate_ADG_model$ETA[[2]]$b
concentrate_ADG_SDb = concentrate_ADG_model$ETA[[2]]$SD.b
plot(concentrate_ADG_b ^2, ylab="Estimated Squared-ORF Effect",
     type="o",cex=.5,col="red",main="concentrate ADG ORF Effects",
     xlab="ORF")
points(concentrate_ADG_b^2,cex=0.5,col="blue")

concentrate_ADG_ORF_effects = concentrate_ADG_model$ETA[[2]]$b %>%
  as.data.frame()
colnames(concentrate_ADG_ORF_effects) = "Effect"
concentrate_ADG_ORF_effects$Effect_sq = (concentrate_ADG_ORF_effects$Effect)^2
concentrate_ADG_ORF_effects = concentrate_ADG_ORF_effects %>%
  arrange(desc(Effect_sq))
full_ADG_ORF_effects$ORF_ID = rownames(full_ADG_ORF_effects)
full_ADG_ORF_effects = full_ADG_ORF_effects %>%
  relocate(ORF_ID, .before = Effect)
concentrate_ADG_ORF_effects$ORF_ID = rownames(concentrate_ADG_ORF_effects)
concentrate_ADG_ORF_effects = concentrate_ADG_ORF_effects %>%
  relocate(ORF_ID, .before = Effect)

##FtG####
concentrate_FtG_ETA<-list( list(~factor(Management_Group) + heterosis +
                           perANS + perHHS + perARS +
                           perSHS + perDSS + perBMS + 
                           perBRS + perBNS + perSGS + 
                           perBVS + perCHS + perCAS + 
                           perGVS + perLMS + perMAS +        
                           perSAS + perSMS + perTAS + 
                           perBVSo + perCHSo + perLMSo +
                           perHH +  perAN +  perSM +  
                           perCH +  perM2 +  perM3 + 
                           perRS +  perBV +  perBR +  
                           perXB +  perRO ,
                         data=concentrate_data, model="FIXED"),
                    list(X = ORFcs, model="BRR")
)

concentrate_FtG_model<-BGLR(y=concentrate_data$FtG,ETA=concentrate_FtG_ETA, 
                            nIter=100000, burnIn=5000, thin = 20, saveAt = "concentrate_FtG_")

#Diagnostics
list.files()
# Residual variance
varE_concentrate_FtG <-scan("concentrate_FtG_varE.dat")
plot(varE_concentrate_FtG, type="o", col=2, cex=.5,
     ylab=expression(sigma[epsilon]^2),
     xlab="Sample",main="Residual Variance");
abline(h=concentrate_FtG_model$varE,col=4,lwd=2);
abline(v=concentrate_FtG_model$burnIn/concentrate_FtG_model$thin,col=4)
effectiveSize(mcmc(varE_concentrate_FtG))

# Metagenomic Variance across chain
lambda_concentrate_FtG <-scan("concentrate_FtG_ETA_2_varB.dat")
plot(lambda_concentrate_FtG, type="o", col=2, cex=.5,
     xlab="Sample",ylab=expression(lambda),
     main="Metagenome Variance");
abline(h=concentrate_FtG_model$ETA[[2]]$varB,col=4,lwd=2);
abline(v=concentrate_FtG_model$burnIn/concentrate_FtG_model$thin,col=4)
effectiveSize(mcmc(lambda_concentrate_FtG))

#microbiome variance
sigma2m_concentrate_FtG <- ncol(ORFcs) * concentrate_FtG_model$ETA[[2]]$varB
# m2
sigma2m_concentrate_FtG / (sigma2m_concentrate_FtG + concentrate_FtG_model$varE) #0.17

#Mean and sd m^2
m2_concentrate_FtG =  (ncol(ORFcs) * lambda_concentrate_FtG) / ( (ncol(ORFcs) * lambda_concentrate_FtG) + varE_concentrate_FtG) 
total_m2$concentrate_FtG = c(round(mean(m2_concentrate_FtG) - 2*sd(m2_concentrate_FtG), 2),
                      round(mean(m2_concentrate_FtG), 2),
                      round(mean(m2_concentrate_FtG) + 2*sd(m2_concentrate_FtG), 2))

#ORF Effects
concentrate_FtG_b <- concentrate_FtG_model$ETA[[2]]$b
concentrate_FtG_SDb = concentrate_FtG_model$ETA[[2]]$SD.b
plot(concentrate_FtG_b ^2, ylab="Estimated Squared-ORF Effect",
     type="o",cex=.5,col="red",main="concentrate FtG ORF Effects",
     xlab="ORF")
points(concentrate_FtG_b^2,cex=0.5,col="blue")

concentrate_FtG_ORF_effects = concentrate_FtG_model$ETA[[2]]$b %>%
  as.data.frame()
colnames(concentrate_FtG_ORF_effects) = "Effect"
concentrate_FtG_ORF_effects$Effect_sq = (concentrate_FtG_ORF_effects$Effect)^2
concentrate_FtG_ORF_effects = concentrate_FtG_ORF_effects %>%
  arrange(desc(Effect_sq))
concentrate_FtG_ORF_effects$ORF_ID = rownames(concentrate_FtG_ORF_effects)
concentrate_FtG_ORF_effects = concentrate_FtG_ORF_effects %>%
  relocate(ORF_ID, .before = Effect)

#Forage data####
#Subset data to forage animals (heifers)
forage_data = data %>%
  filter(., Gender == "Heifer")

#Get ORF of only animals in this diet group
ORFcs = ORF[rownames(ORF) %in% forage_data$AnimalID,]

#Remove any ORF which is NOT present in this diet group
ORFcs = ORFcs[,colSums(ORFcs) != 0]

#Center and scale remaining ORF
ORFcs = scale(ORFcs, center = FALSE, scale = apply(ORFcs, 2, sd, na.rm = TRUE))

#Logic checks
if(!(identical(rownames(ORFcs), as.character(forage_data$AnimalID)))){
  stop("Data order does NOT match!")
}

if(!(identical(nrow(ORFcs), length(forage_data$AnimalID)))){
  stop("Data order does NOT match!")
}

##ADDMI####
forage_ADDMI_ETA<-list( list(~factor(Management_Group) + heterosis +
                             perANS + perHHS + perARS +
                             perSHS + perDSS + perBMS + 
                             perBRS + perBNS + perSGS + 
                             perBVS + perCHS + perCAS + 
                             perGVS + perLMS + perMAS +        
                             perSAS + perSMS + perTAS + 
                             perBVSo + perCHSo + perLMSo +
                             perHH +  perAN +  perSM +  
                             perCH +  perM2 +  perM3 + 
                             perRS +  perBV +  perBR +  
                             perXB +  perRO ,
                           data=forage_data, model="FIXED"),
                      list(X = ORFcs, model="BRR")
)

forage_ADDMI_model<-BGLR(y=forage_data$ADDMI_keep,ETA=forage_ADDMI_ETA, 
                         nIter=100000, burnIn=5000, thin = 20, saveAt = "forage_ADDMI_")

#Diagnostics
list.files()
# Residual variance
varE_forage_ADDMI <-scan("forage_ADDMI_varE.dat")
plot(varE_forage_ADDMI,type="o",col=2,cex=.5,
     ylab=expression(sigma[epsilon]^2),
     xlab="Sample",main="Residual Variance");
abline(h=forage_ADDMI_model$varE,col=4,lwd=2);
abline(v=forage_ADDMI_model$burnIn/forage_ADDMI_model$thin,col=4)
effectiveSize(mcmc(varE_forage_ADDMI))

# Metagenomic Variance across chain
lambda_forage_ADDMI <-scan("forage_ADDMI_ETA_2_varB.dat")
plot(lambda_forage_ADDMI, type="o", col=2, cex=.5,
     xlab="Sample",ylab=expression(lambda),
     main="Metagenome Variance");
abline(h=forage_ADDMI_model$ETA[[2]]$varB,col=4,lwd=2);
abline(v=forage_ADDMI_model$burnIn/forage_ADDMI_model$thin,col=4)
effectiveSize(mcmc(lambda_forage_ADDMI))

#microbiome variance
sigma2m_forage_ADDMI <- ncol(ORFcs) * forage_ADDMI_model$ETA[[2]]$varB
# m2
sigma2m_forage_ADDMI / (sigma2m_forage_ADDMI + forage_ADDMI_model$varE) #0.34

#Mean and sd m^2
m2_forage_ADDMI =  (ncol(ORFcs) * lambda_forage_ADDMI) / ( (ncol(ORFcs) * lambda_forage_ADDMI) + varE_forage_ADDMI) 
total_m2$forage_ADDMI = c(round(mean(m2_forage_ADDMI) - 2*sd(m2_forage_ADDMI), 2),
                        round(mean(m2_forage_ADDMI), 2),
                        round(mean(m2_forage_ADDMI) + 2*sd(m2_forage_ADDMI), 2))

#ORF effects
forage_ADDMI_b <- forage_ADDMI_model$ETA[[2]]$b
forage_ADDMI_SDb = forage_ADDMI_model$ETA[[2]]$SD.b
plot(forage_ADDMI_b ^2, ylab="Estimated Squared-ORF Effect",
     type="o",cex=.5,col="red",main="forage ADDMI ORF Effects",
     xlab="ORF")
points(forage_ADDMI_b^2,cex=0.5,col="blue")

forage_ADDMI_ORF_effects = forage_ADDMI_model$ETA[[2]]$b %>%
  as.data.frame()
colnames(forage_ADDMI_ORF_effects) = "Effect"
forage_ADDMI_ORF_effects$Effect_sq = (forage_ADDMI_ORF_effects$Effect)^2
forage_ADDMI_ORF_effects = forage_ADDMI_ORF_effects %>%
  arrange(desc(Effect_sq))
forage_ADDMI_ORF_effects$ORF_ID = rownames(forage_ADDMI_ORF_effects)
forage_ADDMI_ORF_effects = forage_ADDMI_ORF_effects %>%
  relocate(ORF_ID, .before = Effect)

##ADG####
forage_ADG_ETA<-list( list(~factor(Management_Group) + heterosis +
                           perANS + perHHS + perARS +
                           perSHS + perDSS + perBMS + 
                           perBRS + perBNS + perSGS + 
                           perBVS + perCHS + perCAS + 
                           perGVS + perLMS + perMAS +        
                           perSAS + perSMS + perTAS + 
                           perBVSo + perCHSo + perLMSo +
                           perHH +  perAN +  perSM +  
                           perCH +  perM2 +  perM3 + 
                           perRS +  perBV +  perBR +  
                           perXB +  perRO ,
                         data=forage_data, model="FIXED"),
                    list(X = ORFcs, model="BRR")
)

forage_ADG_model<-BGLR(y=forage_data$ADG_keep,ETA=forage_ADG_ETA, 
                       nIter=100000, burnIn=5000, thin = 20, saveAt = "forage_ADG_")

#Diagnostics
list.files()
# Residual variance
varE_forage_ADG <-scan("forage_ADG_varE.dat")
plot(varE_forage_ADG, type="o", col=2, cex=.5,
     ylab=expression(sigma[epsilon]^2),
     xlab="Sample",main="Residual Variance");
abline(h=forage_ADG_model$varE,col=4,lwd=2);
abline(v=forage_ADG_model$burnIn/forage_ADG_model$thin,col=4)
effectiveSize(mcmc(varE_forage_ADG))

# Metagenomic Variance across chain
lambda_forage_ADG <-scan("forage_ADG_ETA_2_varB.dat")
plot(lambda_forage_ADG, type="o", col=2, cex=.5,
     xlab="Sample",ylab=expression(lambda),
     main="Metagenome Variance");
abline(h=forage_ADG_model$ETA[[2]]$varB,col=4,lwd=2);
abline(v=forage_ADG_model$burnIn/forage_ADG_model$thin,col=4)
effectiveSize(mcmc(lambda_forage_ADG))

#microbiome variance
sigma2m_forage_ADG <- ncol(ORFcs) * forage_ADG_model$ETA[[2]]$varB
# m2
sigma2m_forage_ADG / (sigma2m_forage_ADG + forage_ADG_model$varE) #0.34

#Mean and sd m^2
m2_forage_ADG =  (ncol(ORFcs) * lambda_forage_ADG) / ( (ncol(ORFcs) * lambda_forage_ADG) + varE_forage_ADG) 
total_m2$forage_ADG = c(round(mean(m2_forage_ADG) - 2*sd(m2_forage_ADG), 2),
                      round(mean(m2_forage_ADG), 2),
                      round(mean(m2_forage_ADG) + 2*sd(m2_forage_ADG), 2))

#ORF Effects
forage_ADG_b <- forage_ADG_model$ETA[[2]]$b
forage_ADG_SDb = forage_ADG_model$ETA[[2]]$SD.b
plot(forage_ADG_b ^2, ylab="Estimated Squared-ORF Effect",
     type="o",cex=.5,col="red",main="forage ADG ORF Effects",
     xlab="ORF")
points(forage_ADG_b^2,cex=0.5,col="blue")

forage_ADG_ORF_effects = forage_ADG_model$ETA[[2]]$b %>%
  as.data.frame()
colnames(forage_ADG_ORF_effects) = "Effect"
forage_ADG_ORF_effects$Effect_sq = (forage_ADG_ORF_effects$Effect)^2
forage_ADG_ORF_effects = forage_ADG_ORF_effects %>%
  arrange(desc(Effect_sq))
forage_ADG_ORF_effects$ORF_ID = rownames(forage_ADG_ORF_effects)
forage_ADG_ORF_effects = forage_ADG_ORF_effects %>%
  relocate(ORF_ID, .before = Effect)

##FtG####
forage_FtG_ETA<-list( list(~factor(Management_Group) + heterosis +
                           perANS + perHHS + perARS +
                           perSHS + perDSS + perBMS + 
                           perBRS + perBNS + perSGS + 
                           perBVS + perCHS + perCAS + 
                           perGVS + perLMS + perMAS +        
                           perSAS + perSMS + perTAS + 
                           perBVSo + perCHSo + perLMSo +
                           perHH +  perAN +  perSM +  
                           perCH +  perM2 +  perM3 + 
                           perRS +  perBV +  perBR +  
                           perXB +  perRO ,
                         data=forage_data, model="FIXED"),
                    list(X = ORFcs, model="BRR")
)

forage_FtG_model<-BGLR(y=forage_data$FtG,ETA=forage_FtG_ETA, 
                       nIter=100000, burnIn=5000, thin = 20, saveAt = "forage_FtG_")

#Diagnostics
list.files()
# Residual variance
varE_forage_FtG <-scan("forage_FtG_varE.dat")
plot(varE_forage_FtG, type="o", col=2, cex=.5,
     ylab=expression(sigma[epsilon]^2),
     xlab="Sample",main="Residual Variance");
abline(h=forage_FtG_model$varE,col=4,lwd=2);
abline(v=forage_FtG_model$burnIn/forage_FtG_model$thin,col=4)
effectiveSize(mcmc(varE_forage_FtG))

# Metagenomic Variance across chain
lambda_forage_FtG <-scan("forage_FtG_ETA_2_varB.dat")
plot(lambda_forage_FtG, type="o", col=2, cex=.5,
     xlab="Sample",ylab=expression(lambda),
     main="Metagenome Variance");
abline(h=forage_FtG_model$ETA[[2]]$varB,col=4,lwd=2);
abline(v=forage_FtG_model$burnIn/forage_FtG_model$thin,col=4)
effectiveSize(mcmc(lambda_forage_FtG))

#microbiome variance
sigma2m_forage_FtG <- ncol(ORFcs) * forage_FtG_model$ETA[[2]]$varB
# m2
sigma2m_forage_FtG / (sigma2m_forage_FtG + forage_FtG_model$varE) #0.17

#Mean and sd m^2
m2_forage_FtG =  (ncol(ORFcs) * lambda_forage_FtG) / ( (ncol(ORFcs) * lambda_forage_FtG) + varE_forage_FtG) 
total_m2$forage_FtG = c(round(mean(m2_forage_FtG) - 2*sd(m2_forage_FtG), 2),
                      round(mean(m2_forage_FtG), 2),
                      round(mean(m2_forage_FtG) + 2*sd(m2_forage_FtG), 2))

#ORF Effects
forage_FtG_b <- forage_FtG_model$ETA[[2]]$b
forage_FtG_SDb = forage_FtG_model$ETA[[2]]$SD.b
plot(forage_FtG_b ^2, ylab="Estimated Squared-ORF Effect",
     type="o",cex=.5,col="red",main="forage FtG ORF Effects",
     xlab="ORF")
points(forage_FtG_b^2,cex=0.5,col="blue")

forage_FtG_ORF_effects = forage_FtG_model$ETA[[2]]$b %>%
  as.data.frame()
colnames(forage_FtG_ORF_effects) = "Effect"
forage_FtG_ORF_effects$Effect_sq = (forage_FtG_ORF_effects$Effect)^2
forage_FtG_ORF_effects = forage_FtG_ORF_effects %>%
  arrange(desc(Effect_sq))
forage_FtG_ORF_effects$ORF_ID = rownames(forage_FtG_ORF_effects)
forage_FtG_ORF_effects = forage_FtG_ORF_effects %>%
  relocate(ORF_ID, .before = Effect)

#Write out####
library(openxlsx)
write.xlsx(list("Microbiability" = total_m2, 
                "Full_ADDMI_ORF_Eff" = full_ADDMI_ORF_effects, 
                "Full_ADG_ORF_Eff" = full_ADG_ORF_effects,
                "Full_FtG_ORF_Eff" = full_FtG_ORF_effects,
                "Conc_ADDMI_ORF_Eff" = concentrate_ADDMI_ORF_effects, 
                "Conc_ADG_ORF_Eff" = concentrate_ADG_ORF_effects,
                "Conc_FtG_ORF_Eff" = concentrate_FtG_ORF_effects,
                "Forage_ADDMI_ORF_Eff" = forage_ADDMI_ORF_effects, 
                "Forage_ADG_ORF_Eff" = forage_ADG_ORF_effects,
                "Forage_FtG_ORF_Eff" = forage_FtG_ORF_effects),
                "Rumen_MWAS_Results.xlsx")
