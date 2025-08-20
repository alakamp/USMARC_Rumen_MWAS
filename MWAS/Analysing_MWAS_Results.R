setwd("~/Documents/PhD/Microbiome/Subsetting/MWAS")
library(data.table)
library(dplyr)
library(stringr)
library(openxlsx)

sheets = readxl::excel_sheets("Rumen_MWAS_Results.xlsx")

MWAS_sheets = list()
for(i in sheets){
  MWAS_sheets[[length(MWAS_sheets) + 1]] = read.xlsx("Rumen_MWAS_Results.xlsx", sheet = i)
}
names(MWAS_sheets) = sheets
list2env(MWAS_sheets, envir = globalenv())
rm(i, sheets, MWAS_sheets)

#Significant ORF overlap within set between traits####
##Significance defined by 1% top squared effect
library(ggvenn)

##Full####
Full_ORF = list(
  top_Full_ADDMI = Full_ADDMI_ORF_Eff %>%
    arrange(., desc(Effect_sq)) %>%
    filter(., Effect_sq > quantile(Effect_sq, 0.99)) %>%
    select(., ORF_ID) %>%
    unlist()
  ,
  top_Full_ADG = Full_ADG_ORF_Eff %>%
    arrange(., desc(Effect_sq)) %>%
    filter(., Effect_sq > quantile(Effect_sq, 0.99)) %>%
    select(., ORF_ID) %>%
    unlist()
  ,
  top_Full_FtG = Full_FtG_ORF_Eff %>%
    arrange(., desc(Effect_sq)) %>%
    filter(., Effect_sq > quantile(Effect_sq, 0.99)) %>%
    select(., ORF_ID) %>%
    unlist()
)

# All three
Full_all <- Reduce(intersect, Full_ORF)

# Pairs
ADDMI_ADG_Full <- setdiff(intersect(Full_ORF$top_Full_ADDMI, Full_ORF$top_Full_ADG), Full_all)
ADDMI_FtG_Full <- setdiff(intersect(Full_ORF$top_Full_ADDMI, Full_ORF$top_Full_FtG), Full_all)
ADG_FtG_Full <- setdiff(intersect(Full_ORF$top_Full_ADG, Full_ORF$top_Full_FtG), Full_all)

ggvenn(Full_ORF)

##Conc####
Conc_ORF = list(
  top_Conc_ADDMI = Conc_ADDMI_ORF_Eff %>%
    arrange(., desc(Effect_sq)) %>%
    filter(., Effect_sq > quantile(Effect_sq, 0.99)) %>%
    select(., ORF_ID) %>%
    unlist()
  ,
  top_Conc_ADG = Conc_ADG_ORF_Eff %>%
    arrange(., desc(Effect_sq)) %>%
    filter(., Effect_sq > quantile(Effect_sq, 0.99)) %>%
    select(., ORF_ID) %>%
    unlist()
  ,
  top_Conc_FtG = Conc_FtG_ORF_Eff %>%
    arrange(., desc(Effect_sq)) %>%
    filter(., Effect_sq > quantile(Effect_sq, 0.99)) %>%
    select(., ORF_ID) %>%
    unlist()
)

# All three
Conc_all <- Reduce(intersect, Conc_ORF)

# Pairs
ADDMI_ADG_Conc <- setdiff(intersect(Conc_ORF$top_Conc_ADDMI, Conc_ORF$top_Conc_ADG), Conc_all)
ADDMI_FtG_Conc <- setdiff(intersect(Conc_ORF$top_Conc_ADDMI, Conc_ORF$top_Conc_FtG), Conc_all)
ADG_FtG_Conc <- setdiff(intersect(Conc_ORF$top_Conc_ADG, Conc_ORF$top_Conc_FtG), Conc_all)

ggvenn(Conc_ORF)

##Forage####
Forage_ORF = list(
  top_Forage_ADDMI = Forage_ADDMI_ORF_Eff %>%
    arrange(., desc(Effect_sq)) %>%
    filter(., Effect_sq > quantile(Effect_sq, 0.99)) %>%
    select(., ORF_ID) %>%
    unlist()
  ,
  top_Forage_ADG = Forage_ADG_ORF_Eff %>%
    arrange(., desc(Effect_sq)) %>%
    filter(., Effect_sq > quantile(Effect_sq, 0.99)) %>%
    select(., ORF_ID) %>%
    unlist()
  ,
  top_Forage_FtG = Forage_FtG_ORF_Eff %>%
    arrange(., desc(Effect_sq)) %>%
    filter(., Effect_sq > quantile(Effect_sq, 0.99)) %>%
    select(., ORF_ID) %>%
    unlist()
)

# All three
Forage_all <- Reduce(intersect, Forage_ORF)

# Pairs
ADDMI_ADG_Forage <- setdiff(intersect(Forage_ORF$top_Forage_ADDMI, Forage_ORF$top_Forage_ADG), Forage_all)
ADDMI_FtG_Forage <- setdiff(intersect(Forage_ORF$top_Forage_ADDMI, Forage_ORF$top_Forage_FtG), Forage_all)
ADG_FtG_Forage <- setdiff(intersect(Forage_ORF$top_Forage_ADG, Forage_ORF$top_Forage_FtG), Forage_all)

ggvenn(Forage_ORF)

#Significant ORF overlap within trait between sets####
##Significance defined by 1% top squared effect
library(ggvenn)

##ADDMI####
ADDMI_ORF = list(
top_Full_ADDMI = Full_ADDMI_ORF_Eff %>%
  arrange(., desc(Effect_sq)) %>%
  filter(., Effect_sq > quantile(Effect_sq, 0.99)) %>%
  select(., ORF_ID) %>%
  unlist()
,
top_Conc_ADDMI = Conc_ADDMI_ORF_Eff %>%
  arrange(., desc(Effect_sq)) %>%
  filter(., Effect_sq > quantile(Effect_sq, 0.99)) %>%
  select(., ORF_ID) %>%
  unlist()
,
top_Forage_ADDMI = Forage_ADDMI_ORF_Eff %>%
  arrange(., desc(Effect_sq)) %>%
  filter(., Effect_sq > quantile(Effect_sq, 0.99)) %>%
  select(., ORF_ID) %>%
  unlist()
)

# All three
ADDMI_all <- Reduce(intersect, ADDMI_ORF)

# Pairs
Full_Conc_ADDMI <- setdiff(intersect(ADDMI_ORF$top_Full_ADDMI, ADDMI_ORF$top_Conc_ADDMI), ADDMI_all)
Full_Forage_ADDMI <- setdiff(intersect(ADDMI_ORF$top_Full_ADDMI, ADDMI_ORF$top_Forage_ADDMI), ADDMI_all)
Conc_Forage_ADDMI <- setdiff(intersect(ADDMI_ORF$top_Conc_ADDMI, ADDMI_ORF$top_Forage_ADDMI), ADDMI_all)

ggvenn(ADDMI_ORF)

##ADG####
ADG_ORF = list(
  top_Full_ADG = Full_ADG_ORF_Eff %>%
    arrange(., desc(Effect_sq)) %>%
    filter(., Effect_sq > quantile(Effect_sq, 0.99)) %>%
    select(., ORF_ID) %>%
    unlist()
  ,
  top_Conc_ADG = Conc_ADG_ORF_Eff %>%
    arrange(., desc(Effect_sq)) %>%
    filter(., Effect_sq > quantile(Effect_sq, 0.99)) %>%
    select(., ORF_ID) %>%
    unlist()
  ,
  top_Forage_ADG = Forage_ADG_ORF_Eff %>%
    arrange(., desc(Effect_sq)) %>%
    filter(., Effect_sq > quantile(Effect_sq, 0.99)) %>%
    select(., ORF_ID) %>%
    unlist()
)

# All three
ADG_all <- Reduce(intersect, ADG_ORF)

# Pairs
Full_Conc_ADG <- setdiff(intersect(ADG_ORF$top_Full_ADG, ADG_ORF$top_Conc_ADG), ADG_all)
Full_Forage_ADG <- setdiff(intersect(ADG_ORF$top_Full_ADG, ADG_ORF$top_Forage_ADG), ADG_all)
Conc_Forage_ADG <- setdiff(intersect(ADG_ORF$top_Conc_ADG, ADG_ORF$top_Forage_ADG), ADG_all)

ggvenn(ADG_ORF)

##FtG####
FtG_ORF = list(
  top_Full_FtG = Full_FtG_ORF_Eff %>%
    arrange(., desc(Effect_sq)) %>%
    filter(., Effect_sq > quantile(Effect_sq, 0.99)) %>%
    select(., ORF_ID) %>%
    unlist()
  ,
  top_Conc_FtG = Conc_FtG_ORF_Eff %>%
    arrange(., desc(Effect_sq)) %>%
    filter(., Effect_sq > quantile(Effect_sq, 0.99)) %>%
    select(., ORF_ID) %>%
    unlist()
  ,
  top_Forage_FtG = Forage_FtG_ORF_Eff %>%
    arrange(., desc(Effect_sq)) %>%
    filter(., Effect_sq > quantile(Effect_sq, 0.99)) %>%
    select(., ORF_ID) %>%
    unlist()
)

# All three
FtG_all <- Reduce(intersect, FtG_ORF)

# Pairs
Full_Conc_FtG <- setdiff(intersect(FtG_ORF$top_Full_FtG, FtG_ORF$top_Conc_FtG), FtG_all)
Full_Forage_FtG <- setdiff(intersect(FtG_ORF$top_Full_FtG, FtG_ORF$top_Forage_FtG), FtG_all)
Conc_Forage_FtG <- setdiff(intersect(FtG_ORF$top_Conc_FtG, FtG_ORF$top_Forage_FtG), FtG_all)

ggvenn(FtG_ORF)

#Conversely Related####
##Full####
#Select ORF which have extreme effect size for that trait
Full_ADDMI_ORF_Eff$Effect_Z = scale(Full_ADDMI_ORF_Eff$Effect) %>%
  as.vector()
converse_Full_ADDMI = Full_ADDMI_ORF_Eff %>%
  filter(., abs(Effect_Z) > 3) %>%
  select(., c(ORF_ID, Effect, Effect_Z))
colnames(converse_Full_ADDMI) = c("ORF_ID", "ADDMI_Effect", "ADDMI_Effect_Z")

Full_ADG_ORF_Eff$Effect_Z = scale(Full_ADG_ORF_Eff$Effect) %>%
  as.vector()
converse_Full_ADG = Full_ADG_ORF_Eff %>%
  filter(., abs(Effect_Z) > 3) %>%
  select(., c(ORF_ID, Effect, Effect_Z))
colnames(converse_Full_ADG) = c("ORF_ID", "ADG_Effect", "ADG_Effect_Z")

#Merge the extreme ORF together
converse_Full = merge(converse_Full_ADDMI, converse_Full_ADG)

#Find the ORF which are extreme in different directions for ADDMI and ADG
converse_Full = converse_Full[converse_Full$ADDMI_Effect_Z * converse_Full$ADG_Effect_Z < 0, ]
rm(converse_Full_ADDMI, converse_Full_ADG)

##Conc####
#Select ORF which have extreme effect size for that trait
Conc_ADDMI_ORF_Eff$Effect_Z = scale(Conc_ADDMI_ORF_Eff$Effect) %>%
  as.vector()
converse_Conc_ADDMI = Conc_ADDMI_ORF_Eff %>%
  filter(., abs(Effect_Z) > 3) %>%
  select(., c(ORF_ID, Effect, Effect_Z))
colnames(converse_Conc_ADDMI) = c("ORF_ID", "ADDMI_Effect", "ADDMI_Effect_Z")

Conc_ADG_ORF_Eff$Effect_Z = scale(Conc_ADG_ORF_Eff$Effect) %>%
  as.vector()
converse_Conc_ADG = Conc_ADG_ORF_Eff %>%
  filter(., abs(Effect_Z) > 3) %>%
  select(., c(ORF_ID, Effect, Effect_Z))
colnames(converse_Conc_ADG) = c("ORF_ID", "ADG_Effect", "ADG_Effect_Z")

#Merge the extreme ORF together
converse_Conc = merge(converse_Conc_ADDMI, converse_Conc_ADG)

#Find the ORF which are extreme in different directions for ADDMI and ADG
converse_Conc = converse_Conc[converse_Conc$ADDMI_Effect_Z * converse_Conc$ADG_Effect_Z < 0, ]
rm(converse_Conc_ADDMI, converse_Conc_ADG)

##Forage####
#Select ORF which have extreme effect size for that trait
Forage_ADDMI_ORF_Eff$Effect_Z = scale(Forage_ADDMI_ORF_Eff$Effect) %>%
  as.vector()
converse_Forage_ADDMI = Forage_ADDMI_ORF_Eff %>%
  filter(., abs(Effect_Z) > 3) %>%
  select(., c(ORF_ID, Effect, Effect_Z))
colnames(converse_Forage_ADDMI) = c("ORF_ID", "ADDMI_Effect", "ADDMI_Effect_Z")

Forage_ADG_ORF_Eff$Effect_Z = scale(Forage_ADG_ORF_Eff$Effect) %>%
  as.vector()
converse_Forage_ADG = Forage_ADG_ORF_Eff %>%
  filter(., abs(Effect_Z) > 3) %>%
  select(., c(ORF_ID, Effect, Effect_Z))
colnames(converse_Forage_ADG) = c("ORF_ID", "ADG_Effect", "ADG_Effect_Z")

#Merge the extreme ORF together
converse_Forage = merge(converse_Forage_ADDMI, converse_Forage_ADG)

#Find the ORF which are extreme in different directions for ADDMI and ADG
converse_Forage = converse_Forage[converse_Forage$ADDMI_Effect_Z * converse_Forage$ADG_Effect_Z < 0, ]
rm(converse_Forage_ADDMI, converse_Forage_ADG)
