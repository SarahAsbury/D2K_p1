#D2K 
#Taxa summary table
#Identify taxa that have been pulled out of multiple analysis tools. I
#Then identify overlap in sig taxa across multiple tools. 
#Aldex, RF, β-diversity
#MDD, MDD_GAD


#Updated for PhyseqV3-paper 6jan21

#Updated for Paper_v2 9feb21

# Load Packages -----------------------------------------------------------
#Load packages
source("/Users/sar/Dropbox/Dallas2K/PhyloseqV2/Load Pipeline Packages.R")
library(ALDEx2)
library("DESeq2")
library(FactoMineR)
library(factoextra)
library(viridis)
library(RColorBrewer)
library(mclust)
library(NbClust)
library(ggfortify)
library(rgl)
library(car)
library(factoextra)
library(pca3d)
library("tidyverse")
library("Rmisc")
library("reshape")
library("RColorBrewer")
library("plotly")
library("foreign")
library(DataExplorer)
library(MASS)
library(scales)
library(cowplot)
library(gridGraphics)
library(metagMisc)
library(janitor)
library(FDRestimation)
library(scales)
library(FSA)
library(pgirmess)
library(data.table)
library("GGally")
library(gridExtra)
library(janitor)
library(gtools)

# Clear Environment -------------------------------------------------------
#Clear environment
rm(list = ls(all.names = TRUE))


# Load glom to genus ------------------------------------------------------
load("/Users/sar/Dropbox/Dallas2K/PhyloseqV3/Paper_v2/Physeq_glommed_objects/physeq_ResolveGenus.RData")
dat_glom #Filtered phyloseq object - Counts
dat_rel_glom #Filtered phyloseq object - Rel Abundance
dat_glom_aldex #Contains OTHER taxa so that seq depth is preserved for CLR calculation 


# Functions - Microbiome data toolbox ---------------------------------------------------------------
sp_genus <- function(physeq, #object with reference taxa table
                     taxa #To convert. Must contain species code (sp) assigned by phyloseq (character list)
) 
  #Use species code (e.g Sp1) to extract clean genus-level names from phyloseq object 
  #Extracts family_NA when genus is not available
{
  taxdf <- physeq %>% tax_table %>% data.frame() %>% rownames_to_column(var = "OTU") %>% 
    mutate(sp = gsub("_.*", "", OTU)) %>% 
    mutate(genus_level = ifelse(is.na(Genus),
                                paste(sp, Family, Genus, sep = "_"), 
                                paste(sp, Genus, sep = "_")
    )
    ) %>% select(sp, genus_level)
  
  taxa.in <- taxa %>% data.frame %>% mutate(sp = gsub("_.*", "", taxa))
  taxa.out <- left_join(taxa.in, taxdf, by = "sp")
  return(taxa.out$genus_level)
}

aldex.contglm.extract <- function(aldex_result, #aldex glm output
                                  voi #variable of interest 
)
  #Extracts BH p-value, glm coefficient (slope), and glm x-intercept from aldex results for each taxon
  #Renames columns to human read-able columns
  #Taxa names must be rownames
{
  cont.glm <- aldex_result %>% 
    rename(p = paste0("model.", voi, ".Pr...t.."),
           p.BH = paste0("model.", voi, ".Pr...t...BH"),
           X.Intercept = X.Intercept..Estimate,
           coef = paste0("model.", voi, ".Estimate"))%>% 
    select(p, p.BH, X.Intercept, coef) %>% 
    rownames_to_column(var = "OTU") 
  return(cont.glm)
}





# Import WCNA ------------------------------------------------------------
#Import WCNA module membership
#Brown module
module.edit <- function(df, color){
  df <- df %>% rename(taxa = X) %>% filter(moduleColors == color) %>% rename(moduleColor = moduleColors) %>%
    mutate(across(starts_with("MM"), ~round(.x,digits = 2))) %>% 
    mutate(across(starts_with("p.MM"), ~scientific(.x,digits = 3))) %>%
    mutate(taxa = sp_genus(taxa = taxa, physeq = dat_glom))
  return(df)
}

brown.member <- read.csv("/Users/sar/Dropbox/Dallas2K/PhyloseqV3/Paper_v2/wcna/signed/minmod50/brown/browntaxa_membership.csv") %>% module.edit("brown") 
head(brown.member)


#Import WCNA correlations - all brown module taxa
wcna.alledit <- function(df){
  p.name <- df %>% select(starts_with("p.GS")) %>% colnames
  p.sub <- p.name %>% sub("GS.", "", .)
  sig.name <- p.name %>% sub("p.GS.", "sig.", .)
  
  df.return <- df %>% rename(taxa = X) %>% select(-moduleColors) %>%
    mutate(across(starts_with("GS"), ~round(.x,digits = 2))) %>% 
    rename_at(vars(starts_with("GS")), funs(sub("GS.", "", .))) %>%
    mutate(p = scientific(eval(parse(text = p.name)), digits =3)) %>%
    relocate(p, .before = paste0(p.name)) %>%
    rename_with(.fn = ~paste0(p.sub), .cols = p) %>%
    mutate(across(starts_with("p.GS"), ~p.adjust(.x, method = "BH"), digits = 3)) %>%
    mutate(sig = stars.pval(eval(parse(text = p.name)))) %>% 
    rename_with(.fn = ~paste0(sig.name), .cols = sig) %>%
    mutate(across(starts_with("p.GS"), ~ifelse(.x < 0.001, "<0.001", round(.x, digits = 3)))) %>%
    rename_at(vars(starts_with("p.GS")), funs(sub("p.GS.", "p.BH.", .))) %>%
    mutate(taxa = sp_genus(taxa = taxa, physeq = dat_glom)) %>%
    
    return(df.return)
}

wcna.all.gad <- read.csv("/Users/sar/Dropbox/Dallas2K/PhyloseqV3/Paper_v2/wcna/signed/minmod50/brown/gad7_total/brown_gad7_total_ALLclinicaltaxa.csv") %>% wcna.alledit
wcna.all.phq <- read.csv("/Users/sar/Dropbox/Dallas2K/PhyloseqV3/Paper_v2/wcna/signed/minmod50/brown/phq9_total/brown_phq9_total_ALLclinicaltaxa.csv") %>% wcna.alledit
wcna.all.age <- read.csv("/Users/sar/Dropbox/Dallas2K/PhyloseqV3/Paper_v2/wcna/signed/minmod50/brown/age/brown_age_ALLclinicaltaxa.csv") %>% wcna.alledit

head(wcna.all.gad)


# WCNA - Brown Module Membership and Taxa Significance --------------------
#Supplementary table 
#Module Membership and Taxa Significance table
wcna.mmsigtab <- join_all(list(brown.member, wcna.all.gad, wcna.all.phq, wcna.all.age), by = "taxa", type = "left")
head(wcna.mmsigtab)




# Export Module Membership Taxa Significance Table ------------------------
setwd("/Users/sar/Dropbox/Dallas2K/PhyloseqV3/Paper_v2/wcna/signed/minmod50")
write.csv(wcna.mmsigtab, "brownMM_taxsig_table.csv", row.names = FALSE)



# Import WCNA hub taxa ----------------------------------------------
#Brown module
wcna.edit <- function(df){
  return(df %>% data.frame() %>% rename(taxa = X) %>% mutate(taxa = sp_genus(taxa = taxa, physeq = dat_glom)))
}
wcna.hub <- read.csv("/Users/sar/Dropbox/Dallas2K/PhyloseqV3/Paper_v2/wcna/signed/minmod50/brown/brown_hubtaxa.csv") %>% wcna.edit %>% select(taxa)


# Import RF results --------------------------------------------------------------
rf.taxa20 <- function(df){
  return(df %>% slice_max(n = 20, order_by = IncNodePurity_name) %>% select(-c(X.IncMSE_name, X)) %>% mutate(taxa = sp_genus(taxa = predictors, physeq = dat_glom))
         )
}

rf.rank <- function(df){
  return(df %>% data.frame %>% arrange(desc(IncNodePurity_name)) %>% select(-c(X.IncMSE_name, X)) %>% mutate(taxa = sp_genus(taxa = predictors, physeq = dat_glom)) %>% 
    rownames_to_column(var = "rank")
  )
}

rf.gad.dir <- "/Users/sar/Dropbox/Dallas2K/PhyloseqV3/Paper_v2/randomforest_out/D2K_8feb21_microbiome_gad7_total/imp_aggregate.csv"
rf.phq.dir <- "/Users/sar/Dropbox/Dallas2K/PhyloseqV3/Paper_v2/randomforest_out/D2K_8feb21_microbiome_phq9_total/imp_aggregate.csv"
rf.age.dir <- "/Users/sar/Dropbox/Dallas2K/PhyloseqV3/Paper_v2/randomforest_out/D2K_8feb21_microbiome_age/imp_aggregate.csv"

rf.gad <- read.csv(rf.gad.dir) %>% rf.taxa20
rf.phq <- read.csv(rf.phq.dir) %>% rf.taxa20
rf.age <- read.csv(rf.age.dir) %>% rf.taxa20


rf.gad.rank <- read.csv(rf.gad.dir) %>% rf.rank
rf.phq.rank <- read.csv(rf.phq.dir) %>% rf.rank
rf.age.rank <- read.csv(rf.age.dir) %>% rf.rank


head(rf.gad)
head(rf.age.rank)

# Import Aldex  ------------------------------------------------------------
aldex.p5 <- function(df){
  return(df %>% rename(taxa = OTU) %>% filter(p <= 0.05) %>% mutate(taxa = sp_genus(taxa = taxa, physeq = dat_glom)))
}

aldex.gad <- read.csv("/Users/sar/Dropbox/Dallas2K/PhyloseqV3/Paper_v2/aldex/gad_aldex.csv") %>% data.frame %>% column_to_rownames(var = "X") %>% 
  aldex.contglm.extract(voi = "gad7_total") %>% aldex.p5 
aldex.phq <- read.csv("/Users/sar/Dropbox/Dallas2K/PhyloseqV3/Paper_v2/aldex/phq_aldex.csv") %>% data.frame %>% column_to_rownames(var = "X") %>% 
  aldex.contglm.extract(voi = "phq9_total")%>% aldex.p5
aldex.age <- read.csv("/Users/sar/Dropbox/Dallas2K/PhyloseqV3/Paper_v2/aldex/age_aldex.csv") %>% data.frame %>% column_to_rownames(var = "X") %>%
  aldex.contglm.extract(voi = "age") %>% aldex.p5




# Taxa Consensus Functions ----------------------------------------------------------
convert_numeric <- function(x){
  ifelse(x == 1, "YES", 
         ifelse(x == 0, "", NA))
}

aldex.format <- function(df){
  df <- df %>% select(-X.Intercept) %>% mutate(p = ifelse(p < 0.001, "<0.001", as.character(round(p, digits = 3)))) %>% 
    mutate(p.BH = round(p.BH, digits = 3)) %>% mutate(coef = round(coef, digits = 2)) %>% 
    relocate(coef, .before = p)
  
  return(df)
}

wcna.filter <- function(wcna = wcna.mmsigtab, module.members, module, voi.in){
  print(voi.in)
  #Naming and filtering variables
  MMmodule <- paste0("MM", module)
  quant.50 <- quantile(module.members %>% pull(eval(parse(text = MMmodule)))) %>% data.frame() %>% 
    rename(val = ".") %>% rownames_to_column(var = "quartile") %>% filter(quartile == "50%") %>% pull(val)
  sig.name <- paste0("sig.", voi.in)
  sig.symbols <- c("***", "**", "*")

  #Edit df
  wcna <- wcna %>%
    rename_at((vars(starts_with(paste(voi.in)))), funs(paste0("r"))) %>% 
    mutate(r = abs(r)) %>%
    rename_at(MMmodule, funs(paste0("MM"))) %>%
    rename_at(vars(contains(sig.name)), funs(paste0("sig"))) %>% 
    filter(MM >= quant.50 & r >= 0.2 & sig %in% sig.symbols) 
  
  return(wcna)
                  
}


consensus.taxa <- function(voi #voi must the way it appears in dataframe names 
                           ) 
  #Extracts consensus taxa from WCNA, RF, and Aldex results
  {
  aldex.df <- get(paste0("aldex.", voi))
  aldex.taxa <- aldex.df %>% select(taxa) %>% pull(taxa)
  rf.taxa <- get(paste0("rf.", voi)) %>% select(taxa) %>% pull(taxa)
  rf.rank <- get(paste0("rf.", voi, ".rank")) %>% select(-IncNodePurity_name)
  wcna.taxa <- wcna.mmsigtab %>% wcna.filter(wcna = ., module.members = brown.member, module = "brown", voi.in = paste(voi)) %>% select(taxa) %>% pull(taxa)
  

  summary <- union(wcna.taxa, rf.taxa) %>% union(., aldex.taxa) %>% 
    data.frame %>% rename(taxa = ".") %>%
    mutate(WCNA = ifelse(taxa %in% wcna.taxa, 1, 0)) %>%
    mutate(RandomForest = ifelse(taxa %in% rf.taxa, 1, 0)) %>% 
    mutate(Aldex = ifelse(taxa %in% aldex.taxa, 1, 0)) %>% 
    rowwise %>% mutate(criteria = sum(WCNA, RandomForest, Aldex)) %>% filter(criteria >= 2) %>% select(-criteria) %>%
    mutate(across(c("WCNA", "RandomForest", "Aldex"),convert_numeric )) %>% 
    left_join(aldex.df, by = "taxa") %>% 
    aldex.format %>%
    left_join(rf.rank, by = "taxa") %>% mutate(RandomForest = ifelse(RandomForest == "YES", RandomForest, rank)) %>% select(-rank) #add rank to taxa not in Top20
  return(summary)
}


# Export consensus results ----------------------------------------------------------
setwd("/Users/sar/Dropbox/Dallas2K/PhyloseqV3/Paper_v2/consensus")
write.csv(consensus.taxa(voi = "gad"), "consensus_gad.csv", row.names = FALSE)
write.csv(consensus.taxa(voi = "phq"), "consensus_phq.csv", row.names = FALSE)
write.csv(consensus.taxa(voi = "age"), "consensus_age.csv", row.names = FALSE)






# Import Aldex - Full Table -----------------------------------------------
aldex.full <- function(df){
  return(df %>% data.frame %>% rename(taxa = X) %>% filter(taxa != "other") %>% 
           mutate(taxa = sp_genus(taxa = taxa, physeq = dat_glom)) %>% 
           column_to_rownames(var = "taxa")
  )
}

aldex.gad <- read.csv("/Users/sar/Dropbox/Dallas2K/PhyloseqV3/Paper/aldex/gad_aldex.csv") %>% aldex.full %>% aldex.contglm.extract(voi = "gad7_total") %>% aldex.format() %>% 
  rename(β = coef)
head(aldex.gad)
aldex.phq <- read.csv("/Users/sar/Dropbox/Dallas2K/PhyloseqV3/Paper/aldex/phq_aldex.csv") %>% aldex.full %>% aldex.contglm.extract(voi = "phq9_total") %>% aldex.format() %>% 
  rename(β = coef)
head(aldex.phq)
aldex.age <- read.csv("/Users/sar/Dropbox/Dallas2K/PhyloseqV3/Paper/aldex/age_aldex.csv") %>% aldex.full %>% aldex.contglm.extract(voi = "age") %>% aldex.format() %>% 
  rename(β = coef)
head(aldex.age)


# Export aldex tables -----------------------------------------------------
setwd("/Users/sar/Dropbox/Dallas2K/PhyloseqV3/Paper_v2/aldextables")
write.csv(aldex.gad, "aldexgadtab.csv", row.names = FALSE)
write.csv(aldex.phq, "aldexphqtab.csv", row.names = FALSE)
write.csv(aldex.age, "aldexagetab.csv", row.names = FALSE)


