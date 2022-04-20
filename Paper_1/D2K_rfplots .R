# Load Packages -----------------------------------------------------------
#Load packages
source("/Users/sar/Dropbox/Dallas2K/PhyloseqV2/Load Pipeline Packages.R")
library("tidyverse")
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
library(GGally)
library(pime)
library(parallel)
library(foreach)
library(doParallel)



# Clear Environment -------------------------------------------------------
#Clear environment
rm(list = ls(all.names = TRUE))





# Set wd -------------------------------------------------------------------
setwd("/Users/sar/Dropbox/Dallas2K/PhyloseqV3/Paper_v2/randomforest_out")



# Import aggregated variable importance dfs ---------------------------------------

age.varimp <- read.csv("/Users/sar/Dropbox/Dallas2K/PhyloseqV3/Paper_v2/randomforest_out/D2K_8feb21_microbiome_age/imp_aggregate.csv") %>% select(-X)
gad.varimp <- read.csv("/Users/sar/Dropbox/Dallas2K/PhyloseqV3/Paper_v2/randomforest_out/D2K_8feb21_microbiome_gad7_total/imp_aggregate.csv") %>% select(-X)
phq.varimp <- read.csv("/Users/sar/Dropbox/Dallas2K/PhyloseqV3/Paper_v2/randomforest_out/D2K_8feb21_microbiome_phq9_total/imp_aggregate.csv") %>% select(-X)

head(age.varimp)
# Import phyloseq ---------------------------------------------------------
load("/Users/sar/Dropbox/Dallas2K/PhyloseqV3/Paper_v2/Physeq_glommed_objects/physeq_ResolveGenus.RData")
dat_glom #Filtered phyloseq object - Counts
dat_rel_glom #Filtered phyloseq object - Rel Abundance
dat_glom_aldex #Contains OTHER taxa so that seq depth is preserved for CLR calculation 

# Function: sp_genus -------------------------------------------------------
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



# Function: Var importance lollipops  --------------------------------------------------
varimport_plot <- function(varimp, #dataframe of variable importance
                           metric = "gini", #One of: mse or gini
                           selection_type = "random_top", #One of: random_top, top
                           top = 20 #number of "top" taxa to plot
                           )
  #Lollipop plot of variable importance 
  {
  
  #X and y variables
  if(metric == "mse"){
    y <- "X.IncMSE_name"
    ylab <- "Increased Mean Square Error"
  }
  
  if(metric == "gini"){
    y <- "IncNodePurity_name"
    ylab <- "Increased node purity"
  }
  
  x <- "taxa"
  
  #Variable importance dataframe
  if(selection_type == "random_top")
    #Take top 20 (or specified) taxa, bottom 5 (specified*0.25) taxa, and 15 (specified*0.75) additional random samples for visualization
    {
    toptax <- varimp %>% slice_max(order_by = get(y), n = 20)
    mintax <- varimp %>% slice_min(order_by = get(y), n = round(top*0.25, digits = 0))
    randomtax <- varimp %>% filter(!(taxa %in% toptax$taxa | taxa %in% mintax$taxa)) %>% 
      slice_sample(n = round(top*0.75, digits = 0))
  
    varimp <- varimp %>% filter(taxa %in% toptax$taxa | taxa %in% mintax$taxa | taxa %in% randomtax$taxa)
  }
  
  if(selection_type == "top")
    #Take top 20 taxa (or specified)
  {
    varimp <- varimp %>% slice_max(order_by = y, n = top)
  }
  
  varimp <- varimp %>% arrange(desc(y)) #re-order variables
  
  #Plot
  ggdotchart(varimp, x = paste(x), y = paste(y),
             sorting = "ascending",                        
             ggtheme = theme_pubr(), 
             add = "segment",
             ylab = ylab, 
             xlab = "") + 
    coord_flip() + 
    theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5))
}
  




# Gini Plots --------------------------------------------------------------------
width = 2300
height = 2300
res = 300

tiff("age_varimp.tiff", width = width, height = height, res = res)
varimport_plot(age.varimp %>% mutate(taxa = sp_genus(physeq = dat_glom, taxa = predictors)), 
               "gini")
dev.off()

tiff("gad_varimp.tiff", width = width, height = height, res = res)
varimport_plot(gad.varimp %>% mutate(taxa = sp_genus(physeq = dat_glom, taxa = predictors)), 
               "gini")
dev.off()

tiff("phq_varimp.tiff", width = width, height = height, res = res)
varimport_plot(phq.varimp %>% mutate(taxa = sp_genus(physeq = dat_glom, taxa = predictors)), 
               "gini")
dev.off()
# MSE plots ---------------------------------------------------------------
tiff("age_varimp_mse.tiff", width = width, height = height, res = res)
varimport_plot(age.varimp %>% mutate(taxa = sp_genus(physeq = dat_glom, taxa = predictors)), 
               "mse")
dev.off()

tiff("gad_varimp_mse.tiff", width = width, height = height, res = res)
varimport_plot(gad.varimp %>% mutate(taxa = sp_genus(physeq = dat_glom, taxa = predictors)), 
               "mse")
dev.off()

tiff("phq_varimp_mse.tiff", width = width, height = height, res = res)
varimport_plot(phq.varimp %>% mutate(taxa = sp_genus(physeq = dat_glom, taxa = predictors)), 
               "mse")
dev.off()

