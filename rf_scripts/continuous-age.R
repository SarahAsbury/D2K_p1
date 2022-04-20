# Clear Environment -------------------------------------------------------
#Clear environment
rm(list = ls(all.names = TRUE))

# Load packages -------------------------------------------------------
myPackages <- c("GGally", "e1071", "cowplot","randomForest","tidyverse","nnet","ROCR",
                "Hmisc", "NCmisc", "phyloseq", "forecast", "janitor", "ggpubr", "cowplot", "ViewPipeSteps")
tryCount <- 0    

while( !all(myPackages %in% (.packages())) ){
  
  try(require(GGally))
  try(require(e1071))
  try(require(cowplot))
  try(require(randomForest))
  try(require(tidyverse))
  try(require(nnet))
  try(require(ROCR))
  try(require(Hmisc))
  try(require(NCmisc))
  try(require(phyloseq))
  try(require(forecast))
  try(require(janitor))
  try(require(ggpubr))
  try(require(cowplot))
  library(ViewPipeSteps) 
  
  tryCount <- tryCount + 1
  
  if( !all(myPackages %in% (.packages()))  ){
    cat(paste0("Failure: ", tryCount, "\n"))
    cat("Failed to load: ")
    cat(myPackages[ !myPackages %in% (.packages()) ])
    cat("\n")
  } else {
    print(paste0("Packages loaded successfully!"))
  }
  
  Sys.sleep(5)
  
}





# Load rf functions -------------------------------------------------------
source("/Users/sar/Documents/GitHub/asbury-datatools/RandomForest_Functions.R")



# Resolve to genus  ----------------------------------------------------------
load("/Users/sar/Dropbox/Dallas2K/PhyloseqV3/Paper_v2/Physeq_glommed_objects/physeq_ResolveGenus.RData")
dat_rel_glom #glommed ASVs 




# Session inputs ----------------------------------------------------------
vpred <- "age"
rf.type <- "reg"



# Set wd ------------------------------------------------------------------
setwd("/Users/sar/Dropbox/Dallas2K/PhyloseqV3/Paper_v2/randomforest_out")

# Define rf df ------------------------------------------------------------
dat_clr <- dat_glom_aldex %>% microbiome::transform("clr")
x <- psmelt(dat_clr) %>%
  group_by(StudyID) %>%
  dplyr::select(c(OTU,Sample, Abundance, StudyID, vpred)) %>%
  arrange(Sample, desc(Abundance))%>%
  group_by(Sample) %>%
  tidyr::pivot_wider(values_from = Abundance, names_from = OTU) %>%
  column_to_rownames("Sample") %>% 
  dplyr::select(-StudyID) %>% select(-other) %>% clean_names() %>%
  filter(!is.na(get(vpred)))


# Create 10 test sets -----------------------------------------------------
split.ratio(df = x, rf.type = rf.type)



# Mtry guide --------------------------------------------------------------
mtry.guide.save <- mtry.guide(npred = x %>% select(-vpred) %>% colnames %>% length, 
                              rf.type = rf.type)
mtry.guide.save$range




# Run rf ------------------------------------------------------------------
sa_rfreg(vpred = vpred, 
         custom.name = TRUE,
         date = "8feb21",
         dataframe.name = "D2K",
         predictors.name = "microbiome",
         mtry = mtry.guide.save$range)
  

