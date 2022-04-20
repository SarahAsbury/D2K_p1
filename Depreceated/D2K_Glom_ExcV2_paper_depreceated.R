#Glom to genus
#Exectuable
#Algorithim Updated: 13oct21
#V2 - changed pre-glom wrangling to include clin phen (MDD/GAD/MDD_GAD/NONE)
#V2 - This code now includes OTHER taxa update.

#To do; Next version updates: actually convert names of taxa between relative and count phyloseq objects


#Updated: 9nov21
#Changed outputs for PhyseqV3 
#Removed z-score scaling for clinical variables 

#Updated 3jan22
#No longer remove samples that have control depression/anxiety but have anhedonia


#Updated 7feb22
#Export phyloseq object with all visit codes (i.e before visit = 1 subset)

# Clear Environment -------------------------------------------------------
#Clear environment
rm(list = ls(all.names = TRUE))





# Set date ----------------------------------------------------------------
date <- readline("Input date: ")

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

# Set wd (physeq input) ------------------------------------------------------------------

#Set wd
setwd("/Users/sar/Dropbox/Dallas2K/PhyloseqV3")

# Ensure each of these correspond to the file names of your asv, taxa and map files from your 'data' folder. Run these lines to load your data into R: 
asvfile = paste0(getwd(),'/PhlyoseqInput/seqtab_nochim_transposed_JFUTSouthwestern_v34_SA21feb21.csv')
taxfile = paste0(getwd(), '/PhlyoseqInput/taxa_JFUTSouthwestern_v34_silva132.csv')
mapfile = paste0(getwd(), '/PhlyoseqInput/mapdf_v3.csv')


# Set wd (output) ----------------------------------------------------------
setwd("/Users/sar/Dropbox/Dallas2K/PhyloseqV3/paper")


# Import and inspect ASV --------------------------------------------------
# Convert ASV table to a data frame and designate the first row as the row names:
asvdf = read.csv(asvfile, row.names=1)
asvdf <- asvdf
str(asvdf)

# This will give us the dimensions (or columns and rows) of the data frame. To view the dimensions of the ASV table, run this line: 
dim(asvdf)

# Run this line to view the first 10 columns and rows within the data in the ASV data frame:
asvdf[1:10,1:10]

# The row names need to be removed to properly merge the tables. Set the value of the ASV sequence names (that are the row names) to 'NULL': 
seqs = rownames(asvdf)
rownames(asvdf) = NULL



# Import and Inspect Taxonomy Table ---------------------------------------
# Do the same thing with the taxonomy table. Convert the taxonomy table to a data frame and designate the first row as the row names: 
taxdf = read.csv(taxfile, row.names = 1)
head(taxdf)

# Create a data frame from the taxa table and remove the row names to match the ASV data frame:
all(rownames(taxdf) == seqs)
rownames(taxdf) = NULL

# The taxonomy table and ASV table should have the same number of rows. Run this line to verify this:
dim(taxdf)



# Import and Inspect Map File ---------------------------------------------
# Convert the map file to a data frame and designating the first row as the row names: 
mapdf = read.csv(mapfile, 
                 strip.white = TRUE,
                 na.strings = 'NA',
                 stringsAsFactors = FALSE)

# To verify the dimensions of the map file to ensure all samples have imported properly, run these lines: 
dim(mapdf)
head(mapdf)


#Check that D2K IDs removed for physeqV3 are absent in dataset
qc <- mapdf %>% filter(record == "D2K000196" | record == "D2K000151") %>% nrow()
ifelse(qc == 0, "QC passed - correct df for physeqV3 imported", "QC failed")

# Convert clinical data to categorical variables --------------------------
mapdf <- mapdf %>% mutate(depression = ifelse(phq9_total >= 10, "depressed", ifelse(phq9_total <=9, "control", NA)),
                          anxiety = ifelse(gad7_total >=10, "anxious", ifelse(gad7_total <=9, "control", NA)),
                          anhedonia = ifelse(dars17_total <= 44, "anhedonia", ifelse(dars17_total > 44, "control", NA)))


# Verify File Formatting and Alignment ------------------------------------
# Verify the structure of the files. There should not be row names in the taxonomy file or ASV file. The files should also be loaded and converted to data frames. 
str(asvdf)
str(taxdf)
str(mapdf)

#mapdf wrangle
mapdf <- mapdf %>% dplyr::rename(StudyID = JF_ID)


#D2K000310 JF1785	was removed from mapdf because it failed clincial data QC 
#Identify and remove D2K000310 JF1785
d2k.rem <- data.frame(colnames(asvdf)) %>% filter(!colnames.asvdf. %in% mapdf$StudyID)
asvdf <- asvdf %>% dplyr::select(-c((paste(d2k.rem))))


# If the column names (with subject IDs) line up among the tables, this command will return 'TRUE' and character(0), respectively:
all(colnames(asvdf) %in% mapdf$StudyID)
mapdf$StudyID[!(mapdf$StudyID %in% colnames(asvdf))]
data.frame(colnames(asvdf)) %>% filter(!colnames.asvdf. %in% mapdf$StudyID) #Should return nothing (0 x1 df)



# Create Phyloseq Object --------------------------------------------------
# Verify the row names in the 'mapfile' data frame are the same as the column names in the ASV data frame:
rownames(mapdf) = mapdf$StudyID
all(rownames(mapdf) %in% colnames(asvdf))
all(colnames(asvdf) %in% rownames(mapdf))

# Convert the ASV and taxonomy data frames to matrices. 
tax_mat = as.matrix(taxdf)
asv_mat = as.matrix(asvdf)

# To filter by time point for later analysis, export the ASV matrix and group the time points together that you are interested in comparing and load these new tables into R. For example, if you compare day 0 to day seven, you may choose to HERE

# To verify matrix structure, run these lines:  
head(tax_mat)
head(asv_mat)

# To make a phyloseq object, you must specify the orientation of your ASV table:
dat = phyloseq(otu_table(asv_mat, taxa_are_rows = TRUE), tax_table(tax_mat), sample_data(mapdf))
str(dat)


# Export ASV - sp code ----------------------------------------------------
taxdf.asv <- read.csv(taxfile, row.names = 1)

qc1 <- (dat %>% tax_table() %>% data.frame() %>% replace(is.na(.), "UNKNOWN") == taxdf %>% replace(is.na(.), "UNKNOWN")) %>% 
  as.logical %>% as.data.frame() %>% rename(qc = ".") #physeq vs. taxdf 
qc2 <- (taxdf %>% replace(is.na(.), "UNKNOWN") == taxdf.asv %>% replace(is.na(.), "UNKNOWN")) %>% as.data.frame() #taxdf.asv vs. taxdf
qc3 <- qc2 %>% mutate_all(as.integer) #convert qc2 TRUE/FALSE to integer 

qc <- ifelse(unique(qc1) %>% length == 1 & unique(qc1)[1,1] == TRUE &
             unique(qc2) %>% nrow == 1 & rowSums(unique(qc3)) == ncol(qc3),
             "pass", "failed")



print(paste("Phyloseq taxa table (dat) and ASV sequence taxa df (taxdf.asv) aligned with original taxa df (taxdf) :", qc))
if(qc == "pass"){
  print("Create ASV - sp code")
  sp_asv <- (taxdf.asv %>% rownames) %>% data.frame() %>% rename(asv_seq = ".") %>% mutate(sp = dat %>% tax_table %>% rownames())
}

write.csv(sp_asv, "asv_spCode.csv")

# Add clinical phenotypes (MDD, MDD_GAD, NONE) to Phyloseq object ------------------------------------
#Remove anhedonia from analysis, but use anhedonia as exclusion criteria for NONE 
clin.phen <- sample_data(dat) %>% data.frame() %>% 
  mutate(mood = ifelse(depression == "control" & anxiety == "control" & anhedonia == "control", "NONE", 
                       ifelse(depression == "depressed" & anxiety == "control", "MDD",
                              ifelse(depression == "control" & anxiety == "anxious", "GAD", 
                                     ifelse(depression == "depressed" & anxiety == "anxious", "MDD_GAD", NA)
                              )
                       )
  )
  )
nrow(clin.phen) #n = 382

#Subset samples in Phyloseq object
dat.prep2 <- dat %>% subset_samples(StudyID %in% clin.phen$StudyID) #Only include StudyIDs that made it through previous clin phenotype filtering on mapdf_alpha (QC wad done on just df; not phyloseq )
dat.prep2
#QC that StudyIDs are aligned between mapdf_alpha (which contains mood column) and Phyloseq object 
qc.df <- sample_data(dat.prep2) %>% data.frame() %>% 
  mutate(StudyID_clin.phen = clin.phen$StudyID) %>% 
  mutate(QC = ifelse(StudyID_clin.phen == StudyID, "Yes", "No")) %>% count('QC')
qc.df
qc.pass <- ifelse("No" %in% qc.df$QC, FALSE, TRUE)

#Add mood column to phyloseq object
if (qc.pass == TRUE){
  print("QC passed. Mood column will be added to Phyloseq object (dat.prep3)")
  dat.prep3 <- dat.prep2
  sample_data(dat.prep3)$mood <- clin.phen$mood
} 
if (qc.pass == FALSE){
  print("Warning: QC failed. Mood column not added.")
}
dat.prep3
dat.prep3 %>% sample_data() %>% head()
dat <- dat.prep3 #assign newly generated phyloseq object containig mood in sample data to dat 

# Add taxonomic label to Phyloseq ASVs ---------------------------------------------
#Get taxa names from phyloseq object (df)
df <- psmelt(dat)
head(df)
asvtab <- data.frame(df)
colnames(asvtab)

#Taxa names 
taxacol <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
asvtaxa <- asvtab[,c("OTU", taxacol)] %>% unique()
head(asvtaxa)

#Add taxonomy to ASV name - assign asvtaxa2 df
head(asvtaxa)
asvtaxa2 <- asvtaxa %>% #Create new col tname with OTU and full taxonomic details 
  mutate(tname = paste(OTU, Kingdom, Phylum, Class, Order, Family, Genus, sep = "_"),
         OTU_order = parse_number(OTU) %>% as.numeric()) %>%
  arrange(OTU_order) %>% dplyr::select(-OTU_order)
head(asvtaxa2)

#Add taxa names to phyloseq object
#QC will = 0 if there asvtaxa2$OTU and phyloseq object ASV row names are aligned. 
qc.taxa_names <- (asvtaxa2$OTU == tax_table(dat) %>% row.names()) %>% data.frame() %>% filter(FALSE %in% .) %>% nrow()
qc.taxa_names
if(qc.taxa_names == 0){
  print("QC passed, full ASV taxonomic names will be added to phyloseq object")
  taxa_names(dat) <- asvtaxa2$tname
}

# Remove QC'ed Samples from Pipeline --------------------------------------
#Samples: 
#D2K000240 Reads >10,000 & Low alpha diversity
#D2K000106 Reads ~ 17.5k & low alpha diversity 


dat
dat_pr <- subset_samples(dat, !StudyID %in% c("JF1789", "JF1793")) #Remove QC'ed samples (Low counts and/or α diversity)
dat_pr




# Remove Host Sequences ---------------------------------------------------
dat_pr = subset_taxa(dat_pr,
                     !(Kingdom == 'Eukaryota' |
                         is.na(Phylum) |
                         (!is.na(Family) &
                            Family == 'Mitochondria')))
dat_pr





# Subset visit = 1 --------------------------------------------------------
dat_all <- dat_pr
dat_pr <- subset_samples(dat_pr, visit == 1)
dat_pr



# Convert Counts to Relative Abundance ------------------------------------
# Include the conversion to relative abundance as an anonymous function:
dat_rel = transform_sample_counts(dat_pr,function(x) x/sum(x))
dat_rel
# QC Relative Abundance Sums to 1 -----------------------------------------
# Verify that all relative abundances sum to 1
melt_rel_group = psmelt(dat_rel) %>% group_by(Sample) %>%
  dplyr::summarize(total_Abundance = sum(Abundance))
melt_rel_group
summary(melt_rel_group)















# Melt Relative abundance into df --------------------------------------------
# Write table for dat_rel and 'melt' 
melt_rel = psmelt(dat_rel)
melt_rel = melt_rel %>%
  group_by(StudyID) %>%
  arrange(desc(Abundance))

# Relative Abundance ASV Table as df ---------------------------------------------------------
# Melt ASV table into data frame and inspect 
#Select the first 4 rows of the melted phyloseq object containing relative abundance
rel_abund_plotting <- melt_rel %>%
  dplyr::select(c(OTU,Sample, Abundance, StudyID))
head(rel_abund_plotting)

rel_abund_long <- rel_abund_plotting 
# Check number of unique OTUs/ASVs in your data frame. This should contain as many as original data frame after filtering of host sequences
length(unique(rel_abund_plotting$OTU))

# ASV Table: Taxa as Columns (Wide) ----------------------------------------------
rel_abund_plotting = rel_abund_plotting %>%
  arrange(Sample, desc(Abundance))%>%
  group_by(Sample) %>%
  tidyr::pivot_wider(values_from = Abundance, names_from = OTU) %>%
  replace(is.na(.), 0)
head(rel_abund_plotting)


# ==== Resolve to Genus - prepare df ==== -------------------------------------------------
#Set new dir 
dir <- "/Users/sar/Dropbox/Dallas2K/PhyloseqV3/paper"
newdir <- "ResolveGenus"
dir.create(newdir)
setwd(paste0(dir, "/", newdir))

#Prepare df
rel <- rel_abund_plotting %>% ungroup() %>% select(-Sample)


#Extract taxa frequency and mean relative abundance 
relmean <- rel %>% select(-StudyID) %>% summarise_all(mean) %>% mutate(stat = "mean") %>% 
  pivot_longer(!stat, names_to = "OTU", values_to = "mean") %>% select(-stat)

relprop <- rel %>% select(-StudyID) %>% 
  mutate_all(function(x) (ifelse(x>0, 1, 0))) %>% #converts taxa presence as binary (1 = yes, 0 = no)
  summarise_all(mean) %>%
  mutate(stat = "prop") %>% pivot_longer(!stat, names_to = "OTU", values_to = "prop") %>% select(-stat)

rel.criteria <- relmean %>% left_join(relprop, by = "OTU") %>% dplyr::rename(mean_relabund = mean)
summary(rel.criteria)




# Resolve to Genus: Data exploration  -------------------------------------------------------
#Data exploration
#Density
density(log10(rel.criteria$mean_relabund)) %>% plot(main = "Log 10 Relative Abundance") #density
p1 <- recordPlot()
ggdraw(p1)

density(rel.criteria$prop) %>% plot(main = "Proportion")
p2 <- recordPlot()
ggdraw(p2)


tiff(paste0("meanRA_prop_density_", date, ".tiff"), height = 500, width = 1200, res = 90)
cowplot::plot_grid(p1, p2)
dev.off()

#Scatter
p <- ggplot(rel.criteria, aes(x = prop, y = mean_relabund)) + geom_point() + 
  theme_classic() + scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                                  labels = trans_format("log10", math_format(10^.x))) + 
  annotation_logticks(sides = "l")  

tiff(paste0("meanRA_prop_scatter_", date, ".tiff"), height = 500, width = 500, res = 90)
p
dev.off()

# Resolve to Genus: QC ----------------------------------------------------
#QC
#Identified empty taxa, QCing to ensure they are actually empty in original rel df
#These taxa should be removed from analyses to increase statistical power (too many features = more comparisons = higher multiple comparisons adjustment)
emptytaxa <- rel.criteria %>% filter(prop == 0)
emptytaxa.df <- rel %>% select(emptytaxa$OTU) %>% summarise_all(mean) %>% mutate(dummy = "dummy") %>% pivot_longer(!dummy, names_to = "OTU", values_to = "mean_abund") %>% select(-dummy)
summary(emptytaxa.df)


# Empty taxa removed - physeq - for export only --------------------------------
emptytaxa.prune <- rel.criteria %>% filter(!(OTU %in% emptytaxa$OTU))
dat_empty <- prune_taxa(emptytaxa.prune$OTU, dat_pr)
dat_empty_rel <- prune_taxa(emptytaxa.prune$OTU, dat_rel)

# Step 1.0: Check abund/prop and glom  ---------------------------------
#Create pass and fail (glom) df
pass <- rel.criteria %>% filter(mean_relabund > 0.001 & prop > 0.1 | prop > 0.5) %>% select(OTU) #Taxa that do not need to be glommed (meet criteria)
fail <- rel.criteria %>% filter(!OTU %in% pass$OTU & !OTU %in% emptytaxa.df$OTU) %>% select(OTU) #Taxa that require glom (below abund/proportion thresholds) & with empty taxa removed 

#Pass criteria taxa
dat_g1_pass_rel <- prune_taxa(pass$OTU, dat_rel)
dat_g1_pass_counts <- prune_taxa(pass$OTU, dat_pr)

dat_g1_pass_rel 
dat_g1_pass_counts #use for downstream analyses

setdiff(otu_table(dat_g1_pass_rel) %>% data.frame() %>% row.names(),
        otu_table(dat_g1_pass_counts) %>% data.frame() %>% row.names()) #should return character(0) if the OTUs are the same in both Physeq objects




#Fail criteria taxa
dat_g1_fail_rel <- prune_taxa(fail$OTU, dat_rel)
dat_g1_fail_counts <- prune_taxa(fail$OTU, dat_pr)

dat_g1_fail_rel
dat_g1_fail_counts

setdiff(otu_table(dat_g1_fail_rel) %>% data.frame() %>% row.names(),
        otu_table(dat_g1_fail_counts) %>% data.frame() %>% row.names()) #should return character(0) if the OTUs are the same in both Physeq objects

taxrel.preglom <- as(access(dat_g1_fail_rel, "tax_table"), "matrix")[, 6] %>% as.data.frame()
taxcount.preglom <- as(access(dat_g1_fail_counts, "tax_table"), "matrix")[, 6] %>% as.data.frame()
setdiff(taxrel.preglom, taxcount.preglom)

write.csv(taxrel.preglom, "taxrel.preglom.csv")
write.csv(taxcount.preglom, "taxcount.preglom.csv")

#QC relative abundance sum of fail + pass = 1 (before glom)
melt_rel_group = rbind(psmelt(dat_g1_pass_rel), psmelt(dat_g1_fail_rel)) %>% group_by(Sample) %>% #QC relative abundance; should be 1 because haven't discarded genus yet
  dplyr::summarize(total_Abundance = sum(Abundance))
melt_rel_group
summary(melt_rel_group)

#Glom failed taxa
rank_names(dat_g1_fail_rel)
rank_names(dat_g1_fail_counts)


dat_g1_fail_rel <- tax_glom(dat_g1_fail_rel, NArm = TRUE, taxrank = "Genus")
dat_g1_fail_counts <- tax_glom(dat_g1_fail_counts, NArm = TRUE, taxrank = "Genus")


#List of taxa postglom
taxrel <- as(access(dat_g1_fail_rel, "tax_table"), "matrix")[, 6] %>% as.data.frame() 
taxcount <- as(access(dat_g1_fail_counts, "tax_table"), "matrix")[, 6] %>% as.data.frame() 
diff <- setdiff(taxrel, taxcount)

#Conversion between relative taxa sp# and count taxa sp#
tax.conversion <- taxrel %>% rownames_to_column("Rel.OTU") %>% 
  left_join(taxcount %>% rownames_to_column("Count.OTU"), by = ".")



#Export postglom taxa list 
write.csv(taxrel, "taxrel.postglom.csv")
write.csv(taxcount, "taxcount.postglom.csv")
write.csv(tax.conversion, "taxaconversion_postglom.csv", row.names = FALSE)

#Issue 1.0: OTUs don't match between relative and counts phyloseq object 
setdiff(otu_table(dat_g1_fail_rel) %>% data.frame() %>% row.names(),
        otu_table(dat_g1_fail_counts) %>% data.frame() %>% row.names()) #should return character(0) if the OTUs are the same in both Physeq objects
#Solution 1.0: Will convert rel OTUs to count OTUs prior to pruning taxa in count phyloseq object (Step 2.0; QC below)
duplicated(taxrel$.) %>% count(.) #QC - Should return FALSE - If FALSE, it means that each genus only appears once, therefore can convert SPs between CLR and Rel dfs without risk that they are different taxa
duplicated(taxcount$.) %>% count(.) #QC - Should return FALSE - If FALSE, it means that each genus only appears once, therefore can convert SPs between CLR and Rel dfs without risk that they are different taxa



#Psmelt glommed relative Physeq
physeq <- dat_g1_fail_rel
intersect(sample_variables(physeq), rank_names(physeq)) #Should return: character(o)
rel <- psmelt(physeq) %>% dplyr::select(c(OTU,Sample, Abundance, StudyID)) %>%
  arrange(Sample, desc(Abundance))%>%
  group_by(Sample) %>%
  tidyr::pivot_wider(values_from = Abundance, names_from = OTU) %>%
  replace(is.na(.), 0) %>% ungroup() %>% select(-Sample)

#QC relative abundance sum of fail + pass (after glom)
melt_rel_group = rbind(psmelt(dat_g1_pass_rel), psmelt(dat_g1_fail_rel)) %>% group_by(Sample) %>% #QC relative abundance; should be 1 because haven't discarded genus yet - value less than 1 means that some ASVs were not glommed (no genus to glom to)
  dplyr::summarize(total_Abundance = sum(Abundance))
melt_rel_group
summary(melt_rel_group)


# Step 1.1: Check abund/prop of genus post-glom ----------------------------------------------
#Extract taxa frequency and mean relative abundance 
relmean <- rel %>% select(-StudyID) %>% summarise_all(mean) %>% mutate(stat = "mean") %>% 
  pivot_longer(!stat, names_to = "OTU", values_to = "mean") %>% select(-stat)

relprop <- rel %>% select(-StudyID) %>% 
  mutate_all(function(x) (ifelse(x>0, 1, 0))) %>% #converts taxa presence as binary (1 = yes, 0 = no)
  summarise_all(mean) %>%
  mutate(stat = "prop") %>% pivot_longer(!stat, names_to = "OTU", values_to = "prop") %>% select(-stat)

rel.criteria <- relmean %>% left_join(relprop, by = "OTU") %>% dplyr::rename(mean_relabund = mean)
summary(rel.criteria)

pass <- rel.criteria %>% filter(mean_relabund > 0.0001 & prop > 0.1 | prop > 0.5)
fail <- rel.criteria %>% filter(!OTU %in% pass$OTU)


# Step 1.2: Data exploration of glommed df  -------------------------------------------------------
#Data exploration
#Density
density(log10(rel.criteria$mean_relabund)) %>% plot(main = "Log10 Relative Abundance") #density
p1 <- recordPlot()
cowplot::ggdraw(p1)

density(rel.criteria$prop) %>% plot(main = "Proportion")
p2 <- recordPlot()
cowplot::ggdraw(p2)


tiff("meanRA_prop_density_glom1Failedtaxa.tiff", height = 500, width = 1200, res = 90)
cowplot::plot_grid(p1, p2)
dev.off()

#Scatter
p <- ggplot(rel.criteria, aes(x = prop, y = mean_relabund)) + geom_point() + 
  theme_classic() + scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                                  labels = trans_format("log10", math_format(10^.x))) + 
  annotation_logticks(sides = "l")  

tiff("meanRA_prop_scatter_glom1Failedtaxa.tiff", height = 500, width = 500, res = 90)
p
dev.off()


# Step 2.0: Prune taxa (genus level) ----------------------------------------------------------------
#Pass criteria taxa 
pass.count <- tax.conversion %>% filter(Rel.OTU %in% pass$OTU) #Convert relative abundance OTUs to correct Count OTUs (Solution 1.0 to Issue 1.0)

setdiff(pass.count$Rel.OTU, pass$OTU)
setdiff(pass$OTU, pass.count$Count.OTU)
tax.CountReldifference <- setdiff(pass$OTU, pass.count$Count.OTU)
tax.conversion

dat_g2_pass_rel <- prune_taxa(pass$OTU, dat_g1_fail_rel)
dat_g2_pass_counts <- prune_taxa(pass.count$Count.OTU, 
                                 dat_g1_fail_counts)


# QC Relative abundances sum -----------------------------------------
dat_g2_pass_rel
dat_g1_pass_rel

melt_rel_group = rbind(psmelt(dat_g1_pass_rel), psmelt(dat_g2_pass_rel)) %>% group_by(Sample) %>% #Relative abundance will be less than 1 because discard genus instead of glom
  dplyr::summarize(total_Abundance = sum(Abundance))
melt_rel_group
summary(melt_rel_group)



sink("postglom_abundance.txt")
print("Post-glom relative abundance (ASV + Resolved to Genus)")
print("Summary statistics")
summary(melt_rel_group)
print("Sample where total relative abundance >0.9")
melt_rel_group %>% filter(total_Abundance < 0.9)
sink()

tiff("Postglom_abundance.tiff", res = 180, height = 1000, width = 500)
boxplot(melt_rel_group$total_Abundance)
dev.off()

# Combine Phyloseq objects ------------------------------------------------
dat_glom <- merge_phyloseq(dat_g1_pass_counts, dat_g2_pass_counts)
dat_rel_glom <- merge_phyloseq(dat_g1_pass_rel, dat_g2_pass_rel)
dat_glom_aldex <- merge_phyloseq(dat_g1_pass_counts, #ASVs (filtered)
                                 dat_g1_fail_counts) #genus (non-filtered)


#QC Relative abundance of merged relative abundance Phyloseq object
melt_rel_group_1 = rbind(psmelt(dat_g1_pass_rel), psmelt(dat_g2_pass_rel)) %>% group_by(Sample) %>% dplyr::summarize(total_Abundance = sum(Abundance)) #Rel abundance of unmerged objects combined as dataframes
melt_rel_group_2 = psmelt(dat_rel_glom) %>% group_by(Sample) %>% dplyr::summarize(total_Abundance = sum(Abundance)) #Rel abundance of merged objects

# Add "Other" taxa to account for filtered or lost ASVs/genus -------------------------
abund_diff <- cbind(dat_glom %>% otu_table %>% data.frame() %>% colSums() %>% data.frame(), 
                    dat_pr %>% otu_table %>% data.frame() %>% colSums() %>% data.frame())
colnames(abund_diff) <- c("glom_abund", "total_abund")
abund_diff <- abund_diff %>% mutate(Other = total_abund - glom_abund) %>% mutate(percent_diff = round(Other/total_abund * 100, 2))
other_taxa <- abund_diff %>% dplyr::select(Other) %>% t() #other = other taxa

# Create phyloseq other ---------------------------------------------------
#Create phyloseq object containing the taxa lost during taxa agglomeration or filtering (Taxa = Other)
#mapdf
mapdf.other <- mapdf

#ASV df
asvdf.other <- other_taxa %>% as.matrix()
row.names(asvdf.other) <- "other"

#Taxa df
taxdf.other <- data.frame(c("Other", "Other", "Other", "Other", "Other", "Other")) %>% t() %>% data.frame()
row.names(taxdf.other) <- "other"
colnames(taxdf.other) <-  colnames(taxdf)
taxdf.other <- taxdf.other %>% as.matrix()

#Check structure
head(taxdf.other)
head(asvdf.other)
head(mapdf.other)


#Create phyloseq object 
dat_other <- phyloseq(otu_table(asvdf.other, taxa_are_rows = TRUE), tax_table(taxdf.other), sample_data(mapdf.other))
dat_other %>% otu_table()


# Merge Phyloseq Other ----------------------------------------------------
dat_glom_aldex <- merge_phyloseq(dat_glom, #Filtered ASVs and Genus
                                 dat_other) #Lost taxa





# Export QC tests ---------------------------------------------------------
sink("postglom_merged_QCtests.txt")
print("===== Relative Abundance Tests =====")
print("Tests whether total relative abundance of each sample is equal between filtered/glommedrelative abundance merged and unmerged ASV-level & genus-level phyloseq objects")
print("If it returns TRUE, it means it's the same between the merged and unmerged merged and unmerged ASV-level & genus-level phyloseq objects (i.e. behaving as expected")
print(melt_rel_group_1 == melt_rel_group_2)
print("Summary statistics of filtered/glommed relative abundance ASV-level and genus-level phyloseq objects that have been merged to 1 phyloseq object")
print(summary(melt_rel_group_2))

print("===== Total count test =====")
print("Tests whether total ASV/genus abundance of each sample is equal between filtered/glommedrelative merged and unmerged ASV-level & genus-level phyloseq objects")
print("If it returns TRUE, it means it's the same between the merged and unmerged merged and unmerged ASV-level & genus-level phyloseq objects (i.e. behaving as expected")
print(psmelt(dat_glom) %>% group_by(Sample) %>% dplyr::summarize(total_Abundance = sum(Abundance)) == 
  rbind(psmelt(dat_g1_pass_counts), psmelt(dat_g2_pass_counts)) %>% group_by(Sample) %>% dplyr::summarize(total_Abundance = sum(Abundance)))


print("===== Species Tests =====")
print("Tests whether thet species in the taxa table of the merged phyloseq object is the same as the combination of the two individual ASV-level and genus-level phyloseq object taxa tables")
print("If it returns 0, then it means all of the species are the same (i.e. behaving as expected)")
print("Count phyloseq object")
print(setdiff(tax_table(dat_glom) %>% data.frame, 
        rbind(tax_table(dat_g1_pass_counts) %>% data.frame(), tax_table(dat_g2_pass_counts) %>% data.frame())))
print("Relative abundance phyloseq object")
print(setdiff(tax_table(dat_rel_glom) %>% data.frame, 
        rbind(tax_table(dat_g1_pass_rel) %>% data.frame(), tax_table(dat_g2_pass_rel) %>% data.frame())))


print("==== Species Test - Unmerged ====")
print("Tests whether if any of the species in the unmerged ASV-level and genus-level taxa are shared - this could cause issues with merging")
print("If it returns 0, then it means none of the species are the same (i.e. behaving as expected)")
print("Count ASV/genus phyloseq")
print(rbind(tax_table(dat_g1_pass_counts) %>% data.frame(), tax_table(dat_g2_pass_counts)) %>% rownames_to_column(("name")) %>% filter(duplicated(.)))
print("Relative abundance ASV/genus phyloseq objects ")
print(rbind(tax_table(dat_g1_pass_rel) %>% data.frame(), tax_table(dat_g2_pass_rel)) %>% rownames_to_column(("name")) %>% filter(duplicated(.)))

print("==== Total Abundance - dat_glom_aldex ====")
print("Total abundance of original physeq (dat_pr) should be equivalent to total abundance of dat_glom_aldex, which has Other taxa incorporoated")
print("Other taxa is the taxa lost during taxa agglomeration or prevalence filtering")
print("If this test returns 0, then there are no differences in total abundance between dat_pr and dat_glom_aldex (working as intended")
qc.1 <- dat_pr %>% otu_table %>% data.frame() %>% colSums()
qc.2 <- dat_glom_aldex %>% otu_table %>% data.frame() %>% colSums()
setdiff(qc.1, qc.2)


sink()



# Print merged phyloseq object names  -------------------------------------
print("Glommed phyloseq object - counts")
print("dat_glom")
print(dat_glom)
print("Glommed phyloseq object - relative abundance")
print("dat_rel_glom")
print(dat_rel_glom)
print("Glommed phyloseq object - unfiltered genus-level - for Aldex (Counts)")
print("Contains OTHER taxa group for accurate CLR calculation and/or Aldex comparisons")
print("dat_glom_aldex")
print(dat_glom_aldex)


# Export Phyloseq objects -------------------------------------------------
setwd("/Users/sar/Dropbox/Dallas2K/PhyloseqV3/paper/Physeq_glommed_objects")
save(dat_glom, dat_rel_glom, dat_glom_aldex, 
     dat_pr, dat_empty, dat_empty_rel, dat_all, file = "physeq_ResolveGenus.RData")
setwd("/Users/sar/Dropbox/Dallas2K/PhyloseqV3/paper/ResolveGenus")

