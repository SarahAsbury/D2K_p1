#D2K
#Add α diveristy and β diveristy to analysis of glommed phyloseq object

#Updated: 20oct21 
#Load glommed phyloseq object instead of executing glom. Remove redundant code deal with in new glomming algorithim.

#Updated 9nov21
#Update inputs for physeq3 dfs 

#Updated 20jan21
#Updated inputs for Paper (n = 205) phyloseq objects and dataframes
#Adjusted figures for paper 

#Updated 9feb21 
#Change to paper_v2 dataset


# Set wd ------------------------------------------------------------------
#Clear environment
rm(list = ls(all.names = TRUE))

#Set wd
main.wd <- "/Users/sar/Dropbox/Dallas2K/PhyloseqV3/Paper_v2/diversity"
setwd(main.wd)

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



# Load colour palette -----------------------------------------------------
load.colors_contrast <- c("#114f89ff", "#c11319ff", "#41721eff", "#e8b91eff", "#5c1f7cff", "#140909ff",
                     "#D991BA", "#65743A", "#b1f8f2", "#32E875", "#90A9B7", "#820263")
pie(rep(1, length(load.colors_contrast)), col = load.colors_contrast)

colors_contrast <- load.colors_contrast[c(4:1)]
pie(rep(1, length(colors_contrast)), col = colors_contrast)

# Load glommed data--------------------------------------------------------
load("/Users/sar/Dropbox/Dallas2K/PhyloseqV3/Paper_v2/Physeq_glommed_objects/physeq_ResolveGenus.RData")
dat_glom #Filtered phyloseq object - Counts
dat_rel_glom #Filtered phyloseq object - Rel Abundance
dat_glom_aldex #Contains OTHER taxa so that seq depth is preserved for CLR calculation 

# ===== Alpha Diversity ===== ---------------------------------------------
# Alpha diversity - calculate ---------------------------------------------------------
#Set physeq object
dat <- dat_glom_aldex

#Calculate alpha diversity
alpha_summary = estimate_richness(dat, measures=c("Shannon", "Chao1", "InvSimpson", "Observed"))
str(alpha_summary)

# Merge phyloseq object map df and alpha summary tables
mapdf_alpha = cbind(sample_data(dat) %>% data.frame %>% mutate(mood = ifelse(depression == "control" & anxiety == "control" & anhedonia == "anhedonia",
                                                              "NONE", mood)),  #Anhedonic non-depressed non-anxious patients added as controls as well
                    alpha_summary)

str(mapdf_alpha) # Verify structure
summary(mapdf_alpha) # Check the table to make sure the sample IDs are included.

#Add total reads to mapdf_alpha
total_reads = sample_sums(dat)
mapdf_alpha$TotalReads = total_reads # Add total reads column 
mapdf_alpha %>% head# Re-Assess

# Alpha diversity - graphs ------------------------------------------------
theme_alpha <- ggplot2::theme(
  plot.title = ggplot2::element_text(color="black", size=11, face="bold", hjust = 0.75),
  axis.title.x = ggplot2::element_text(color="black", size=9, face="bold"),
  axis.text.x = ggplot2::element_text(angle = 90, size = 6, face="plain", hjust =1), 
  axis.title.y = ggplot2::element_text(color="black", size=9, face="bold"),
  legend.title = ggplot2::element_text(colour="black", size=9, face="bold"))


p1 <-  mapdf_alpha %>% drop_na(mood) %>% 
  ggpubr::ggboxplot(x = "mood", y = "Shannon", color = "mood", fill = "mood", alpha = 0.3) +
  labs(title="",
       x = "", 
       y = 'Shannon Index') + theme_alpha
p1 <- p1 %>% set_palette(rev(colors_contrast[1:4]))

p2 <- mapdf_alpha %>% drop_na(mood) %>%
  ggpubr::ggboxplot(x = "mood", y = "Chao1", color = "mood", fill = "mood", alpha = 0.3) +
  labs(title="",
       x = "", 
       y = 'Chao1') +
  theme_alpha 
p2 <- p2 %>% set_palette(rev(colors_contrast[1:4]))



p3 <- mapdf_alpha %>% drop_na(mood) %>%
  ggpubr::ggboxplot(x = "mood", y = "InvSimpson", color = "mood", fill = "mood", alpha = 0.3) +
  labs(title="",
       x = "", 
       y = 'InvSimpson') +
  theme_alpha 
p3 <- p3 %>% set_palette(rev(colors_contrast[1:4]))



p4 <- mapdf_alpha %>% drop_na(mood) %>%
  ggpubr::ggboxplot(x = "mood", y = "Observed", color = "mood", fill = "mood", alpha = 0.3) +
  labs(title="",
       x = "", 
       y = 'Observed') +
  theme_alpha 
p4 <- p4 %>% set_palette(rev(colors_contrast[1:4]))


p.alpha <- cowplot::plot_grid(p1 + theme(legend.position = "none"), p2 + theme(legend.position = "none"), 
                   p3 + theme(legend.position = "none"), p4 + theme(legend.position = "none"), ncol = 2)
p.alpha




graphG2 = ggscatter(mapdf_alpha, x = "TotalReads", y = "Shannon",
                    add = "reg.line", conf.int = TRUE, 
                    cor.coef = TRUE, cor.method = "spearman",
                    xlab = "Total Reads", ylab = "Shannon Index",
                    title = "Relationship Between Total Reads and Shannon Index")



# α-diversity Kruskal Wallis - Shannon ----------------------------------------------
# Shannon Index
# Mood
kruskal <- kruskal.test(Shannon ~ mood, mapdf_alpha)


kruskal_mood <- data.frame(kruskal$data.name,
                      kruskal$p.value,
                      kruskal$statistic,
                      kruskal$parameter)
kruskal.multi <- kruskalmc(Shannon ~ mood, mapdf_alpha)
kruskal.multi
dunn <- dunnTest(Shannon ~ mood, data = mapdf_alpha, method="bh")
dunn <- dunn$res


#Age
#mapdf_alpha <- mapdf_alpha %>% mutate(age_categorical = )
#kruskal <- kruskal.test(Shannon ~ age, mapdf_alpha)
#kruskal_age <- data.frame(kruskal$data.name,
                      #kruskal$p.value,
                      #kruskal$statistic,
                      #kruskal$parameter)


# Alpha diveristy Other indices -----------------------------------------------------------
kruskal.test(Chao1 ~ mood, mapdf_alpha)
kruskal.test(InvSimpson ~ mood, mapdf_alpha)
kruskal.test(Observed ~ mood, mapdf_alpha)

kruskalmc(Chao1 ~ mood, mapdf_alpha)
kruskalmc(InvSimpson ~ mood, mapdf_alpha)
kruskalmc(Observed ~ mood, mapdf_alpha)
# Export alpha diversity results ------------------------------------------
setwd(main.wd)
dir.create("alpha diversity")
setwd(paste0(main.wd,"/alpha diversity"))

write.csv(kruskal_mood, "kruskal.csv", row.names = FALSE)
write.csv(kruskal.multi, "kruskal_multiple_comparisons.csv")
write.csv(dunn, "dunn.csv", row.names = FALSE)
tiff("alpha_plots.tiff", width = 1800, height = 1800, res = 300)
p.alpha
dev.off()
tiff("alpha_seqdepth.tiff", width = 1000, height = 400, res = 100)
graphG2
dev.off()
tiff("shannon_plot.tiff", width = 1500, height = 1500, res = 300)
p1 + theme(plot.margin = unit(c(0,2,0,1), "cm"))
dev.off()

setwd(main.wd)

# Beta-Diversity Prep: Convert to CLR -------------------------------------
#List of all taxa except other
taxa <- dat_glom %>% tax_table() %>% data.frame() %>% 
  rownames_to_column(var = "OTU") %>% dplyr::select(OTU)

#Convert to CLR
dat_clr <- microbiome::transform(dat_glom_aldex, "clr") %>% 
  prune_taxa(taxa$OTU, .)
otu_table(dat_clr)[1:5, 1:5]

# Beta Diveristy - Aitchisons ---------------------------------------------
#https://www.nicholas-ollberding.com/post/introduction-to-the-statistical-analysis-of-microbiome-data-in-r/
dat_input <- dat_clr %>% subset_samples(!(is.na(mood))) #181 samples
ord_clr <- phyloseq::ordinate(dat_input, "RDA")

#Scree plot
phyloseq::plot_scree(ord_clr) + 
  geom_bar(stat="identity", fill = "blue") +
  labs(x = "\nAxis", y = "Proportion of Variance\n")

#Assess proportion of variance explained
head(ord_clr$CA$eig)                                                  
sapply(ord_clr$CA$eig[1:5], function(x) x / sum(ord_clr$CA$eig))


#Scale axes and plot ordination
clr1 <- ord_clr$CA$eig[1] / sum(ord_clr$CA$eig)
clr2 <- ord_clr$CA$eig[2] / sum(ord_clr$CA$eig)
p_aitch <- phyloseq::plot_ordination(dat_glom, ord_clr, type="samples", color="mood") +
  scale_colour_manual(values=colors_contrast) +
  geom_point(size = 2) +
  coord_fixed(clr2 / clr1) +
  theme_classic()
p_aitch

#Generate distance matrix
distance <- phyloseq::distance(dat_input, method = "euclidean") 
dim(distance)
class(distance)


#PERMANOVA: Mood
distance <- phyloseq::distance(dat_input, method = "euclidean") 
mood <- phyloseq::sample_data(dat_input)$mood
permanova.aitch <- vegan::adonis2(distance ~ mood,
                    permutations=9999, paralell = 8)
permanova.aitch

#Dispersion test
dispr <- vegan::betadisper(distance, phyloseq::sample_data(dat_input)$mood)
dispr
plot(dispr, main = "Ordination Centroids and Dispersion Labeled: Aitchison Distance", sub = "")
boxplot(dispr, main = "", xlab = "")
dispr.test <- permutest(dispr)
dispr.test


#PERMANOVA: Mood + Age
dat_input_perm <- dat_input %>% subset_samples(!(is.na(age)))
distance <- phyloseq::distance(dat_input_perm, method = "euclidean") 
meta <- dat_input_perm %>% sample_data %>% data.frame() %>% 
  dplyr::select(mood, age)

permanova.aitch.age_effect <- vegan::adonis2(distance ~ mood + age,
                                  permutations=9999, paralell = 8,
                                  data = meta)
permanova.aitch.age_interact <- vegan::adonis2(distance ~ mood + age + mood:age,
                                             permutations=9999, paralell = 8,
                                             data = meta)
permanova.aitch.age_effect
permanova.aitch.age_interact


#PERMANOVA: Mood + Sex
dat_input_perm <- dat_input %>% subset_samples(!(is.na(dsex)))
distance <- phyloseq::distance(dat_input_perm, method = "euclidean") 
meta <- dat_input_perm %>% sample_data %>% data.frame() %>% 
  dplyr::select(mood, dsex)

permanova.aitch.sex_effect <- vegan::adonis2(distance ~ mood + dsex,
                                             permutations=9999, paralell = 8,
                                             data = meta)
permanova.aitch.sex_interact <- vegan::adonis2(distance ~ mood + dsex + mood:dsex,
                                               permutations=9999, paralell = 8,
                                               data = meta)

permanova.aitch.sex_effect
permanova.aitch.sex_interact

#PERMANOVA: Mood + Race
dat_input_perm <- dat_input %>% subset_samples(!(is.na(race)))
distance <- phyloseq::distance(dat_input_perm, method = "euclidean") 
meta <- dat_input_perm %>% sample_data %>% data.frame() %>% 
  dplyr::select(mood, race)

permanova.aitch.race_effect <- vegan::adonis2(distance ~ mood + race,
                                             permutations=9999, paralell = 8,
                                             data = meta)
permanova.aitch.race_interact <- vegan::adonis2(distance ~ mood + race + mood:race,
                                               permutations=9999, paralell = 8,
                                               data = meta)
permanova.aitch.race_effect
permanova.aitch.race_interact



# Beta Diveristy - Aitchison's - Remove GAD ---------------------------------------
dat_input <- dat_clr %>% subset_samples(!(is.na(mood)) & mood!= "GAD") #176 samples
ord_clr <- phyloseq::ordinate(dat_input, "RDA")

#Scree plot
phyloseq::plot_scree(ord_clr) + 
  geom_bar(stat="identity", fill = "blue") +
  labs(x = "\nAxis", y = "Proportion of Variance\n")

#Assess proportion of variance explained
head(ord_clr$CA$eig)                                                  
sapply(ord_clr$CA$eig[1:5], function(x) x / sum(ord_clr$CA$eig))


#Scale axes and plot ordination
clr1 <- ord_clr$CA$eig[1] / sum(ord_clr$CA$eig)
clr2 <- ord_clr$CA$eig[2] / sum(ord_clr$CA$eig)
p_aitch_mdd <- phyloseq::plot_ordination(dat_glom, ord_clr, type="samples", color="mood") +
  scale_colour_manual(values=colors_contrast[-1]) +
  geom_point(size = 2) +
  coord_fixed(clr2 / clr1) +
  theme_classic()
p_aitch_mdd

#Generate distance matrix
distance <- phyloseq::distance(dat_input, method = "euclidean") 
dim(distance)
class(distance)


#PERMANOVA: Mood
distance <- phyloseq::distance(dat_input, method = "euclidean") 
mood <- phyloseq::sample_data(dat_input)$mood
permanova.aitch.mdd <- vegan::adonis2(distance ~ mood,
                                  permutations=9999, paralell = 8)
permanova.aitch.mdd

#Dispersion test
dispr.mdd <- vegan::betadisper(distance, phyloseq::sample_data(dat_input)$mood)
dispr.mdd
plot(dispr.mdd, main = "Ordination Centroids and Dispersion Labeled: Aitchison Distance", sub = "")
boxplot(dispr.mdd, main = "", xlab = "")
dispr.mdd.test <- permutest(dispr.mdd)


# Beta Diversity - Combine results ----------------------------------------
#Adonis PERMANOVA continuous variables: https://www.researchgate.net/post/How_ADONIS_is_calculated_for_continues_variable#view=5be7816bf0fb6228f959acda
  #I.e combines it with a linear regression
#Results combined
#bdres.fun
#Creates a dataframe with β-diversity results from adonis 
#x = list containing the name of adonis objects to extract results from   
bdres.fun <- function(x){
  df <- data.frame()
  for (i in x){
    x <- eval(parse(text = paste0(i, "$aov.tab"))) %>% data.frame() %>% 
      mutate(formula = eval(parse(text = paste(i, "$call$formula"))) %>% as_label, 
                                                                               bdiverse = eval(parse(text = paste(i, "$call$method"))) %>% as_label,
                                                                               result_name = paste(i)) %>%
      rownames_to_column(var = "result_var")
    
    df <- rbind(df, x)
  }
  df <- df %>% mutate(formula_vars = gsub(".*~ ", "", formula)) %>% select(-formula)
  print(df)
}




bdres.fun.aitch <- function(x){
  df <- data.frame()
  for (i in x){
    x <- eval(parse(text = paste0(i))) %>% data.frame() %>% 
      mutate(formula_vars = eval(parse(text = paste0("attributes(", i, ")$heading[2]"))) %>%
               gsub(".*distance ~ ", "", .) %>% gsub(", .*", "", .),
             bdiverse = "aitch",
             result_name = i
             ) %>%
      rownames_to_column(var = "result_var")
      
    df <- rbind(df, x)
    
  }
  print(df)
}


#Beta diversity results - extract from adonis2 
bdiverse.results.2 <- bdres.fun.aitch(x = c("permanova.aitch",
                                    "permanova.aitch.age_effect", "permanova.aitch.age_interact",
                                    "permanova.aitch.sex_effect", "permanova.aitch.race_interact",
                                    "permanova.aitch.race_effect", "permanova.aitch.race_interact",
                                    "permanova.aitch.mdd"))
bdiverse.results <- bdiverse.results.2 %>%  mutate(formula_vars = ifelse(grepl("race_interact", result_name) == TRUE,
                               "mood + race + mood:race", formula_vars))

# Beta-Diversity: Exports -------------------------------------------------
setwd(main.wd)
dir.create("b_diverse")
setwd(paste0(main.wd, "/b_diverse"))

write.csv(bdiverse.results, "bdiverse_results.csv", row.names = FALSE)
write.csv(bdiverse.results %>% filter(!(is.na(Pr..F.)) & Pr..F. < 0.1), "bdiverse_sigresults.csv", row.names = FALSE)
write.csv(dispr.test$tab, "bdiverse_dispersion.csv")
write.csv(dispr.mdd.test$tab, "bdiverse_nogad_dispersion.csv")

#Beta diversity PCA plots:
pheight = 1250
pwidth = 1250
res = 180

tiff("aitch_beta_plot.tiff", width = pheight, height = pwidth, res = res)
p_aitch
dev.off()

tiff("aitch_mdd_beta_plot.tiff", width = pheight, height = pwidth, res = res)
p_aitch_mdd
dev.off()
