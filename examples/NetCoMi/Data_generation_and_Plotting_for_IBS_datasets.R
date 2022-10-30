library(phyloseq)
library(NetCoMi)
library(plyr)
library(dplyr)
library("RColorBrewer")
library(stringr)

################################################################################

#Load all 12 IBS datasets phyloseq obects (except labus as it does not have species level information & except ringel as it does not have healthy / non-IBS classification)
agp_phyloseq <- readRDS("./input_data/phyloseq_objects/physeq_agp.rds")
fukui_phyloseq  <- readRDS("./input_data/phyloseq_objects/physeq_fukui.rds")
hugerth_phyloseq  <- readRDS("./input_data/phyloseq_objects/physeq_hugerth.rds")
labus_phyloseq <- NULL
liu_phyloseq  <- readRDS("./input_data/phyloseq_objects/physeq_liu.rds")
lopresti_phyloseq  <- readRDS("./input_data/phyloseq_objects/physeq_lopresti.rds")
mars_phyloseq  <- readRDS("./input_data/phyloseq_objects/physeq_mars.rds")
nagel_phyloseq  <- readRDS("./input_data/phyloseq_objects/physeq_nagel.rds")
pozuelo_phyloseq  <- readRDS("./input_data/phyloseq_objects/physeq_pozuelo.rds")
ringel_phyloseq  <- readRDS("./input_data/phyloseq_objects/physeq_ringel.rds")
zeber_phyloseq  <- readRDS("./input_data/phyloseq_objects/physeq_zeber.rds")
zhu_phyloseq  <- readRDS("./input_data/phyloseq_objects/physeq_zhu.rds")
zhuang_phyloseq  <- readRDS("./input_data/phyloseq_objects/physeq_zhuang.rds")

################################################################################

#Load 12 IBS datasets bacdive information from Python output
agp_bacdive <- read.table("./output_data/agp/NetCoMi/NetCoMi.tsv")
fukui_bacdive  <- read.table("./output_data/fukui/NetCoMi/NetCoMi.tsv")
hugerth_bacdive  <- read.table("./output_data/hugerth/NetCoMi/NetCoMi.tsv")
labus_bacdive <- NULL
liu_bacdive  <- read.table("./output_data/liu/NetCoMi/NetCoMi.tsv")
lopresti_bacdive  <- read.table("./output_data/lopresti/NetCoMi/NetCoMi.tsv")
mars_bacdive  <- read.table("./output_data/mars/NetCoMi/NetCoMi.tsv")
nagel_bacdive  <- read.table("./output_data/nagel/NetCoMi/NetCoMi.tsv")
pozuelo_bacdive  <- read.table("./output_data/pozuelo/NetCoMi/NetCoMi.tsv")
ringel_bacdive  <- read.table("./output_data/ringel/NetCoMi/NetCoMi.tsv")
zeber_bacdive  <- read.table("./output_data/zeber/NetCoMi/NetCoMi.tsv")
zhu_bacdive  <- read.table("./output_data/zhu/NetCoMi/NetCoMi.tsv")
zhuang_bacdive  <- read.table("./output_data/zhuang/NetCoMi/NetCoMi.tsv")

################################################################################

# Agglomerate to species level
agp_agglom_data_species <- tax_glom(agp_phyloseq, taxrank = "Species")
fukui_agglom_data_species <- tax_glom(fukui_phyloseq, taxrank = "Species")
hugerth_agglom_data_species <- tax_glom(hugerth_phyloseq, taxrank = "Species")
labus_agglom_data_species <- NULL
liu_agglom_data_species <- tax_glom(liu_phyloseq, taxrank = "Species")
lopresti_agglom_data_species <- tax_glom(lopresti_phyloseq, taxrank = "Species")
mars_agglom_data_species <- tax_glom(mars_phyloseq, taxrank = "Species")
nagel_agglom_data_species <- tax_glom(nagel_phyloseq, taxrank = "Species")
pozuelo_agglom_data_species <- tax_glom(pozuelo_phyloseq, taxrank = "Species")
ringel_agglom_data_species <- tax_glom(ringel_phyloseq, taxrank = "Species")
zeber_agglom_data_species <- tax_glom(zeber_phyloseq, taxrank = "Species")
zhu_agglom_data_species <- tax_glom(zhu_phyloseq, taxrank = "Species")
zhuang_agglom_data_species <- tax_glom(zhuang_phyloseq, taxrank = "Species")

################################################################################

#Subset to only contain agglomeration level data
agp_ASV_input <- row.names(agp_agglom_data_species@tax_table)
agp_bacdive_data <- agp_bacdive[row.names(agp_bacdive) %in% agp_ASV_input, ]

fukui_ASV_input <- row.names(fukui_agglom_data_species@tax_table)
fukui_bacdive_data <- fukui_bacdive[row.names(fukui_bacdive) %in% fukui_ASV_input, ]

hugerth_ASV_input <- row.names(hugerth_agglom_data_species@tax_table)
hugerth_bacdive_data <- hugerth_bacdive[row.names(hugerth_bacdive) %in% hugerth_ASV_input, ]

labus_ASV_input <- NULL
labus_bacdive_data <- NULL

liu_ASV_input <- row.names(liu_agglom_data_species@tax_table)
liu_bacdive_data <- liu_bacdive[row.names(liu_bacdive) %in% liu_ASV_input, ]

lopresti_ASV_input <- row.names(lopresti_agglom_data_species@tax_table)
lopresti_bacdive_data <- lopresti_bacdive[row.names(lopresti_bacdive) %in% lopresti_ASV_input, ]

mars_ASV_input <- row.names(mars_agglom_data_species@tax_table)
mars_bacdive_data <- mars_bacdive[row.names(mars_bacdive) %in% mars_ASV_input, ]

nagel_ASV_input <- row.names(nagel_agglom_data_species@tax_table)
nagel_bacdive_data <- nagel_bacdive[row.names(nagel_bacdive) %in% nagel_ASV_input, ]

pozuelo_ASV_input <- row.names(pozuelo_agglom_data_species@tax_table)
pozuelo_bacdive_data <- pozuelo_bacdive[row.names(pozuelo_bacdive) %in% pozuelo_ASV_input, ]

ringel_ASV_input <- row.names(ringel_agglom_data_species@tax_table)
ringel_bacdive_data <- ringel_bacdive[row.names(ringel_bacdive) %in% ringel_ASV_input, ]

zeber_ASV_input <- row.names(zeber_agglom_data_species@tax_table)
zeber_bacdive_data <- zeber_bacdive[row.names(zeber_bacdive) %in% zeber_ASV_input, ]

zhu_ASV_input <- row.names(zhu_agglom_data_species@tax_table)
zhu_bacdive_data <- zhu_bacdive[row.names(zhu_bacdive) %in% zhu_ASV_input, ]

zhuang_ASV_input <- row.names(zhuang_agglom_data_species@tax_table)
zhuang_bacdive_data <- zhuang_bacdive[row.names(zhuang_bacdive) %in% zhuang_ASV_input, ]

################################################################################

# Rename taxonomic table and make Rank7 (species) unique
agp_species_renamed <- renameTaxa(agp_agglom_data_species, 
                                   pat = "<name>", 
                                   substPat = "<name>_<subst_name>(<subst_R>)", numDupli = "Species")

# Rename taxonomic table and make Rank7 (species) unique
fukui_species_renamed <- renameTaxa(fukui_agglom_data_species, 
                                   pat = "<name>", 
                                   substPat = "<name>_<subst_name>(<subst_R>)", numDupli = "Species")


# Rename taxonomic table and make Rank7 (species) unique
hugerth_species_renamed <- renameTaxa(hugerth_agglom_data_species, 
                                   pat = "<name>", 
                                   substPat = "<name>_<subst_name>(<subst_R>)", numDupli = "Species")


# Rename taxonomic table and make Rank7 (species) unique
labus_species_renamed <- NULL

# Rename taxonomic table and make Rank7 (species) unique
liu_species_renamed <- renameTaxa(liu_agglom_data_species, 
                                   pat = "<name>", 
                                   substPat = "<name>_<subst_name>(<subst_R>)", numDupli = "Species")

# Rename taxonomic table and make Rank7 (species) unique
lopresti_species_renamed <- renameTaxa(lopresti_agglom_data_species, 
                                   pat = "<name>", 
                                   substPat = "<name>_<subst_name>(<subst_R>)", numDupli = "Species")


# Rename taxonomic table and make Rank7 (species) unique
mars_species_renamed <- renameTaxa(mars_agglom_data_species, 
                                   pat = "<name>", 
                                   substPat = "<name>_<subst_name>(<subst_R>)", numDupli = "Species")


# Rename taxonomic table and make Rank7 (species) unique
nagel_species_renamed <- renameTaxa(nagel_agglom_data_species, 
                                   pat = "<name>", 
                                   substPat = "<name>_<subst_name>(<subst_R>)", numDupli = "Species")


# Rename taxonomic table and make Rank7 (species) unique
pozuelo_species_renamed <- renameTaxa(pozuelo_agglom_data_species, 
                                   pat = "<name>", 
                                   substPat = "<name>_<subst_name>(<subst_R>)", numDupli = "Species")


# Rename taxonomic table and make Rank7 (species) unique
ringel_species_renamed <- renameTaxa(ringel_agglom_data_species, 
                                   pat = "<name>", 
                                   substPat = "<name>_<subst_name>(<subst_R>)", numDupli = "Species")


# Rename taxonomic table and make Rank7 (species) unique
zeber_species_renamed <- renameTaxa(zeber_agglom_data_species, 
                                   pat = "<name>", 
                                   substPat = "<name>_<subst_name>(<subst_R>)", numDupli = "Species")


# Rename taxonomic table and make Rank7 (species) unique
zhu_species_renamed <- renameTaxa(zhu_agglom_data_species, 
                                   pat = "<name>", 
                                   substPat = "<name>_<subst_name>(<subst_R>)", numDupli = "Species")


# Rename taxonomic table and make Rank7 (species) unique
zhuang_species_renamed <- renameTaxa(zhuang_agglom_data_species, 
                                   pat = "<name>", 
                                   substPat = "<name>_<subst_name>(<subst_R>)", numDupli = "Species")


################################################################################

#Combine Genus and Species column to new Species column
update_phyloseq <- function(data_species_renamed, bacdive_data) {
  combined <- as.data.frame(data_species_renamed@tax_table)
  combined$concat_species <- str_c(combined$Genus, " ", combined$Species)
  combined$Species <- NULL
  colnames(combined)[which(names(combined) == "concat_species")] <- "Species"
  taxtab <- as(combined, "matrix")
  tax_table(data_species_renamed) <- taxtab
  
  return(data_species_renamed)
}

agp_phyloseq <- update_phyloseq(agp_species_renamed, agp_bacdive_data)
fukui_phyloseq <- update_phyloseq(fukui_species_renamed, fukui_bacdive_data)
hugerth_phyloseq <- update_phyloseq(hugerth_species_renamed, hugerth_bacdive_data)
labus_phyloseq <- NULL
liu_phyloseq <- update_phyloseq(liu_species_renamed, liu_bacdive_data)
lopresti_phyloseq <- update_phyloseq(lopresti_species_renamed, lopresti_bacdive_data)
mars_phyloseq <- update_phyloseq(mars_species_renamed, mars_bacdive_data)
nagel_phyloseq <- update_phyloseq(nagel_species_renamed, nagel_bacdive_data)
ringel_phyloseq <- update_phyloseq(ringel_species_renamed, ringel_bacdive_data)
pozuelo_phyloseq <- update_phyloseq(pozuelo_species_renamed, pozuelo_bacdive_data)
zeber_phyloseq <- update_phyloseq(zeber_species_renamed, zeber_bacdive_data)
zhu_phyloseq <- update_phyloseq(zhu_species_renamed, zhu_bacdive_data)
zhuang_phyloseq <- update_phyloseq(zhuang_species_renamed, zhuang_bacdive_data)

################################################################################

#Construct network
networks <- function(data_species_renamed, bacdive_data, taxa) {
  
  colnames(data_species_renamed@otu_table@.Data) <- data_species_renamed@tax_table@.Data[, "Species"]
  countMat <- data_species_renamed@otu_table@.Data
  head(countMat)
  
  #split into healthy and IBS
  group_vec <- phyloseq::get_variable(data_species_renamed, "host_disease")
  
  # Network construction NetCoMi command
  net_single3 <- netConstruct(data = countMat,
                              group = group_vec,
                              filtTax = "highestFreq",
                              filtTaxPar = list(highestFreq = taxa),
                              jointPrepro = TRUE,
                              verbose = 3,
                              seed = 112233)
  
  return(net_single3)
}


analyze_network <- function(net_single3) {
  props_single3 <- netAnalyze(net_single3, 
                      centrLCC = TRUE,
                      avDissIgnoreInf = TRUE,
                      sPathNorm = FALSE,
                      clustMethod = "cluster_fast_greedy",
                      hubPar = "eigenvector",
                      weightDeg = FALSE, 
                      normDeg = TRUE)
  return(props_single3)
}

################################################################################

naming <- c("agp", "fukui", "hugerth", "liu", "lopresti", "mars", "nagel", "pozuelo", "ringel", "zeber", "zhu", "zhuang")
datasets.species <- list(agp_phyloseq, fukui_phyloseq, hugerth_phyloseq, liu_phyloseq, lopresti_phyloseq, mars_phyloseq, nagel_phyloseq, pozuelo_phyloseq, ringel_phyloseq, zeber_phyloseq, zhu_phyloseq, zhuang_phyloseq)
names(datasets.species) <- naming
datasets.species

# Saving on object in RData format
save(datasets.species, file = "./input_data/RObjects/datasets.species.RData")

#Reload object
#datasets.species <- load("./input_data/RObjects/datasets.species.RData")

################################################################################

bacdive.data <- list(agp_bacdive_data, fukui_bacdive_data, hugerth_bacdive_data, liu_bacdive_data, lopresti_bacdive_data, mars_bacdive_data, nagel_bacdive_data, pozuelo_bacdive_data, ringel_bacdive_data, zeber_bacdive_data, zhu_bacdive_data, zhuang_bacdive_data)
names(bacdive.data) <- naming
#bacdive.data$labus

# Saving on object in RData format
save(bacdive.data, file = "./input_data/RObjects/bacdive.data.RData")

#Reload object
#bacdive_data <- load("./input_data/RObjects/bacdive.data.RData")

################################################################################

# Number of taxa that should be used for each network scaled to the ln-samplesize
ntaxa.species <- c(agp = 170,
                   fukui = 115,
                   hugerth = 125,
                   liu = 115,
                   loPresti = 90,
                   mars = 100,
                   nagel = 80,
                   pozuelo = 135,
                   ringel = 105,
                   zeber = 110,
                   zhu = 80,
                   zhuang = 80)


networks.species <- mapply(networks, datasets.species[-which(names(datasets.species) == "ringel")], 
                           bacdive.data[-which(names(bacdive.data) == "ringel")], 
                           taxa = ntaxa.species[-which(names(ntaxa.species) == "ringel")], SIMPLIFY = FALSE)

probs.species <- lapply(networks.species, analyze_network)

#save object
save(networks.species, probs.species, file = "./input_data/RObjects/network_object")

#load object
network_object <- load("./input_data/RObjects/network_object")
network_object


################################################################################
plot_networks <- function(net_single3, props_single3, data_species_renamed, bacdive_data, plot_column, title_text) {
  # Compute layout
  graph3 <- igraph::graph_from_adjacency_matrix(net_single3$adjaMat1)
  set.seed(123456)
  lay_fr <- igraph::layout_with_fr(graph3)
  
  # Row names of the layout matrix must match the node names
  rownames(lay_fr) <- rownames(net_single3$adjaMat1)
  
  #Add Bacdive info to plot
  taxtab <- as(tax_table(data_species_renamed), "matrix")
  c_vec <- names(taxtab[, "Species"])
  
  # Rename all not included in Bacdive species from "no" to "not_in_bacdive"
  # in order to distinguish "unknown" in the plot from those species which are not even found in Bacdive.
  
  bacdive_data[, plot_column][bacdive_data$included_in_Bacdive == "no"] <- "not_in_Bacdive"  # rename all values in <column of interest> where "included_in_Bacdive" == "no"
  bacdive_info <- bacdive_data[, plot_column]
  
  #concatenate Genus and species columns into one
  concat_df <- as.data.frame(data_species_renamed@tax_table)
  concat_df$combined <- str_c(concat_df$Genus, " ", concat_df$Species)
  
  names(bacdive_info) <- taxtab[, "Species"]
  bacdive_info 
  table(bacdive_info)
  
  num_col <- length(unique(as.data.frame(table(bacdive_info))$bacdive_info))
  color <- rainbow(num_col)
  
  plot(props_single3,
       layout = "spring",
       sameLayout = TRUE,
       repulsion = 0.7,
       shortenLabels = "none",
       charToRm = "g__",
       labelScale = TRUE,
       rmSingles = TRUE,
       nodeSize = "clr",
       nodeSizeSpread = 4,
       nodeColor = "feature", 
       featVecCol = bacdive_info, 
       colorVec =  color,
       posCol = "#009900", 
       negCol = "red",
       groupNames = c("Healthy", "IBS"),
       edgeTranspLow = 0,
       edgeTranspHigh = 40,
       cexNodes = 2.5,
       cexLabels = 4.5,
       title1 = title_text, 
       showTitle = TRUE,
       cexTitle = 1.0, 
       highlightHubs = TRUE)
  
  # Colors used in the legend should be equally transparent as in the plot
  phylcol_transp <- NetCoMi:::colToTransp(color, 60)
  
  legend("bottom", xpd=TRUE, title = "estimated association:", legend = c("+","-"), 
         col = c("#009900","red"), cex = 0.7,  lty = 1, lwd = 3, box.lty=2, box.lwd=4, box.col="black",
         bty = "n", horiz = TRUE)
  
  legend(1.05, 1.3, cex = 0.6, pt.cex = 1,xpd=TRUE, 
         legend=unique(as.data.frame(table(bacdive_info))$bacdive_info), col = color, bty = "n", pch = 16) 
  
  title(main = title_text,
        cex.main = 0.9, font.main= 4, col.main= "black",
  )
}

################################################################################

#Plotting motility information

#Agp
pdf("./output_data/agp/NetCoMi/Motility_network.pdf", paper="a4r", height = 8.27, width =  11.69)
plot_networks(networks.species$agp, probs.species$agp, agp_phyloseq,
              agp_bacdive_data, plot_column = "Morphology_cell_morphology_motility",
              title_text = "Agp network displaying motility information") 
dev.off()

#Fukui
pdf("./output_data/fukui/NetCoMi/Motility_network.pdf", paper="a4r", height = 8.27, width =  11.69)
plot_networks(networks.species$fukui, probs.species$fukui, fukui_phyloseq,
              fukui_bacdive_data, plot_column = "Morphology_cell_morphology_motility",
              title_text = "Fukui network displaying motility information") 
dev.off()

#Hugerth
pdf("./output_data/hugerth/NetCoMi/Motility_network.pdf", paper="a4r", height = 8.27, width =  11.69)
plot_networks(networks.species$hugerth, probs.species$hugerth, hugerth_phyloseq,
              hugerth_bacdive_data, plot_column = "Morphology_cell_morphology_motility",
              title_text = "Hugerth network displaying motility information") 
dev.off()

#Liu
pdf("./output_data/liu/NetCoMi/Motility_network.pdf", paper="a4r", height = 8.27, width =  11.69)
plot_networks(networks.species$liu, probs.species$liu, liu_phyloseq,
              liu_bacdive_data, plot_column = "Morphology_cell_morphology_motility",
              title_text = "Liu network displaying motility information") 
dev.off()

#Lopresti
pdf("./output_data/lopresti/NetCoMi/Motility_network.pdf", paper="a4r", height = 8.27, width =  11.69)
plot_networks(networks.species$lopresti, probs.species$lopresti, lopresti_phyloseq,
              lopresti_bacdive_data, plot_column = "Morphology_cell_morphology_motility",
              title_text = "Lopresti network displaying motility information") 
dev.off()

#Mars
pdf("./output_data/mars/NetCoMi/Motility_network.pdf", paper="a4r", height = 8.27, width =  11.69)
plot_networks(networks.species$mars, probs.species$mars, mars_phyloseq,
              mars_bacdive_data, plot_column = "Morphology_cell_morphology_motility",
              title_text = "Mars network displaying motility information") 
dev.off()

#Nagel
pdf("./output_data/nagel/NetCoMi/Motility_network.pdf", paper="a4r", height = 8.27, width =  11.69)
plot_networks(networks.species$nagel, probs.species$nagel, nagel_phyloseq,
              nagel_bacdive_data, plot_column = "Morphology_cell_morphology_motility",
              title_text = "Nagel network displaying motility information") 
dev.off()

#Pozuelo
pdf("./output_data/pozuelo/NetCoMi/Motility_network.pdf", paper="a4r", height = 8.27, width =  11.69)
plot_networks(networks.species$pozuelo, probs.species$pozuelo, pozuelo_phyloseq,
              pozuelo_bacdive_data, plot_column = "Morphology_cell_morphology_motility",
              title_text = "Pozuelo network displaying motility information") 
dev.off()

#Zeber
pdf("./output_data/zeber/NetCoMi/Motility_network.pdf", paper="a4r", height = 8.27, width =  11.69)
plot_networks(networks.species$zeber, probs.species$zeber, zeber_phyloseq,
              zeber_bacdive_data, plot_column = "Morphology_cell_morphology_motility",
              title_text = "Zeber network displaying motility information") 
dev.off()

#Zhu
pdf("./output_data/zhu/NetCoMi/Motility_network.pdf", paper="a4r", height = 8.27, width =  11.69)
plot_networks(networks.species$zhu, probs.species$zhu, zhu_phyloseq,
              zhu_bacdive_data, plot_column = "Morphology_cell_morphology_motility",
              title_text = "Zhu network displaying motility information") 
dev.off()

#Zhuang
pdf("./output_data/zhuang/NetCoMi/Motility_network.pdf", paper="a4r", height = 8.27, width =  11.69)
plot_networks(networks.species$zhuang, probs.species$zhuang, zhuang_phyloseq,
              zhuang_bacdive_data, plot_column = "Morphology_cell_morphology_motility",
              title_text = "Zhuang network displaying motility information") 
dev.off()

################################################################################

#Plotting oxygen tolerance information

#Agp
pdf("./output_data/agp/NetCoMi/Oxygen_tolerance_network.pdf", paper="a4r", height = 8.27, width =  11.69)
plot_networks(networks.species$agp, probs.species$agp, agp_phyloseq,
              agp_bacdive_data, plot_column = "Physiology_and_metabolism_oxygen_tolerance_oxygen_tolerance",
              title_text = "Agp network displaying oxygen tolerance information") 
dev.off()

#Fukui
pdf("./output_data/fukui/NetCoMi/Oxygen_tolerance_network.pdf", paper="a4r", height = 8.27, width =  11.69)
plot_networks(networks.species$fukui, probs.species$fukui, fukui_phyloseq,
              fukui_bacdive_data, plot_column = "Physiology_and_metabolism_oxygen_tolerance_oxygen_tolerance",
              title_text = "Fukui network displaying oxygen tolerance information") 
dev.off()

#Hugerth
pdf("./output_data/hugerth/NetCoMi/Oxygen_tolerance_network.pdf", paper="a4r", height = 8.27, width =  11.69)
plot_networks(networks.species$hugerth, probs.species$hugerth, hugerth_phyloseq,
              hugerth_bacdive_data, plot_column = "Physiology_and_metabolism_oxygen_tolerance_oxygen_tolerance",
              title_text = "Hugerth network displaying oxygen tolerance information") 
dev.off()

#Liu
pdf("./output_data/liu/NetCoMi/Oxygen_tolerance_network.pdf", paper="a4r", height = 8.27, width =  11.69)
plot_networks(networks.species$liu, probs.species$liu, liu_phyloseq,
              liu_bacdive_data, plot_column = "Physiology_and_metabolism_oxygen_tolerance_oxygen_tolerance",
              title_text = "Liu network displaying oxygen tolerance information") 
dev.off()

#Lopresti
pdf("./output_data/lopresti/NetCoMi/Oxygen_tolerance_network.pdf", paper="a4r", height = 8.27, width =  11.69)
plot_networks(networks.species$lopresti, probs.species$lopresti, lopresti_phyloseq,
              lopresti_bacdive_data, plot_column = "Physiology_and_metabolism_oxygen_tolerance_oxygen_tolerance",
              title_text = "Lopresti network displaying oxygen tolerance information") 
dev.off()

#Mars
pdf("./output_data/mars/NetCoMi/Oxygen_tolerance_network.pdf", paper="a4r", height = 8.27, width =  11.69)
plot_networks(networks.species$mars, probs.species$mars, mars_phyloseq,
              mars_bacdive_data, plot_column = "Physiology_and_metabolism_oxygen_tolerance_oxygen_tolerance",
              title_text = "Mars network displaying oxygen tolerance information") 
dev.off()

#Nagel
pdf("./output_data/nagel/NetCoMi/Oxygen_tolerance_network.pdf", paper="a4r", height = 8.27, width =  11.69)
plot_networks(networks.species$nagel, probs.species$nagel, nagel_phyloseq,
              nagel_bacdive_data, plot_column = "Physiology_and_metabolism_oxygen_tolerance_oxygen_tolerance",
              title_text = "Nagel network displaying oxygen tolerance information") 
dev.off()

#Pozuelo
pdf("./output_data/pozuelo/NetCoMi/Oxygen_tolerance_network.pdf", paper="a4r", height = 8.27, width =  11.69)
plot_networks(networks.species$pozuelo, probs.species$pozuelo, pozuelo_phyloseq,
              pozuelo_bacdive_data, plot_column = "Physiology_and_metabolism_oxygen_tolerance_oxygen_tolerance",
              title_text = "Pozuelo network displaying oxygen tolerance information") 
dev.off()

#Zeber
pdf("./output_data/zeber/NetCoMi/Oxygen_tolerance_network.pdf", paper="a4r", height = 8.27, width =  11.69)
plot_networks(networks.species$zeber, probs.species$zeber, zeber_phyloseq,
              zeber_bacdive_data, plot_column = "Physiology_and_metabolism_oxygen_tolerance_oxygen_tolerance",
              title_text = "Zeber network displaying oxygen tolerance information") 
dev.off()

#Zhu
pdf("./output_data/zhu/NetCoMi/Oxygen_tolerance_network.pdf", paper="a4r", height = 8.27, width =  11.69)
plot_networks(networks.species$zhu, probs.species$zhu, zhu_phyloseq,
              zhu_bacdive_data, plot_column = "Physiology_and_metabolism_oxygen_tolerance_oxygen_tolerance",
              title_text = "Zhu network displaying oxygen tolerance information") 
dev.off()

#Zhuang
pdf("./output_data/zhuang/NetCoMi/Oxygen_tolerance_network.pdf", paper="a4r", height = 8.27, width =  11.69)
plot_networks(networks.species$zhuang, probs.species$zhuang, zhuang_phyloseq,
              zhuang_bacdive_data, plot_column = "Physiology_and_metabolism_oxygen_tolerance_oxygen_tolerance",
              title_text = "Zhuang network displaying oxygen tolerance information") 
dev.off()

################################################################################

#Plotting gram stain information

#Agp
pdf("./output_data/agp/NetCoMi/Gram_stain_network.pdf", paper="a4r", height = 8.27, width =  11.69)
plot_networks(networks.species$agp, probs.species$agp, agp_phyloseq,
              agp_bacdive_data, plot_column = "Morphology_cell_morphology_gram_stain",
              title_text = "Agp network displaying gram stain information") 
dev.off()

#Fukui
pdf("./output_data/fukui/NetCoMi/Gram_stain_network.pdf", paper="a4r", height = 8.27, width =  11.69)
plot_networks(networks.species$fukui, probs.species$fukui, fukui_phyloseq,
              fukui_bacdive_data, plot_column = "Morphology_cell_morphology_gram_stain",
              title_text = "Fukui network displaying gram stain information") 
dev.off()

#Hugerth
pdf("./output_data/hugerth/NetCoMi/Gram_stain_network.pdf", paper="a4r", height = 8.27, width =  11.69)
plot_networks(networks.species$hugerth, probs.species$hugerth, hugerth_phyloseq,
              hugerth_bacdive_data, plot_column = "Morphology_cell_morphology_gram_stain",
              title_text = "Hugerth network displaying gram stain information") 
dev.off()

#Liu
pdf("./output_data/liu/NetCoMi/Gram_stain_network.pdf", paper="a4r", height = 8.27, width =  11.69)
plot_networks(networks.species$liu, probs.species$liu, liu_phyloseq,
              liu_bacdive_data, plot_column = "Morphology_cell_morphology_gram_stain",
              title_text = "Liu network displaying gram stain information") 
dev.off()

#Lopresti
pdf("./output_data/lopresti/NetCoMi/Gram_stain_network.pdf", paper="a4r", height = 8.27, width =  11.69)
plot_networks(networks.species$lopresti, probs.species$lopresti, lopresti_phyloseq,
              lopresti_bacdive_data, plot_column = "Morphology_cell_morphology_gram_stain",
              title_text = "Lopresti network displaying gram stain information") 
dev.off()

#Mars
pdf("./output_data/mars/NetCoMi/Gram_stain_network.pdf", paper="a4r", height = 8.27, width =  11.69)
plot_networks(networks.species$mars, probs.species$mars, mars_phyloseq,
              mars_bacdive_data, plot_column = "Morphology_cell_morphology_gram_stain",
              title_text = "Mars network displaying gram stain information") 
dev.off()

#Nagel
pdf("./output_data/nagel/NetCoMi/Gram_stain_network.pdf", paper="a4r", height = 8.27, width =  11.69)
plot_networks(networks.species$nagel, probs.species$nagel, nagel_phyloseq,
              nagel_bacdive_data, plot_column = "Morphology_cell_morphology_gram_stain",
              title_text = "Nagel network displaying gram stain information") 
dev.off()

#Pozuelo
pdf("./output_data/pozuelo/NetCoMi/Gram_stain_network.pdf", paper="a4r", height = 8.27, width =  11.69)
plot_networks(networks.species$pozuelo, probs.species$pozuelo, pozuelo_phyloseq,
              pozuelo_bacdive_data, plot_column = "Morphology_cell_morphology_gram_stain",
              title_text = "Pozuelo network displaying gram stain information") 
dev.off()

#Zeber
pdf("./output_data/zeber/NetCoMi/Gram_stain_network.pdf", paper="a4r", height = 8.27, width =  11.69)
plot_networks(networks.species$zeber, probs.species$zeber, zeber_phyloseq,
              zeber_bacdive_data, plot_column = "Morphology_cell_morphology_gram_stain",
              title_text = "Zeber network displaying gram stain information") 
dev.off()

#Zhu
pdf("./output_data/zhu/NetCoMi/Gram_stain_network.pdf", paper="a4r", height = 8.27, width =  11.69)
plot_networks(networks.species$zhu, probs.species$zhu, zhu_phyloseq,
              zhu_bacdive_data, plot_column = "Morphology_cell_morphology_gram_stain",
              title_text = "Zhu network displaying gram stain information") 
dev.off()

#Zhuang
pdf("./output_data/zhuang/NetCoMi/Gram_stain_network.pdf", paper="a4r", height = 8.27, width =  11.69)
plot_networks(networks.species$zhuang, probs.species$zhuang, zhuang_phyloseq,
              zhuang_bacdive_data, plot_column = "Morphology_cell_morphology_gram_stain",
              title_text = "Zhuang network displaying gram stain information") 
dev.off()

################################################################################

#Plotting antibiotic resistence information

#Agp
pdf("./output_data/agp/NetCoMi/Antibiotic_resistance_network.pdf", paper="a4r", height = 8.27, width =  11.69)
plot_networks(networks.species$agp, probs.species$agp, agp_phyloseq,
              agp_bacdive_data, plot_column = "Physiology_and_metabolism_antibiotic_resistance_is_resistant",
              title_text = "Agp network displaying antibiotic resistance information") 
dev.off()

#Fukui
pdf("./output_data/fukui/NetCoMi/Antibiotic_resistance_network.pdf", paper="a4r", height = 8.27, width =  11.69)
plot_networks(networks.species$fukui, probs.species$fukui, fukui_phyloseq,
              fukui_bacdive_data, plot_column = "Physiology_and_metabolism_antibiotic_resistance_is_resistant",
              title_text = "Fukui network displaying antibiotic resistance information") 
dev.off()

#Hugerth
pdf("./output_data/hugerth/NetCoMi/Antibiotic_resistance_network.pdf", paper="a4r", height = 8.27, width =  11.69)
plot_networks(networks.species$hugerth, probs.species$hugerth, hugerth_phyloseq,
              hugerth_bacdive_data, plot_column = "Physiology_and_metabolism_antibiotic_resistance_is_resistant",
              title_text = "Hugerth network displaying antibiotic resistance information") 
dev.off()

#Liu
pdf("./output_data/liu/NetCoMi/Antibiotic_resistance_network.pdf", paper="a4r", height = 8.27, width =  11.69)
plot_networks(networks.species$liu, probs.species$liu, liu_phyloseq,
              liu_bacdive_data, plot_column = "Physiology_and_metabolism_antibiotic_resistance_is_resistant",
              title_text = "Liu network displaying antibiotic resistance information") 
dev.off()

#Lopresti
pdf("./output_data/lopresti/NetCoMi/Antibiotic_resistance_network.pdf", paper="a4r", height = 8.27, width =  11.69)
plot_networks(networks.species$lopresti, probs.species$lopresti, lopresti_phyloseq,
              lopresti_bacdive_data, plot_column = "Physiology_and_metabolism_antibiotic_resistance_is_resistant",
              title_text = "Lopresti network displaying antibiotic resistance information") 
dev.off()

#Mars
pdf("./output_data/mars/NetCoMi/Antibiotic_resistance_network.pdf", paper="a4r", height = 8.27, width =  11.69)
plot_networks(networks.species$mars, probs.species$mars, mars_phyloseq,
              mars_bacdive_data, plot_column = "Physiology_and_metabolism_antibiotic_resistance_is_resistant",
              title_text = "Mars network displaying antibiotic resistance information") 
dev.off()

#Nagel
pdf("./output_data/nagel/NetCoMi/Antibiotic_resistance_network.pdf", paper="a4r", height = 8.27, width =  11.69)
plot_networks(networks.species$nagel, probs.species$nagel, nagel_phyloseq,
              nagel_bacdive_data, plot_column = "Physiology_and_metabolism_antibiotic_resistance_is_resistant",
              title_text = "Nagel network displaying antibiotic resistance information") 
dev.off()

#Pozuelo
pdf("./output_data/pozuelo/NetCoMi/Antibiotic_resistance_network.pdf", paper="a4r", height = 8.27, width =  11.69)
plot_networks(networks.species$pozuelo, probs.species$pozuelo, pozuelo_phyloseq,
              pozuelo_bacdive_data, plot_column = "Physiology_and_metabolism_antibiotic_resistance_is_resistant",
              title_text = "Pozuelo network displaying antibiotic resistance information") 
dev.off()

#Zeber
pdf("./output_data/zeber/NetCoMi/Antibiotic_resistance_network.pdf", paper="a4r", height = 8.27, width =  11.69)
plot_networks(networks.species$zeber, probs.species$zeber, zeber_phyloseq,
              zeber_bacdive_data, plot_column = "Physiology_and_metabolism_antibiotic_resistance_is_resistant",
              title_text = "Zeber network displaying antibiotic resistance information") 
dev.off()

#Zhu
pdf("./output_data/zhu/NetCoMi/Antibiotic_resistance_network.pdf", paper="a4r", height = 8.27, width =  11.69)
plot_networks(networks.species$zhu, probs.species$zhu, zhu_phyloseq,
              zhu_bacdive_data, plot_column = "Physiology_and_metabolism_antibiotic_resistance_is_resistant",
              title_text = "Zhu network displaying antibiotic resistance information") 
dev.off()

#Zhuang
pdf("./output_data/zhuang/NetCoMi/Antibiotic_resistance_network.pdf", paper="a4r", height = 8.27, width =  11.69)
plot_networks(networks.species$zhuang, probs.species$zhuang, zhuang_phyloseq,
              zhuang_bacdive_data, plot_column = "Physiology_and_metabolism_antibiotic_resistance_is_resistant",
              title_text = "Zhuang network displaying antibiotic resistance information") 
dev.off()

################################################################################




