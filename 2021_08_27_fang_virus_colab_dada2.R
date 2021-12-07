# Title and author information --------------------------------------------
#!/usr/bin/R

#################################################
#                                               #
#       2021_08_27_fang_virus_colab_dada2.R     #
#                                               #
#################################################

#Title: Effects of vaccination on the porcine  gut microbiome
#
#Copyright (C) 2021-2022  Christopher A. Gaulke
#author contact: chris.gaulke@gmail.com
#
#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#
#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#For a copy of the GNU General Public License see
#<http://www.gnu.org/licenses/>.

#Purpose: The purpose of this script is to generate ASV profiles for
# fang virus collaboration


# SET ENVIRONMENT ---------------------------------------------------------

library(dada2)
library(ggplot2)
library(vegan)
library(cowplot)
library(lme4)
library(glmmTMB)
library(qvalue)
library(reshape2)

options("stringsAsFactors" = F)


# FUNCTIONS ---------------------------------------------------------------


###
#        Function phylotype_analysis             #
###

phylotype_analysis <- function(obj, tax) {
  #obj: microbiome object with at least 1 slot (data)
  #tax: a tax object (named list taxa as names values in the list are seq ids)
  obj.out <- NULL
  for (h in 1:length(tax)) {
    df <- NULL
    #print(h)#debugging
    for (i in 1:length(tax[[h]])) {
      print(i)#debugging
      v1       <- obj$data[, unlist(tax[[h]][[i]])]
      v2       <- names(tax[[h]])[i]
      if (is.null(dim(v1))) {
        df[[v2]] <- v1
      } else{
        df[[v2]] <- rowSums(v1)
      }
    }
    obj.out[[names(tax)[h]]] <- as.data.frame(df)
  }
  return(obj.out)
}

make_taxa_df <- function(tax){

  kingdom.df <- replicate(length(unique(tax[, 2])), c())
  names(kingdom.df) <- unique(tax[, 2])
  phylum.df  <- replicate(length(unique(tax[, 3])), c())
  names(phylum.df) <- unique(tax[, 3])
  class.df   <- replicate(length(unique(tax[, 4])), c())
  names(class.df) <- unique(tax[, 4])
  order.df   <- replicate(length(unique(tax[, 5])), c())
  names(order.df) <- unique(tax[, 5])
  family.df  <- replicate(length(unique(tax[, 6])), c())
  names(family.df) <- unique(tax[, 6])
  genus.df   <- replicate(length(unique(tax[, 7])), c())
  names(genus.df) <- unique(tax[, 7])

  for (i in 1:nrow(tax)) {
    kingdom.df[[tax[i, 2]]] <-
      c(kingdom.df[[tax[i, 2]]], tax[i, 1])
    phylum.df[[tax[i, 3]]]  <-
      c(phylum.df[[tax[i, 3]]], tax[i, 1])
    class.df[[tax[i, 4]]]   <-
      c(class.df[[tax[i, 4]]], tax[i, 1])
    order.df[[tax[i, 5]]]   <-
      c(order.df[[tax[i, 5]]], tax[i, 1])
    family.df[[tax[i, 6]]]  <-
      c(family.df[[tax[i, 6]]], tax[i, 1])
    genus.df[[tax[i, 7]]]   <-
      c(genus.df[[tax[i, 7]]], tax[i, 1])
  }

  tax.obj <- NULL
  tax.obj$kingdom <- kingdom.df
  tax.obj$phylum  <- phylum.df
  tax.obj$class   <- class.df
  tax.obj$order   <- order.df
  tax.obj$family  <- family.df
  tax.obj$genus   <- genus.df

  return(tax.obj)
}

# IMPORT DATA -------------------------------------------------------------

path <- "~/unsynced_projects/raw_data/2021_08_13_yuan_pig_vaccine_16S_raw/"
filt.path <- "/Users/cgaulke/Documents/research/fang_virus_colab/data/filtered_data/" #filtered file directory make sure to update
#git ignore to ignore this
list.files(path)

fnFs <- sort(list.files(path, pattern="_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq.gz", full.names = TRUE))


# ANALYSIS: QC ------------------------------------------------------------

#make sure the lengths are the same
length(fnFs) == length(fnRs)

#get sample names, in our case we have to leave some extra on the end
sample.names <- sapply(strsplit(basename(fnFs), "_R1"), `[`, 1)

#preview
plotQualityProfile(fnFs[1])
plotQualityProfile(fnRs[1])

#aggregate all data together now
fnFs_qual.plot <- plotQualityProfile(fnFs,aggregate = T)
fnRs_qual.plot <- plotQualityProfile(fnRs,aggregate = T)

#set up for filtering
filtFs <- file.path(filt.path, "filtered",
                    paste0(sample.names, "_F_filt.fastq.gz"))

filtRs <- file.path(filt.path, "filtered",
                    paste0(sample.names, "_R_filt.fastq.gz"))

#make these named
names(filtFs) <- sample.names
names(filtRs) <- sample.names

#filter and trim
#by looking at the data we see a rapid drop in quality around 200bp. Since the
#average drops below ~25 around 230 we will truncate at 220 for the reverse. The
#forward looks better (this is usual) so we will truncate around 260

#note the original tutorial uses generic variable names

#stopped here

filter.out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(280,250),
                            maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                            compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
#take a look
View(filter.out)

colMeans(filter.out)
mean(1-(filter.out[,2]/filter.out[,1]))

#lets look at these numbers in a little more detail
fivenum(1-(filter.out[,2]/filter.out[,1]))
hist(1-(filter.out[,2]/filter.out[,1]))

# since some libraries look like there is a higher level of filtration
# than others lets take a closer look at this

sort(1-(filter.out[,2]/filter.out[,1]))

#let's keep this in mind moving forward

# ANALYSIS: ERROR ---------------------------------------------------------

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)

dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)


# ANALYSIS: MERGE AND FILTER -----------------------------------------------
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

seqtab <- makeSequenceTable(mergers)
dim(seqtab)

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus",
                                    multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

#how much data was wrapped up in chimeras
sum(seqtab.nochim)/sum(seqtab)


# ANALYSIS: TRACK READS ---------------------------------------------------

getN <- function(x) sum(getUniques(x))
track <- cbind(filter.out, sapply(dadaFs, getN),
               sapply(dadaRs, getN),
               sapply(mergers, getN),
               rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)


# ANALYSIS: ADD TAX -------------------------------------------------------

#***Here is where you will need to go and download the silva databases.
#***Be sure to get the right ones (the names are the same as the ones below)
#***These files can be downloaded here:https://zenodo.org/record/4587955#.YSlzKC1h1hA

taxa <- assignTaxonomy(seqtab.nochim,
                       "/Users/cgaulke/unsynced_projects/db/silva_dada2/silva_nr99_v138.1_train_set.fa",
                       multithread=TRUE)

taxa <- addSpecies(taxa, "/Users/cgaulke/unsynced_projects/db/silva_dada2/silva_species_assignment_v138.1.fa")

taxa.print <- taxa
rownames(taxa.print) <- NULL
head(taxa.print)


# IMPORT: METADATA --------------------------------------------------------
Master_metadata.df <- read.table("/Users/cgaulke/Documents/research/fang_virus_colab/data/metadata.txt",
                                sep = "\t",
                                header = T,
                                row.names = 1)

add_metadata.df <- read.table("/Users/cgaulke/Documents/research/fang_virus_colab/data/additional_metadata.txt",
                              sep = "\t",
                              header = T,
                              row.names = NULL
                              )

add_metadata.df$combined_id <- paste0(add_metadata.df$ID, "_", add_metadata.df$date)

serum_metadata.df <- add_metadata.df[which(add_metadata.df$Sample_type == "Serum"),c(1:5,16:29)]
lung_metadata.df <- add_metadata.df[which(add_metadata.df$Sample_type == "Lung"),c(1:9, 29)]
pbmc_metadata.df <- add_metadata.df[which(add_metadata.df$Sample_type == "PBMC"),c(1:5,10:15,29)]

fecal_add_metadata <- read.table("/Users/cgaulke/Documents/research/fang_virus_colab/data/additionaLfecal_metadata.txt",
                                 sep = "\t",
                                 header = T,
                                 row.names = 1)

fecal_add_metadata <- fecal_add_metadata[which(fecal_add_metadata$date == "45d"),]
fecal_add_metadata <- fecal_add_metadata[1:18,] #remove the all NA samples

fecal_add_metadata_cells <- na.omit(fecal_add_metadata)


# DATA: MAKE SEQ DICTIONARY -----------------------------------------------

asv.dict <- data.frame(names = paste0("asv", seq(1:ncol(seqtab.nochim))),
                       seq   = colnames(seqtab.nochim)
                       )
rownames(asv.dict) <- asv.dict$names

colnames(seqtab.nochim) <- asv.dict$names

#update taxa.print
rownames(taxa.print) <- asv.dict$names

# ANALYSIS: NORMALIZE -----------------------------------------------------

seqtab_nochim.relabd <- sweep(seqtab.nochim,
                              MARGIN = 1,
                              rowSums(seqtab.nochim),
                              FUN = '/'
)


rarecurve(x = seqtab.nochim,
          step = 100,
          label = F
)
abline(v = 1000)
abline(v = 5000) # 5000 seems to hit the curves after the hinge

#get rid of samples with less than 5000 reads and ASVs with 0 abundance
seqtab_nochim.rare <- seqtab.nochim[which(rowSums(seqtab.nochim) > 4999), ]
seqtab_nochim.rare <- seqtab_nochim.rare[,which(colSums(seqtab_nochim.rare) > 0)]

#rarefy
seqtab_nochim.rare <- rrarefy(seqtab_nochim.rare, sample = 5000)

seqtab_nochim.rare <- seqtab_nochim.rare[,which(colSums(seqtab_nochim.rare) > 0)]

rare_metadata.df <- Master_metadata.df[which(rownames(Master_metadata.df) %in% rownames(seqtab_nochim.rare)),]


# ANALYSIS: PHYLOTYPE ANALYSIS --------------------------------------------

#we need to add an ASV identifier for the tax function to work properly
#we also need to restrict this analysis to only ASV IDs present in the data
#the will be used in pig.obj$data (below)

tax <- taxa.print[which(rownames(taxa.print) %in% colnames(seqtab_nochim.rare)),]
tax <- data.frame(asv = rownames(tax), tax)

#make a data frame of taxonomy
pig.tax <- make_taxa_df(tax = tax)

# aggregate phylotype counts for this we will need to provide an pseudoobject with a
# slot called data
pig.obj <- NULL
pig.obj$data <- seqtab_nochim.rare

#classifying taxonomy as also called phyotyping
pig_phylotype <- phylotype_analysis(pig.obj, tax = pig.tax)


# ANALYSIS: PCA ASV ALL-----------------------------------------------------------

#make PCA object
seqtab_nochim.pca <- prcomp(seqtab_nochim.rare,scale. = F, center = T)

#check that rownames are the same
all(rownames(seqtab_nochim.pca$x) == rownames(rare_metadata.df))

#make a PCA data frame
seqtab_nochim_pca.df <- as.data.frame(seqtab_nochim.pca$x[,1:5])

#add metadata for plotting
seqtab_nochim_pca.df$Vaccination_status <- rare_metadata.df$Vaccination_status
seqtab_nochim_pca.df$Sample_type <- rare_metadata.df$Sample_type
seqtab_nochim_pca.df$date <- rare_metadata.df$date
seqtab_nochim_pca.df$group <- rare_metadata.df$group
seqtab_nochim_pca.df$ID <- rare_metadata.df$ID

seqtab_nochim_pca.plot <- ggplot(data = seqtab_nochim_pca.df,
                                 aes(x = PC1,
                                     y = PC2,
                                     color = Sample_type))

pdf("figures/fang_all_samples_PCA.pdf", width = 4, height = 4)
seqtab_nochim_pca.plot +
  geom_point(size = 3, alpha = .7)+
  scale_fill_brewer("",  type = "qual",
                      palette = 2,
                      direction = 1,
                      aesthetics = "colour",
                      labels = c("Fecal", "Nasal", "Placenta")
  )+
  theme(
    text = element_text(size = 16, colour = "black"),
    panel.grid.major = element_line(colour = "grey99"),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(colour = "black"),
    aspect.ratio = 1,
  )+
  xlab("PC1 (25%)")+
  ylab("PC2 (15%)")
dev.off()

# ANALYSIS: ADONIS --------------------------------------------------------

set.seed(731)

seqtab_nochim.adonis <- adonis2(seqtab_nochim.rare ~  Sample_type + Vaccination_status,
                               permutations = 5000, by = "terms",
                               data = rare_metadata.df)

seqtab_nochim.adonis


# ANALYSIS: FECAL ANALYSIS ------------------------------------------------

#Get names for fecal samples
fecal.names <- rownames(rare_metadata.df[which(rare_metadata.df$Sample_type == "FS"), ])

#filter ASV table
fecal_seqtab_nochim.relabd <- seqtab_nochim.relabd[fecal.names,]
fecal_seqtab_nochim.relabd <- fecal_seqtab_nochim.relabd[,which(colSums(fecal_seqtab_nochim.relabd) > 0 )]

fecal_seqtab_nochin.rare <- seqtab_nochim.rare[fecal.names,]
fecal_seqtab_nochin.rare <- fecal_seqtab_nochin.rare[,which(colSums(fecal_seqtab_nochin.rare) > 0 )]

#filter metadata
fecal_master_metadata.df   <- rare_metadata.df[fecal.names,]
fecal_master_metadata.df$vaccine <- ifelse(fecal_master_metadata.df$group == "G3" & fecal_master_metadata.df$date != "0d", 1,0)
fecal_master_metadata.df$virus <- ifelse(fecal_master_metadata.df$group != "G5" & fecal_master_metadata.df$date == "45d", 1,0)
fecal_master_metadata.df$combined <- paste0(fecal_master_metadata.df$virus, "_", fecal_master_metadata.df$vaccine)

#filter genus
fecal_genus.df <- pig_phylotype$genus[fecal.names,]
fecal_genus.df <- fecal_genus.df[,which(colSums(fecal_genus.df) >0)]


# ANALYSIS: FECAL ALPHA DIVERSITY -----------------------------------------

fecal.richness <- specnumber(x = fecal_seqtab_nochin.rare,MARGIN = 1 )
fecal.shannon  <- diversity(fecal_seqtab_nochin.rare,MARGIN = 1)

# The most complicated - and powerful - of all models covered here are those
# built using the glmmTMB package, which can handle everything from the relatively
# simple gaussian models with a single random effect to more complex zero inflated
# hurdle and other advanced modeling. Here we will us gaussian glmms to quantify
# associations between alpha diversity, time, vaccination, and infection


fecal_richness.glmmtmb <- glmmTMB(fecal.richness~ date + vaccine*virus + (1|ID),
                                  data = fecal_master_metadata.df
                                  )
fecal_shannon.glmmtmb <- glmmTMB(fecal.shannon~ date + vaccine*virus + (1|ID),
                                 data = fecal_master_metadata.df
                                 )

#for both models it appears as though

#richness
fecal_richness.df <- data.frame(richness = fecal.richness,
                                group = fecal_master_metadata.df$combined)

fecal_richness.plot <- ggplot(fecal_richness.df, aes(x = group,
                              y = richness,
                              fill = group))+
  geom_boxplot()+
  theme(
    text = element_text(size = 16, colour = "black"),
    panel.grid.major = element_line(colour = "grey99"),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(colour = "black"),
    axis.text.x = element_blank(),
    axis.ticks.x =  element_blank(),
    aspect.ratio = 1,
    plot.title = element_text(size = 12, hjust = 0.5)
  )+
  scale_fill_brewer("",
                    type = "qual",
                    palette = 2,
                    direction = 1,
                    aesthetics = "fill",
                    labels= c("Control", "Vaccine", "Virus", "Virus+Vaccine")
  )+
  xlab("")+
  ylab("Richness")

#shannon
fecal_shannon.df <- data.frame(shannon = fecal.shannon,
                                group = fecal_master_metadata.df$combined)

fecal_shannon.plot <- ggplot(fecal_shannon.df, aes(x = group,
                                                     y = shannon,
                                                     fill = group))+
  geom_boxplot()+
  theme(
    text = element_text(size = 16, colour = "black"),
    panel.grid.major = element_line(colour = "grey99"),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(colour = "black"),
    axis.text.x = element_blank(),
    axis.ticks.x =  element_blank(),
    aspect.ratio = 1,
    plot.title = element_text(size = 12, hjust = 0.5)
  )+
  scale_fill_brewer(
    type = "qual",
    palette = 2,
    direction = 1,
    aesthetics = "fill",
    labels= c("Control", "Vaccine", "Virus", "Virus+Vaccine")
  )+
  xlab("")+
  ylab("Shannon")

legend <- get_legend(
  fecal_richness.plot
)

pdf("figures/fecal_alpha_diversity.pdf",height = 4)
plot_grid(fecal_richness.plot + theme(legend.position = "none"),
          fecal_shannon.plot + theme(legend.position = "none"),
          legend, ncol = 3,rel_widths = c(2,2,1),
          labels = c("A","B"),
          label_size = 20
            )
dev.off()

# ANALYSIS: FECAL PCA ASV-----------------------------------------------------

#make PCA object
fecal_seqtab_nochim.pca <- prcomp(fecal_seqtab_nochin.rare,scale. = F, center = T)

#get portion of variance explained by each axis
summary(fecal_seqtab_nochim.pca)$importance[,1:5]

#plot variances in skree plot
plot(fecal_seqtab_nochim.pca)

#check that rownames are the same
all(rownames(fecal_seqtab_nochim.pca$x) == rownames(fecal_master_metadata.df))

#make a PCA data frame
fecal_seqtab_nochim_pca.df <- as.data.frame(fecal_seqtab_nochim.pca$x[,1:5])

#add metadata for plotting
fecal_seqtab_nochim_pca.df$Vaccination_status <- fecal_master_metadata.df$Vaccination_status
fecal_seqtab_nochim_pca.df$date <- fecal_master_metadata.df$date
fecal_seqtab_nochim_pca.df$group <- fecal_master_metadata.df$group
fecal_seqtab_nochim_pca.df$ID <- fecal_master_metadata.df$ID

fecal_seqtab_nochim_pca.plot <- ggplot(data = fecal_seqtab_nochim_pca.df,
                                 aes(x = PC1,
                                     y = PC2,
                                     color = group,
                                     shape = date))

#ugly duckling example

fecal_seqtab_nochim_pca.plot +
  geom_point(size = 5) #+
  # scale_color_brewer("Group", type = "qual",
  #                                        palette = 2,
  #                                        direction = 1,
  #                                        aesthetics = "colour",
  #                                        labels= c("Vir+Vax", "Vir", "Control")
  # )+
  # stat_ellipse(inherit.aes = FALSE,
  #              type = "t",
  #              level = .3, #since we don't care about the stat, we can lower confidence interval
  #              geom = "polygon",
  #              alpha = .3,
  #              data = fecal_seqtab_nochim_pca.df,
  #              aes(x = PC1,
  #                  y = PC2,
  #                  color = group,
  #                  fill = group)
  # )+
  # scale_fill_brewer(type = "qual",
  #                   palette = 2,
  #                   direction = 1,
  #                   aesthetics = "fill"
  # )+
  # theme(
  #   text = element_text(size = 20, colour = "black"),
  #   panel.grid.major = element_line(colour = "grey99"),
  #   panel.grid.minor = element_blank(),
  #   panel.background = element_blank(),
  #   axis.line = element_line(colour = "black"),
  #   axis.text = element_text(colour = "black"),
  #   aspect.ratio = 1,
  #   plot.title = element_text(hjust = 0.5)
  # )+
  # guides(fill = "none") +ggtitle("All")
  #


# ANALYSIS: FECAL PCA ASV PANEL PLOT --------------------------------------

#all
fecal_seqtab_nochim_pca.plotall <- fecal_seqtab_nochim_pca.plot +
  geom_point(size = 5)+
  scale_color_brewer("Group", type = "qual",
                     palette = 2,
                     direction = 1,
                     aesthetics = "colour",
                     labels= c("Vir+Vax", "Vir", "Control")
  )+
  # stat_ellipse(inherit.aes = FALSE,
  #              type = "t",
  #              level = .3, #since we don't care about the stat, we can lower confidence interval
  #              geom = "polygon",
  #              alpha = .3,
  #              data = fecal_seqtab_nochim_pca.df,
  #              aes(x = PC1,
  #                  y = PC2,
  #                  color = group,
  #                  fill = group)
  # )+
  scale_fill_brewer(type = "qual",
                    palette = 2,
                    direction = 1,
                    aesthetics = "fill"
  ) +
  theme(
    text = element_text(size = 20, colour = "black"),
    panel.grid.major = element_line(colour = "grey99"),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(colour = "black"),
    aspect.ratio = 1,
    plot.title = element_text(hjust = 0.5)
  )+
  guides(fill = "none") +ggtitle("All")


#day 0
fecal_seqtab_nochim_pca.plot0d <- fecal_seqtab_nochim_pca.plot +
  geom_point(data = ~ subset(., date == "0d") ,size = 5)+
  scale_color_brewer("Group",
                     type = "qual",
                     palette = 2,
                     direction = 1,
                     aesthetics = "colour",
                     labels= c("Vir+Vax", "Vir", "Control")
  )+
  stat_ellipse(inherit.aes = FALSE,
                type = "t",
                level = .3, #since we don't care about the stat, we can lower confidence interval
                geom = "polygon",
                alpha = .3,
                data = subset(fecal_seqtab_nochim_pca.df,date == "0d"),
                aes(x = PC1,
                    y = PC2,
                    color = group,
                    fill = group)
  )+
  scale_fill_brewer(type = "qual",
                    palette = 2,
                    direction = 1,
                    aesthetics = "fill"
  ) +
  theme(
    text = element_text(size = 20, colour = "black"),
    panel.grid.major = element_line(colour = "grey99"),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(colour = "black"),
    aspect.ratio = 1,
    plot.title = element_text(hjust = 0.5)
  )+
  guides(shape = "none", fill = "none") +  ggtitle("Day 0")


# day 35
fecal_seqtab_nochim_pca.plot35d <- fecal_seqtab_nochim_pca.plot +
  geom_point(data = ~ subset(., date == "35d") ,size = 5) +
  scale_color_brewer("Group", type = "qual",
                     palette = 2,
                     direction = 1,
                     aesthetics = "colour",
                     labels= c("Vir+Vax", "Vir", "Control")
  ) +
  stat_ellipse(inherit.aes = FALSE,
                  type = "t",
                  level = .3, #since we don't care about the stat, we can lower confidence interval
                  geom = "polygon",
                  alpha = .3,
                  data = subset(fecal_seqtab_nochim_pca.df,date == "35d"),
                  aes(x = PC1,
                      y = PC2,
                      color = group,
                      fill = group)
  )+
  scale_fill_brewer(type = "qual",
                    palette = 2,
                    direction = 1,
                    aesthetics = "fill"
  ) +
  theme(
    text = element_text(size = 20, colour = "black"),
    panel.grid.major = element_line(colour = "grey99"),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(colour = "black"),
    aspect.ratio = 1,
    plot.title = element_text(hjust = 0.5)
  )+
  guides(shape = "none", fill = "none")+ ggtitle("Day 35")



# day 45
fecal_seqtab_nochim_pca.plot45d <- fecal_seqtab_nochim_pca.plot +
  geom_point(data = ~ subset(., date == "45d") ,size = 5) +
  scale_color_brewer("Group", type = "qual",
                     palette = 2,
                     direction = 1,
                     aesthetics = "colour",
                     labels= c("Vir+Vax", "Vir", "Control")
  )+
  stat_ellipse(inherit.aes = FALSE,
               type = "t",
               level = .3, #since we don't care about the stat, we can lower confidence interval
               geom = "polygon",
               alpha = .3,
               data = subset(fecal_seqtab_nochim_pca.df,date == "45d"),
               aes(x = PC1,
                   y = PC2,
                   color = group,
                   fill = group)
  )+
  scale_fill_brewer(type = "qual",
                    palette = 2,
                    direction = 1,
                    aesthetics = "fill"
  ) +
  theme(
    text = element_text(size = 20, colour = "black"),
    panel.grid.major = element_line(colour = "grey99"),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(colour = "black"),
    aspect.ratio = 1,
    plot.title = element_text(hjust = 0.5)
  )+
  guides(shape = "none", fill = "none")+ ggtitle("Day 45")

#This is where we put it all together in one plot
fecal_asv_by_day.plot <- plot_grid(fecal_seqtab_nochim_pca.plotall,
                                   fecal_seqtab_nochim_pca.plot0d,
                                   fecal_seqtab_nochim_pca.plot35d,
                                   fecal_seqtab_nochim_pca.plot45d,
                                   labels = "AUTO", label_size = 30

)

pdf("figures/fecal_asv_pca_plots.pdf", width = 14, height = 14)
fecal_asv_by_day.plot
dev.off()

# ANALYSIS: FECAL ADONIS ASV --------------------------------------------------------

set.seed(731)

fecal_seqtab_nochim.adonis <- adonis(fecal_seqtab_nochim.relabd ~  group * date +ID ,
                               permutations = 5000,
                               data = fecal_master_metadata.df)

fecal_seqtab_nochim.adonis

fecal_seqtab_nochim2.adonis <- adonis(fecal_seqtab_nochim.relabd ~  virus * vaccine + date  ,
                                      permutations = 5000,
                                      data = fecal_master_metadata.df)
fecal_seqtab_nochim2.adonis


#notice that this will completely change the outputs and our interpretation
fecal_seqtab_nochim2.adonis <- adonis(fecal_seqtab_nochim.relabd ~  date + virus * vaccine  ,
                               permutations = 5000,
                               data = fecal_master_metadata.df)
fecal_seqtab_nochim2.adonis


#we can kinda settle this using a slightly different function
x <- adonis2(fecal_seqtab_nochim.relabd ~  date + virus * vaccine  ,
            permutations = 5000,
            data = fecal_master_metadata.df, by = "margin")

y <- adonis2(fecal_seqtab_nochim.relabd ~   virus * vaccine + date ,
             permutations = 5000,
             data = fecal_master_metadata.df, by = "margin")

x
y

# ANALYSIS: FECAL PCA GENUS-----------------------------------------------------

#make PCA object
fecal_genus.pca <- prcomp(fecal_genus.df,scale. = F, center = T)

#get variance %
summary(fecal_genus.pca)

#plot variances in skree plot
plot(fecal_genus.pca)

#check that rownames are the same
all(rownames(fecal_genus.pca$x) == rownames(fecal_master_metadata.df))

#make a PCA data frame
fecal_genus_pca.df <- as.data.frame(fecal_genus.pca$x[,1:5])

#add metadata for plotting
fecal_genus_pca.df$Vaccination_status <- fecal_master_metadata.df$Vaccination_status
fecal_genus_pca.df$date <- fecal_master_metadata.df$date
fecal_genus_pca.df$group <- fecal_master_metadata.df$group
fecal_genus_pca.df$ID <- fecal_master_metadata.df$ID

fecal_genus_pca.plot <- ggplot(data = fecal_genus_pca.df,
                               aes(x = PC1,
                                   y = PC2,
                                   color = group,
                                   shape = date))



# ANALYSIS: FECAL PCA GENUS PANEL PLOTS -----------------------------------

#all
fecal_genus_pca.plotall <- fecal_genus_pca.plot +
  geom_point(size = 5) +
  scale_color_brewer("Group", type = "qual",
                     palette = 2,
                     direction = 1,
                     aesthetics = "colour",
                     labels= c("Vir+Vax", "Vir", "Control")
  )+
  # stat_ellipse(inherit.aes = FALSE,
  #              type = "t",
  #              level = .3, #since we don't care about the stat, we can lower confidence interval
  #              geom = "polygon",
  #              alpha = .3,
  #              data = fecal_genus_pca.df,
  #              aes(x = PC1,
  #                  y = PC2,
  #                  color = group,
  #                  fill = group)
  # )+
  # scale_fill_brewer(type = "qual",
  #                   palette = 2,
  #                   direction = 1,
  #                   aesthetics = "fill"
  # ) +
  theme(
    text = element_text(size = 20, colour = "black"),
    panel.grid.major = element_line(colour = "grey99"),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(colour = "black"),
    aspect.ratio = 1,
    plot.title = element_text(hjust = 0.5)
  )+
  guides(fill = "none") +ggtitle("All")


# start
fecal_genus_pca.plot0d <- fecal_genus_pca.plot +
  geom_point(data = ~ subset(., date == "0d") ,size = 5)+
  scale_color_brewer("Group", type = "qual",
                     palette = 2,
                     direction = 1,
                     aesthetics = "colour",
                     labels= c("Vir+Vax", "Vir", "Control")
  )+
  stat_ellipse(inherit.aes = FALSE,
               type = "t",
               level = .3, #since we don't care about the stat, we can lower confidence interval
               geom = "polygon",
               alpha = .3,
               data = subset(fecal_genus_pca.df,date == "0d"),
               aes(x = PC1,
                   y = PC2,
                   color = group,
                   fill = group)
  )+
  scale_fill_brewer(type = "qual",
                    palette = 2,
                    direction = 1,
                    aesthetics = "fill"
  ) +
  theme(
    text = element_text(size = 20, colour = "black"),
    panel.grid.major = element_line(colour = "grey99"),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(colour = "black"),
    aspect.ratio = 1,
    plot.title = element_text(hjust = 0.5)
  )+
  guides(shape = "none", fill = "none") +  ggtitle("Day 0")

#35 days

fecal_genus_pca.plot35d <- fecal_genus_pca.plot +
  geom_point(data = ~ subset(., date == "35d") ,size = 5) +
  scale_color_brewer("Group", type = "qual",
                     palette = 2,
                     direction = 1,
                     aesthetics = "colour",
                     labels= c("Vir+Vax", "Vir", "Control")
  )+
  stat_ellipse(inherit.aes = FALSE,
               type = "t",
               level = .3, #since we don't care about the stat, we can lower confidence interval
               geom = "polygon",
               alpha = .3,
               data = subset(fecal_genus_pca.df,date == "35d"),
               aes(x = PC1,
                   y = PC2,
                   color = group,
                   fill = group)
  )+
  scale_fill_brewer(type = "qual",
                    palette = 2,
                    direction = 1,
                    aesthetics = "fill"
  ) +
  theme(
    text = element_text(size = 20, colour = "black"),
    panel.grid.major = element_line(colour = "grey99"),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(colour = "black"),
    aspect.ratio = 1,
    plot.title = element_text(hjust = 0.5)
  )+
  guides(shape = "none", fill = "none")+ ggtitle("Day 35")

#45 days

fecal_genus_pca.plot45d <- fecal_genus_pca.plot +
  geom_point(data = ~ subset(., date == "45d") ,size = 5) +
  scale_color_brewer("Group", type = "qual",
                     palette = 2,
                     direction = 1,
                     aesthetics = "colour",
                     labels= c("Vir+Vax", "Vir", "Control")
  )+
  stat_ellipse(inherit.aes = FALSE,
               type = "t",
               level = .3, #since we don't care about the stat, we can lower confidence interval
               geom = "polygon",
               alpha = .3,
               data = subset(fecal_genus_pca.df,date == "45d"),
               aes(x = PC1,
                   y = PC2,
                   color = group,
                   fill = group)
               )+
  scale_fill_brewer(type = "qual",
                     palette = 2,
                     direction = 1,
                     aesthetics = "fill"
  ) +
  theme(
    text = element_text(size = 20, colour = "black"),
    panel.grid.major = element_line(colour = "grey99"),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(colour = "black"),
    aspect.ratio = 1,
    plot.title = element_text(hjust = 0.5)
  )+
  guides(shape = "none", fill = "none")+ ggtitle("Day 45")

# If you look at the structure it appears like the vaccinated and the the mock infected
# cluster together while the shame vaccinated:infected cluster more distinctly at 45d post
# infection. This makes sense but wont be detectable by the approach I have used here... or at least not easily detectable

fecal_genus_by_day.plot <- plot_grid(fecal_genus_pca.plotall,
                                     fecal_genus_pca.plot0d,
                                     fecal_genus_pca.plot35d,
                                     fecal_genus_pca.plot45d,
                                     labels = "AUTO",
                                     label_size = 30
                                     )


pdf("figures/fecal_genus_pca_plots.pdf", width = 14, height = 14)

fecal_genus_by_day.plot

dev.off()

# something strange is going on here... Why are only G5 and G3 moving at day 45?
# Could the metadata be off, or was something done to these animals that is different
# from G4 ?


# ANALYSIS: FECAL ADONIS GENUS --------------------------------------------------------

set.seed(731)

fecal_genus.adonis <- adonis(fecal_genus.df ~ date *group + ID ,
                             permutations = 5000,
                             data = fecal_master_metadata.df)

fecal_genus.adonis

#spliced differently
fecal_genus2.adonis <- adonis2(fecal_genus.df ~ date + virus* vaccine  ,
                              permutations = 5000,
                              data = fecal_master_metadata.df)
fecal_genus2.adonis

# here there appears to be a significant effect of data on the animals. This
# could be due to a common environmental parameter being changed in the animals
# at these time points. What might have been changing?

# another interesting point is that G3 and G5 cluster more together while G4 is
# separate


# ANALYSIS: FECAL GENUS STATS -------------------------------------------

# In the interest of data reduction we will conduct the statistical analyses on
# genus level taxonomic classification.

# For teaching purposes I will go over various test moving from most simple and
# less appropriate to move complicated but more appropriate.

#basic linear model
fit0 <- lm(fecal_genus.df$Clostridium.sensu.stricto.1 ~ date + virus* vaccine,data = fecal_master_metadata.df )
summary(fit0)

#but our data is not normally distributed

hist(fecal_genus.df$Clostridium.sensu.stricto.1)

#generalized linear models are sometimes used when normality cannot be assumed
fit1 <- glm(fecal_genus.df$Clostridium.sensu.stricto.1 ~ date + virus* vaccine,data = fecal_master_metadata.df )
summary(fit1)

#But we know from the lit and prior experience that these data are best fit by a neg bionomial distribution
fit2 <- glm.nb(fecal_genus.df$Clostridium.sensu.stricto.1 ~ date + virus* vaccine,data = fecal_master_metadata.df )
summary(fit2)

# But still we have the other issue of have repeated sampling of individuals across time
# to account for this variation we can use a mixed model which accounts for the
# variation in response of each individual

fit3 <- glmmTMB(fecal_genus.df[,1] ~ date + virus * vaccine + (1|ID),family=nbinom1 ,data = fecal_master_metadata.df )
summary(fit3)

# But now we have a problem because we only have the fixed effects p-value and
# not the overall model. This presents an issue with FDR control because family-wise
# error has already been accounted for. To solve this problem we can use a NULL
# model which only includes an implicit intercept and the random (mixed) effect

fit3.null <- glmmTMB(fecal_genus.df[,1] ~ factor(date) +  (1|ID), family=nbinom1,data = fecal_master_metadata.df )
summary(fit3.null)

#if the anova is significant adding those extra terms allows significantly more
# variation to be explained by the model
anova(fit3.null, fit3)

# it should be noted that here, while the null model above is appropriate, a null
# that includes date would also be OK (and potentially better). Regardless now the
# anova will be the thing we will correct on

#so to do all the stats we can use a loop and we will cover that later

fecal_nbglm <- NULL

for(i in 1:ncol(fecal_genus.df)){
#apply(fecal_genus.df,2,function(i){

  index  <- names(fecal_genus.df)[[i]]
  counts <- fecal_genus.df[[i]]

  #Null model
  fecal_nbglm[[index]][["null"]] <- NULL
  fecal_nbglm[[index]][["ALT"]]  <- NULL
  fecal_nbglm[[index]][["ANOVA"]] <- NULL

  try(
  fecal_nbglm[[index]][["null"]] <- glmmTMB(counts ~ date + (1|ID) ,
                                             family=nbinom1,
                                             data = fecal_master_metadata.df )
  )
  #Alternative model
  try(
  fecal_nbglm[[index]][["ALT"]] <- glmmTMB(counts ~ date + virus * vaccine + (1|ID) ,
                                            family=nbinom1,
                                            data = fecal_master_metadata.df )
  )

  # ANOVA test
  try(
  fecal_nbglm[[index]][["ANOVA"]] <- anova(fecal_nbglm[[index]][["null"]],
                                           fecal_nbglm[[index]][["ALT"]])
  )

}

# Note that many of these test will likely run into convergence errors which will
# be caught with try() and a NULL value will be passed to the ANOVA element of the
# Taxa list which will in turn result in NA for ANOVA p value.


#Build a data frame with the false discovery rate controlled values

fecal_nbglm.df <- data.frame(aov.pval= unlist(lapply(fecal_nbglm, function(x){x$ANOVA$`Pr(>Chisq)`[2]})),
                             qval= qvalue::qvalue(unlist(lapply(fecal_nbglm, function(x){x$ANOVA$`Pr(>Chisq)`[2]})))$qvalue )

sig_fecal_genera.names <- rownames(fecal_nbglm.df[which(fecal_nbglm.df$qval < 0.2),])

#make some plots from these taxa
sig_fecal_genera.df <- fecal_genus.df[,sig_fecal_genera.names]

sig_fecal_genera.df$status <- paste0(fecal_master_metadata.df$virus, "_", fecal_master_metadata.df$vaccine)

plots <- NULL
my_legend <- NULL
for(i in seq_along(sig_fecal_genera.names)){

  g <- ggplot(sig_fecal_genera.df, aes_string(x = "status",
                                y = sig_fecal_genera.names[i]
                                )) +
  geom_boxplot(aes(fill = factor(status))) +
  scale_fill_brewer("Group",
                    type = "qual",
                    palette = 2,
                    direction = 1,
                    aesthetics = "fill",
                    labels= c("Control", "Vaccine", "Virus", "Virus+Vaccine")
                    )+
  theme(
    text = element_text(size = 12, colour = "black"),
    panel.grid.major = element_line(colour = "grey99"),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(colour = "black"),
    axis.text.x = element_blank(),
    axis.ticks.x =  element_blank(),
    aspect.ratio = 1,
    #plot.title = element_text(size = 12, hjust = 0.5)
  ) +
  xlab("") #+
  #ylab("")+
  #ggtitle(sig_fecal_genera.names[i])

  if(i == length(sig_fecal_genera.names)){
    my_legend <- get_legend(g)
  }
  plots[[i]] <- g + theme(legend.position = "none")

}

pdf("figures/sig_genera_feces.pdf", width = 10, height = 10)
plot_grid(plots[[1]],
          plots[[2]],
          plots[[3]],
          plots[[4]],
          plots[[5]],
          plots[[6]],
          plots[[7]],
          plots[[8]],
          legend)
dev.off()


# ANALYSIS: FECAL HOST CORRS ----------------------------------------------

fecal_corrs.genus <- fecal_genus.df[which(rownames(fecal_genus.df) %in% rownames(fecal_add_metadata)),]
#fecal_corrs.genus <- fecal_corrs.genus[,which(colSums(fecal_corrs.genus) > 0)]
fecal_corrs.genus <- fecal_corrs.genus[,apply(fecal_corrs.genus,2, function(x){sum(x>0) > length(x)*.33})]

fecal.pvals <- NULL
fecal.est <- NULL
fecal.stat <- NULL
fecal.host <- NULL
fecal.taxa <- NULL

#just be aware that this is an inefficient way of running this loop but with
# such limited data it is fine.

for(i in 8:(ncol(fecal_add_metadata))){
  for(j in colnames(fecal_corrs.genus)){
    pval <- cor.test(fecal_add_metadata[,i], fecal_corrs.genus[,j], method = "s")$p.value
    est  <- cor.test(fecal_add_metadata[,i], fecal_corrs.genus[,j], method = "s")$estimate
    stat <- cor.test(fecal_add_metadata[,i], fecal_corrs.genus[,j], method = "s")$statistic
    host <- colnames(fecal_add_metadata)[i]
    taxa <- j
    fecal.pvals <- c(fecal.pvals, pval)
    fecal.est   <- c(fecal.est, est)
    fecal.stat  <- c(fecal.stat, stat)
    fecal.host  <- c(fecal.host, host)
    fecal.taxa  <- c(fecal.taxa,taxa)
  }
}

fecal_host.corrs <- data.frame(
  pvals = fecal.pvals,
  est  = fecal.est,
  stat = fecal.stat,
  host = fecal.host,
  taxa = fecal.taxa,
  qval = qvalue::qvalue(fecal.pvals)$qvalues
)

# nothing makes qvalue cutoffs, but there are some interesting correlations
# to potentially keep an eye on. For example, Turicibacter, Ruminococcus, and
# Prevotella_9

# ANALYSIS: NASAL ANALYSIS ------------------------------------------------

# Get names for nasal samples use rarefied names to be consistent between normalized
# data frames

nasal.names <- rownames(rare_metadata.df[which(rare_metadata.df$Sample_type == "NS"), ])

#make relative abundance data frame
nasal_seqtab_nochim.relabd <- seqtab_nochim.relabd[nasal.names,]
nasal_seqtab_nochim.relabd <- nasal_seqtab_nochim.relabd[,which(colSums(nasal_seqtab_nochim.relabd) > 0 )]

# make rarefied data frame

nasal_seqtab_nochim.rare <- seqtab_nochim.rare[nasal.names,]
nasal_seqtab_nochim.rare <- nasal_seqtab_nochim.rare[,which(colSums(nasal_seqtab_nochim.rare) > 0 )]

# make genus data frame

#filter genus
nasal_genus.df <- pig_phylotype$genus[nasal.names,]
nasal_genus.df <- nasal_genus.df[,which(colSums(nasal_genus.df) >0)]

#update metadata
nasal_master_metadata.df <- rare_metadata.df[nasal.names, ]

nasal_master_metadata.df$vaccine <-
  ifelse(nasal_master_metadata.df$group %in% c("GC", "GB"), 1,0)
nasal_master_metadata.df$virus <-
  ifelse(nasal_master_metadata.df$group %in% c("GC", "GD") & nasal_master_metadata.df$date %in% c("5d", "22d"), 1,0)
nasal_master_metadata.df$combined <- paste0(nasal_master_metadata.df$virus, "_", nasal_master_metadata.df$vaccine)

nasal.5d <- rownames(nasal_master_metadata.df[which(nasal_master_metadata.df$date == "5d"), ])
nasal.22d <- rownames(nasal_master_metadata.df[which(nasal_master_metadata.df$date == "22d"), ])

nasal_ids.5d  <- paste0(nasal_master_metadata.df[which(nasal_master_metadata.df$date == "5d"),  "ID"], "_5d")
nasal_ids.22d <- paste0(nasal_master_metadata.df[which(nasal_master_metadata.df$date == "22d"), "ID"], "_22d")

# ANALYSIS: NASAL ALPHA-DIVERSITY --------------------------------------------

nasal.richness <- specnumber(x = nasal_seqtab_nochim.rare,MARGIN = 1 )
nasal.shannon  <- diversity(nasal_seqtab_nochim.rare,MARGIN = 1)

# Here we will us gaussian glmms to quantify associations between alpha
# diversity, time, vaccination, and infection. It is likely that we will be under
# powered so trends in the data will be more important to look at


nasal_richness.glmmtmb <- glmmTMB(nasal.richness~ date + vaccine*virus + (1|ID),
                                  data = nasal_master_metadata.df
)

nasal_shannon.glmmtmb <- glmmTMB(nasal.shannon~ date + vaccine*virus + (1|ID),
                                 data = nasal_master_metadata.df
)

summary(nasal_richness.glmmtmb)
summary(nasal_shannon.glmmtmb)


#richness
nasal_richness.df <- data.frame(richness = nasal.richness,
                                group = nasal_master_metadata.df$combined)

nasal_richness.plot <- ggplot(nasal_richness.df, aes(x = group,
                                                     y = richness,
                                                     fill = group))+
  geom_boxplot()+
  theme(
    text = element_text(size = 16, colour = "black"),
    panel.grid.major = element_line(colour = "grey99"),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(colour = "black"),
    axis.text.x = element_blank(),
    axis.ticks.x =  element_blank(),
    aspect.ratio = 1,
    plot.title = element_text(size = 12, hjust = 0.5)
  )+
  scale_fill_brewer("",
                    type = "qual",
                    palette = 2,
                    direction = 1,
                    aesthetics = "fill",
                    labels= c("Control", "Vaccine", "Virus", "Virus+Vaccine")
  )+
  xlab("")+
  ylab("Richness")

#shannon
nasal_shannon.df <- data.frame(shannon = nasal.shannon,
                               group = nasal_master_metadata.df$combined)

nasal_shannon.plot <- ggplot(nasal_shannon.df, aes(x = group,
                                                   y = shannon,
                                                   fill = group))+
  geom_boxplot()+
  theme(
    text = element_text(size = 16, colour = "black"),
    panel.grid.major = element_line(colour = "grey99"),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(colour = "black"),
    axis.text.x = element_blank(),
    axis.ticks.x =  element_blank(),
    aspect.ratio = 1,
    plot.title = element_text(size = 12, hjust = 0.5)
  )+
  scale_fill_brewer(
    type = "qual",
    palette = 2,
    direction = 1,
    aesthetics = "fill",
    labels= c("Control", "Vaccine", "Virus", "Virus+Vaccine")
  )+
  xlab("")+
  ylab("Shannon")

legend <- get_legend(
  nasal_richness.plot
)

pdf("figures/nasal_alpha_diversity.pdf",height = 4, width = 8)
plot_grid(nasal_richness.plot + theme(legend.position = "none"),
          nasal_shannon.plot + theme(legend.position = "none"),
          legend, ncol = 3,rel_widths = c(2,2,1),
          labels = c("A","B"),
          label_size = 20
)
dev.off()


# ANALYSIS: NASAL PCA PANEL PLOTS --------------------------------------------

#make PCA object
nasal_seqtab_nochim.pca <- prcomp(nasal_seqtab_nochim.rare,scale. = F, center = T)

#get variance %
summary(nasal_seqtab_nochim.pca)

#plot variances in skree plot
plot(nasal_seqtab_nochim.pca)

#check that rownames are the same
all(rownames(nasal_seqtab_nochim.pca$x) == rownames(nasal_master_metadata.df))

#make a PCA data frame
nasal_seqtab_nochim_pca.df <- as.data.frame(nasal_seqtab_nochim.pca$x[,1:5])

#add metadata for plotting
nasal_seqtab_nochim_pca.df$status <- nasal_master_metadata.df$combined
nasal_seqtab_nochim_pca.df$date <- nasal_master_metadata.df$date
nasal_seqtab_nochim_pca.df$group <- nasal_master_metadata.df$group
nasal_seqtab_nochim_pca.df$ID <- nasal_master_metadata.df$ID

nasal_seqtab_nochim_pca.plot <- ggplot(data = nasal_seqtab_nochim_pca.df,
                                       aes(x = PC2,
                                           y = PC3,
                                           color = status,
                                           shape = date
                                           )
                                           )
#all
nasal_seqtab_nochim_pca.plotall <- nasal_seqtab_nochim_pca.plot +
  geom_point(size = 5)+
  scale_color_brewer("Group", type = "qual",
                   palette = 2,
                   direction = 1,
                   aesthetics = "colour",
                   labels= c("Control", "Vaccine", "Virus", "Virus+Vaccine")
  )+
  scale_fill_brewer(type = "qual",
                  palette = 2,
                  direction = 1,
                  aesthetics = "fill"
  ) +
  theme(
    text = element_text(size = 20, colour = "black"),
    panel.grid.major = element_line(colour = "grey99"),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(colour = "black"),
    aspect.ratio = 1,
    plot.title = element_text(hjust = 0.5)
  )+
  guides(fill = "none") +ggtitle("All")


#0d
nasal_seqtab_nochim_pca.plot0d <- nasal_seqtab_nochim_pca.plot +
  geom_point(data = ~ subset(., date == "0d") ,size = 5)+
  scale_color_brewer("Group", type = "qual",
                     palette = 2,
                     direction = 1,
                     aesthetics = "colour",
                     labels= c("Control", "Vaccine")
  )+
  scale_fill_brewer(type = "qual",
                    palette = 2,
                    direction = 1,
                    aesthetics = "fill"
  ) +
  theme(
    text = element_text(size = 20, colour = "black"),
    panel.grid.major = element_line(colour = "grey99"),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(colour = "black"),
    aspect.ratio = 1,
    plot.title = element_text(hjust = 0.5)
  )+
  guides(fill = "none", shape = "none") +ggtitle("Day 0")

#Day 5
nasal_seqtab_nochim_pca.plot5d <-
  nasal_seqtab_nochim_pca.plot +
  geom_point(data = ~ subset(., date == "5d") ,size = 5)+
  scale_color_brewer("Group", type = "qual",
                     palette = 2,
                     direction = 1,
                     aesthetics = "colour",
                     labels= c("Control", "Vaccine","Virus", "Virus+Vaccine")
  )+
  scale_fill_brewer(type = "qual",
                    palette = 2,
                    direction = 1,
                    aesthetics = "fill"
  ) +
  theme(
    text = element_text(size = 20, colour = "black"),
    panel.grid.major = element_line(colour = "grey99"),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(colour = "black"),
    aspect.ratio = 1,
    plot.title = element_text(hjust = 0.5)
  )+
  guides(fill = "none", shape = "none") +ggtitle("Day 5")

#Day 22
nasal_seqtab_nochim_pca.plot22d <- nasal_seqtab_nochim_pca.plot +
  geom_point(data = ~ subset(., date == "22d") ,size = 5)+
  scale_color_brewer("Group", type = "qual",
                     palette = 2,
                     direction = 1,
                     aesthetics = "colour",
                     labels= c("Control", "Vaccine","Virus", "Virus+Vaccine")
  )+
  scale_fill_brewer(type = "qual",
                    palette = 2,
                    direction = 1,
                    aesthetics = "fill"
  ) +
  theme(
    text = element_text(size = 20, colour = "black"),
    panel.grid.major = element_line(colour = "grey99"),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(colour = "black"),
    aspect.ratio = 1,
    plot.title = element_text(hjust = 0.5)
  )+
  guides(fill = "none", shape = "none") +ggtitle("Day 22")


nasal.legend <- get_legend(nasal_seqtab_nochim_pca.plotall)


nasal_asv_by_day.plot <-
plot_grid(nasal_seqtab_nochim_pca.plotall + theme(legend.position = "none"),
          nasal_seqtab_nochim_pca.plot0d  + theme(legend.position = "none"),
          nasal_seqtab_nochim_pca.plot5d  + theme(legend.position = "none"),
          nasal_seqtab_nochim_pca.plot22d + theme(legend.position = "none"),
          nasal.legend,
          labels = c("A", "B", "C", "D", ''),
          label_size = 30, nrow = 2
)


pdf("figures/nasal_asv_pca_plots.pdf", width = 14, height = 14)
nasal_asv_by_day.plot
dev.off()

# ANALYSIS: NASAL ADONIS --------------------------------------------------------

set.seed(731)

nasal_seqtab_nochim.adonis <- adonis2(nasal_seqtab_nochim.rare ~  group * date ,
                                     permutations = 5000,
                                     by = "margin",
                                     data = nasal_master_metadata.df)

nasal_seqtab_nochim.adonis

nasal_seqtab_nochim.adonis2 <- adonis2(nasal_seqtab_nochim.rare ~  date + virus*vaccine  ,
                                     permutations = 5000,
                                     by = "margin",
                                     data = nasal_master_metadata.df)
nasal_seqtab_nochim.adonis2


# ANALYSIS: NASAL MICROBIOME:PBMC CORRELATIONS ---------------------------

# for this analysis it would be best to have a reduced set of taxa to correlate
# with these metadata because many of these taxa might have very low prevalence

nasal_genus_corr.df <- nasal_genus.df[,which(colMeans(nasal_genus.df) > 10)]

#serum_metadata.df
#lung_metadata.df
#pbmc_metadata.df

#nasal_ids.22d nasal_ids.5d nasal.22d nasal.5d

#filter and reorder
pbmc_metadata.df <- pbmc_metadata.df[which(pbmc_metadata.df$combined_id %in% nasal_ids.22d),]
rownames(pbmc_metadata.df) <- pbmc_metadata.df$combined_id
pbmc_metadata.df <- pbmc_metadata.df[nasal_ids.22d,]

pbmc_genus.df <- nasal_genus_corr.df[which(rownames(nasal_genus_corr.df) %in% nasal.22d),]
pbmc_genus.df <- pbmc_genus.df[,which(colSums(pbmc_genus.df) > 0)]

#any sig corrs ?
pbmc.pvals <- NULL

for(i in 6:(ncol(pbmc_metadata.df)-1)){
  for(j in colnames(pbmc_genus.df)){
    #x <- cor.test(pbmc_metadata.df[,i], pbmc_genus.df[,j], method = "s")$p.value
    x <- cor.test(pbmc_metadata.df[,i], pbmc_genus.df[,j], method = "p")$p.value
    pbmc.pvals <- c(pbmc.pvals, x)
  }
}

hist(qvalue::qvalue(pbmc.pvals)$qvalues)


# ANALYSIS: NASAL MICROBIOME:LUNG CORRELATIONS ---------------------------

lung_metadata.df <- lung_metadata.df[which(lung_metadata.df$combined_id %in% nasal_ids.5d),]
rownames(lung_metadata.df) <- lung_metadata.df$combined_id
lung_metadata.df <- lung_metadata.df[nasal_ids.5d,]

lung_genus.df <- nasal_genus_corr.df[which(rownames(nasal_genus_corr.df) %in% nasal.5d),]
lung_genus.df <- lung_genus.df[,which(colSums(lung_genus.df) > 0)]

#any sig corrs ?

#make a series of vectors that we can later bind together
lung.pvals <- NULL
lung.est <- NULL
lung.stat <- NULL
lung.host <- NULL
lung.taxa <- NULL

for(i in 6:(ncol(lung_metadata.df)-1)){
  for(j in colnames(lung_genus.df)){
    pval <- cor.test(lung_metadata.df[,i], lung_genus.df[,j], method = "s")$p.value
    est  <- cor.test(lung_metadata.df[,i], lung_genus.df[,j], method = "s")$estimate
    stat <- cor.test(lung_metadata.df[,i], lung_genus.df[,j], method = "s")$statistic
    host <- colnames(lung_metadata.df)[i]
    taxa <- j
    lung.pvals <- c(lung.pvals, pval)
    lung.est   <- c(lung.est, est)
    lung.stat  <- c(lung.stat, stat)
    lung.host  <- c(lung.host, host)
    lung.taxa  <- c(lung.taxa,taxa)
  }
}


lung.corrs <- data.frame(
  pvals = lung.pvals,
  est  = lung.est,
  stat = lung.stat,
  host = lung.host,
  taxa = lung.taxa,
  qval = qvalue::qvalue(lung.pvals)$qvalues
)

# There are a few correlations that pass fdr control of 0.2, which is impressive
# because we had so few samples

# ANALYSIS: NASAL MICROBIOME:SERUM CORRELATIONS ---------------------------

serum_metadata.df <- serum_metadata.df[which(serum_metadata.df$combined_id %in% c(nasal_ids.5d, nasal_ids.22d)),]
rownames(serum_metadata.df) <- serum_metadata.df$combined_id
serum_metadata.df <- serum_metadata.df[c(nasal_ids.5d, nasal_ids.22d),]

serum_genus.df <- nasal_genus_corr.df[which(rownames(nasal_genus_corr.df) %in% c(nasal.5d,nasal.22d)),]
serum_genus.df <- serum_genus.df[c(nasal.5d,nasal.22d), ]
serum_genus.df <- serum_genus.df[,which(colSums(serum_genus.df) > 0)]


serum_metadata.df5d <- serum_metadata.df[1:6,]
serum_genus.df5d <- serum_genus.df[1:6,]

serum_metadata.df22d <- serum_metadata.df[7:14,]
serum_genus.df22d <- serum_genus.df[7:14,]

#any sig corrs ?

#all
serum.pvals <- NULL
serum.est <- NULL
serum.stat <- NULL
serum.host <- NULL
serum.taxa <- NULL

for(i in 6:(ncol(serum_metadata.df)-1)){
  for(j in colnames(serum_genus.df)){
    pval <- cor.test(serum_metadata.df[,i], serum_genus.df[,j], method = "s")$p.value
    est  <- cor.test(serum_metadata.df[,i], serum_genus.df[,j], method = "s")$estimate
    stat <- cor.test(serum_metadata.df[,i], serum_genus.df[,j], method = "s")$statistic
    host <- colnames(serum_metadata.df)[i]
    taxa <- j
    serum.pvals <- c(serum.pvals, pval)
    serum.est   <- c(serum.est, est)
    serum.stat  <- c(serum.stat, stat)
    serum.host  <- c(serum.host, host)
    serum.taxa  <- c(serum.taxa,taxa)
  }
}

serum.corrs <- data.frame(
  pvals = serum.pvals,
  est  = serum.est,
  stat = serum.stat,
  host = serum.host,
  taxa = serum.taxa,
  qval = qvalue::qvalue(serum.pvals)$qvalues
)

#nothing here but what if we split?

#5d
#make a series of vectors that we can later bind together
serum.pvals <- NULL
serum.est <- NULL
serum.stat <- NULL
serum.host <- NULL
serum.taxa <- NULL

for(i in 6:(ncol(serum_metadata.df5d)-1)){
  for(j in colnames(serum_genus.df5d)){
    pval <- cor.test(serum_metadata.df5d[,i], serum_genus.df5d[,j], method = "s")$p.value
    est  <- cor.test(serum_metadata.df5d[,i], serum_genus.df5d[,j], method = "s")$estimate
    stat <- cor.test(serum_metadata.df5d[,i], serum_genus.df5d[,j], method = "s")$statistic
    host <- colnames(serum_metadata.df5d)[i]
    taxa <- j
    serum.pvals <- c(serum.pvals, pval)
    serum.est   <- c(serum.est, est)
    serum.stat  <- c(serum.stat, stat)
    serum.host  <- c(serum.host, host)
    serum.taxa  <- c(serum.taxa,taxa)
  }
}

serum5d.corrs <- data.frame(
  pvals = serum.pvals,
  est  = serum.est,
  stat = serum.stat,
  host = serum.host,
  taxa = serum.taxa,
  qval = qvalue::qvalue(serum.pvals)$qvalues
)


#now 22d

serum.pvals <- NULL
serum.est <- NULL
serum.stat <- NULL
serum.host <- NULL
serum.taxa <- NULL

for(i in 6:(ncol(serum_metadata.df22d)-1)){
  for(j in colnames(serum_genus.df22d)){
    pval <- cor.test(serum_metadata.df22d[,i], serum_genus.df22d[,j], method = "s")$p.value
    est  <- cor.test(serum_metadata.df22d[,i], serum_genus.df22d[,j], method = "s")$estimate
    stat <- cor.test(serum_metadata.df22d[,i], serum_genus.df22d[,j], method = "s")$statistic
    host <- colnames(serum_metadata.df22d)[i]
    taxa <- j
    serum.pvals <- c(serum.pvals, pval)
    serum.est   <- c(serum.est, est)
    serum.stat  <- c(serum.stat, stat)
    serum.host  <- c(serum.host, host)
    serum.taxa  <- c(serum.taxa,taxa)
  }
}

serum22d.corrs <- data.frame(
  pvals = serum.pvals,
  est  = serum.est,
  stat = serum.stat,
  host = serum.host,
  taxa = serum.taxa,
  qval = qvalue::qvalue(serum.pvals)$qvalues
)

#probably not powered for this, but there are some interesting associations here
#a larger sampling would potentially


# ANALYSIS: NASAL MICROBIOME CORR PLOTS -----------------------------------

lacto_infa.plot <- ggplot(data = data.frame(x1 = lung_metadata.df$ifnapcr,
                                            y1 = lung_genus.df$Lactobacillus),
                          aes(x = x1, y = y1)
)+
  geom_point(size = 3)+
  geom_smooth(method = "lm",se = FALSE, color = "black")+
  annotate("text", x=4, y=100, label= "-0.93", size =10)+
  xlab("IFNa (pcr)")+
  ylab("Lactobacillus")+
  theme(
    text = element_text(size = 16, colour = "black"),
    panel.grid.major = element_line(colour = "grey99"),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(colour = "black"),
    aspect.ratio = 1,
  )

# terrisporobacter
terrisporobacter_infa.plot <- ggplot(data = data.frame(x1 = lung_metadata.df$ifnapcr,
                                            y1 = lung_genus.df$Terrisporobacter),
                          aes(x = x1, y = y1)
)+
  geom_point(size = 3)+
  geom_smooth(method = "lm",se = FALSE, color = "black")+
  annotate("text", x=4, y=300, label= "-0.94", size =10)+
  xlab("IFNa (pcr)")+
  ylab("Terrisporobacter")+
  theme(
    text = element_text(size = 16, colour = "black"),
    panel.grid.major = element_line(colour = "grey99"),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(colour = "black"),
    aspect.ratio = 1,
  )


#Turicibacter

turicibacter_infa.plot <- ggplot(data = data.frame(x1 = lung_metadata.df$ifnapcr,
                                                       y1 = lung_genus.df$Turicibacter),
                                     aes(x = x1, y = y1)
)+
  geom_point(size = 3)+
  geom_smooth(method = "lm",se = FALSE, color = "black")+
  annotate("text", x=4, y=300, label= "-0.94", size =10)+
  xlab("IFNa (pcr)")+
  ylab("Turicibacter")+
  theme(
    text = element_text(size = 16, colour = "black"),
    panel.grid.major = element_line(colour = "grey99"),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(colour = "black"),
    aspect.ratio = 1,
  )

#Romboutsia

romboutsia_infa.plot <- ggplot(data = data.frame(x1 = lung_metadata.df$ifnapcr,
                                                   y1 = lung_genus.df$Romboutsia),
                                 aes(x = x1, y = y1)
)+
  geom_point(size = 3)+
  geom_smooth(method = "lm",se = FALSE, color = "black")+
  annotate("text", x=4, y=300, label= "-0.94", size =10)+
  xlab("IFNa (pcr)")+
  ylab("Romboutsia")+
  theme(
    text = element_text(size = 16, colour = "black"),
    panel.grid.major = element_line(colour = "grey99"),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(colour = "black"),
    aspect.ratio = 1,
  )

#Cellulosilyticum


cellulosilyticum_infa.plot <- ggplot(data = data.frame(x1 = lung_metadata.df$ifnapcr,
                                                 y1 = lung_genus.df$Cellulosilyticum),
                               aes(x = x1, y = y1)
)+
  geom_point(size = 3)+
  geom_smooth(method = "lm",se = FALSE, color = "black")+
  annotate("text", x=4, y=150, label= "-0.93", size =10)+
  xlab("IFNa (pcr)")+
  ylab("Cellulosilyticum")+
  theme(
    text = element_text(size = 16, colour = "black"),
    panel.grid.major = element_line(colour = "grey99"),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(colour = "black"),
    aspect.ratio = 1,
  )




pdf("figures/nasal_genus_host_corr_plots.pdf", width = 14, height = 14)
plot_grid(lacto_infa.plot,
          terrisporobacter_infa.plot,
          turicibacter_infa.plot,
          romboutsia_infa.plot,
          cellulosilyticum_infa.plot,
          labels = c("A", "B", "C", "D", 'E'),
          label_size = 30, nrow = 2)
dev.off()



# ANALYSIS: PLACENTA ANALYSIS ---------------------------------------------

# Get names for placenta samples use rarefied names to be consistent between normalized
# data frames

placenta.names <- rownames(rare_metadata.df[which(rare_metadata.df$Sample_type == "Pla"), ])

#make relative abundance data frame
placenta_seqtab_nochim.relabd <- seqtab_nochim.relabd[placenta.names,]
placenta_seqtab_nochim.relabd <- placenta_seqtab_nochim.relabd[,which(colSums(placenta_seqtab_nochim.relabd) > 0 )]

# make rarefied data frame

placenta_seqtab_nochim.rare <- seqtab_nochim.rare[placenta.names,]
placenta_seqtab_nochim.rare <- placenta_seqtab_nochim.rare[,which(colSums(placenta_seqtab_nochim.rare) > 0 )]

# make genus data frame

#filter genus
placenta_genus.df <- pig_phylotype$genus[placenta.names,]
placenta_genus.df <- placenta_genus.df[,which(colSums(placenta_genus.df) >0)]

#update metadata
placenta_master_metadata.df <- rare_metadata.df[placenta.names, ]


placenta_master_metadata.df$virus <-
  ifelse(placenta_master_metadata.df$group == "GD" , 1,0)

placenta.5d <- rownames(placenta_master_metadata.df[which(placenta_master_metadata.df$date == "5d"), ])
placenta.22d <- rownames(placenta_master_metadata.df[which(placenta_master_metadata.df$date == "22d"), ])

placenta_ids.5d  <- paste0(placenta_master_metadata.df[which(placenta_master_metadata.df$date == "5d"),  "ID"], "_5d")
placenta_ids.22d <- paste0(placenta_master_metadata.df[which(placenta_master_metadata.df$date == "22d"), "ID"], "_22d")

# ANALYSIS: PLACENTA ALPHA-DIVERSITY --------------------------------------------

placenta.richness <- specnumber(x = placenta_seqtab_nochim.rare,MARGIN = 1 )
placenta.shannon  <- diversity(placenta_seqtab_nochim.rare,MARGIN = 1)

#use lm because the pig IDS are different
placenta_richness.lm <- glm(placenta.richness~ date + virus ,
                                     data = placenta_master_metadata.df
)

placenta_shannon.lm <- lm(placenta.shannon~ date + virus,
                                    data = placenta_master_metadata.df
)

summary(placenta_richness.lm)
summary(placenta_shannon.lm)


#richness
placenta_richness.df <- data.frame(richness = placenta.richness,
                                   group = placenta_master_metadata.df$group,
                                   date  = placenta_master_metadata.df$date
                                  )

placenta_richness.plot <-
  ggplot(placenta_richness.df, aes(x = group,
                                   y = richness,
                                   color = group
                                   )
         )+
  geom_point(size = 7)+
  theme(
    text = element_text(size = 16, colour = "black"),
    panel.grid.major = element_line(colour = "grey99"),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(colour = "black"),
    axis.text.x = element_blank(),
    axis.ticks.x =  element_blank(),
    aspect.ratio = 1,
    plot.title = element_text(size = 12, hjust = 0.5),
    strip.background = element_blank()
  )+
  scale_color_brewer("",
                    type = "qual",
                    palette = 2,
                    direction = 1,
                    aesthetics = "color",
                    labels= c("Control",  "Virus")
  )+
  xlab("")+
  ylab("Richness")+
  facet_grid(.~factor(date, levels = c("5d", "22d")),
             labeller = as_labeller(  c("5d" = "Day 5", "22d" = "Day 22")))

#shannon
placenta_shannon.df <- data.frame(shannon = placenta.shannon,
                                  group = placenta_master_metadata.df$group,
                                  date = placenta_master_metadata.df$date)

placenta_shannon.plot <- ggplot(placenta_shannon.df, aes(x = group,
                                                         y = shannon,
                                                         color = group))+

  geom_point(size = 7)+
  theme(
    text = element_text(size = 16, colour = "black"),
    panel.grid.major = element_line(colour = "grey99"),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(colour = "black"),
    axis.text.x = element_blank(),
    axis.ticks.x =  element_blank(),
    aspect.ratio = 1,
    plot.title = element_text(size = 12, hjust = 0.5),
    strip.background = element_blank()
  )+
  scale_color_brewer("",
                     type = "qual",
                     palette = 2,
                     direction = 1,
                     aesthetics = "color",
                     labels= c("Control",  "Virus")
  )+
  xlab("")+
  ylab("Shannon")+
  facet_grid(.~factor(date, levels = c("5d", "22d")),
             labeller = as_labeller(  c("5d" = "Day 5", "22d" = "Day 22")))


pdf("figures/placenta_alpha_diversity.pdf",height = 7, width = 7)
plot_grid(placenta_richness.plot ,
          placenta_shannon.plot,
          labels = c("A","B"),
          label_size = 20,
          nrow =2
)
dev.off()


# ANALYSIS: PLACENTA PCA --------------------------------------------------

#make PCA object
placenta_seqtab_nochim.pca <- prcomp(placenta_seqtab_nochim.rare,scale. = F, center = T)

#get variance %
summary(placenta_seqtab_nochim.pca)

#plot variances in skree plot
plot(placenta_seqtab_nochim.pca)

#check that rownames are the same
all(rownames(placenta_seqtab_nochim.pca$x) == rownames(placenta_master_metadata.df))

#make a PCA data frame
placenta_seqtab_nochim_pca.df <- as.data.frame(placenta_seqtab_nochim.pca$x[,1:5])

#add metadata for plotting
placenta_seqtab_nochim_pca.df$date <- placenta_master_metadata.df$date
placenta_seqtab_nochim_pca.df$group <- placenta_master_metadata.df$group

placenta_seqtab_nochim_pca.plot <- ggplot(data = placenta_seqtab_nochim_pca.df,
                                          aes(x = PC1,
                                              y = PC2,
                                              color = group,
                                              shape = date
                                          )
)

placenta_pca.plot <-  placenta_seqtab_nochim_pca.plot +
  geom_point(size = 7)+
  scale_shape_discrete("Date")+
  scale_color_brewer("Group", type = "qual",
                     palette = 2,
                     direction = 1,
                     aesthetics = "colour",
                     labels= c("Control",  "Virus")
  )+
  theme(
    text = element_text(size = 20, colour = "black"),
    panel.grid.major = element_line(colour = "grey99"),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(colour = "black"),
    aspect.ratio = 1,
    plot.title = element_text(hjust = 0.5)
  )

pdf("figures/placenta_pca_plot.pdf", width = 7, height = 7)
placenta_pca.plot
dev.off()

# ANALYSIS: PLACENTA ADONIS --------------------------------------------------------

set.seed(731)

placenta_seqtab_nochim.adonis2 <- adonis2(placenta_seqtab_nochim.rare ~  date + virus,
                                       permutations = 5000,
                                       by = "margin",
                                       data = placenta_master_metadata.df)
placenta_seqtab_nochim.adonis2

# ANALYSIS: PLACENTA MICROBIOME:SERUM CORRELATIONS ---------------------------

#make placenta correlation data frame
placenta_genus_corr.df <- placenta_genus.df[,which(colMeans(placenta_genus.df) > 10)]

#reset the changes to the metadata
serum_metadata.df <- add_metadata.df[which(add_metadata.df$Sample_type == "Serum"),c(1:5,16:29)]
pbmc_metadata.df <- add_metadata.df[which(add_metadata.df$Sample_type == "PBMC"),c(1:5,10:15,29)]

#subset
serum_metadata.df <- serum_metadata.df[which(serum_metadata.df$combined_id %in% c(placenta_ids.5d, placenta_ids.22d)),]
rownames(serum_metadata.df) <- serum_metadata.df$combined_id
serum_metadata.df <- serum_metadata.df[c(placenta_ids.5d, placenta_ids.22d),]

#subset genus data frame
serum_genus.df <- placenta_genus_corr.df[which(rownames(placenta_genus_corr.df) %in% c(placenta.5d,placenta.22d)),]
serum_genus.df <- serum_genus.df[c(placenta.5d,placenta.22d), ]
serum_genus.df <- serum_genus.df[,which(colSums(serum_genus.df) > 0)]


#any sig corrs ?

#all
serum.pvals <- NULL
serum.est <- NULL
serum.stat <- NULL
serum.host <- NULL
serum.taxa <- NULL

for(i in 6:(ncol(serum_metadata.df)-1)){
  for(j in colnames(serum_genus.df)){
    pval <- cor.test(serum_metadata.df[,i], serum_genus.df[,j], method = "s")$p.value
    est  <- cor.test(serum_metadata.df[,i], serum_genus.df[,j], method = "s")$estimate
    stat <- cor.test(serum_metadata.df[,i], serum_genus.df[,j], method = "s")$statistic
    host <- colnames(serum_metadata.df)[i]
    taxa <- j
    serum.pvals <- c(serum.pvals, pval)
    serum.est   <- c(serum.est, est)
    serum.stat  <- c(serum.stat, stat)
    serum.host  <- c(serum.host, host)
    serum.taxa  <- c(serum.taxa,taxa)
  }
}

placenta_serum.corrs <- data.frame(
  pvals = serum.pvals,
  est  = serum.est,
  stat = serum.stat,
  host = serum.host,
  taxa = serum.taxa,
  qval = qvalue::qvalue(serum.pvals)$qvalues
)

#plot of sig
romboutsia_ifng.plot <- ggplot(data = data.frame(x1 = serum_metadata.df$IFN.g,
                                            y1 = serum_genus.df$Romboutsia),
                          aes(x = x1, y = y1)
)+
  geom_point(size = 3)+
  geom_smooth(method = "lm",se = FALSE, color = "black")+
  annotate("text", x=7, y=300, label= "-0.95", size =6)+
  xlab("IFNg")+
  ylab("Romboutsia")+
  theme(
    text = element_text(size = 16, colour = "black"),
    panel.grid.major = element_line(colour = "grey99"),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(colour = "black"),
    aspect.ratio = 1,
  )

# ANALYSIS: PLACENTA MICROBIOME:PBMC CORRELATIONS ---------------------------

placenta_genus_corr.df <- placenta_genus.df[,which(colMeans(placenta_genus.df) > 10)]

#filter and reorder
pbmc_metadata.df <- pbmc_metadata.df[which(pbmc_metadata.df$combined_id %in% c(placenta_ids.5d, placenta_ids.22d)),]
rownames(pbmc_metadata.df) <- pbmc_metadata.df$combined_id
pbmc_metadata.df <- pbmc_metadata.df[c(placenta_ids.5d, placenta_ids.22d),]
pbmc_metadata.df <- na.omit(pbmc_metadata.df)


pbmc_genus.df <- placenta_genus_corr.df[which(rownames(placenta_genus_corr.df) %in% c("Pla_GA_120_GTGCCATAATCG_L001",
                                                                                      "Pla_GA_122_ACTCTTACTTAG_L001",
                                                                                      "Pla_GD_143_GTTCATTAAACT_L001",
                                                                                      "Pla_GD_144_TTCGATGCCGCA_L001"
                                                                                      )
                                                                                      ),
                                        ]

pbmc_genus.df <- pbmc_genus.df[,which(colSums(pbmc_genus.df) > 0)]

#any sig corrs ?
pbmc.pvals <- NULL
pbmc.est <- NULL
pbmc.stat <- NULL
pbmc.host <- NULL
pbmc.taxa <- NULL

for(i in 6:(ncol(pbmc_metadata.df)-1)){
  for(j in colnames(pbmc_genus.df)){
    pval <- cor.test(pbmc_metadata.df[,i], pbmc_genus.df[,j], method = "s")$p.value
    est  <- cor.test(pbmc_metadata.df[,i], pbmc_genus.df[,j], method = "s")$estimate
    stat <- cor.test(pbmc_metadata.df[,i], pbmc_genus.df[,j], method = "s")$statistic
    host <- colnames(pbmc_metadata.df)[i]
    taxa <- j
    pbmc.pvals <- c(pbmc.pvals, pval)
    pbmc.est   <- c(pbmc.est, est)
    pbmc.stat  <- c(pbmc.stat, stat)
    pbmc.host  <- c(pbmc.host, host)
    pbmc.taxa  <- c(pbmc.taxa,taxa)
  }
}

placenta_pbmc.corrs <- data.frame(
  pvals = pbmc.pvals,
  est  = pbmc.est,
  stat = pbmc.stat,
  host = pbmc.host,
  taxa = pbmc.taxa,
  qval = qvalue::qvalue(pbmc.pvals)$qvalues
)

#Oscillospira
oscillospira_cd8.plot <- ggplot(data = data.frame(x1 = pbmc_metadata.df$cd8illinois,
                                                 y1 = pbmc_genus.df$Oscillospira),
                               aes(x = x1, y = y1)
)+
  geom_point(size = 3)+
  geom_smooth(method = "lm",se = FALSE, color = "black")+
  annotate("text", x=1, y=20, label= "1.00", size =6)+
  xlab("CD8(Illinois)")+
  ylab("Oscillospira")+
  theme(
    text = element_text(size = 16, colour = "black"),
    panel.grid.major = element_line(colour = "grey99"),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(colour = "black"),
    aspect.ratio = 1,
  )

#Succinivibrio
succinivibrio_cd8.plot <- ggplot(data = data.frame(x1 = pbmc_metadata.df$cd8illinois,
                                                  y1 = pbmc_genus.df$Succinivibrio),
                                aes(x = x1, y = y1)
)+
  geom_point(size = 3)+
  geom_smooth(method = "lm",se = FALSE, color = "black")+
  annotate("text", x=1, y=70, label= "1.00", size =6)+
  xlab("CD8(Illinois)")+
  ylab("Succinivibrio")+
  theme(
    text = element_text(size = 16, colour = "black"),
    panel.grid.major = element_line(colour = "grey99"),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(colour = "black"),
    aspect.ratio = 1,
  )


#dgA.11.gut.group
dga_cd8.plot <- ggplot(data = data.frame(x1 = pbmc_metadata.df$cd8illinois,
                                                   y1 = pbmc_genus.df$dgA.11.gut.group),
                                 aes(x = x1, y = y1)
)+
  geom_point(size = 3)+
  geom_smooth(method = "lm",se = FALSE, color = "black")+
  annotate("text", x=1, y=120, label= "1.00", size =6)+
  xlab("CD8(Illinois)")+
  ylab("dgA.11.gut.group")+
  theme(
    text = element_text(size = 16, colour = "black"),
    panel.grid.major = element_line(colour = "grey99"),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(colour = "black"),
    aspect.ratio = 1,
  )

#Parabacteroides
parabacteroides_cd8.plot <- ggplot(data = data.frame(x1 = pbmc_metadata.df$cd8illinois,
                                         y1 = pbmc_genus.df$Parabacteroides),
                       aes(x = x1, y = y1)
)+
  geom_point(size = 3)+
  geom_smooth(method = "lm",se = FALSE, color = "black")+
  annotate("text", x=1, y=35, label= "1.00", size =6)+
  xlab("CD8(Illinois)")+
  ylab("Parabacteroides")+
  theme(
    text = element_text(size = 16, colour = "black"),
    panel.grid.major = element_line(colour = "grey99"),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(colour = "black"),
    aspect.ratio = 1,
  )


# ANALYSIS: PLACENTA MICROBIOME:HOST CORR PANEL PLOTS ---------------------

pdf("figures/placenta_corrs.pdf")
plot_grid(oscillospira_cd8.plot,
          succinivibrio_cd8.plot,
          dga_cd8.plot,
          parabacteroides_cd8.plot,
          romboutsia_ifng.plot,
          labels = "AUTO",label_size = 20)
dev.off()



# ANALYSIS: COMMUNITY BARPLOTS --------------------------------------------


fecal_bar_genus.df <- sweep(fecal_genus.df,
                            MARGIN  =1,
                            rowSums(fecal_genus.df),
                            "/"
                            )

nasal_bar_genus.df <- sweep(nasal_genus.df,
                            MARGIN  =1,
                            rowSums(nasal_genus.df),
                            "/"
)

placenta_bar_genus.df <- sweep(placenta_genus.df,
                            MARGIN  =1,
                            rowSums(placenta_genus.df),
                            "/"
)

#grab the top ten taxa
fecal.keeps    <- names(sort(colMeans(fecal_bar_genus.df), decreasing = T))[1:10]
nasal.keeps    <- names(sort(colMeans(nasal_bar_genus.df), decreasing = T))[1:10]
placenta.keeps <- names(sort(colMeans(placenta_bar_genus.df), decreasing = T))[1:10]

#filter data frame so only includes these top ten taxa
fecal_bar_genus.df    <- fecal_bar_genus.df[,fecal.keeps]
nasal_bar_genus.df    <- nasal_bar_genus.df[,nasal.keeps]
placenta_bar_genus.df <- placenta_bar_genus.df[,placenta.keeps]

#add an other category so everything adds to 1


fecal_bar_genus.df$other    <- 1 - rowSums(fecal_bar_genus.df)
nasal_bar_genus.df$other    <- 1 - rowSums(nasal_bar_genus.df)
placenta_bar_genus.df$other <- 1 - rowSums(placenta_bar_genus.df)


#add some metadata
fecal_bar_genus.df$id <- fecal_master_metadata.df$ID
fecal_bar_genus.df$virus <- fecal_master_metadata.df$virus
fecal_bar_genus.df$vaccine <- fecal_master_metadata.df$vaccine
fecal_bar_genus.df$date <- fecal_master_metadata.df$date
fecal_bar_genus.df$group <- fecal_master_metadata.df$group

#need to changes the genus df into long format for ggplot to play nice with it
fecal_bar_genus.melt <- melt(fecal_bar_genus.df,
                             id.vars = c("id", "virus", "vaccine", "group","date"))

fecal_bar_genus.plot <- ggplot(fecal_bar_genus.melt,
                                aes(x = factor(id),
                                    y = value,
                                    fill = variable))


fecal_bar_genus.plot <-
  fecal_bar_genus.plot+
  geom_col()+
  ylab("Proportional Genus Abundance")+
  xlab("")+
  scale_y_continuous(expand = c(0,0))+
  scale_fill_brewer("Genus", type = "qual",
                    palette = 3,
                    direction = 1,
                    aesthetics = "fill"
  ) +
  facet_grid(.~date,
             scales = "free_x",
             labeller = as_labeller(c("0d"  = "Day 0",
                                      "35d" = "Day 35",
                                      "45d" = "Day 45")
                                      )
             )+
  theme(
    text = element_text(size = 16, colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(colour = "black"),
    axis.text.x = element_blank(),
    axis.ticks.x =element_blank(),
    #aspect.ratio = 1,
    strip.background = element_blank()
  )

pdf("figures/fecal_genus_abudance_barplots.pdf", width = 10 )
fecal_bar_genus.plot
dev.off()

#nasal


#add some metadata
nasal_bar_genus.df$id <- nasal_master_metadata.df$ID
nasal_bar_genus.df$virus <- nasal_master_metadata.df$virus
nasal_bar_genus.df$vaccine <- nasal_master_metadata.df$vaccine
nasal_bar_genus.df$date <- nasal_master_metadata.df$date
nasal_bar_genus.df$group <- nasal_master_metadata.df$group

#need to changes the genus df into long format for ggplot to play nice with it
nasal_bar_genus.melt <- melt(nasal_bar_genus.df,
                             id.vars = c("id", "virus", "vaccine", "group","date"))

nasal_bar_genus.plot <- ggplot(nasal_bar_genus.melt,
                               aes(x = factor(id),
                                   y = value,
                                   fill = variable))


nasal_bar_genus.plot <-
  nasal_bar_genus.plot+
  geom_col()+
  ylab("Proportional Genus Abundance")+
  xlab("")+
  scale_y_continuous(expand = c(0,0))+
  scale_fill_brewer("Genus", type = "qual",
                    palette = 3,
                    direction = 1,
                    aesthetics = "fill"
  ) +
  facet_grid(.~factor(date,levels = c("0d", "5d", "22d")),
             scales = "free_x",
             labeller = as_labeller(c("0d"  = "Day 0",
                                      "5d"  = "Day 5",
                                      "22d" = "Day 22")
             )
  )+
  theme(
    text = element_text(size = 16, colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(colour = "black"),
    axis.text.x = element_blank(),
    axis.ticks.x =element_blank(),
    #aspect.ratio = 1,
    strip.background = element_blank()
  )

pdf("figures/nasal_genus_abudance_barplots.pdf", width = 10 )
nasal_bar_genus.plot
dev.off()

#placenta

#add some metadata
placenta_bar_genus.df$id <- placenta_master_metadata.df$ID
placenta_bar_genus.df$virus <- placenta_master_metadata.df$virus
placenta_bar_genus.df$date <- placenta_master_metadata.df$date
placenta_bar_genus.df$group <- placenta_master_metadata.df$group

#need to changes the genus df into long format for ggplot to play nice with it
placenta_bar_genus.melt <- melt(placenta_bar_genus.df,
                                id.vars = c("id", "virus", "group","date"))

placenta_bar_genus.plot <- ggplot(placenta_bar_genus.melt,
                                  aes(x = factor(id),
                                      y = value,
                                      fill = variable))


placenta_bar_genus.plot <-
  placenta_bar_genus.plot+
  geom_col()+
  ylab("Proportional Genus Abundance")+
  xlab("")+
  scale_y_continuous(expand = c(0,0))+
  scale_fill_brewer("Genus", type = "qual",
                    palette = 3,
                    direction = 1,
                    aesthetics = "fill"
  ) +
  facet_grid(.~factor(date, levels = c("5d", "22d")),
             scales = "free_x",
             labeller = as_labeller(c("5d"  = "Day 5",
                                      "22d" = "Day 22"
                                      )
             )
  )+
  theme(
    text = element_text(size = 16, colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(colour = "black"),
    axis.text.x = element_blank(),
    axis.ticks.x =element_blank(),
    #aspect.ratio = 1,
    strip.background = element_blank()
  )

pdf("figures/placenta_genus_abudance_barplots.pdf", width = 10 )
placenta_bar_genus.plot
dev.off()
