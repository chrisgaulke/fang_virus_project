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

rownames(pig_phylotype$genus) == rownames(rare_metadata.df)

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
                                     color = group,
                                     shape = Sample_type))

seqtab_nochim_pca.plot +
  geom_point(size = 3)+
  scale_fill_brewer(  type = "qual",
                      palette = 2,
                      direction = 1,
                      aesthetics = "colour"
  )

# ANALYSIS: ADONIS --------------------------------------------------------

set.seed(731)

seqtab_nochim.adonis <- adonis(seqtab_nochim.rare ~ group + Sample_type + Vaccination_status,
                               permutations = 5000,
                               data = rare_metadata.df)

seqtab_nochim.adonis


# ANALYSIS: FECAL ANALYSIS ------------------------------------------------

#Get names for fecal samples

fecal.names <- rownames(rare_metadata.df[which(rare_metadata.df$Sample_type == "FS"), ])

#filter ASV table
fecal_seqtab_nochim.relabd <- seqtab_nochim.relabd[fecal.names,]
fecal_seqtab_nochim.relabd <- fecal_seqtab_nochim.relabd[,which(colSums(fecal_seqtab_nochim.relabd) > 0 )]

#filter metadata
fecal_master_metadata.df   <- rare_metadata.df[fecal.names,]
fecal_master_metadata.df$vaccine <- ifelse(fecal_master_metadata.df$group == "G3" & fecal_master_metadata.df$date != "0d", 1,0)
fecal_master_metadata.df$virus <- ifelse(fecal_master_metadata.df$group != "G5" & fecal_master_metadata.df$date == "45d", 1,0)

#filter genus
fecal_genus.df <- pig_phylotype$genus[fecal.names,]
fecal_genus.df <- fecal_genus.df[,which(colSums(fecal_genus.df) >0)]


# ANALYSIS: FECAL PCA ASV-----------------------------------------------------

#make PCA object
fecal_seqtab_nochim.pca <- prcomp(fecal_seqtab_nochim.relabd,scale. = F, center = T)

#get varaince %
summary(fecal_seqtab_nochim.pca)

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
  geom_point(size = 5) +
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
               data = fecal_seqtab_nochim_pca.df,
               aes(x = PC1,
                   y = PC2,
                   color = group,
                   fill = group)
  )+
  scale_fill_brewer(type = "qual",
                    palette = 2,
                    direction = 1,
                    aesthetics = "fill"
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
  )+
  guides(fill = "none") +ggtitle("All")



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
                                   labels = "AUTO"
)
fecal_asv_by_day.plot

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



#all
fecal_genus_pca.plotall <- fecal_genus_pca.plot +
  geom_point(size = 5) +
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
               data = fecal_genus_pca.df,
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
                                     labels = "AUTO"
                                     )
fecal_genus_by_day.plot

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
fecal_genus2.adonis <- adonis(fecal_genus.df ~ date + virus* vaccine  ,
                              permutations = 5000,
                              data = fecal_master_metadata.df)
fecal_genus2.adonis

# here there appears to be a significant effect of data on the animals. This
# could be due to a common environmental parameter being changed in the animals
# at these time points. What might have been changing?

# another interesting point is that G3 and G5 cluster more together while G4 is
# separate


# ANALYSIS: FECAL GENUS STATS -------------------------------------------

#basic linear model
fit0 <- lm(fecal_genus.df$Clostridium.sensu.stricto.1 ~ date + virus* vaccine,data = fecal_master_metadata.df )
summary(fit0)

#but our data is not normally distributed

hist(fecal_genus.df$Clostridium.sensu.stricto.1)
hist(fecal_genus.df$Lactobacillus)
hist(fecal_genus.df$Blautia)

#generalized linear models are sometimes used when normality cannot be assumed
fit1 <- glm(fecal_genus.df$Clostridium.sensu.stricto.1 ~ date + virus* vaccine,data = fecal_master_metadata.df )
summary(fit1)

#But we know from the lit and prior experience that these data are best fit by a neg bionomial distribution
fit2 <- glm.nb(fecal_genus.df$Clostridium.sensu.stricto.1 ~ date + virus* vaccine,data = fecal_master_metadata.df )
summary(fit2)

# But still we have the other issue of have repeated sampling of individuals across time
# to account for this variation we can use a mixed model which accounts for the
# variation in response of each individual

library(lme4)
fit3 <- glmer.nb(fecal_genus.df$Clostridium.sensu.stricto.1 ~ date + virus * vaccine + (1|ID) ,data = fecal_master_metadata.df )
summary(fit3)

# But now we have a problem because we only have the fixed effects p-value and
# not the overall model. This presents an issue with FDR control because family-wise
# error has already been accounted for. To solve this problem we can use a NULL
# model which only includes an implicit intercept and the random (mixed) effect

fit3.null <- glmer.nb(fecal_genus.df$Clostridium.sensu.stricto.1 ~ (1|ID) ,data = fecal_master_metadata.df )
summary(fit3.null)

#if the anova is significant adding those extra terms allows significantly more
# variation to be explained by the model
anova(fit3.null, fit3)

# it should be noted that here, while the null model above is appropriate, a null
# that includes date would also be OK (and potentially better). Regardless now the
# anova will be the thing we will correct on

#so to do all the stats we can use a loop and we will cover that later


# ANALYSIS: nasal ANALYSIS ------------------------------------------------

#Get names for nasal samples

nasal.names <- rownames(Master_metadata.df[which(Master_metadata.df$Sample_type == "NS"), ])

nasal_seqtab_nochim.relabd <- seqtab_nochim.relabd[nasal.names,]
nasal_master_metadata.df   <- Master_metadata.df[nasal.names,]


# ANALYSIS: nasal PCA -----------------------------------------------------


#make PCA object
nasal_seqtab_nochim.pca <- prcomp(nasal_seqtab_nochim.relabd,scale. = F, center = T)

#get varaince %
summary(nasal_seqtab_nochim.pca)

#plot variances in skree plot
plot(nasal_seqtab_nochim.pca)

#check that rownames are the same
all(rownames(nasal_seqtab_nochim.pca$x) == rownames(nasal_master_metadata.df))

#make a PCA data frame
nasal_seqtab_nochim_pca.df <- as.data.frame(nasal_seqtab_nochim.pca$x[,1:5])

#add metadata for plotting
nasal_seqtab_nochim_pca.df$Vaccination_status <- nasal_master_metadata.df$Vaccination_status
nasal_seqtab_nochim_pca.df$date <- nasal_master_metadata.df$date
nasal_seqtab_nochim_pca.df$group <- nasal_master_metadata.df$group
nasal_seqtab_nochim_pca.df$ID <- nasal_master_metadata.df$ID

nasal_seqtab_nochim_pca.plot <- ggplot(data = nasal_seqtab_nochim_pca.df,
                                       aes(x = PC2,
                                           y = PC3,
                                           color = group,
                                           shape = date))

nasal_seqtab_nochim_pca.plot +
  geom_point(size = 3)



nasal_seqtab_nochim_pca.plot +
  geom_point(data = ~ subset(., date == "0d") ,size = 3)

nasal_seqtab_nochim_pca.plot +
  geom_point(data = ~ subset(., date == "5d") ,size = 3)

nasal_seqtab_nochim_pca.plot +
  geom_point(data = ~ subset(., date == "22d") ,size = 3)


# ANALYSIS: nasal ADONIS --------------------------------------------------------

set.seed(731)

nasal_seqtab_nochim.adonis <- adonis(nasal_seqtab_nochim.relabd ~  group * date +ID ,
                                     permutations = 5000,
                                     data = nasal_master_metadata.df)

nasal_seqtab_nochim.adonis
