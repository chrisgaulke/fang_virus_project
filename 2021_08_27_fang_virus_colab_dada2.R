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

## Install Bioconductor

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install(version = "3.13")

## Install dada2
#BiocManager::install("dada2")
library(dada2)
library(ggplot2)
library(vegan)


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


# ANALYSIS: NORMALIZE -----------------------------------------------------

seqtab_nochim.relabd <- sweep(seqtab.nochim,
                              MARGIN = 1,
                              rowSums(seqtab.nochim),
                              FUN = '/'
)


# ANALYSIS: PCA -----------------------------------------------------------

#make PCA object
seqtab_nochim.pca <- prcomp(seqtab_nochim.relabd,scale. = F, center = T)

#check that rownames are the same
all(rownames(seqtab_nochim.pca$x) == rownames(Master_metadata.df))

#make a PCA data frame
seqtab_nochim_pca.df <- as.data.frame(seqtab_nochim.pca$x[,1:5])

#add metadata for plotting
seqtab_nochim_pca.df$Vaccination_status <- Master_metadata.df$Vaccination_status
seqtab_nochim_pca.df$Sample_type <- Master_metadata.df$Sample_type
seqtab_nochim_pca.df$date <- Master_metadata.df$date
seqtab_nochim_pca.df$group <- Master_metadata.df$group
seqtab_nochim_pca.df$ID <- Master_metadata.df$ID

seqtab_nochim_pca.plot <- ggplot(data = seqtab_nochim_pca.df,
                                 aes(x = PC2,
                                     y = PC3,
                                     color = group,
                                     shape = Sample_type))

seqtab_nochim_pca.plot +
  geom_point(size = 3)


# ANALYSIS: ADONIS --------------------------------------------------------

set.seed(731)

seqtab_nochim.adonis <- adonis(seqtab_nochim.relabd ~ group + Sample_type + Vaccination_status,
                               permutations = 5000,
                               data = Master_metadata.df)

seqtab_nochim.adonis


# ANALYSIS: FECAL ANALYSIS ------------------------------------------------


#Get names for fecal samples

fecal.names <- rownames(Master_metadata.df[which(Master_metadata.df$Sample_type == "FS"), ])

fecal_seqtab_nochim.relabd <- seqtab_nochim.relabd[fecal.names,]
fecal_master_metadata.df   <- Master_metadata.df[fecal.names,]


# ANALYSIS: FECAL PCA -----------------------------------------------------


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

fecal_seqtab_nochim_pca.plot +
  geom_point(size = 3)+
  scale_color_brewer(palette = "Set1") # look into better colors

fecal_seqtab_nochim_pca.plot +
  geom_point(data = ~ subset(., date == "35d") ,size = 3)


fecal_seqtab_nochim_pca.plot +
  geom_point(data = ~ subset(., date == "45d") ,size = 3)


# ANALYSIS: FECAL ADONIS --------------------------------------------------------

set.seed(731)

fecal_seqtab_nochim.adonis <- adonis(fecal_seqtab_nochim.relabd ~  group * date +ID ,
                               permutations = 5000,
                               data = fecal_master_metadata.df)

fecal_seqtab_nochim.adonis


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
