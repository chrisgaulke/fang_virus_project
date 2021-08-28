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
