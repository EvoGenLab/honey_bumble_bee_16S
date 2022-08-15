####---- DADA2 pipeline ---####


#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install()

library("BiocManager")

# didn't work
biocLite(suppressUpdates = FALSE)
#biocLite("ShortRead", suppressUpdates = FALSE)

#BiocManager::install("devtools") #THIS WORKED
library("devtools")
#devtools::install_github("benjjneb/dada2")

BiocManager::install("biocLite")
install.packages("biocLite")
#biocLite("devtools")
#devtools::install_github("benjjneb/dada2") # use the ref="v1.10" (e.g.) argument to get specific versions
#asks which packages I want to update? Updated all...

library(dada2)
packageVersion("dada2")
## received warning:no function found corresponding to methods exports from ‘GenomicAlignments’ for: ‘concatenateObjects’
# https://benjjneb.github.io/dada2/tutorial.html

#### ---- set working directory ----####
setwd("~/Bee Data/BeeMicro11.2019")
path <- "~/Bee Data/BeeMicro11.2019"
list.files(path)

#HC's computer
setwd("/home/student/Desktop/BeeMicrobiome/BeeMicro11.2019")
path <-"/home/student/Desktop/BeeMicrobiome/BeeMicro11.2019"

  #### ---- Sample names #### 
# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
sample.names

####---INSPECTION OF QUALITY PROFILES----####
# when looking at quality profile; want Quality score to be above Q30
## For our primers you need trunLen to come to 460 bp for reads to merge well later
# amplicon is ~385 bp and need additional 50 bp to pair end.  This puts you between 440-450 bp needed.
# the longer the bp, the more sequences you get. Try and find happy medium.

plotQualityProfile(fnFs[1:4])
#numbers are number of different samples looking at
#truncate at 260 ish

plotQualityProfile(fnRs[1:4])
#truncate at try 195


####----FILTERING AND TRIMMING----####
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

## Try different run parameters to see how it impacts output 
#maxEE = Max number of expected errors allowed in read

## For our primers you need trunLen to come to 460 bp for reads to merge well later
out1 <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(260,195),
                      maxN=0, maxEE=c(2,2), trimLeft = 19, trimRight = 23, 
                      truncQ=2, rm.phix=TRUE,
                      compress=TRUE, multithread=TRUE) 
head(out1)
out1
str(out1)
## average reads.out per sample
mean(out1[,2])
#on average = 2064.797 reads

#out2 <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(270,200),
                      #maxN=0, maxEE=c(2,2), trimLeft = 19, trimRight = 23, 
                      #truncQ=2, rm.phix=TRUE,
                      #compress=TRUE, multithread=TRUE)
#head(out2)
#mean(out2[,2])
#average = 2018.273
#going with out1

###--Learn Error rate----####
#ONCE you have decided on the best parameters, run the out command last.  
#going to go with out1 for now
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

plotErrors(errF, nominalQ=TRUE)
#Warning messages: Transformation introduced infinite values in continious y axis- MK 4.15.19

####----dereplicate? ----####
derepFs <- derepFastq(filtFs, verbose=TRUE)

derepRs <- derepFastq(filtRs, verbose=TRUE)

# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

###---sample inferences---####
#filter and trip sequence data
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)

dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
dadaFs[[1]]
#5 seq variants from 73 inputs 11.15.19
#26 sequence variants were inferred from 489 input unique sequences

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
seqtab <- makeSequenceTable(mergers)

dim(seqtab)
#128 samples, 2165 sequences
table(nchar(getSequences(seqtab)))

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
#128 287
sum(seqtab.nochim)/sum(seqtab)
#0.8647883
#chimeras account for ~14%

####----Track what reads made it through dada2 pipeline ---####
getN <- function(x) sum(getUniques(x))

track <- cbind(out1, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
#write.csv(track, "dada2output11.8.19.csv")

#in excel file divide nochim sequences by input.  
#33% of the data went through - kinda low

####----Assign taxonomy ----####
getwd()
## Assign taxonomy is up three directories so that I can use these files for multiple projects
taxa <- assignTaxonomy(seqtab.nochim, "../rdp_train_set_16.fa.gz", multithread=TRUE)

taxa <- addSpecies(taxa, "../rdp_species_assignment_16.fa.gz")

#write.csv(taxa, "taxa_trial.csv")
#write.csv(seqtab, "seqtab.csv")
#write.csv(seqtab.nochim, "seqtab.nochim.csv")


#  inspect the taxonomic assignments:
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

####----Install phyloseq---####
#https://joey711.github.io/phyloseq/install.html



#use to insall phyloseq
BiocManager::install("phyloseq")


library(Biostrings)
library(ggplot2)

#source('http://bioconductor.org/biocLite.R')
#biocLite('phyloseq')


library(phyloseq)
packageVersion("phyloseq")


#OTU = otu_table(seqtab.nochim, taxa_are_rows = TRUE)

#try changin seq tab from integer to character?
#class(taxa)
#class(seqtab.nochim)

#conver to matrix
#m.seqtab.nochim <- as.matrix(seqtab.nochim)
#m.taxa <- as.matrix(taxa)

## combine feature table and taxonomy table in same order
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = FALSE), 
               tax_table(taxa))
ps
#287 taxa 

## rename ASVs to numbers
new.names <- paste0("ASV", seq(ntaxa(ps))) # Define new names ASV1, ASV2, ...
seqs <- taxa_names(ps) # Store sequences
names(seqs) <- new.names # Make map from ASV1 to full sequence
taxa_names(ps) <- new.names # Rename to human-friendly format

install.packages("seqRFLP")

library(seqRFLP)
## should work, but, below can also convert .csv file to fasta
#seq_data <- dataframe2fas(seqs, file = 'IMLS_DNAsequences.fasta')


## convert feature table to matrix
site_species <-as(otu_table(ps), "matrix")

## need to change this to match mapping file later
rownames(site_species)

## transpose to make a species by site matrix

species_site <- t(site_species)

# taxon table #This will give you large sequences for each ASV
tax <- as(tax_table(ps), "matrix")
seqs
tax
####----Write csv of tax table---####
getwd()
setwd("/home/student/Desktop/BeeMicrobiome/")

## Write this file out and look at it. Determine what to do with ASVs in Neg Controls
## select whole worksheet and filter by reads in neg controls, 
## MAKE SURE you have selected all columns/rows before you filter

## Remove ASVs if they are in all negative controls or in one negative control, but almost all samples
## pretty subjective, no guidelines really on what to do here
## now we are using decontam to deal with this!
## Once done making neg control decisions, remove negative controls from file 
## and rename _final.csv
write.csv(species_site, "Bee_featuretable.csv")

write.csv(tax, "Bee_taxonomy.csv")
write.csv(seqs, 'Bee_DNAsequences.csv')


#library(seqRFLP) 
## can also convert .csv file to fasta

seq_data <- read.csv("Bee_DNAsequences.csv", header = T)
seq_data <- dataframe2fas(seq_data, file = "Bee_DNAsequences2.fasta")

