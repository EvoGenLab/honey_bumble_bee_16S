#Bee Pre-process 11/15/19

library(phyloseq)
library(ape)
library(vegan)
library(ggplot2)
library(devtools)
#devtools::install_github("benjjneb/decontam")
library(decontam)

####----Phyloseq Object----####
path
##Set working directory to bring in files##
setwd("/home/student/Desktop/BeeMicrobiome/")

#read in feature table; DNA sequences with sample IDs
feature_tab <- read.csv("Bee_featuretable.csv", header = T, row.names = 1)
#Make compatible for phyloseq format
asv_tab = otu_table(feature_tab, taxa_are_rows = TRUE)
dim(asv_tab)

#read in meta data file:
meta_data <- read.csv("Mastersheet_11.15.19.csv", header=T, row.names = 3)
# dataframe is expected for sample_data
#md <- subset(meta_data, rownames =129)
#rm(md)

#md <- meta_data[-c(129, 130),]

#class(md)
# make compatible for phyloseq
mf <- sample_data(meta_data)

#Error in make.names(col.names, unique = TRUE) : 
#  invalid multibyte string 13
#use this code if you can't get it to import and want to check what may be causing it
#x = read.csv("Mastersheet_11.15.19.csv", check.names = F)
#iconv(names(x), to = "ASCII", sub = "")
#iconv(names(x), to = "ASCII", sub = "")
#iconv(x, from = "", to = "", sub = NA, mark = TRUE, toRaw = FALSE)
#iconv(x, from = "\n", to= "", sub=NA)


#read in taxonomy file
taxonomy <- read.csv("Bee_taxonomy.csv", row.names=1)
# Needs to be a matrix
taxonomy <- as.matrix(taxonomy)
# Make compatible for phyloseq
taxonomy_final = tax_table(taxonomy)


library(Biostrings)
DNAseqs <- readDNAStringSet("Bee_DNAsequences2.fasta")

# You can also add a phylogenetic tree here, if you have one
tree = read.tree("tree.nwk")
#to make tree:
#open up terminal
#follow textedit file detailing download of qiime 
# and activating environment in terminal
#not sure what we do next?


#inspect sample names:
#sample_names(ps)
#sn <- sample_names(ps)
#sn
#write.csv(sn, "samplenames.csv")
#missing osmia145 and chrysidadae163

#gplot::venn(list(metadata = rownames(meta_data), featuretable = colnames(feature_tab)))

#View(meta_data$feature_tab)

# Merge it all together
#Make phyloseq object
ps <- merge_phyloseq(asv_tab, taxonomy_final, mf, DNAseqs, tree)
#componenet sample names do not match
#sam_data is empty - part of meta data?
ps
#287 taxa by 128 samples
head(sample_data(ps))
str(ps)

sum(sample_sums(ps))
#bee data: 211,132
#Pilot data:1,329,079
#IMLS: 2,128740

df.bee.samp <- as(sample_data(ps), "data.frame")
m.bee.samp <- as(sample_data(ps), "matrix")
## I am missing most of my negative controls.  At first it was a naming issue, but i corrected it
#still not working. 


####----Filter----####
## Are there ASVs that are no longer present in subsetted data.
sum(taxa_sums(ps) == 0)

## getting rid of ASvs that are in negative controls only
psFT = filter_taxa(ps, function(x) sum(x) !=0, TRUE)

psFT
#287 taxa
sort(sample_sums(psFT))

## filter singletons
psFS <- filter_taxa(psFT, function (x) {sum(x > 0) > 1}, prune=TRUE)
sort(sample_sums(psFS))
psFS
#down to 144 taxa

bee_data <- psFS

sum(sample_sums(bee_data))
#Total sequences 202046
sort(sample_sums(bee_data))

min(sample_sums(bee_data))
#3
max(sample_sums(bee_data))
#10929
max(sample_sums(bee_data))/min(sample_sums(bee_data))
#3643
mean(sample_sums(bee_data))
#1578.484

#phyloseq object is called bee_data from now on

#remove samples with sample sum below 500
#figure out how to write code. 

####----Alpha Diversity PreProcess----####
#Alpha Diversity pre process, Plot Tree, Phylogenetic Diversity and Species Richness
#sample_data and otu_table are generic names given to files in phyloseq object
#we turn them into a data.frame and matrix

#Making data frame of sample data 
DF.bee <- as(sample_data(bee_data), "data.frame")

#making matrix of otu table. transposed (flipping col and rows). 
t_otu <-t(as(otu_table(bee_data), "matrix"))

install.packages("picante")
library(picante)

tree$tip.label

prunedTree <- prune.sample(t_otu,tree)
plot(tree)
## pd estimates phylogentic diversity (Faith's Phylogenetic Diversity), higher Faith's PD value = more phylogenetically diverse community 
PD <- pd(t_otu, prunedTree, include.root = F)
#Returns a dataframe. SR is species richness
#ses.pd() compares observed PD to values expected under randomizations.  Takes awhile to run.

#need to have both alpha and df having the same column info
PD$SampleID <- row.names(PD)

SUMseqs <- as.data.frame(sample_sums(bee_data))
SUMseqs$SampleID <- row.names(SUMseqs)
SUMseqs

DF.bee$SampleID <- row.names(DF.bee)
DF.bee
#now merge to get sequence counts and SR and PD in mapping file
meta_SRPD <- merge(DF.bee, PD, by = "SampleID")
meta_SRPDseqs<- merge(meta_SRPD, SUMseqs, by = "SampleID")
#figured this out! Needed to add SampleID column to DFele dataframe. 
meta_SRPDseqs
#This is a dataframe NOT a phyloseq object.  Will be used in phyloseq object in future code.
#this is my new meta data file
#after I add rarefaction, I will write it to a csv.

write.csv(meta_SRPDseqs, "Meta_PDSR.csv")

#trials:
#option - filter by sample sums AFTER making phyloseq object. 
#meta2 <- subset_samples(meta_SRPDseqs, (sample_sums(bee_data)> 500))
#meta2 <- meta_SRPDseqs$`sample_sums(bee_data)`>500
#meta2 <-subset(meta_SRPDseqs, sample_sums(bee_data) >500)


####----Rarefaction----####
#### 11.20.19 - wait until we filter out low sample sums ###
## Because there is a 22.94x difference then rarefying
eleRare = rarefy_even_depth(ele_data, replace=FALSE, rngseed = 711)
sort(sample_sums(eleRare))
#all samples have 1738 sequence

dfEleR <- as(sample_data(eleRare), "data.frame")

t_otuR <-t(as(otu_table(eleRare), "matrix"))

tree$tip.label

prunedTree <- prune.sample(t_otuR,tree)
## pd estimates phylogentic diversity (Faith's Phylogenetic Diversity), higher Faith's PD value = more phylogenetically diverse community 
PD2 <- pd(t_otuR, prunedTree, include.root = F)

colnames(PD2)[colnames(PD2)=="PD"] <- "PD_rare"
colnames(PD2)[colnames(PD2)=="SR"] <- "SR_rare"

#need to have both alpha and df having the same column info
PD2$SampleID <- row.names(PD2)

#now merge to get SR and PD of rarefied in mapping file
meta_rare <- merge(meta_SRPDseqs, PD2, by = "SampleID")

#This is my new meta data file to use in future phyloseq objects for analysis
write.csv(meta_rare, "Meta_EleRare.csv")



####----RAREFACTION PLOTS----####

## Load packages
library(ggplot2)
library(phyloseq)
library(reshape2)
library(ape)
library(gridExtra)

scripts <- c("graphical_methods.R",
             "tree_methods.R",
             "plot_merged_trees.R",
             "specificity_methods.R",
             "ternary_plot.R",
             "richness.R",
             "edgePCA.R",
             "copy_number_correction.R",
             "import_frogs.R",
             "prevalence.R",
             "compute_niche.R")
urls <- paste0("https://raw.githubusercontent.com/mahendra-mariadassou/phyloseq-extended/master/R/", scripts)

for (url in urls) {
  source(url)
}


p <- ggrare(bee_data, step = 50, color="Species", se = FALSE)
# took out: color = "TYPE2"
p + xlab("Sequence count") + ylab("Bacterial ASV richness") +
  theme(text = element_text(size = 15)) + labs(color='Bee Species') +xlim(0,1000) 

####----Write CSV files----####
#11.20.19 did not do this yet
#Write csv files of final taxonomy and meta data files...

## OTU TABLE
## convert feature table to matrix
## ele_data is your phyloseq object
species_site <-as(otu_table(ele_data),"matrix")

write.csv(species_site, "Ele_feature_table_final2.csv")

## TAXONOMY TABLE
# taxon table 
tax <- as(tax_table(ele_data), "matrix")

write.csv(tax, "Ele_taxonomy_final.csv")

## SEQUENCE DATA
DNAseqs <- taxa_names(ele_data) # Store sequences
library(seqRFLP)
## should work, but might not, below can also convert .csv file to fasta
DNAseq_data2 <- dataframe2fas(DNAseqs, file = 'Ele_DNAsequences_final.fasta')

write.csv(DNAseqs, 'Ele_DNAsequences_final.csv')

## can also convert .csv file to fasta

#seq_data <- read.csv("SalAMP_feature_DNAsequences.csv", header = T)
#seq_data <- dataframe2fas(seq_data, file = "SalAMP_DNAsequences.fasta")

## SAMPLE DATA
## don't know if this will work
meta_final <- as(sample_data(ele_data), "data.frame")

#Didn't write csv below because already made this after rarefying. 
write.csv(meta_rare, "Ele_Meta_Rare_final.csv")
#write.csv(meta_final, 'Ele_meta_final.csv')


