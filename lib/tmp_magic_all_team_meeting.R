# Quick MAGIC figure for the all-team meeting
# Suzie Hoops

library(vegan)
library(ggplot2)


# Loading temporary files



# NEED TO MAKE A QUALITY CONTROL R SCRIPT!!! Where rarefied table is generated




# Load the rarefied data & aitchison distances
load('data/magic/magic_sample_dataframe.Rdata') # magic_c
load('data/magic/magic_sample_shortinfo.Rdata') # other_info
load('data/magic/aitchison_dists_magic.Rdata')  # aitch_all
meta <- read.delim('~/Dropbox/Magic analysis/data/map/suzie_clean_metadata_short_february2023.txt', sep="\t", header=T)
rownames(meta) <- meta$Sample_ID
## Rarefy and reformat to match metadata
magic_c <- magic_c[rowSums(magic_c) > 500000,] # remove 263 samples below 500k reads
rarefied <- rrarefy(magic_c, 1000000)          # rarefy to 1mil reads

meta <- meta[rownames(rarefied),]

# Alpha diversity over time (age) - box plot w/ jitter and lines
meta$alpha_shannon <- diversity(magic_c, index="shannon")



# 3D PCoA plot of all samples (Aitchison distance)

# Interaction between intrapartum abx and not in first month


# Sankey plot of samples at each age?