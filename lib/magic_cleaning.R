# Cleaning up the magic data
# Suzie Hoops
# Nov. 11 2022
# currently available runs from magic: 046, 055, 068 (1-4), 070, 076-1 
library(stringr)
library(dplyr)


##### Read in Files & Clean Names #####
# Run 046
r46 <- read.csv("data/magic/run046-taxatable-burst-capitalist.txt", row.names=1, header=1, sep="\t")
colnames(r46) <- gsub("[.]S.*[.]fa", "", colnames(r46))
r46 <- as.data.frame(t(r46))
r46$Sample_ID <- rownames(r46)

# Run 055
r55 <- read.csv("data/magic/run055-taxatable-burst-capitalist.txt", row.names=1, header=1, sep="\t")
colnames(r55) <- paste0("magic.", str_pad(gsub("X", "", gsub("[.]S.*[.]fa", "", colnames(r55))), 4, side="left", pad="0"))
r55 <- as.data.frame(t(r55))
r55$Sample_ID <- rownames(r55)

# Run 068-1
r681 <- read.csv("data/magic/run068-1-taxatable-burst-captialist.txt", row.names=1, header=1, sep="\t")
colnames(r681) <- gsub("[.]S.*[.]fa", "", colnames(r681))
r681 <- as.data.frame(t(r681))
r681$Sample_ID <- rownames(r681)

# Run 068-2
r682 <- read.csv("data/magic/run068-2-taxatable-burst-capitalist.txt", row.names=1, header=1, sep="\t")
colnames(r682) <- gsub("[.]S.*[.]fa", "", colnames(r682))
r682 <- as.data.frame(t(r682))
r682$Sample_ID <- rownames(r682)

# Run 068-3
r683 <- read.csv("data/magic/run068-3-taxatable-burst-capitalist.txt", row.names=1, header=1, sep="\t")
colnames(r683) <- gsub("[.]S.*[.]fa", "", colnames(r683))
r683 <- as.data.frame(t(r683))
r683$Sample_ID <- rownames(r683)

# Run 068-4
r684 <- read.csv("data/magic/run068-4-taxatable-burst-capitalist.txt", row.names=1, header=1, sep="\t")
colnames(r684) <- gsub("[.]S.*[.]fa", "", colnames(r684))
r684 <- as.data.frame(t(r684))
r684$Sample_ID <- rownames(r684)

# Run 070
r70 <- read.csv("data/magic/run070-taxatable-burst-capitalist.txt", row.names=1, header=1, sep="\t")
colnames(r70) <- gsub("[.]S.*[.]fa", "", colnames(r70))
r70 <- as.data.frame(t(r70))
r70$Sample_ID <- rownames(r70)

# Run 076-1
r761 <- read.csv("data/magic/run076-1-taxatable-burst-capitalist.txt", row.names=1, header=1, sep="\t")
colnames(r761) <- gsub("[.]S.*[.]fa", "", colnames(r761))
r761 <- as.data.frame(t(r761))
r761$Sample_ID <- rownames(r761)


##### Combine into one matrix & Clean #####
magic_c <- rbind.fill(r46, r55, r681, r682, r683, r684, r70, r761)
# Remove duplicates & blanks & .H or ATCC or .G
magic_c <- magic_c[-grep("[.][H|G]", magic_c$Sample_ID),]         # removes 2
magic_c <- magic_c[-grep("blank", tolower(magic_c$Sample_ID)),]   # removes 1
magic_c <- magic_c[-grep("control", tolower(magic_c$Sample_ID)),] # removes 1
magic_c <- magic_c[-grep("ATCC", magic_c$Sample_ID),]             # removes 3
magic_c <- magic_c[-which(duplicated(magic_c$Sample_ID)),]        # removes 11
rownames(magic_c) <- magic_c$Sample_ID
magic_c$Sample_ID <- NULL
magic_c <- t(magic_c)
magic_c[is.na(magic_c)] = 0
magic_c <- as.data.frame(magic_c)


##### Map File #####
# Read in table
meta_magic <- read.table("data/magic/clean_metafile.txt", sep="\t", header=T)
rownames(meta_magic) <- meta_magic$Sample_ID
sum(colnames(magic_c) %in% meta_magic$Sample_ID) # 1750 of 1757 match
# Only keep samples matching metadata
magic_c <- magic_c[,(colnames(magic_c) %in% meta_magic$Sample_ID)]


##### Normalize and Reduce Sizes #####
# Check out count totals per smaple... range is crazy (1k --> 6m)
head(sort(colSums(magic_c), decreasing=F),100)
tail(sort(colSums(magic_c), decreasing=F),20)
# Drop everything with less than 140k (drops 70ish samples)
magic_c <- magic_c[,(colSums(magic_c) >= 140000)] # 1750 to 1679 samples
# Drop non-existent taxa in table
magic_c <- magic_c[(rowSums(magic_c) > 0),] # drops ~200 taxa non-existent
magic_c <- magic_c[(rowSums(magic_c) > 1),] # drops ~1410 singletons
# Create a normalized counts table
magic_n <- sweep(magic_c, 2, colSums(magic_c), FUN="/")
# Drop samples from metadata (match order to magic_c)
meta_magic <- meta_magic[colnames(magic_c),]


##### Write Out: Counts, Norm, Map #####
write.table(magic_c, "data/magic/mycleaned_magic_counts.txt", sep="\t", row.names=T, col.names=T)
write.table(magic_n, "data/magic/mycleaned_magic_norm.txt", sep="\t", row.names=T, col.names=T)
write.table(meta_magic, "data/magic/mycleaned_magic_metadata.txt", sep="\t", row.names=T, col.names=T)


