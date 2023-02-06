# Cleaning up the magic data
# Suzie Hoops
# Nov. 11 2022, updated February 2023
# currently available runs from magic: 046, 055, 068 (1-4), 070, 076-1 (rest of 076 and 077 in process)
library(stringr)
library(dplyr)
library(plyr)
library(vegan)

# Optional/Temporary : load datasets (these were created here!)
# load("data/magic/magic_sample_dataframe.Rdata") # loads magic_c
# load("data/magic/magic_sample_shortinfo.Rdata") # loads other_info
# load("data/magic/aitchison_dists_magic.Rdata")  # loads aitch_all

# NOTE: BLANK/Control Counts (29 total)
# 046 - 1 Blank, 1 Positive.Control
# 055 - none
# 068 - none
# 069 - 3 blanks: BLANK_PPE, Blank__PE, Blank_MagAE (skipping 069)
# 070 - none
# 076 - 11 ATCC (3 pool1, 3 pool2, 5 pool3)
# 077 - 3 ATCC, 1 PositiveControl (pool1), 4 ATCC (pool 2), 5 ATCC (pool3)


##### Read in Files & Clean Names #####
# INITIALIZE: Dataframe keeping track of Sequencing & Extraction
seq_name_df <- data.frame(Sample_ID=c(),Sequencing_Project=c(), Filename=c(), Extraction=c())

# Run 046 
r46 <- read.csv("data/magic/run046-taxatable-burst-capitalist.txt", row.names=1, header=1, sep="\t")
r46_filenames <- colnames(r46)
colnames(r46) <- gsub("[.]S.*[.]fa", "", colnames(r46))
r46 <- as.data.frame(t(r46))
r46$Sample_ID <- rownames(r46)
seq_name_df <- rbind(seq_name_df,data.frame(Sample_ID=r46$Sample_ID,
                                            Sequencing_Project="Knights_Project_046",
                                            Filename=r46_filenames,
                                            Extraction="MagAttract"))

# Run 055
r55 <- read.csv("data/magic/run055-taxatable-burst-capitalist.txt", row.names=1, header=1, sep="\t")
r55_filenames <- colnames(r55)
colnames(r55) <- paste0("magic.", str_pad(gsub("X", "", gsub("[.]S.*[.]fa", "", colnames(r55))), 4, side="left", pad="0"))
r55 <- as.data.frame(t(r55))
r55$Sample_ID <- rownames(r55)
seq_name_df <- rbind(seq_name_df,data.frame(Sample_ID=r55$Sample_ID,
                                            Sequencing_Project="Knights_Project_055",
                                            Filename=r55_filenames,
                                            Extraction="PowerSoil"))

# Run 068-1
r681 <- read.csv("data/magic/run068-1-taxatable-burst-captialist.txt", row.names=1, header=1, sep="\t")
r681_filenames <- colnames(r681)
colnames(r681) <- gsub("[.]S.*[.]fa", "", colnames(r681))
r681 <- as.data.frame(t(r681))
r681$Sample_ID <- rownames(r681)
seq_name_df <- rbind(seq_name_df,data.frame(Sample_ID=r681$Sample_ID,
                                            Sequencing_Project="Knights_Project_068_Pool1",
                                            Filename=r681_filenames,
                                            Extraction="PowerSoil"))

# Run 068-2
r682 <- read.csv("data/magic/run068-2-taxatable-burst-capitalist.txt", row.names=1, header=1, sep="\t")
r682_filenames <- colnames(r682)
colnames(r682) <- gsub("[.]S.*[.]fa", "", colnames(r682))
r682 <- as.data.frame(t(r682))
r682$Sample_ID <- rownames(r682)
seq_name_df <- rbind(seq_name_df,data.frame(Sample_ID=r682$Sample_ID,
                                            Sequencing_Project="Knights_Project_068_Pool2",
                                            Filename=r682_filenames,
                                            Extraction="PowerSoil"))

# Run 068-3
r683 <- read.csv("data/magic/run068-3-taxatable-burst-capitalist.txt", row.names=1, header=1, sep="\t")
r683_filenames <- colnames(r683)
colnames(r683) <- gsub("[.]S.*[.]fa", "", colnames(r683))
r683 <- as.data.frame(t(r683))
r683$Sample_ID <- rownames(r683)
seq_name_df <- rbind(seq_name_df,data.frame(Sample_ID=r683$Sample_ID,
                                            Sequencing_Project="Knights_Project_068_Pool3",
                                            Filename=r683_filenames,
                                            Extraction="PowerSoil"))

# Run 068-4
r684 <- read.csv("data/magic/run068-4-taxatable-burst-capitalist.txt", row.names=1, header=1, sep="\t")
r684_filenames <- colnames(r684)
colnames(r684) <- gsub("[.]S.*[.]fa", "", colnames(r684))
r684 <- as.data.frame(t(r684))
r684$Sample_ID <- rownames(r684)
seq_name_df <- rbind(seq_name_df,data.frame(Sample_ID=r684$Sample_ID,
                                            Sequencing_Project="Knights_Project_068_Pool4",
                                            Filename=r684_filenames,
                                            Extraction="PowerSoil"))

# Run 069
r69 <- read.csv("data/magic/run069-taxatable-burst-capitalist.txt", row.names=1, header=1, sep="\t")
r69_filenames <- colnames(r69)
r69_sampids <- gsub("[.][P|M].*[.]fa", "", colnames(r69))
r69_ext <- gsub("[.]S.*[.]fa", "", r69_filenames)
r69_ext <- factor(gsub(".*[.]", "", r69_ext))
levels(r69_ext) <- c("MagAttract","PowerSoil","PowerSoilPro")
r69 <- as.data.frame(t(r69))
r69$Sample_ID <- r69_sampids
seq_name_df <- rbind(seq_name_df,data.frame(Sample_ID=r69$Sample_ID,
                                            Sequencing_Project="Knights_Project_069",
                                            Filename=r69_filenames,
                                            Extraction=r69_ext))

# Run 070
r70 <- read.csv("data/magic/run070-taxatable-burst-capitalist.txt", row.names=1, header=1, sep="\t")
r70_filenames <- colnames(r70)
colnames(r70) <- gsub("[.]S.*[.]fa", "", colnames(r70))
r70 <- as.data.frame(t(r70))
r70$Sample_ID <- rownames(r70)
seq_name_df <- rbind(seq_name_df,data.frame(Sample_ID=r70$Sample_ID,
                                            Sequencing_Project="Knights_Project_070",
                                            Filename=r70_filenames,
                                            Extraction="PowerSoil"))

# Run 076-1
r761 <- read.csv("data/magic/run076-1-taxatable-burst-capitalist.txt", row.names=1, header=1, sep="\t")
r761_filenames <- colnames(r761)
colnames(r761) <- gsub("[.]S.*[.]fa", "", colnames(r761))
r761 <- as.data.frame(t(r761))
r761$Sample_ID <- rownames(r761)
seq_name_df <- rbind(seq_name_df,data.frame(Sample_ID=r761$Sample_ID,
                                            Sequencing_Project="Knights_Project_076_Pool1",
                                            Filename=r761_filenames,
                                            Extraction="PowerSoil"))

# Run 076-2
r762 <- read.csv("data/magic/run076-2-taxatable-burst-capitalist.txt", row.names=1, header=1, sep="\t")
r762_filenames <- colnames(r762)
colnames(r762) <- gsub("[.]S.*[.]fa", "", colnames(r762))
r762 <- as.data.frame(t(r762))
r762$Sample_ID <- rownames(r762)
seq_name_df <- rbind(seq_name_df,data.frame(Sample_ID=r762$Sample_ID,
                                            Sequencing_Project="Knights_Project_076_Pool2",
                                            Filename=r762_filenames,
                                            Extraction="PowerSoil"))

# Run 076-3
r763 <- read.csv("data/magic/run076-3-taxatable-burst-capitalist.txt", row.names=1, header=1, sep="\t")
r763_filenames <- colnames(r763)
colnames(r763) <- gsub("[.]S.*[.]fa", "", colnames(r763))
r763 <- as.data.frame(t(r763))
r763$Sample_ID <- rownames(r763)
seq_name_df <- rbind(seq_name_df,data.frame(Sample_ID=r763$Sample_ID,
                                            Sequencing_Project="Knights_Project_076_Pool3",
                                            Filename=r763_filenames,
                                            Extraction="PowerSoil"))

# Run 077-1
r771 <- read.csv("data/magic/run077-1-taxatable-burst-capitalist.txt", row.names=1, header=1, sep="\t")
r771_filenames <- colnames(r771)
colnames(r771) <- gsub("[.]S.*[.]fa", "", colnames(r771))
r771 <- as.data.frame(t(r771))
r771$Sample_ID <- rownames(r771)
seq_name_df <- rbind(seq_name_df,data.frame(Sample_ID=r771$Sample_ID,
                                            Sequencing_Project="Knights_Project_077_Pool1",
                                            Filename=r771_filenames,
                                            Extraction="PowerSoilPro"))

# Run 077-2
r772 <- read.csv("data/magic/run077-2-taxatable-burst-capitalist.txt", row.names=1, header=1, sep="\t")
r772_filenames <- colnames(r772)
colnames(r772) <- gsub("[.]S.*[.]fa", "", colnames(r772))
r772 <- as.data.frame(t(r772))
r772$Sample_ID <- rownames(r772)
seq_name_df <- rbind(seq_name_df,data.frame(Sample_ID=r772$Sample_ID,
                                            Sequencing_Project="Knights_Project_077_Pool2",
                                            Filename=r772_filenames,
                                            Extraction="PowerSoilPro"))

# Run 077-3
r773 <- read.csv("data/magic/run077-3-taxatable-burst-capitalist.txt", row.names=1, header=1, sep="\t")
r773_filenames <- colnames(r773)
colnames(r773) <- gsub("[.]S.*[.]fa", "", colnames(r773))
r773 <- as.data.frame(t(r773))
r773$Sample_ID <- gsub("[.]PS", "", rownames(r773))
r773_ext <- rep("PowerSoilPro",nrow(r773))
r773_ext[grepl("PS",rownames(r773))] <- "PowerSoil"
seq_name_df <- rbind(seq_name_df,data.frame(Sample_ID=r773$Sample_ID,
                                            Sequencing_Project="Knights_Project_077_Pool3",
                                            Filename=r773_filenames,
                                            Extraction=r773_ext))



##### COMBINE & CLEAN! #####
# Combine into one matrix
magic_c <- rbind.fill(r46, r55, r681, r682, r683, r684, r69, r70, r761, r762, r763, r771, r772, r773)
identical(magic_c$Sample_ID, seq_name_df$Sample_ID)
magic_c$Sequencing_Project <- seq_name_df$Sequencing_Project
magic_c$Filename <- seq_name_df$Filename
magic_c$Extraction <- seq_name_df$Extraction
rm(seq_name_df, r46, r55, r681, r682, r683, r684, r69, r70, r761, r762, r763, r771, r772, r773)
## fill with zeros where NA values are
sample_names <- magic_c$Sample_ID
seq_proj <- magic_c$Sequencing_Project
fnames <- magic_c$Filename
extr_type <- magic_c$Extraction
magic_c <- magic_c[,-which(colnames(magic_c) %in% c("Sample_ID","Sequencing_Project","Filename","Extraction"))]
magic_c <- as.matrix(magic_c)
magic_c[is.na(magic_c)] <- 0
magic_c <- as.data.frame(magic_c)
magic_c$Sample_ID <- sample_names
magic_c$Sequencing_Project <- seq_proj
magic_c$Filename <- fnames
magic_c$Extraction <- extr_type


# Find blanks & controls (saving positive controls elsewhere)
blank_idx <- grep("blank", tolower(magic_c$Sample_ID))                     # 4 samples found
control_idx <- grep("control", tolower(magic_c$Sample_ID))                 # 1 positive control sample
control_idx <- c(control_idx, grep("positiv", tolower(magic_c$Sample_ID))) # 1 positive control sample
control_idx <- c(control_idx, grep("atcc", tolower(magic_c$Sample_ID)))    # 23 ATCC control samples
## save controls elsewhere
controls <- magic_c[control_idx,]
idx <- sapply(controls, class)=="numeric"
controls <- controls[,-which(colSums(controls[,idx]) == 0)]
controls <- controls[,c((ncol(controls)-3):ncol(controls),1:(ncol(controls)-4))]
write.table(controls, "data/magic/positive_control_samples.txt", sep="\t", row.names=F)
## remove blanks and controls from master list
magic_c <- magic_c[-c(control_idx,blank_idx),]


# Find duplicates (sample magic identifier or .H and .G suffixes)
magic_c$Sample_ID <- gsub(".H03", "", magic_c$Sample_ID)         # cleaning within project naming convention ()
magic_c$Sample_ID <- gsub(".G03", "", magic_c$Sample_ID)         # same but other name
duplicates <- unique(magic_c$Sample_ID[duplicated(magic_c$Sample_ID)])  # 88 samples duplicated, 26 from project 069
sum(magic_c$Sample_ID %in% duplicates)                                  # 195 samples total
## SELECTION CRITERIA:
##   (1) Want powersoil extraction, except for project 077
##   (2) Otherwise, figure out which of each sample is closest on average to other samples
tmp <- as.matrix(magic_c[,-which(colnames(magic_c) %in% c("Sequencing_Project","Sample_ID","Filename","Extraction"))])
tmp <- sweep(tmp, 1, rowSums(tmp), FUN="/")
aitch_all <- vegdist(tmp+0.000001, method="aitchison") # takes a few minutes, saved in next step
## applying selection criteria
remove_dup_idx <- c()
for (samp in duplicates) {
  # get set of duplicates for this sample
  idx <- which(samp == magic_c$Sample_ID)
  # criterion 1: powersoil extr. except run077 (powersoil pro)
  if ("Knights_Project_077_Pool3" %in% magic_c$Sequencing_Project[idx]) {
    good_dup <- idx[which("PowerSoilPro" %in% magic_c$Extraction[idx])]
  } else if (sum("PowerSoil" == magic_c$Extraction[idx]) == 1) {
    good_dup <- idx[which("PowerSoil" == magic_c$Extraction[idx])]
  }
  # criterion 2: closest fit to other data
  else { # (is.null(good_dup) | length(good_dup) > 1)
    good_dup <- idx[ which.min( rowSums(as.matrix(aitch_all)[idx,]) ) ]
  }
  # add to removal list
  bad_dup <- idx[which(!(idx == good_dup))]
  remove_dup_idx <- c(remove_dup_idx, bad_dup)
  # reset good_dup and bad_dup
  #good_dup <- NULL; bad_dup <- NULL;
}
## save the removed duplicates
dup_table <- magic_c
dup_table$Total_Counts <- rowSums(dup_table[,-which(colnames(dup_table) %in% c("Sample_ID","Sequencing_Project","Filename","Extraction"))])
dup_table$Remove <- 0; dup_table$Remove[remove_dup_idx] <- 1;
dup_table <- dup_table[dup_table$Sample_ID %in% duplicates,c("Sample_ID","Sequencing_Project","Filename","Extraction","Total_Counts","Remove")]
dup_table <- dup_table[order(dup_table$Sample_ID),]
write.table(dup_table, "data/magic/magic_jan2023_duplicates.txt", row.names=F, sep="\t")
## remove duplicates from master magic_c:
magic_c <- magic_c[-remove_dup_idx,]
sum(duplicated(magic_c$Sample_ID))

# Clean Up: magic_c becomes just the samples & otus
## cleaning aitch_all as well
aitch_all <- as.matrix(aitch_all)[-remove_dup_idx,-remove_dup_idx]
rownames(aitch_all) <- magic_c$Sample_ID; colnames(aitch_all) <- magic_c$Sample_ID;
aitch_all <- aitch_all[order(rownames(aitch_all)),order(colnames(aitch_all))]
aitch_all <- as.dist(aitch_all)
## cleaning magic_c
magic_c <- magic_c[order(magic_c$Sample_ID),]
rownames(magic_c) <- 1:nrow(magic_c)
sum(is.na(magic_c))
other_info <- magic_c[,c("Sample_ID","Sequencing_Project","Filename","Extraction")]
rownames(magic_c) <- magic_c$Sample_ID
magic_c <- magic_c[,-which(colnames(magic_c) %in% c("Sample_ID","Sequencing_Project","Filename","Extraction"))]
magic_c <- magic_c[,colSums(magic_c) > 0]
## note: magic_c is now 3442 samples x 13337 OTUs
## saving large data files for easy access later
identical(rownames(as.matrix(aitch_all)), rownames(magic_c))
save(aitch_all, file="data/magic/aitchison_dists_magic.Rdata")
save(magic_c, file="data/magic/magic_sample_dataframe.Rdata")
save(other_info,file="data/magic/magic_sample_shortinfo.Rdata")


##### Map File #####
# Read in table
meta_magic <- read.table("data/magic/clean_metafile.txt", sep="\t", header=T) # current metadata: 4864 x 654
rownames(meta_magic) <- meta_magic$Sample_ID
sum(rownames(magic_c) %in% meta_magic$Sample_ID) # only 3078 of 3442 matching...
# Only keep samples matching metadata
magic_c <- magic_c[,(colnames(magic_c) %in% meta_magic$Sample_ID)]
seq_name_df <- seq_name_df[(seq_name_df$Sample_ID %in% meta_magic$Sample_ID),]



##### Normalize and Reduce Sizes #####
# Check out count totals per smaple... range is crazy (1k --> 6m)
head(sort(colSums(magic_c), decreasing=F),100)
tail(sort(colSums(magic_c), decreasing=F),20)
# Drop everything with less than 140k (drops 70ish samples)
low_yield <- which(colSums(magic_c) < 140000)
magic_c <- magic_c[,-low_yield] # 1750 to 1679 samples
seq_name_df <- seq_name_df[-low_yield,]
# Drop non-existent taxa in table
magic_c <- magic_c[(rowSums(magic_c) > 0),] # drops ~200 taxa non-existent
magic_c <- magic_c[(rowSums(magic_c) > 1),] # drops ~1410 singletons
# Create a normalized counts table
magic_n <- sweep(magic_c, 2, colSums(magic_c), FUN="/")
# Drop samples from metadata (match order to magic_c)
meta_magic <- meta_magic[colnames(magic_c),]
# Append sequencing project names to metadata
meta_magic$Sequencing_Run <- seq_name_df$Sequencing_Project


##### Write Out: Counts, Norm, Map #####
write.table(magic_c, "data/magic/mycleaned_magic_counts.txt", sep="\t", row.names=T, col.names=T)
write.table(magic_n, "data/magic/mycleaned_magic_norm.txt", sep="\t", row.names=T, col.names=T)
write.table(meta_magic, "data/magic/mycleaned_magic_metadata.txt", sep="\t", row.names=T, col.names=T)


