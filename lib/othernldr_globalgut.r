# othernldr_humangut.r
# Testing Variety of Other NLDR Methods on HMP data
# Knights Lab - University of Minnesota
# February 2020
# usage : source('othernldr_humangut.r')

##### Set Up #####
current_dir <- "/project/flatiron2/suzie/detrending/humangut/"
library(Rtsne)
library(diffusionMap)
library(lle)
source('/project/flatiron2/suzie/detrending/fake/lgd_source.r')

##### Load Data #####
# Cleaned datasets
dat_hg <- read.delim(paste0(current_dir, "processed_data/421_clean_otus.txt"),
                     header=T, sep="\t", quote="\"", row.names=1, comment.char="", as.is=T)
meta_hg <- read.delim(paste0(current_dir, "processed_data/clean_map.txt"),
                      header=T, sep="\t", quote="\"", comment.char="", as.is=T)
pop_cols <- c("red", "blue", "green") # malawians, usa, venezuelans
# Bring in unweighted UniFrac matrix
d_hg <- read.delim(paste0(current_dir, "processed_data/clean_uw_unifrac_dist.txt"),
                   header=T, row.names=1, sep="\t", as.is=T)
# Convert to distance object
d_hg <- as.dist(as.matrix(d_hg))
# Split into per population
malawi <- as.dist(as.matrix(d_hg)[meta_hg$geo_loc_name == "Malawi", meta_hg$geo_loc_name == "Malawi"])
usa <- as.dist(as.matrix(d_hg)[meta_hg$geo_loc_name == "USA", meta_hg$geo_loc_name == "USA"])
venezuela <- as.dist(as.matrix(d_hg)[meta_hg$geo_loc_name == "Venezuela", meta_hg$geo_loc_name == "Venezuela"])

##### PCA #####
pc_cmd <- cmdscale(d_hg)
pc_df <- as.data.frame(pc_cmd)
colnames(pc_df) <- paste0("Dim", 1:ncol(pc_df))
pc_df <- pc_df[meta_hg$SampleID,]
pc_df <-cbind(pc_df, Population=factor(meta_hg$geo_loc_name), Age=meta_hg$age)
# Populations
ggplot(pc_df, aes(x=Dim1, y=Dim2, color=Population)) +
    geom_point(size = 4, alpha=0.8) +
    scale_color_manual(values = pop_cols) +
    labs(title = "PCA: Human Microbiome Project", x = "Dimension 1", y = "Dimension 2") +
    theme_classic() + NULL
# Age
ggplot(pc_df, aes(x=Dim1, y=Dim2, color=Age)) +
  geom_point(size = 4, alpha=0.8) +
  scale_color_viridis_c() +
  labs(title = "PCA: Human Microbiome Project", x = "Dimension 1", y = "Dimension 2") +
  theme_classic() + NULL
# Split by population
pc_m <- as.data.frame(cmdscale(malawi))
pc_u <- as.data.frame(cmdscale(usa))
pc_v <- as.data.frame(cmdscale(venezuela))
colnames(pc_m) <- paste0("Dim", 1:ncol(pc_m))
colnames(pc_u) <- paste0("Dim", 1:ncol(pc_u))
colnames(pc_v) <- paste0("Dim", 1:ncol(pc_v))
pc_m$Age <- meta_hg$age[meta_hg$geo_loc_name == "Malawi"]
pc_u$Age <- meta_hg$age[meta_hg$geo_loc_name == "USA"]
pc_v$Age <- meta_hg$age[meta_hg$geo_loc_name == "Venezuela"]
ggplot(pc_m, aes(x=Dim1, y=Dim2, color=Age)) +
  geom_point(size = 4, alpha=0.8) +
  scale_color_viridis_c() +
  labs(title = "PCA: HMP Malawi", x = "Dimension 1", y = "Dimension 2") +
  theme_classic() + NULL
ggplot(pc_u, aes(x=Dim1, y=Dim2, color=Age)) +
  geom_point(size = 4, alpha=0.8) +
  scale_color_viridis_c() +
  labs(title = "PCA: HMP USA", x = "Dimension 1", y = "Dimension 2") +
  theme_classic() + NULL
ggplot(pc_v, aes(x=Dim1, y=Dim2, color=Age)) +
  geom_point(size = 4, alpha=0.8) +
  scale_color_viridis_c() +
  labs(title = "PCA: HMP Venzuela", x = "Dimension 1", y = "Dimension 2") +
  theme_classic() + NULL

##### t-SNE #####
tsne <- Rtsne(d_hg)
tsne_df <- as.data.frame(tsne$Y)
rownames(tsne_df) <- rownames(as.matrix(d_hg))
colnames(tsne_df) <- paste0("Dim", 1:ncol(tsne_df))
tsne_df <- tsne_df[meta_hg$SampleID,]
tsne_df <- cbind(tsne_df, Population=factor(meta_hg$geo_loc_name), Age=meta_hg$age)
ggplot(tsne_df, aes(x=Dim1, y=Dim2, color=Population)) +
  geom_point(size = 4, alpha=0.8) +
  scale_color_manual(values = pop_cols) +
  labs(title = "t-SNE: Human Microbiome Project", x = "Dimension 1", y = "Dimension 2") +
  theme_classic() + NULL
ggplot(tsne_df, aes(x=Dim1, y=Dim2, color=Age)) +
  geom_point(size = 4, alpha=0.8) +
  scale_color_viridis_c() +
  labs(title = "t-SNE: Human Microbiome Project", x = "Dimension 1", y = "Dimension 2") +
  theme_classic() + NULL

##### Diff Map #####
dmap <- diffuse(d_hg, neigen = 3) # default epsilon is median dist to 0.01*n nearest neighbor
dmap_df <- as.data.frame(dmap$X)
rownames(dmap_df) <- rownames(as.matrix(d_hg))
colnames(dmap_df) <- paste0("Dim", 1:ncol(dmap_df))
dmap_df <- dmap_df[meta_hg$SampleID,]
dmap_df <- cbind(dmap_df, Population=meta_hg$geo_loc_name, Age=meta_hg$age)
ggplot(dmap_df, aes(x=Dim1, y=Dim2, color=Population)) +
  geom_point(size = 4, alpha=0.8) +
  scale_color_manual(values = pop_cols) +
  labs(title = "DiffMap: Human Microbiome Project", x = "Dimension 1", y = "Dimension 2") +
  theme_classic() + NULL
ggplot(dmap_df, aes(x=Dim1, y=Dim2, color=Age)) +
  geom_point(size = 4, alpha=0.8) +
  scale_color_viridis_c() +
  labs(title = "DiffMap: Human Microbiome Project", x = "Dimension 1", y = "Dimension 2") +
  theme_classic() + NULL


##### LGD #####
# PCA of transformation
lgd_hg <- lg.dist(d_hg, neighborhood.size = 15, weighted=TRUE)
lgd_pc <- cmdscale(lgd_hg)
lgd_pc <- as.data.frame(lgd_pc)
colnames(lgd_pc) <- paste0("Dim", 1:ncol(lgd_pc))
lgd_pc <- lgd_pc[meta_hg$SampleID,]
lgd_pc <-cbind(lgd_pc, Population=factor(meta_hg$geo_loc_name), Age=meta_hg$age)
ggplot(lgd_pc, aes(x=Dim1, y=Dim2, color=Population)) +
  geom_point(size = 4, alpha=0.8) +
  scale_color_manual(values = pop_cols) +
  labs(title = "LGD & PCA: Human Microbiome Project", x = "Dimension 1", y = "Dimension 2") +
  theme_classic() + NULL
ggplot(lgd_pc, aes(x=Dim1, y=Dim2, color=Age)) +
  geom_point(size = 4, alpha=0.8) +
  scale_color_viridis_c() +
  labs(title = "PCA: Human Microbiome Project", x = "Dimension 1", y = "Dimension 2") +
  theme_classic() + NULL
# Within populations
lgd_m <- lg.dist(malawi, neighborhood.size = 7, weighted=TRUE)
lgd_u <- lg.dist(usa, neighborhood.size = 7, weighted=TRUE)
lgd_v <- lg.dist(venezuela, neighborhood.size = 5, weighted=TRUE)
lgd_pcm <- as.data.frame(cmdscale(malawi))
lgd_pcu <- as.data.frame(cmdscale(usa))
lgd_pcv <- as.data.frame(cmdscale(venezuela))
colnames(lgd_pcm) <- paste0("Dim", 1:ncol(lgd_pcm))
colnames(lgd_pcu) <- paste0("Dim", 1:ncol(lgd_pcu))
colnames(lgd_pcv) <- paste0("Dim", 1:ncol(lgd_pcv))
lgd_pcm$Age <- meta_hg$age[meta_hg$geo_loc_name == "Malawi"]
lgd_pcu$Age <- meta_hg$age[meta_hg$geo_loc_name == "USA"]
lgd_pcv$Age <- meta_hg$age[meta_hg$geo_loc_name == "Venezuela"]
ggplot(lgd_pcm, aes(x=Dim1, y=Dim2, color=Age)) +
  geom_point(size = 4, alpha=0.8) +
  scale_color_viridis_c() +
  labs(title = "LGD & PCA: HMP Malawi", x = "Dimension 1", y = "Dimension 2") +
  theme_classic() + NULL
ggplot(lgd_pcu, aes(x=Dim1, y=Dim2, color=Age)) +
  geom_point(size = 4, alpha=0.8) +
  scale_color_viridis_c() +
  labs(title = "LGD & PCA: HMP USA", x = "Dimension 1", y = "Dimension 2") +
  theme_classic() + NULL
ggplot(lgd_pcv, aes(x=Dim1, y=Dim2, color=Age)) +
  geom_point(size = 4, alpha=0.8) +
  scale_color_viridis_c() +
  labs(title = "LGD & PCA: HMP Venezuela", x = "Dimension 1", y = "Dimension 2") +
  theme_classic() + NULL


