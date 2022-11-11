# Written Preliminary Exam Figures
library(vegan)
library(ggplot2)
library(GGally)
library(gridExtra)
library(grid)
library(cowplot)
library(gridGraphics)
set.seed(25)

##### LOADING DATASETS #####
# Immigration Microbiome Project
meta_imp <- read.table("data/imp/map.txt", sep="\t", header=T)
rownames(meta_imp) <- meta_imp$SampleID
rare_imp <- read.delim("data/imp/final_otu.txt", header=T, row=1) # OTU table rarefied to ~10mil counts
## refine to no first gen samples (these were all over the place)
meta_imp <- meta_imp[!(grepl("1st", meta_imp$Sample.Group)),]     # this omits 137 Hmong1st and 236 Karen1st (278 total left)
meta_imp <- meta_imp[meta_imp$SampleID %in% colnames(rare_imp),]  # ignore 9 samples not matching otu table (269 total left)
rare_imp <- rare_imp[,meta_imp$SampleID]
rare_imp <- rare_imp[rowSums(rare_imp) > 0,]                      # restrict to non-empty OTU rows (6803 remian)

# Soil dataset (88 soils)
meta_soil <- read.table("data/soil/clean_map.txt", header=T, sep="\t")
rownames(meta_soil) <- meta_soil$SampleID
soil_n <- read.table("data/soil/44766_clean_otus_norm.txt", row=1, header=T, sep="\t")
#identical(meta_soil$SampleID, colnames(soil_n))

# Guerrero Negro dataset
meta_gn <- 
  
# Cecum dataset

# MAGIC dataset
magic_c <- read.table("data/magic/mycleaned_magic_counts.txt", row=1, header=T, sep="\t")
magic_n <- read.table("data/magic/mycleaned_magic_norm.txt", row=1, header=T, sep="\t")
meta_magic <- read.table("data/magic/mycleaned_magic_metadata.txt", row=1, header=T, sep="\t")

##### FIGURE 1: Distances Comparison #####
# Comparing: B-C, Euclidean, Jaccard, uUniFrac, Aitchison, rAitchison

# Compute distances:
## Bray-Curtis
bray <- vegdist(t(rare_imp), method="bray")
## Euclidean Distance
euclid <- vegdist(t(rare_imp), method="euclidean")
## Jaccard
jacc <- vegdist(t(rare_imp), method="jaccard", binary = T)
## Unweighted UniFrac (load from IMP analyses GitHub)
u_unifrac <- read.delim('data/imp/unweighted_unifrac_dm.txt', sep='\t', row=1)
u_unifrac <- as.dist(u_unifrac[colnames(rare_imp),colnames(rare_imp)])
## Aitchison Distance (Euclidean distance in CLR transformed data, requires non-zero data)
aitch <- vegdist(t(rare_imp), method="aitchison", pseudocount=1)
## Robust Aitchison
deicode <- read.table("data/imp/deicode_out_imp/distance-matrix.tsv", sep="\t", header=T, row=1)
robust_a <- as.dist(deicode)

# Pairwise Scatterplots Visual
## color by within group or outside group distance
group_dist <- dist(as.numeric(as.factor(meta_imp$Sample.Group)))
group_dist <- as.factor(c(group_dist) < 1) # T/F list of w/in group or not
levels(group_dist) <- c("out", "in")
## create data frame of distances for pairwise graphing
dist_df <- data.frame(BrayCurtis=c(bray), Euclidean=c(euclid), Jaccard=c(jacc), 
                      uUniFrac=c(u_unifrac), Aitchison=c(aitch), rAitchison=c(robust_a),
                      Color_Within_Group=c(group_dist))
## CUSTOM FUNCTION for upper triangle in ggpairs
custom_corr <- function (data, mapping, ...) {
  x <- eval_data_col(data, mapping$x)
  y <- eval_data_col(data, mapping$y)
  my_corr <- cor(x, y, method="pearson")
  adjx <- mean(x)
  adjy <- mean(y)
  p <- ggplot(data) +
    geom_blank(mapping) +
    geom_smooth(mapping, method="lm", linetype="dashed") +
    annotation_custom(grid::textGrob(paste0("corr: ", round(my_corr, 2)), gp=gpar(fontsize=16)),
                      xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)
  p
}
## ggplot for publication
pairwise_plt <- ggpairs(dist_df, columns=1:6, ggplot2::aes(col=Color_Within_Group, fill=Color_Within_Group, alpha=0.4),
                        lower = list(continuous = wrap("points", alpha=0.1, pch=16)),
                        upper = list(continuous = custom_corr) ) + 
  scale_color_manual(values=c("#9d02d7", "#ffb14e")) +
  scale_fill_manual(values=c("#9d02d7", "#ffb14e")) +
  theme_bw() +
  theme(title = element_text(face="bold"),
        panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        text = element_text(family="Helvetica", size=16))
pairwise_plt # exported as 1000x1000



##### FIGURE 2: Ordination Comparison (unw UniFrac) #####
# Custom colors and shapes to match the IMP paper
imp_col <- c("#e38cbd", "#7e7e7e", "#d65ca0", "#4ea99b")
names(imp_col) <- c("Hmong2nd", "Control", "HmongThai", "KarenThai")
imp_shp <- c(17, 17, 16, 16)
names(imp_shp) <- c("Hmong2nd", "Control", "HmongThai", "KarenThai")
# Custom 2D plot for IMP
plot_2d <- function (dims, meta_df, title_text="", eig=NULL, legendpos="none") {
  df_2d <- as.data.frame(dims)
  colnames(df_2d) <- paste0("Dim", 1:ncol(df_2d))
  df_2d$imp_group <- meta_df$Sample.Group
  df_2d$sampleid <- rownames(df_2d)
  p <- ggplot(df_2d, aes(x=Dim1, y=Dim2)) +
    geom_point(size=4, alpha=0.7, aes(shape=imp_group, col=imp_group)) +
    scale_color_manual(values = imp_col[unique(df_2d$imp_group)]) +
    scale_shape_manual(values = imp_shp[unique(df_2d$imp_group)]) +
    theme_bw() +
    theme(legend.position = legendpos,
          legend.box.background = element_rect(colour = "black"),
          legend.spacing.y = unit(0.1, 'cm'),
          legend.margin = margin(0,0,0,0,'cm'),
          legend.title = element_blank(),
          panel.grid = element_blank(),
          text = element_text(family="Helvetica", size=16))
  if (is.null(eig)) {
    p <- p + labs(title=title_text, x="Dimension 1", y="Dimension 2", col="", shape="")
  } else {
    eig <- eig/sum(eig)
    var <- round(eig*100, 1)
    p <- p + labs(title=title_text, x=paste0("Dimension 1 [",var[1],"%]"), y=paste0("Dimension 2 [",var[2],"%]"), col="", shape="")
  }
  return(p)
}
## PCoA
pc <- cmdscale(u_unifrac, k=3, eig=T)
tmp <- pc$points
tmp[,1] <- tmp[,1]*-1
tmp[,2] <- tmp[,2]*-1
pc_plt <- plot_2d(tmp, meta_imp, "PCoA", pc$eig)
## NMDS
nmds <- metaMDS(u_unifrac, k=3, try=20, trymax=50, maxit=500)
tmp <- nmds$points
tmp[,1] <- tmp[,1]*-1
tmp[,2] <- tmp[,2]*-1
nmds_plt <- plot_2d(tmp, meta_imp, "NMDS")
## t-SNE
tsne <- Rtsne(u_unifrac, is_distance=T, dims=3)
tmp <- tsne$Y
tmp[,2] <- tmp[,2]*-1
tsne_plt <- plot_2d(tmp, meta_imp, "t-SNE")
## Constrained (RDA)
### compute rda (try using all samples though)
rda_imp <- rda(t(rare_imp) ~ Age + BMI + Ethnicity, data=meta_imp)
### my plot of rda output
rda_plt <- ~{
  par(mar=c(2,2,1,0.5), oma=c(1,1,1,0))
  plot(rda_imp, type="n")
  mtext(side=3, line=0.5, adj=0, cex=1.5, "RDA (constrained)")
  points(rda_imp, pch=imp_shp[meta_imp$Sample.Group], col=alpha(imp_col[meta_imp$Sample.Group], 0.4), cex=2)
  text(rda_imp, dis="cn", cex=1)
  legend(-2250, -250, legend=names(imp_col), col=imp_col, pch=imp_shp, cex=1, y.intersp=1, bty="n")
}
## Plot them all together
plot_grid(pc_plt, nmds_plt, tsne_plt, rda_plt, ncol = 2, align="hv", axis="lrtb") # exported as 800x800
par(mar=c(5,4,4,2), oma=c(3,3,3,3))



##### FIGURE 3: Examples of Arches (hypothesis) #####




