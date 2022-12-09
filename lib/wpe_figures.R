# Written Preliminary Exam Figures
# November 2022
library(vegan)
library(Rtsne)
library(ggplot2)
library(GGally)
library(gridExtra)
library(grid)
library(cowplot)
library(gridGraphics)
source('lib/lgd_source.r')
set.seed(25)

##### HELPER FUNCTIONS #####
# In PCoA plots, get the percent variance
calc.perc.var <- function (eigen, dimension) {
  percents <- round((eigen/sum(eigen))*100, 1) # rounds percentages to one decimal place
  return(percents[dimension])
}



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
meta_gn <- read.table("data/guerrero_negro/clean_map.txt", header=T, sep="\t")
rownames(meta_gn) <- meta_gn$SampleID
meta_gn$depth_mm <- c("0-1","1-2","2-3","3-4","4-5","5-6","6-10","10-22","22-34")[as.factor(meta_gn$depth)]
gn_n <- read.table("data/guerrero_negro/47908_clean_otus_norm.txt", row=1, header=T, sep="\t")
#identical(rownames(meta_gn), colnames(gn_n))

# Cecum dataset
cecum_d <- readRDS("data/cecum/cecum_bray_dist.rds") # Bray-Curtis distances
meta_cecum <- read.table("data/cecum/cecum_meta.csv", sep=",", row=1, header=T)
#identical(rownames(meta_cecum), rownames(as.matrix(cecum_d)))
## get rid of 7 samples removed for the paper
remove <- c("WPC2FHBD15B06C","WPC2FHBD18B07C","WPC2FHBD15B05C","WPC2FD15B04C",
            "WPC2FD15B05C","WPC2FD08B08C","WPC2FD01B10C")
remove_idx <- which(meta_cecum$SeqID %in% remove)
cecum_d <- as.dist(as.matrix(cecum_d)[-remove_idx,-remove_idx])
meta_cecum <- meta_cecum[-remove_idx,]

# MAGIC dataset
magic_c <- read.table("data/magic/mycleaned_magic_counts.txt", row=1, header=T, sep="\t")
magic_n <- read.table("data/magic/mycleaned_magic_norm.txt", row=1, header=T, sep="\t")
meta_magic <- read.table("data/magic/mycleaned_magic_metadata.txt", row=1, header=T, sep="\t")
#identical(meta_magic$Sample_ID,colnames(magic_c))
#identical(meta_magic$Sample_ID,colnames(magic_n))
## refine dataset
meta_magic$age_grp <- cut(meta_magic$Timeline_Weeks, breaks = c(2,5,13,25,37,49,61,73,100)) # excluding first 2 weeks
remove_idx <- which(is.na(meta_magic$age_grp) | is.na(meta_magic$currentfeed_bf))           # refine to age set and breastfeeding info
meta_magic <- meta_magic[-remove_idx,]
magic_c <- magic_c[,-remove_idx]
magic_n <- magic_n[,-remove_idx]

# Simulated Dataset
sim_n <- read.table("data/sim_gradient/fake_rel_abun_g1000_sd50_n50.txt", row=1, header=T, sep="\t")
meta_sim <- data.frame(SampleID=rownames(sim_n), noise=c(rep("n",50),rep("y",50)),
                       gradient=rep(round(seq(1,1000,length.out=50)),2))


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
    annotation_custom(grid::textGrob(paste0("corr: ", round(my_corr, 2)), gp=gpar(fontsize=12)),
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
## SAVE PNG
### note: width of page w/ 1in margins is 6.5in*96dpi = 624px
###       height of page w/ 1in margins is 9in*96dpi = 864px
# png("wpe_figures/figure1_pairwiseDists.png",width=624, height=624, units="px", res=96)
# pairwise_plt
# dev.off()



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
          legend.box.background = element_rect(colour = "white"),
          legend.spacing.y = unit(0.1, 'cm'),
          legend.margin = margin(0,0,0,0,'cm'),
          legend.title = element_blank(),
          panel.grid = element_blank(),
          legend.key.size = unit(0.5, "lines"),
          text = element_text(family="Helvetica", size=16),
          axis.title = element_text(family="Helvetica",size=12),
          legend.text = element_text(family="Helvetica",size=10))
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
tsne_plt <- plot_2d(tmp, meta_imp, "t-SNE", legendpos=c(0.78,0.2))

## Constrained (RDA)
### compute rda (try using all samples though)
rda_imp <- rda(t(rare_imp) ~ Age + BMI + Ethnicity, data=meta_imp)
### my plot of rda output
rda_plt <- ~{
  par(mar=c(2,2,1,0.5), oma=c(1,1,1,0), mgp=c(3,0.5,0))
  plot(rda_imp, type="n")
  mtext(side=3, line=0.5, adj=0, cex=1.5, "RDA (constrained)")
  points(rda_imp, pch=imp_shp[meta_imp$Sample.Group], col=alpha(imp_col[meta_imp$Sample.Group], 0.4), cex=2)
  text(rda_imp, dis="cn", cex=0.9)
  #legend(-1700, -600, legend=names(imp_col), col=imp_col, pch=imp_shp, cex=1, y.intersp=1, bty="n")
}

## Plot them all together (exported as 900x800)
plot_grid(pc_plt, nmds_plt, tsne_plt, rda_plt, ncol = 2, align="hv", axis="lrtb")
## SAVE PNG
### note: width of page w/ 1in margins is 6.5in*96dpi = 624px
###       height of page w/ 1in margins is 9in*96dpi = 864px
# png("wpe_figures/figure2_ordination.png",width=624, height=554, units="px", res=96)
# plot_grid(pc_plt, nmds_plt, tsne_plt, rda_plt, ncol = 2, align="hv", axis="lrtb")
# dev.off()



##### FIGURE 3: Examples of Arches (hypothesis) #####
dev.off() # resets par
par(mgp=c(2, 0.5, 0))

# 88 soils
## colors & shapes from publication
soil_col1 <- c("black", "#002dd5", "#1951f4", "#6e8df7", "#aabcf9", "black") # colors (borders)
soil_bg1 <- c("black", "#002dd5", "#1951f4", "#6e8df7", "#aabcf9", "white") # fills (background)
soil_shp1 <- c(22, 25, 24, 21, 23, 23) # shapes
ord <- cut(meta_soil$ph, breaks = c(0,4:8,14)) # pH groups (<4, 4-5, 5-6, 6-7, 7-8, >8)
soil_col <- soil_col1[ord]; soil_bg <- soil_bg1[ord]; soil_shp <- soil_shp1[ord];
## calculate distances (Jaccard) and PCoA
soil_dat <- vegdist(t(soil_n), method="jaccard")
soil_pc <- cmdscale(soil_dat, k=2, eig=T)
soil_pc$points[,2] <- soil_pc$points[,2] * -1
## plot PCoA (exported as 600x600)
## SAVE PNG
### note: width of page w/ 1in margins is 6.5in*96dpi = 624px
###       height of page w/ 1in margins is 9in*96dpi = 864px
# png("wpe_figures/figure3_arch_cases/figure3_soil_pcoajaccard.png",width=500, height=500, units="px", res=96)
# par(mgp=c(2, 0.5, 0))
plot(soil_pc$points,
     xlim=range(soil_pc$points)+c(-0.1, 0.1), ylim=range(soil_pc$points)+c(-0.1, 0.1),
     col=soil_col, bg=soil_bg, pch=soil_shp, cex=2, lwd=3,
     xlab=paste0("PC1 [", calc.perc.var(soil_pc$eig,1), "%]"),
     ylab=paste0("PC2 [", calc.perc.var(soil_pc$eig,2), "%]"), cex.lab=2)
title("88 Soils", adj=0, cex.main=2)
legend("topright", inset=0.02, legend=c("pH", "<4","4-5","5-6","6-7","7-8",">8"),
       col=c("white",soil_col1), pt.bg=c("white",soil_bg1), pch=c(16,soil_shp1),
       cex=1, y.intersp=0.8, bty="o")
# dev.off()
cor(soil_pc$points[,1], meta_soil$ph, method="spearman") # 0.92 PC1 with gradient
cor.test(soil_pc$points[,1], meta_soil$ph, method="spearman") # p < 0.001
## statistics done in publication
anosim(soil_dat, soil_col, permutations=999) # ANOSIM with pH groups: R=0.429, p=0.001
mantel(soil_dat, vegdist(meta_soil$ph, method="euclidean"), method="spearman") # Mantel test: R=0.733, p=0.001


# Cecum dataset
## colors and shapes from publication
cecum_col1 <- c("#fd2beb","#8e0013","#fa4d3f","#f99314","#00ae03","#57f89f","#0e2dfc","#6033f7")
cecum_shp1 <- c(16, 17)
ord_c <- as.factor(meta_cecum$Collection)
ord_s <- as.factor(meta_cecum$System)
cecum_col = cecum_col1[ord_c]; cecum_shp = cecum_shp1[ord_s]
## calculate PCoA, already have Bray-Curtis distances
cecum_pc <- cmdscale(cecum_d, k=2, eig=T)
cecum_pc$points[,2] <- cecum_pc$points[,2] * -1 # flipping y axis to make traditional arch shape
## plot PCoA (exported as 600x600)
## SAVE PNG
### note: width of page w/ 1in margins is 6.5in*96dpi = 624px
###       height of page w/ 1in margins is 9in*96dpi = 864px
# png("wpe_figures/figure3_arch_cases/figure3_cecum_pcoabraycurtis.png",width=500, height=500, units="px", res=96)
# par(mgp=c(2, 0.5, 0))
plot(cecum_pc$points,
     xlim=range(cecum_pc$points)+c(-0.1,0.1), ylim=range(cecum_pc$points)+c(-0.1,0.1),
     col=alpha(cecum_col,0.8), pch=cecum_shp, cex=2,
     xlab=paste0("PC1 [", calc.perc.var(cecum_pc$eig,1), "%]"),
     ylab=paste0("PC2 [", calc.perc.var(cecum_pc$eig,2), "%]"), cex.lab=2)
title("Turkey Cecum", adj=0, cex.main=2)
legend("bottomright", inset=0.04, legend=c("Day","1","4","8","10","15","18","22","29","","System","Pen","Hatch"),
       col=c("white",cecum_col1,"white","white","black","black"), pch=c(rep(16,9),16,16,16,17), 
       cex=1, y.intersp=0.8, x.intersp=0.5, text.width=0.15, bty="o", ncol=3)
# dev.off()
cor(cecum_pc$points[,1], as.numeric(gsub("D", "", meta_cecum$Collection)), method="spearman") # 0.736 PC1 with gradient
cor.test(cecum_pc$points[,1], as.numeric(gsub("D", "", meta_cecum$Collection)), method="spearman") # p < 0.001
## statistics from the publication
adonis2(cecum_d ~ Collection, meta_cecum, permutations=999)
mantel(as.dist(as.matrix(cecum_d)[meta_cecum$System=="Hatch_Brood",meta_cecum$System=="Hatch_Brood"]),
       as.dist(as.matrix(cecum_d)[meta_cecum$System=="Pen",meta_cecum$System=="Pen"]), method="spearman") # r = 0.473, p = 0.001


# Guerrero Negro
## colors and shapes from publication
depth_cols1 <- c("#F2182F", "#FF6A4F", "#FFA96A", "#FFDE95", "#FCF7C2",
                 "#D6F0F6", "#88D9E8", "#23AED1", "#0076B4")
ord_d <- as.factor(meta_gn$depth)
depth_cols <- depth_cols1[ord_d]
## calculate distances (Jaccard since no tree available), PCoA
gn_d <- vegdist(t(gn_n), method="jaccard")
gn_pc <- cmdscale(gn_d, k=2, eig=T)
## plot PCoA (exported as 600x600)
## SAVE PNG
### note: width of page w/ 1in margins is 6.5in*96dpi = 624px
###       height of page w/ 1in margins is 9in*96dpi = 864px
# png("wpe_figures/figure3_arch_cases/figure3_microbialmat_pcoajaccard.png", width=500, height=500, units="px", res=96)
# par(mgp=c(2, 0.5, 0))
plot(gn_pc$points,
     xlim=range(gn_pc$points)+c(-0.15,0.15), ylim=range(gn_pc$points)+c(-0.15,0.15),
     col=alpha(depth_cols,0.8), pch=16, cex=3,
     xlab=paste0("PC1 [", calc.perc.var(gn_pc$eig,1), "%]"),
     ylab=paste0("PC2 [", calc.perc.var(gn_pc$eig,2), "%]"), cex.lab=2)
title("Microbial Mat", adj=0, cex.main=2)
legend("bottomleft", inset=0.04, legend=c("Depth",unique(meta_gn$depth_mm)),
       col=c("white",depth_cols1), pch=16, cex=1, y.intersp=0.8, bty="o", ncol=3)
# dev.off()
cor(gn_pc$points[,1], meta_gn$end_depth, method="spearman") # 0.941 PC1 with gradient
cor.test(gn_pc$points[,1], meta_gn$end_depth, method="spearman") # p < 0.001


# MAGIC infant gut study
## colors for age categories
magic_cols1 <- c("#ffb14e","#fa7b67","#cd34b5","#9d02d7","#0000ff","#009aa7","#40c557","#036003")
magic_shp1 <- c(17,16)
ord_c <- as.factor(meta_magic$age_grp)
ord_s <- as.factor(meta_magic$currentfeed_bf)
magic_cols <- magic_cols1[ord_c]; magic_shp <- magic_shp1[ord_s];
## compute distances (Bray-Curtis) & PCoA
magic_d <- vegdist(t(magic_n), method="bray")
magic_pc <- cmdscale(magic_d, k=2, eig=T)
## plot PCoA (exported as 600x600)
## SAVE PNG
### note: width of page w/ 1in margins is 6.5in*96dpi = 624px
###       height of page w/ 1in margins is 9in*96dpi = 864px
# png("wpe_figures/figure3_arch_cases/figure3_infantgut_pcoabraycurtis.png", width=500, height=500, units="px", res=96)
# par(mgp=c(2, 0.5, 0))
plot(magic_pc$points,
     xlim=range(magic_pc$points)+c(-0.1,0.1), ylim=range(magic_pc$points)+c(0,0.1),
     col=alpha(magic_cols,0.6), pch=magic_shp, cex=1.2,
     xlab=paste0("PC1 [", calc.perc.var(magic_pc$eig,1), "%]"),
     ylab=paste0("PC2 [", calc.perc.var(magic_pc$eig,2), "%]"), cex.lab=2)
title("Infant Gut", adj=0, cex.main=2)
legend("topright", inset=0.02,
       legend=c("Age","1mo","3mo","6mo","9mo","1yr","15mo","18mo","2yr","","BF","Yes","No"),
       col=c("white",magic_cols1,"white","white","black","black"), pch=c(rep(16,12),17),
       cex=0.85, y.intersp=0.8, bty="o")
# dev.off()
cor(-1*magic_pc$points[,1], meta_magic$age, method="spearman") # 0.478 PC1 with age gradient
cor.test(-1*magic_pc$points[,1], meta_magic$age, method="spearman") # p < 0.001
## Plot again but with jaccard (just to see - looks extra strange (2 arches?? idk))
# magic_d2 <- vegdist(t(magic_n),method="jaccard")
# magic_pc2 <- cmdscale(magic_d2, k=2, eig=T)
# plot(magic_pc2$points,
#      xlim=range(magic_pc2$points)+c(-0.1,0.1), ylim=range(magic_pc2$points)+c(0,0.1),
#      col=alpha(magic_cols,0.7), pch=magic_shp, cex=1,
#      xlab=paste0("PC1 [", calc.perc.var(magic_pc2$eig,1), "%]"),
#      ylab=paste0("PC2 [", calc.perc.var(magic_pc2$eig,2), "%]"), cex.lab=2)
# title("Infant Gut", adj=0, cex.main=2)
# legend("topright", inset=0.04,
#        legend=c("Age","1mo","3mo","6mo","9mo","1yr","15mo","18mo","2yr","","BF","Yes","No"),
#        col=c("white",magic_cols1,"white","white","black","black"), pch=c(rep(16,12),17),
#        cex=1, y.intersp=0.5, bty="o")



##### FIGURE 4: Algortihm Overview #####
dev.off() # reset the par settings
par(mar=c(3,3,1,1), mgp=c(1, 0.5, 0))

# Small simulated dataset (pick 20 samples)
## note: grab from older gradient with 100 samples makes a neater plot for small example
mini_sim <-read.table("data/old_gradient/fake_rel_abun_long_n100.txt", row=1, header=T, sep="\t")
mini_sim <- mini_sim[c(1,10,20,30,40,50,60,70,80,90,100),]

# PCoA plot (horseshoe shape) - exported as 450x450
mini_d <- vegdist(mini_sim, method="euclidean")
mini_pc <- cmdscale(mini_d, k=2, eig=T)
mini_cols <- viridis::viridis(11, alpha=0.8)
mini_pc$points[,2] <- mini_pc$points[,2] * -1
## SAVE PNG
### note: width of page w/ 1in margins is 6.5in*96dpi = 624px
###       height of page w/ 1in margins is 9in*96dpi = 864px
# png("wpe_figures/figure4_algorithm/figure4_alg1_start.png",width=300, height=300, units="px", res=96)
# par(mar=c(2,2,1,1), mgp=c(0.5, 0.5, 0))
plot(mini_pc$points, xaxt='n', yaxt='n',
     xlim=range(mini_pc$points)+c(-0.02,0.02), ylim=range(mini_pc$points)+c(-0.01,0.02),
     col=mini_cols, pch=16, cex=6,
     xlab="PC 1", ylab="PC 2", cex.lab=1.5)
# dev.off()

# Radius circle - exported as 450x450
library(plotrix)
## SAVE PNG
### note: width of page w/ 1in margins is 6.5in*96dpi = 624px
###       height of page w/ 1in margins is 9in*96dpi = 864px
# png("wpe_figures/figure4_algorithm/figure4_alg2_radiuscircle.png",width=300, height=300, units="px", res=96)
# par(mar=c(2,2,1,1), mgp=c(0.5, 0.5, 0))
plot("n", xaxt='n', yaxt='n',
     xlim=range(mini_pc$points)+c(-0.02,0.02), ylim=range(mini_pc$points)+c(-0.01,0.02),
     xlab="PC 1", ylab="PC 2", cex.lab=1.5)
points(mini_pc$points[4,1], mini_pc$points[4,2], pch=16, cex=2)
points(mini_pc$points, col=mini_cols, pch=16, cex=6)
draw.circle(mini_pc$points[4,1], mini_pc$points[4,2], 0.09, lwd=3, lty="dashed")
# dev.off()

# Radius circle + connections
## SAVE PNG
### note: width of page w/ 1in margins is 6.5in*96dpi = 624px
###       height of page w/ 1in margins is 9in*96dpi = 864px
# png("wpe_figures/figure4_algorithm/figure4_alg3_radiuscirclelines.png",width=300, height=300, units="px", res=96)
# par(mar=c(2,2,1,1), mgp=c(0.5, 0.5, 0))
plot("n", xaxt='n', yaxt='n',
     xlim=range(mini_pc$points)+c(-0.02,0.02), ylim=range(mini_pc$points)+c(-0.01,0.02),
     xlab="PC 1", ylab="PC 2", cex.lab=1.5)
points(mini_pc$points[4,1], mini_pc$points[4,2], pch=16, cex=2)
draw.circle(mini_pc$points[4,1], mini_pc$points[4,2], 0.09, lwd=3, lty="dashed")
segments(mini_pc$points[4,1],mini_pc$points[4,2],mini_pc$points[5,1],mini_pc$points[5,2], lwd=7, col="black")
segments(mini_pc$points[4,1],mini_pc$points[4,2],mini_pc$points[3,1],mini_pc$points[3,2], lwd=7, col="black")
points(mini_pc$points, col=mini_cols, pch=16, cex=6)
# dev.off()

# Neighborhood connections - exported as 450x450
## get connections for same distance as radius displays
mini_seg <- which(as.matrix(mini_d) <= 0.19, arr.ind=T)
mini_seg <- mini_seg[-which(mini_seg[,1] == mini_seg[,2]),]
rownames(mini_seg) <- 1:nrow(mini_seg)
mini_seg <- as.data.frame(mini_seg)
mini_seg$startx <- mini_pc$points[mini_seg$row,1]; mini_seg$starty <- mini_pc$points[mini_seg$row,2]
mini_seg$endx <- mini_pc$points[mini_seg$col,1]; mini_seg$endy <- mini_pc$points[mini_seg$col,2]
## SAVE PNG
### note: width of page w/ 1in margins is 6.5in*96dpi = 624px
###       height of page w/ 1in margins is 9in*96dpi = 864px
# png("wpe_figures/figure4_algorithm/figure4_alg4_neighbors.png",width=300, height=300, units="px", res=96)
# par(mar=c(2,2,1,1), mgp=c(0.5, 0.5, 0))
plot("n", xaxt='n', yaxt='n',
     xlim=range(mini_pc$points)+c(-0.02,0.02), ylim=range(mini_pc$points)+c(-0.01,0.02),
     xlab="PC 1", ylab="PC 2", cex.lab=1.5)
segments(mini_seg$startx, mini_seg$starty, mini_seg$endx, mini_seg$endy, lwd=7, col="black")
segments(mini_pc$points[9,1],mini_pc$points[9,2],mini_pc$points[11,1],mini_pc$points[11,2], lwd=7, col="black")
points(mini_pc$points, col=mini_cols, pch=16, cex=6)
# dev.off()

# Corrected with LGD - exported as 450x450
mini_lgd <- lg.dist(mini_d, neighborhood.radius=0.15)
mini_pc_lgd <- cmdscale(mini_lgd, k=2, eig=T)
## SAVE PNG
### note: width of page w/ 1in margins is 6.5in*96dpi = 624px
###       height of page w/ 1in margins is 9in*96dpi = 864px
# png("wpe_figures/figure4_algorithm/figure4_alg5_lgdfixed.png", width=300, height=300, units="px", res=96)
# par(mar=c(2,2,1,1), mgp=c(0.5, 0.5, 0))
plot(mini_pc_lgd$points, xaxt='n', yaxt='n',
     xlim=range(mini_pc_lgd$points)+c(-0.02,0.02), ylim=range(mini_pc_lgd$points)+c(-0.02,0.02),
     col=mini_cols, pch=16, cex=6,
     xlab="PC 1", ylab="PC 2", cex.lab=1.5)
# dev.off()



##### FIGURE 5: Simulated dataset #####
dev.off() # reset par values for plot
par(mgp=c(2, 0.5, 0))

# Creation of the simulated dataset (coenoclines)
## creating and plotting simulation dataset:
## create_simulation_data.r -g 1000 -s 50 -n 50 -o data/sim_gradient
## - makes coenoclines plot (in sim_gradient)
## - pie charts also generated along gradient (in sim_gradient)
## custom color bar to go under the x-axis gradient - exported as 750x200
## SAVE PNG
### note: width of page w/ 1in margins is 6.5in*96dpi = 624px
###       height of page w/ 1in margins is 9in*96dpi = 864px
# png("wpe_figures/figure5_simulated_gradient/figure5_sim_gradient_image.png", width=624, height=100, units="px", res=96)
# par(mgp=c(1, 0.5, 0), mar=c(2,1,1,1))
image(seq(1,1000,length.out=50),1,matrix(1:50,ncol=1),col=viridis::viridis(50),axes=FALSE,xlab="",ylab="")
axis(1)
# dev.off()

# Density of distances - exported as 450x450
sim_d <- vegdist(sim_n, method="euclidean")
## SAVE PNG
### note: width of page w/ 1in margins is 6.5in*96dpi = 624px
###       height of page w/ 1in margins is 9in*96dpi = 864px
# png("wpe_figures/figure5_simulated_gradient/figure5_sim2_ogDistDensiity.png", width=310, height=310, units="px", res=96)
# par(mgp=c(2, 0.5, 0), mar=c(3,3,3,1))
plot(density(sim_d), main="Original Dist. Density",
     xlab="Euclidean distances", ylab="Frequency", lwd=4, col="blue")
# dev.off()

# PCoA : Large dataset w/noise - exported as 600x600
sim_pc <- cmdscale(sim_d, k=2, eig=T)
sim_cols <- rep(viridis::viridis(50, alpha=0.8),2)
sim_pc$points[,1] <- sim_pc$points[,1]*-1
## SAVE PNG
### note: width of page w/ 1in margins is 6.5in*96dpi = 624px
###       height of page w/ 1in margins is 9in*96dpi = 864px
#png("wpe_figures/figure5_simulated_gradient/figure5_sim2_ogDistPCoA.png", width=310, height=310, units="px", res=96)
#par(mgp=c(2, 0.5, 0), mar=c(3,3,3,1))
plot(sim_pc$points,
     xlim=range(sim_pc$points)+c(-0.05,0.05), ylim=range(sim_pc$points)+c(-0.05,0.05),
     col=sim_cols, pch=16, cex=2,
     xlab=paste0("PC1 [",calc.perc.var(sim_pc$eig,1),"%]"),
     ylab=paste0("PC2 [",calc.perc.var(sim_pc$eig,2),"%]"), cex.lab=1)
title("Original Distances", adj=0, cex.main=1.5)
legend(x=0.15, y=0.22, legend=c("1",rep(NA,8),"1000"), fill=viridis::viridis(10),
       border=NA, bty="n", y.intersp=0.5, x.intersp=0.2, pt.cex=2, cex=0.7, ncol=1)
#dev.off()
## CHECK CORR (PEARSON)
### PC 1 and gradient
gradient_nums <- rep(round(seq(0, 1000, length.out=(50+2)))[-c(1,50+2)],2) # took this from create_simulation_data.r
cor(sim_pc$points[,1], gradient_nums, method="pearson") # correlation 0.855 for OG PC1 & gradient

# Dist Density : LGD applied - exported as 450x450
sim_lgd <- lg.dist(sim_d, neighborhood.radius=0.15)
## SAVE PNG
### note: width of page w/ 1in margins is 6.5in*96dpi = 624px
###       height of page w/ 1in margins is 9in*96dpi = 864px
# png("wpe_figures/figure5_simulated_gradient/figure5_sim3_lgdDistDensiity.png", width=310, height=310, units="px", res=96)
# par(mgp=c(2, 0.5, 0), mar=c(3,3,3,1))
plot(density(sim_lgd), main="Adjusted Dist. Density",
     xlab="Euclidean distances", ylab="Frequency", lwd=4, col="purple")
# dev.off()

# PCoA : LGD applied - exported as 600x600
sim_pc_lgd <- cmdscale(sim_lgd, k=2, eig=T)
## SAVE PNG
### note: width of page w/ 1in margins is 6.5in*96dpi = 624px
###       height of page w/ 1in margins is 9in*96dpi = 864px
# png("wpe_figures/figure5_simulated_gradient/figure5_sim3_lgdDistPCoA.png", width=310, height=310, units="px", res=96)
# par(mgp=c(2, 0.5, 0), mar=c(3,3,3,1))
plot(sim_pc_lgd$points,
     xlim=range(sim_pc_lgd$points)+c(-0.05,0.05), ylim=range(sim_pc_lgd$points)+c(-0.05,0.05),
     col=sim_cols, pch=16, cex=2,
     xlab=paste0("PC1 [",calc.perc.var(sim_pc_lgd$eig,1),"%]"),
     ylab=paste0("PC2 [",calc.perc.var(sim_pc_lgd$eig,2),"%]"), cex.lab=1)
title("Adjusted Distances", adj=0, cex.main=1.5)
legend(x=1.1, y=2, legend=c("1",rep(NA,8),"1000"), fill=viridis::viridis(10),
       border=NA, bty="n", y.intersp=0.5, x.intersp=0.2, pt.cex=2, cex=0.7, ncol=1)
# dev.off()
## CHECK CORR
### PC 1 and gradient
cor(sim_pc_lgd$points[,1], gradient_nums, method="pearson") # correlation 0.999 LGD PC1 w/ gradient
# PC 2 and noise (not equal variance between groups...)
tmp_lgd_dat <- data.frame(pc2=sim_pc_lgd$points[,2], noise=as.factor(c(rep(0,50),rep(1,50))))
var(tmp_lgd_dat[1:50,"pc2"])
var(tmp_lgd_dat[51:100,"pc2"])
leveneTest(pc2 ~ noise, tmp_lgd_dat) # Levene's test of variance, null hyp is homogenous variance between groups



##### FIGURE 6: LGD applied to Case Studies #####
# Note: colors and original distances created as part of Figure 3!
dev.off()  # reset the par settings
par(mgp=c(2, 0.5, 0))

# 88 soils
## recalculate distances
soil_lgd <- lg.dist(soil_dat, neighborhood.radius = 0.9)  # LGD distances
soil_pc_lgd <- cmdscale(soil_lgd, k=2, eig=T)             # new PCoA
## plot PCoA (exported as 600x600)
## SAVE PNG
### note: width of page w/ 1in margins is 6.5in*96dpi = 624px
###       height of page w/ 1in margins is 9in*96dpi = 864px
# png("wpe_figures/figure6_resolved_cases/figure6_soil_pcoalgd.png", width=500, height=500, units="px", res=96)
# par(mgp=c(2, 0.5, 0))
plot(soil_pc_lgd$points, #xaxt='n', yaxt='n',
     xlim=range(soil_pc_lgd$points)+c(-0.1, 0.1), ylim=range(soil_pc_lgd$points)+c(-0.1, 0.1),
     col=soil_col, bg=soil_bg, pch=soil_shp, cex=2, lwd=3,
     xlab=paste0("PC1 [", calc.perc.var(soil_pc_lgd$eig,1), "%]"),
     ylab=paste0("PC2 [", calc.perc.var(soil_pc_lgd$eig,2), "%]"), cex.lab=2)
title("88 Soils", adj=0, cex.main=2)
legend("topright", inset=0.02, legend=c("ph", "<4","4-5","5-6","6-7","7-8",">8"),
       col=c("white",soil_col1), pt.bg=c("white",soil_bg1), pch=c(16,soil_shp1),
       cex=1, y.intersp=0.8, bty="o")
# dev.off()
cor(soil_pc_lgd$points[,1], meta_soil$ph, method="spearman") # 0.929 PC1 with gradient
cor.test(soil_pc_lgd$points[,1], meta_soil$ph, method="spearman") # p < 0.001
## statistics done in publication
anosim(soil_lgd, soil_col, permutations=999) # ANOSIM with pH groups: R=0.411, p=0.001
mantel(soil_lgd, vegdist(meta_soil$ph, method="euclidean"), method="spearman") # Mantel test: R=0.704, p=0.001
## PC2 variation linked to something?
corr_meta_soil <- apply(meta_soil[,unlist(lapply(meta_soil, is.numeric))], 2, FUN=function (x) {cor(x, soil_pc_lgd$points[,2], method="pearson")})
head(corr_meta_soil[order(abs(corr_meta_soil), decreasing=T)])
cor.test(meta_soil$annual_season_temp, soil_pc_lgd$points[,2], method="pearson")
cor.test(meta_soil$latitude, soil_pc_lgd$points[,2], method="pearson")
cor(meta_soil$annual_season_temp, meta_soil$latitude, method="spearman")
cor.test(meta_soil$annual_season_temp, meta_soil$latitude, method="spearman")


# Cecum dataset
## recalculate distances
cecum_lgd <- lg.dist(cecum_d, neighborhood.radius = 0.93)   # LGD distances
cecum_pc_lgd <- cmdscale(cecum_lgd, k=2, eig=T)            # new PCoA
## plot PCoA (exported as 700x700)
## SAVE PNG
### note: width of page w/ 1in margins is 6.5in*96dpi = 624px
###       height of page w/ 1in margins is 9in*96dpi = 864px
# png("wpe_figures/figure6_resolved_cases/figure6_cecum_pcoalgd.png", width=500, height=500, units="px", res=96)
# par(mgp=c(2, 0.5, 0))
plot(cecum_pc_lgd$points,
     xlim=range(cecum_pc_lgd$points)+c(-0.1,0.1), ylim=range(cecum_pc_lgd$points)+c(-0.05,0.25),
     col=alpha(cecum_col,0.8), pch=cecum_shp, cex=2,
     xlab=paste0("PC1 [", calc.perc.var(cecum_pc_lgd$eig,1), "%]"),
     ylab=paste0("PC2 [", calc.perc.var(cecum_pc_lgd$eig,2), "%]"), cex.lab=2)
title("Turkey Cecum", adj=0, cex.main=2)
legend("bottomright", inset=0.04, legend=c("Day","1","4","8","10","15","18","22","29","","System","Pen","Hatch"),
       col=c("white",cecum_col1,"white","white","black","black"), pch=c(rep(16,9),16,16,16,17), 
       cex=1, y.intersp=0.8, x.intersp=0.5, text.width=0.26, bty="o", ncol=3)
# dev.off()
cor(cecum_pc_lgd$points[,1], as.numeric(gsub("D","",meta_cecum$Collection)), method="pearson") # 0.71 PC1 with gradient
cor.test(cecum_pc_lgd$points[,1], as.numeric(gsub("D","",meta_cecum$Collection)), method="spearman") # p < 0.001
## statistics from the publication
adonis2(cecum_lgd ~ Collection, meta_cecum, permutations=999)
mantel(as.dist(as.matrix(cecum_lgd)[meta_cecum$System=="Hatch_Brood",meta_cecum$System=="Hatch_Brood"]),
       as.dist(as.matrix(cecum_lgd)[meta_cecum$System=="Pen",meta_cecum$System=="Pen"]), method="spearman") # r = 0.473, p = 0.001



# Guerrero Negro
## recalculate distances
gn_lgd <- lg.dist(gn_d, neighborhood.radius = 0.74)
gn_pc_lgd <- cmdscale(gn_lgd, k=2, eig=T)
## plot PCoA (exported as 700x700)
## SAVE PNG
### note: width of page w/ 1in margins is 6.5in*96dpi = 624px
###       height of page w/ 1in margins is 9in*96dpi = 864px
# png("wpe_figures/figure6_resolved_cases/figure6_microbialmat_pcoalgd.png", width=500, height=500, units="px", res=96)
# par(mgp=c(2, 0.5, 0))
plot(gn_pc_lgd$points,
     xlim=range(gn_pc_lgd$points)+c(0,0.2), ylim=range(gn_pc_lgd$points),
     col=alpha(depth_cols,0.8), pch=16, cex=3,
     xlab=paste0("PC1 [", calc.perc.var(gn_pc_lgd$eig,1), "%]"),
     ylab=paste0("PC2 [", calc.perc.var(gn_pc_lgd$eig,2), "%]"), cex.lab=2)
title("Microbial Mat", adj=0, cex.main=2)
legend("bottomleft", inset=0.04, legend=c("Depth",unique(meta_gn$depth_mm)),
       col=c("white",depth_cols1), pch=16, cex=1, y.intersp=0.8, bty="o", ncol=2)
# dev.off()
cor(gn_pc_lgd$points[,1], meta_gn$end_depth, method="spearman")  # 0.979 PC1 with gradient
cor.test(gn_pc_lgd$points[,1], meta_gn$end_depth, method="spearman") # p < 0.001


# MAGIC
## recalculate distances
magic_lgd <- lg.dist(magic_d, neighborhood.radius = 0.93)
magic_pc_lgd <- cmdscale(magic_lgd, k=2, eig=T)
## plot PCoA (exported as 700x700)
## SAVE PNG
### note: width of page w/ 1in margins is 6.5in*96dpi = 624px
###       height of page w/ 1in margins is 9in*96dpi = 864px
# png("wpe_figures/figure6_resolved_cases/figure6_infant_pcoalgd.png", width=500, height=500, units="px", res=96)
# par(mgp=c(2, 0.5, 0))
plot(magic_pc_lgd$points,
     xlim=range(magic_pc_lgd$points)+c(0,0.1), ylim=range(magic_pc_lgd$points)+c(0,0.1),
     col=alpha(magic_cols,0.7), pch=magic_shp, cex=1,
     xlab=paste0("PC1 [", calc.perc.var(magic_pc_lgd$eig,1), "%]"),
     ylab=paste0("PC2 [", calc.perc.var(magic_pc_lgd$eig,2), "%]"), cex.lab=2)
title("Infant Gut", adj=0, cex.main=2)
legend("topright", inset=0.02,
       legend=c("Age","1mo","3mo","6mo","9mo","1yr","15mo","18mo","2yr","","BF","Yes","No"),
       col=c("white",magic_cols1,"white","white","black","black"), pch=c(rep(16,12),17),
       cex=0.85, y.intersp=0.8, x.intersp=0.5, bty="o", ncol=2)
# dev.off()
cor(-1*magic_pc_lgd$points[,1], meta_magic$age, method="spearman") # 0.638 PC1 with gradient
cor.test(magic_pc_lgd$points[,1], meta_magic$age, method="spearman") # p < 0.001
# ## making 20 plots to judge best radius value
# for (r in c(1.00,0.97,0.95,0.92,0.90,0.87,0.84,0.80,0.77,0.74,0.71,0.68,0.64,0.60,0.55,0.50)){
#   ## recalculate distances
#   magic_lgd <- lg.dist(magic_d, neighborhood.radius = r)
#   magic_pc_lgd <- cmdscale(magic_lgd, k=2, eig=T)
#   ## plot PCoA (exported as 700x700)
#   plot(magic_pc_lgd$points,
#        xlim=range(magic_pc_lgd$points), ylim=range(magic_pc_lgd$points),
#        col=alpha(magic_cols,0.7), pch=magic_shp, cex=1,
#        xlab=paste0("PC1 [", calc.perc.var(magic_pc_lgd$eig,1), "%]"),
#        ylab=paste0("PC2 [", calc.perc.var(magic_pc_lgd$eig,2), "%]"), cex.lab=2)
#   title(paste0("r = ", r))
#   print(paste0("Done with r = ", r, "!"))
# }



##### Supplemental: RDA applied to case studies #####
# 88 soils
## compute rda (traditional)
rda_soil <- dbrda(t(soil_n) ~ ph, data=meta_soil)
## my plot of rda output (export as 700x700)
plot(rda_soil, type="n")
mtext(side=3, line=0.5, adj=0, cex=1.5, "88 Soils RDA")
points(rda_soil, pch=soil_shp, col=soil_col, bg=soil_bg, cex=2)
text(rda_soil, dis="cn", cex=1)
legend("topright", inset=0.04, legend=c("ph", "<4","4-5","5-6","6-7","7-8",">8"),
       col=c("white",soil_col1), pt.bg=c("white",soil_bg1), pch=c(16,soil_shp1),
       cex=1, y.intersp=0.5, bty="o")
## problems: pulls out one sample for some reason??
# ## compute rda (distance-based)
# rda_soil <- dbrda(soil_dat ~ ph, data=meta_soil)
# ## my plot of rda output (export as 700x700)
# plot(rda_soil, type="n")
# mtext(side=3, line=0.5, adj=0, cex=1.5, "88 Soils db-RDA")
# points(rda_soil, pch=soil_shp, col=soil_col, bg=soil_bg, cex=2)
# text(rda_soil, dis="cn", cex=1)
# legend("topright", inset=0.04, legend=c("ph", "<4","4-5","5-6","6-7","7-8",">8"),
#        col=c("white",soil_col1), pt.bg=c("white",soil_bg1), pch=c(16,soil_shp1),
#        cex=1, y.intersp=0.5, bty="o")
## problems: same arch with db approach


# Cecum
## compute rda (distance-based)
rda_cecum <- dbrda(cecum_d ~ Collection + System, data=meta_cecum)
## my plot of rda output (export as 700x700)
plot(rda_cecum, type="n")
mtext(side=3, line=0.5, adj=0, cex=1.5, "Turkey Cecum db-RDA")
points(rda_cecum, pch=cecum_shp, col=cecum_col, cex=2)
text(rda_cecum, dis="cn", cex=1)
legend("bottomleft", inset=0.04, legend=c("Day","1","4","8","10","15","18","22","29","","System","Pen","Hatch"),
       col=c("white",cecum_col1,"white","white","black","black"), pch=c(rep(16,9),16,16,16,17), 
       cex=1, y.intersp=0.5, bty="o", ncol=2)
## problems: same arch with db approach


# Guerrero Negro
## compute rda (traditional)
rda_gn <- rda(t(gn_n) ~ start_depth, data=meta_gn)
## my plot of rda output (export as 700x700)
plot(rda_gn, type="n")
mtext(side=3, line=0.5, adj=0, cex=1.5, "Microbial Mat RDA")
points(rda_gn, pch=16, col=depth_cols, cex=2)
text(rda_gn, dis="cn", cex=1)
legend("bottomright", inset=0.04, legend=c("Depth",unique(meta_gn$depth_mm)),
       col=c("white",depth_cols1), pch=16, cex=1, y.intersp=0.5, bty="o", ncol=2)
## problems: same arch even with traditional approach
# ## compute rda (distance-based)
# rda_gn <- dbrda(gn_d ~ start_depth, data=meta_gn)
# ## my plot of rda output (export as 700x700)
# plot(rda_gn, type="n")
# mtext(side=3, line=0.5, adj=0, cex=1.5, "Microbial Mat db-RDA")
# points(rda_gn, pch=16, col=depth_cols, cex=2)
# text(rda_gn, dis="cn", cex=1)
# legend("bottom", inset=0.04, legend=c("Depth",unique(meta_gn$depth_mm)),
#        col=c("white",depth_cols1), pch=16, cex=1, y.intersp=0.5, bty="o", ncol=2)
## problems: same arch with db approach

# MAGIC
## compute rda (traditional)
rda_magic <- rda(t(magic_n) ~ age + currentfeed_bf, data=meta_magic)
## my plot of rda output (export as 700x700)
plot(rda_magic, type="n")
mtext(side=3, line=0.5, adj=0, cex=1.5, "Infant gut RDA")
points(rda_magic, pch=magic_shp, col=alpha(magic_cols,0.8), cex=1)
text(rda_magic, dis="cn", cex=1)
legend("topright", inset=0.02,
       legend=c("Age","1mo","3mo","6mo","9mo","1yr","15mo","18mo","2yr","","BF","Yes","No"),
       col=c("white",magic_cols1,"white","white","black","black"), pch=c(rep(16,12),17),
       cex=1, y.intersp=0.5, x.intersp=0.2, bty="o", ncol=2)
## problems: pretty much useless...
# ## compute rda (distance-based)
# rda_magic <- dbrda(magic_d ~ age + currentfeed_bf, data=meta_magic)
# ## my plot of rda output (export as 700x700)
# plot(rda_magic, type="n")
# mtext(side=3, line=0.5, adj=0, cex=1.5, "Infant gut db-RDA")
# points(rda_magic, pch=magic_shp, col=magic_cols, cex=2)
# text(rda_magic, dis="cn", cex=1)
# legend("topright", inset=0.02,
#        legend=c("Age","1mo","3mo","6mo","9mo","1yr","15mo","18mo","2yr","","BF","Yes","No"),
#        col=c("white",magic_cols1,"white","white","black","black"), pch=c(rep(16,12),17),
#        cex=1, y.intersp=0.5, x.intersp=0.2, bty="o", ncol=2)
## problems: also pretty much useless...




