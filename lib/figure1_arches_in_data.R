# Figure 1 & Global Gut stuff
library(ggplot2)
library(ggpubr)
library(tools)
source("lib/lgd_source.r")
set.seed(125)

# Load the global gut dataset
dat_hg <- read.delim("421_clean_otus.txt", header=T, sep="\t", quote="\"", row.names=1, comment.char="", as.is=T)
meta_hg <- read.delim("data/global_gut/clean_map.txt",header=T, sep="\t", quote="\"", comment.char="", as.is=T)
rownames(meta_hg) <- meta_hg$Sample_ID
pop_cols <- c("red", "blue", "green"); names(pop_cols) <- c("Malawi", "USA", "Venezuela");

# Load the unifrac distances
d_hg <- read.delim("clean_uw_unifrac_dist.txt", header=T, row.names=1, sep="\t", as.is=T)
d_hg <- as.dist(as.matrix(d_hg))
identical(rownames(as.matrix(d_hg)), meta_hg$Sample_ID)


# Helper functions for plotting & conducting LGD
mypcoaplot <- function (mydf, myvar, colorby, mytitle="") {
  if (colorby == "geo_loc_name") {
    mycolors <- pop_cols[unique(mydf$geo_loc_name)]
    mydf$colorcol <- factor(mydf$geo_loc_name)
  } else {
    mydf$colorcol <- mydf[,colorby]
  }
  out <- ggplot(mydf, aes(x=X, y=Y, color=colorcol)) +
    geom_point(alpha=0.8, size=3) +
    labs(title=mytitle,
         x=paste0("PC 1 [",myvar[1],"%]"),
         y=paste0("PC 2 [",myvar[2],"%]")) +
    coord_fixed(ratio = 1) +
    theme_classic()
  if (colorby == "geo_loc_name") {
    out <- out + 
      scale_color_manual(values=mycolors) +
      guides(color = guide_legend("Population")) + NULL
  } else if (colorby == "Sample.Group") {
    out <- ggplot(mydf, aes(x=X, y=Y, color=colorcol, shape=colorcol)) +
      geom_point(alpha=0.8, size=3) +
      labs(title=mytitle,
           x=paste0("PC 1 [",myvar[1],"%]"),
           y=paste0("PC 2 [",myvar[2],"%]")) +
      coord_fixed(ratio = 1) +
      theme_classic() +
      scale_color_manual(name="Group",values=g_col[unique(mydf$colorcol)]) +
      scale_shape_manual(name="Group",values=g_shp[unique(mydf$colorcol)]) +
      NULL
  } else {
    out <- out +
      scale_color_gradient(low = "#ffc6c4", high = "#672044") +
      guides(color = guide_legend(toTitleCase(gsub("_", " ", colorby)))) + NULL
  }
  return(out)
}

mydensityplot <- function (myd, mytitle="") {
  out <- ggplot(data.frame(weight=c(myd)), aes(x=weight)) +
    geom_density(color="black", size=2) +
    labs(title=mytitle) +
    theme_classic() + NULL
  return(out)
}
  
plot4 <- function (d, colorby="geo_loc_name", meta) {
  tmp_meta <- meta[rownames(as.matrix(d)),]
  ## original data and plots
  pc_o <- cmdscale(d, k=3, eig=T)
  var_o <- round((pc_o$eig / sum(pc_o$eig)) * 100, 1)
  df_o <- cbind(data.frame(X=pc_o$points[,1], Y=pc_o$points[,2], Z=pc_o$points[,3]), tmp_meta)
  ogp <- mypcoaplot(df_o, var_o, colorby, "Original PCoA (unweighted UniFrac)")
  dp <- mydensityplot(d, "Density of unw. UniFrac Distances")
  ## lgd data and plots
  lgd <- lg.dist(d)
  pc_lg <- cmdscale(lgd, k=3, eig=T)
  var_lg <- round((pc_lg$eig / sum(pc_lg$eig)) * 100, 1)
  df_lg <- cbind(data.frame(X=pc_lg$points[,1], Y=pc_lg$points[,2], Z=pc_lg$points[,3]), tmp_meta)
  lgp <- mypcoaplot(df_lg, var_lg, colorby, "New PCoA (LGD-adjusted)")
  dlgp <- mydensityplot(lgd, "Density of LGD-adjusted Distances")
  ## arrange all plots together
  ggarrange(ogp, dp, lgp, dlgp, ncol=2, nrow=2, labels=c("A","B","C","D"))
}


# Original : all populations, unw UniFrac
plot4(d_hg, "geo_loc_name", meta_hg) # LGD chose to stick with largest radius, no change!


# PCoA of U.S. Only (not very interesting it would seem)
d.us <- as.dist(as.matrix(d_hg)[meta_hg$geo_loc_name=="USA", meta_hg$geo_loc_name=="USA"])
plot4(d.us, "age", meta_hg) # LGD still chose to stick with largest radius, no change induced!


# PCoA of Venezuelan and Malawian
d.nonus <- as.dist(as.matrix(d_hg)[meta_hg$geo_loc_name!="USA", meta_hg$geo_loc_name!="USA"])
plot4(d.nonus, "geo_loc_name", meta_hg) # LGD still chose to stick with largest radius, no change induced!


# PCoA of Venezuela only
d.vz <- as.dist(as.matrix(d_hg)[meta_hg$geo_loc_name=="Venezuela", meta_hg$geo_loc_name=="Venezuela"])
plot4(d.vz, "age", meta_hg) # LGD still chose to stick with largest radius, no change induced!


# PCoA of Malawi only
d.ma <- as.dist(as.matrix(d_hg)[meta_hg$geo_loc_name=="Malawi", meta_hg$geo_loc_name=="Malawi"])
plot4(d.ma, "age", meta_hg) # LGD did pick a different radius! Seems to exaggerate the difference between some outliers and the rest of the group.


# IMP study data
## loading and prepping data
imp <- read.delim('data/imp/unweighted_unifrac_dm.txt', sep="\t", row=1, header=T)
meta_imp <- read.delim('data/imp/map.txt', sep="\t", header=T)
rownames(meta_imp) <- meta_imp$SampleID
meta_imp <- meta_imp[rownames(imp),]
exclude <- grep("1st", meta_imp$Sample.Group)
exclude <- c(exclude, which(is.na(meta_imp$Sample.Group)))
exclude <- c(exclude, which(meta_imp$SampleID %in% c("T.CS.023","T.CS.054","CS.227")))
imp <- imp[-exclude,-exclude]
meta_imp <- meta_imp[-exclude,]
identical(rownames(imp), meta_imp$SampleID)
meta_imp$Sample.Group <- as.factor(meta_imp$Sample.Group)
## colors and shapes and plotting
g_col <- c("#e38cbd", "#7e7e7e", "#d65ca0", "#4ea99b")
names(g_col) <- c("Hmong2nd", "Control", "HmongThai", "KarenThai")
g_shp <- c(17, 17, 16, 16)
names(g_shp) <- c("Hmong2nd", "Control", "HmongThai", "KarenThai")
# plot4(as.dist(imp), "Sample.Group", meta_imp)
## let's try with a less aggressive radius
imp_lgd <- lg.dist(as.dist(imp), 0.6)
imp_lgd_pc <- cmdscale(imp_lgd, k=3, eig=T)
colnames(imp_lgd_pc$points) <- c("PC1", "PC2", "PC3")
df <- cbind(imp_lgd_pc$points, meta_imp)
imp_var <- (imp_lgd_pc$eig / sum(imp_lgd_pc$eig)) * 100
ggplot(df, aes(x=PC1, y=PC2, color=Sample.Group, shape=Sample.Group)) +
  geom_point(alpha=0.8, size=3) +
  labs(title="IMP LGD-adjusted",
       x=paste0("PC 1 [",imp_var[1],"%]"),
       y=paste0("PC 2 [",imp_var[2],"%]")) +
  coord_fixed(ratio = 1) +
  theme_classic() +
  scale_color_manual(name="Group",values=g_col[unique(df$Sample.Group)]) +
  scale_shape_manual(name="Group",values=g_shp[unique(df$Sample.Group)]) +
  geom_text(stat="identity", aes(label=SampleID))




