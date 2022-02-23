# mca_simulation.r
# Multiple Correspondence Analysis on Simulation Data
# Knights Lab - University of Minnesota
# August 2019
# usage : mca_simulation.r -i input_file.txt

##### Set Up #####
#library(optparse)
library(ggplot2)
library(FactoMineR)
#library(factoextra)


##### Parse Command Line #####
# option_list <- list(make_option(c("-i", "--input_table"), type="character",
#                                 help="Path to input file. Expects otu table in tab-delimited txt file format.") )
# opts <- parse_args(OptionParser(option_list=option_list), args=commandArgs(trailing=T))
dat <- read.delim("/project/flatiron2/suzie/detrending/fake/fake_rel_abun.txt", sep="\t", header=T)
colorfunc <- colorRampPalette(c("white", "navy"))


##### Correspondence Analysis #####
dat <- round(dat*100000) # use "raw counts"
ca_out <- CA(dat, graph=F)
summary(ca_out) # get Chi Square statistic here and put below
pchisq(447147417, df = (nrow(dat) - 1) * (ncol(dat) - 1), lower.tail = F) #no significant relationship found b/w the rows and columns
ca_points <- as.data.frame(ca_out$row$coord); colnames(ca_points) <- gsub(" ", "", colnames(ca_points));
ca_points$color <- rep(1:(nrow(ca_points)/2),2)
ca_points$noise <- c(rep("no",(nrow(ca_points)/2)), rep("yes",(nrow(ca_points)/2)))
ggplot(ca_points, aes(x=Dim1, y=Dim2, color=factor(color), shape=factor(noise), size=factor(noise))) +
  geom_point(alpha = 0.5) + 
  scale_color_viridis_d(nrow(ca_points)/2) +
  scale_shape_manual(values=c(1, 18)) +
  scale_size_manual(values=c(15, 10)) +
  labs(title = "Multiple Correspondence Analysis",
       x = "Dimension 1",
       y = "Dimension 2",
       shape = "noise") + 
  guides(color="none", size="none", shape = guide_legend(override.aes = list(size=6))) +
  theme_classic() + NULL
fviz_screeplot(ca_out, addlabels = T, ylim = c(0,10)) # scree plot has no elbow...


