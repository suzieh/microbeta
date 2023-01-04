# create_simulation_data.r
# Create simulation data for testing Nonlinear Dimensionality Reduction Techniques
# Knights Lab - University of Minnesota
# June 2019
# usage : create_simulation_data.r -o /project/flatiron2/suzie/detrending/poster

##### Set Up #####
suppressMessages(require(optparse))
suppressMessages(require(gtools))
suppressMessages(require(reshape2))
suppressMessages(require(ggplot2))


##### Parse Command Line #####
option_list <- list(make_option(c("-g", "--gradient"), type="integer", default=1000,
                                help="Length of gradient, starting at 0. [default: 1000]"),
                    make_option(c("-s", "--standard_dev"), type="integer", default=75,
                                help="Standard deviation of each normal OTU distribution. [default: 75]"),
                    make_option(c("-n", "--num_samples"), type="integer", default=50,
                                help="Number of samples (will be doubled with noise). [default: 50]"),
                    make_option(c("-d", "--dirichlet_val"), type="double", default=0.2,
                                help="Dirichlet confidence factor, higher is closer to real data . [default: 0.2]"),
                    make_option(c("-o", "--out_dir"), type="character",
                                help="File path for output directory. [required]"),
                    make_option(c("-v", "--visuals"), type="logical", default=TRUE,
                                help="Logical request to produce visuals in output directory. [default: TRUE]"),
                    make_option(c("-c", "--color_viridis"), type="logical", default=TRUE,
                                help="Logical request to use viridis colors. FALSE uses distinct colors. [default: TRUE]"))
opts <- parse_args(OptionParser(option_list=option_list), args=commandArgs(trailing=T))
g <- opts$gradient
sd <- opts$standard_dev
n <- opts$num_samples
conf <- opts$dirichlet_val
v <- opts$visuals
out_dir <- opts$out_dir
if (!dir.exists(out_dir)) {
  stop("Output directory does not exist. Could not run script.\nUse create_simulation_data.r -- help to see help menu.")
} else if (substr(out_dir, nchar(out_dir), nchar(out_dir)) != "/") {
  out_dir <- paste0(out_dir, "/")
}


##### Coenoclines #####
# Create otu set
m <- 1000 # max gradient length - hard coded for simplicity
xvals <- seq(-(4*sd), m+(4*sd), by=sd/4)[-1]  # ensures 0 to m in figure is flat
x <- sapply(xvals, function (x) {dnorm(seq(0,m,1), x, sd)*100})  # create null distributions
x[abs(x) < 0.001] <- 0  # remove otu entirely at tails
x <- sweep(x, 1, rowSums(x), "/")  # normalize rows (i.e. possible sample combinations)
x <- x[,-which(colSums(x) == 0)] # remove buffer otus


##### Choose Samples #####
inds <- round(seq(0, g, length.out=(n+2)))[-c(1,n+2)] # evenly selects samples along distribution
dat <- as.data.frame(x[inds,])
dat <- dat[,-which(colSums(dat) == 0)]
rownames(dat) <- paste0("sample.", rownames(dat))
colnames(dat) <- paste0("otu.", 1:ncol(dat))


##### Add Noise #####
# Dirichlet per "actual" sample
# note: does not randomly add other otus not in actual sample,
#       may want to add this feature to increase noise later
# [conf] confidence factor is in options - higher means noise is closer to actual sample
noise <- as.data.frame(t(apply(dat, 1, function (x) colMeans(rdirichlet(100, x*conf))))) # use each sample as its own prior
colnames(noise) <- colnames(dat)
rownames(noise) <- paste0("sample.n",seq(1,nrow(noise),1))
# par(mar=c(5.1, 4.1, 4.1, 4.1))
# plot(as.matrix(dat)*100, col = viridis)
# plot(as.matrix(noise)*100, col = viridis)
# par(mar=c(5.1, 4.1, 4.1, 2.1))
dat <- rbind(dat, noise)


##### Visualize Dataset #####
# Set number of coenoclines visualized
nv <- ncol(dat)/15   # every nv-th otu coenocline

# Adjust colors
if(opts$color_viridis) {
  colors <- viridis::viridis((ncol(x) - 4) / nv)
} else {
  colors <- c("#000000", "#800000", "#3cb44b", "#911eb4", "#e6194B",
              "#469990", "#f032e6", "#f58231", "#808000", "#000075") }
# Determine "len" to append to file names
len <- ifelse(g > m/2, "long", "medium")
len <- ifelse(g < m/4, "short", len)
# Coenoclines plot function
plot_coenoclines <- function (df, len) {
  m_df <- df
  m_df[m_df == 0] <- NA
  png(paste0(out_dir, "coenoclines_", len, ".png"), width=624, height=416, res=96) # bg = "transparent")
  par(mgp=c(0.8, 0, 0))
  matplot(m_df, type="l", col=colors, lty=1, lwd=15, xlim=c(0,m),
          ylab="relative abundances", xlab="gradient",
          cex.lab=1.8, xaxt="n", yaxt="n")
  dev.off()
  return(TRUE)
}
# Sample pie chart plot function
plot_sample_pie <- function (df, samp) {
  df <- suppressMessages(melt(df))
  p <- ggplot(df, aes(x="", y=value, fill=variable)) +
    geom_bar(stat="identity", width=1) +
    coord_polar("y", start=0) +
    scale_fill_viridis_d(nrow(df)) +
    theme_void() +
    theme(legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "transparent",colour = NA),
          plot.background = element_rect(fill = "transparent",colour = NA)) + NULL
  p
  ggsave(paste0(out_dir, "pie_", samp, ".png"), width=48, height=48, units="px", device="png", bg="transparent", dpi=96)
  return(ifelse(is.null(p), FALSE, TRUE))
}
# Print out visuals
if (v) {
  # set up transparent background
  par(bg="transparent")
  # coenocline plot (size: 600x400)
  rx <- x[,colSums(dat) > 0] # restrict to otus in present sampling
  rx[g:nrow(rx),] <- NA # cut off otus after current gradient end
  plot_coenoclines(rx[,seq(2, ncol(rx)-2, by=nv)], len)   # by: print every nth coenocline
  # sample pie charts (size: 400x400)
  plot_sample_pie(dat[1,], "start") # first sample pie
  plot_sample_pie(dat[nrow(dat)/4,], "mid") # midpoint sample pie
  plot_sample_pie(dat[nrow(dat)/2,], "end") # end sample pie (recall second half is noise)
}


##### Write to file #####
# Refine data to only present otus
dat <- dat[,colSums(dat) > 0] # restrict to otus in present sampling
# Write to provided output directory
write.table(dat, paste0(out_dir, "fake_rel_abun_", len, ".txt"), sep="\t", row.names=T, col.names=T)
# Write out log of command used to specified directory
fileConn<-file(paste0(out_dir, "commands_", len, ".log"))
writeLines(c(paste0("-g --gradient\t\t\t", g),
             paste0("-s --standard-dev\t\t", sd),
             paste0("-n --num_samples\t\t", n),
             paste0("-o --out_dir\t\t\t", out_dir),
             paste0("-v --visuals\t\t\t", v),
             paste0("-c --color_viridis\t\t", opts$color_viridis)), fileConn)
close(fileConn)







##### Additional Stuff #####
# # Example of a sample and corresponding noise - taxa summary plot
# sample_example <- dat[c(1, n+1),] # created for stacked bar plotting purposes
# sample_example <- sample_example[,colSums(sample_example) > 0] # remove otus not present
# sample_example$name <- c("sample", "noise")
# sample_example <- melt(sample_example, id.vars = "name")
# sample_example$name <- factor(sample_example$name, levels = c("sample","noise"))
# ggplot(sample_example, aes(x = name, y = value, fill = variable)) +
#   geom_bar(stat = "identity") +
#   labs(title = "OTU Distribution in Sample", x = "", y = "Relative Abundance") + 
#   scale_fill_manual(values = colors) +
#   theme_classic() +
#   theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1),
#         legend.position = 'none')

