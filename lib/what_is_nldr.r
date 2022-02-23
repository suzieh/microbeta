# what_is_nldr.r
# Scratch Area for Creating Visuals to Explain NLDR
# Knights Lab - University of Minnesota
# June 2019
# usage : source('what_is_nldr.r')

##### Set Up #####
library(plotly)
library(viridis)
library(lle)
library(ggplot2)


##### Create Multi-Dimensional Data #####
radius <- 80
center_x <- 100
center_y <- 100
get_2d_point <- function (radius, center_x, center_y, r_angle, r_radius) {
  alpha <- 2 * pi * r_angle # random angle
  r <- radius * sqrt(r_radius) # random radius
  x_coord <- r * cos(alpha) + center_x
  y_coord <- r * sin(alpha) + center_y
  return(c(x_coord, y_coord))
}
mdd <- t(sapply(seq(0.01,0.99,by=0.002), function (x) {get_2d_point(radius, center_x, center_y, x, x)}))
mdd <- as.data.frame(mdd)
colnames(mdd) <- c("x", "z")
mdd$y <- sample(-50:50, nrow(mdd),replace=TRUE)
mdd$cols <- seq(1, nrow(mdd), 1)


##### Multidimensional Plot #####
p3 <- plot_ly(mdd, x=~x, y=~y, z=~z,
             marker = list(color = ~cols, colorscale = 'Viridis', showscale = FALSE, size = 6),
             hoverinfo = 'text',
             text = ~paste('Sample ', cols)) %>%
      add_markers()
p3


##### NLDR (principal component analysis) of MDD #####
pca_out <- prcomp(mdd, center=T, scale.=F)


##### NLDR Plot (2D) #####
pca_out <- as.data.frame(pca_out$x)
pca_out$cols <- mdd$cols
p2 <- ggplot(pca_out, aes(x=PC1, y=PC2, color=cols)) +
  geom_point(size = 5, pch=20, alpha = 0.3) +
  scale_color_viridis() +
  theme_classic() +
  theme(legend.position = 'none')
p2



