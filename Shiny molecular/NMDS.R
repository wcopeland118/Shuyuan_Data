#  building an NMDS function for Molecular biology

{ #libraries
  library(tidyverse)
  library(lubridate)
  library(zoo)
  library(RColorBrewer)
  library(hms)
  library(vegan)
  library(dendextend)
} #libraries

Eco_data <- read_csv("NMDS_data1.csv")

x <- Eco_data[,3:11] 

BC_matrix <- as.matrix(vegdist(x, method="bray", binary=FALSE, diag=FALSE, upper=T,
        na.rm = FALSE))

NumCols <- length(colnames(Eco_data))
test1 <- print(BC_matrix, diag = NULL, upper = NULL,
      digits = getOption("digits"), justify = "none",
      right = TRUE)


Habitat <- Eco_data[,2]
NMDS <- metaMDS(x, distance = "bray", k = 2)

plot_df <- tibble(x = NMDS$points[,1], y = NMDS$points[,2], 
                 Habitat) %>% 
  rownames_to_column(var = "n")

p <- plot_df %>% 
  ggplot() +
  geom_point(mapping = aes(x = x, y = y, colour = Habitat, 
                           label = n)) + # add label to aes to show in plotly tooltip
  labs(
    title = "NMDS plot
",
    x = "Axis 1",
    y = "Axis 2"
  ) +
  theme(plot.title = element_text(color = "steelblue4", size = 14, face = "bold"))
ggplotly(p, tooltip = "n")

##  clustering diagram
# Euclidean distance
dist <- dist(data[ , c(4:8)] , diag=TRUE)

# Hierarchical Clustering with hclust
dend <- vegdist(x, method="bray", binary=FALSE, diag=FALSE, upper=T,
                     na.rm = FALSE) %>% 
  hclust %>% 
  as.dendrogram
  set("branches_k_color", k=3) %>% set("branches_lwd", c(1.5,1,1.5)) %>%
  set("branches_lty", c(1,1,3,1,1,2)) %>%
  set("labels_colors") %>% set("labels_cex", c(.5,0.8)) %>% 
  set("nodes_pch", 19) %>% set("nodes_col", c("orange", "black", "plum", NA))
# plot the dend in usual "base" plotting engine:

plot(dend)
ggd1 <- as.ggdend(dend)

ggplot(ggd1, horiz = TRUE) # reproducing the above plot in ggplot2 :)

dend %>% highlight_branches_col %>% plot(main = "Coloring branches")
dend %>% highlight_branches %>% plot(main = "Emphasizing color\n and line-width")
