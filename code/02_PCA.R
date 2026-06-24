##---------------------------------------------------------------------------------------------------------------------
## Principal Component Analysis of the SSCI, meanFrac, and ENL
## Title: Code Script for Diversity/Structure Relationships across Tropical Biodiversity Experiments
## Author: Magnus Onyiriagwu Supervised by: Clara Zemp
##---------------------------------------------------------------------------------------------------------------------

## load libraries 
library(data.table)
library(FactoMineR)
library(factoextra)
library(ggplot2)
library(cowplot)
library(MetBrewer)
library(dplyr)


#######################################################################

#--- Import the dataset 
rediv <- fread("data/datatable/rediv_fsi.csv")
rediv.ref <- fread("data/datatable/RedivRef.csv")

## Explanatory variables
envr.preds <- fread("data/auxdata/envr_preds.csv", drop = "V1") 

# remove missing plots
rediv <- df.rediv[complete.cases(rediv)]

# run the PCA
pca.scale <- PCA(scale(rediv[, .(ENL, SSCI, meanFrac)]), graph = F)

## Eigenvalues 
pca.ind <- get_pca_ind(pca.scale)
## variance of each PC
eig_vals <- pca.scale$svd$vs^2
## Proportion of variance
var_prop <- round(eig_vals/sum(eig_vals) * 100, 2)


## Piloting 

# Set color
pal1 <- rev(met.brewer("VanGogh3", 6))
pal2 <- rev(met.brewer("OKeeffe2", 5))
pal <- c(pal1[1], pal1[3], pal1[6], pal2[3], pal2[1])

## -- PCA Biplot
(pcBiplot <- fviz_pca_biplot(pca.scale, fill.ind = df.rediv$Clim.class, col.ind = "white",
                             # shape.ind = df.rediv$veg.type, habillage = df.rediv$Clim.class, 
                             palette = pal, geom = "point", 
                             label = "var", col.var = "black", labelsize = 4, 
                             repel = TRUE, arrowsize = 0.1,  arrow.color = "black",
                             pointshape = 21, pointsize = 2, alpha.ci = 0.5, 
                             addEllipses = T, ellipse.color = df.rediv$Clim.class, alpha.ellipse = 0.5) + 
    
    labs(fill = "Climate", title = "Experiments") + 
    xlab(paste0("PC1 (", var_prop[1], "%)")) +
    ylab(paste0("PC1 (", var_prop[2], "%)")) +                   
    scale_color_manual(values = pal, guide = "none") +
    theme_bw() +
    theme(#legend.position = "none", 
          legend.position = c(.9, .85),
          legend.background = element_blank(),
          legend.key.size = unit(0.02, "cm"),
          legend.key.width = unit(0.05,"cm"),
          legend.text = element_text(face = "bold", size = 12),
          legend.title = element_text(face = "bold", size = 13),
          plot.title = element_text(face = "bold", size = 13),
          axis.title = element_text(size = 14, face = "bold"),
          axis.text = element_text(size = 14, face = "bold"), 
          panel.border = element_blank()) + 
     coord_cartesian(xlim = c(-4,7), ylim = c(-5.5,7)))

# Export panel
save_plot(pcBiplot, filename = paste0(dir, "output/Panels/pcBiplot.png"), dpi = 600)

