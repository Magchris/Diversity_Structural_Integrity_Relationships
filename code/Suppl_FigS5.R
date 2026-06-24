## Set up the path to the working directory
dir <- selectDirectory()


library(data.table)
library(FactoMineR)
library(factoextra)
library(ggplot2)
library(cowplot)
library(MetBrewer)
library(rstudioapi)
library(dplyr)

## Environmental variables
envr.preds <- fread("data/envr_preds.csv", drop = "V1")

# calculate Carbon/Nitrogen ratio
envr.preds[, C_N := 10 * C / N]

## Predictor Correlation
cor.pred <- corrplot::corrplot(cor(as.matrix(envr.preds[, .(prec, C_N, temp, WA, cec, pH, Elev, Age, Nspace)]), 
                                   method = "spearman"), method = 'circle', type = 'full', order = "hclust")

# run the PCA
pca.envr <- PCA(scale(envr.preds[, .(C_N, cec, pH)]), graph = F)


## Eigenvalues 
ind.envr <- get_pca_ind(pca.envr)
envr.preds$pc1 <- ind.envr$coord[,1]
envr.preds$pc2 <- ind.envr$coord[,2]
envr.preds$pc3 <- ind.envr$coord[,3]

## variance of each PC
eig_vals <- pca.envr$svd$vs^2
## Proportion of variance
var_prop <- round(eig_vals/sum(eig_vals) * 100, 2)


# plot 
(pc.envr <- fviz_pca_var(
  pca.envr, col.var = "contrib",      # color by contribution
  gradient.cols = c("grey70", "blue", "red"), repel = TRUE) +
  labs(title = "")+
    theme(axis.text = element_text(face = "bold"), 
          axis.title = element_text(face = "bold")))


## contribution
(pc.envr.contr <- fviz_contrib(pca.envr, choice = "var", axes = 1, fill = "gray30") + 
    labs(title = "") +
    theme(axis.text = element_text(face = "bold"), 
          axis.title = element_text(face = "bold")))
fviz_contrib(pca.envr, choice = "var", axes = 2)


## merge and export
(pc <- cowplot::plot_grid(pc.envr, pc.envr.contr,
                          rel_widths = c(.6,.4), rel_heights = c(.6,.4)))
save_plot(pc, filename = paste0(dir, "/scripts/result/chapt_two/Panels/pc_envr
                                      1.png"), dpi = 600)


