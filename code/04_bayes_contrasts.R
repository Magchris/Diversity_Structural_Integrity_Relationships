##---------------------------------------------------------------------------------------------------------------------
## Estimating Bayesian Contrasts
## Title: Code Script for Diversity/Structure Relationships across Tropical Biodiversity Experiments
## Author: Magnus Onyiriagwu Supervised by: Clara Zemp
##---------------------------------------------------------------------------------------------------------------------

## load required libraries
library(brms)
library(data.table)
library(ggplot2)

## Import structural complexity and context datasets
rediv <- fread("data/rediv_fsi.csv")

## load RData 
load("bayes_model.RData")

##  equal species richness prediction grid 
newdat <- data.frame(div.lev = sort(unique(rediv$div.lev)))

# log-transform as in the model
newdat$SR.log2 <- log2(newdat$div.lev)

## obtain posterior predictions
post <- posterior_epred(fit_fsi, newdata = newdat, re_formula = ~ 1 + SR.log2)

## calculate posterior difference between richness levels
delta <- post[, -1] - post[, -ncol(post)]

colnames(delta) <- paste0(
  newdat$div.lev[-1], "-",
  newdat$div.lev[-length(newdat$div.lev)]
)

## summarize posterior gains 
delta_summary <- data.frame(
  contrast = colnames(delta),
  mean = apply(delta, 2, mean),
  lwr = apply(delta, 2, quantile, 0.025),
  upr = apply(delta, 2, quantile, 0.975)
)
delta_summary 

## calculate marginal gain per added species
richness_steps <- diff(newdat$div.lev)
delta_per_species <- sweep(delta, 2, richness_steps,"/")

# summarize
gain_summary <- data.frame(
  contrast = colnames(delta_per_species),
  mean = apply(delta_per_species, 2, mean),
  lwr = apply(delta_per_species, 2, quantile, 0.025),
  upr = apply(delta_per_species, 2, quantile, 0.975)
)
gain_summary$contrast <- factor(gain_summary$contrast,levels = unique(gain_summary$contrast))

## estimate the richness level at which most of the effect is achieved
post_rel <- sweep(post, 1, post[, ncol(post)], "/")
apply(post_rel, 2, mean)

## visualization 
(figS2 <- ggplot(gain_summary,aes(x = contrast, y = mean)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lwr, ymax = upr),width = 0.1) +
  labs(y = expression(bold(paste("SS"["int"]," gains per increase in diversity levels")))) +
    theme_bw())

cowplot::save_plot(figS2, filename = "output/Panels/figS2.png", dpi = 300)
