##---------------------------------------------------------------------------------------------------------------------
## Bayesian model
## Title: Code Script for Diversity/Structure Relationships across Tropical Biodiversity Experiments
## Author: Magnus Onyiriagwu Supervised by: Clara Zemp
##---------------------------------------------------------------------------------------------------------------------

## load required libraries
pkgs <- c("data.table", "dplyr", "MetBrewer", "ggplot2", "cowplot", "brms", "cmdstanr")
librarian::shelf(pkgs)


#######################################################################
## Import structural complexity and context datasets
rediv <- fread("data/rediv_fsi.csv")
envr.preds <- fread("data/envr_preds.csv", drop = "V1")

## standardize predictors 
envr.preds[, `:=` (Age_sc = Age, WA_sc = WA, FI_sc = FI, Nspace_sc=Nspace)]
subset <- c("Age_sc", "WA_sc", "FI_sc", "Nspace_sc")
for (c in subset) envr.preds[[c]] <- scale(envr.preds[[c]])[,1]

## merge context variables to structural data
rediv <- merge.data.table(rediv, envr.preds, by = "Site")
rediv$Clim.class <- factor(rediv$Clim.class, levels = c("Af", "Am", "Cfa", "Aw", "Bwh"))
rediv$Site <- factor(rediv$Site, levels = c("Sabah","EBee","Sardinillia","AguaSalud","BEF-China",
                                              "CADE","MataDiv","UADY","BrazilDry","IDENT"))
rediv[, comp := factor(comp)]


#######################################################################
## Modelling

# Base model 
fit_fsi <- brm(
  formula = FSI_mean | se(FSI_se, sigma = TRUE) ~ SR.log2 + (1 + SR.log2 | Site),
  data = rediv,
  family = gaussian(),
  prior = c(
    prior(normal(0, 2), class = "Intercept"),
    prior(normal(0, 1), class = "b", coef = "SR.log2"),
    prior(student_t(3, 0, 1), class = "sd", group = "Site"),
    prior(student_t(3, 0, 1), class = "sigma")
  ),
  chains = 4, iter = 6000, warmup = 3000, cores = 4, seed = 123, backend = "cmdstanr")

##########################################################################
## Base model x interaction effect 


# Stand age 
fit_fsi_age <- brm(
  formula = FSI_sc ~ SR.log2 * Age_sc + (1 + SR.log2 | Site) + (1|comp),
  data = rediv,
  family = gaussian(),
  seed = 123, chains = 4, iter = 6000, warmup = 3000,
  cores = 4, backend = "cmdstanr",
  control = list(adapt_delta = 0.99, max_treedepth = 12)
)
summary(fit_fsi_age)

## Water availability
fit_fsi_wa <- brm(
  formula = FSI_sc ~ SR.log2 * WA_sc + (1 + SR.log2 | Site) + (1|comp),
  data = rediv,
  family = gaussian(),
  seed = 123, chains = 4, iter = 6000, warmup = 3000,
  cores = 4, backend = "cmdstanr",
  control = list(adapt_delta = 0.99, max_treedepth = 12)
)
summary(fit_fsi_wa)

## Soil fertility
fit_fsi_fi <- brm(
  formula = FSI_sc ~ SR.log2 * FI + (1 + SR.log2 | Site) + (1|comp),
  data = rediv,
  family = gaussian(),
  seed = 123, chains = 4, iter = 6000, warmup = 3000,
  cores = 4, backend = "cmdstanr",
  control = list(adapt_delta = 0.99, max_treedepth = 12)
)

summary(fit_fsi_fi)


## Growing space
fit_fsi_space <- brm(
  formula = FSI_sc ~ SR.log2 * Nspace_sc + (1 + SR.log2 | Site) + (1|comp),
  data = rediv,
  family = gaussian(),
  seed = 123, chains = 4, iter = 6000, warmup = 3000,
  cores = 4, backend = "cmdstanr",
  control = list(adapt_delta = 0.99, max_treedepth = 12)
)

summary(fit_fsi_space)

########################################################################
# Extract the model coefficients
########################################################################

## --- Base model 

# coefficients
fsi_site_coefs <- coef(fit_fsi)$Site

# random slopes
fsi_site_slopes <- data.table(
  Site = rownames(fsi_site_coefs[, , "SR.log2"]),
  Slope = fsi_site_coefs[, "Estimate", "SR.log2"],
  Slope_se = fsi_site_coefs[, "Est.Error", "SR.log2"],
  Q2.5 = fsi_site_coefs[, "Q2.5", "SR.log2"],
  Q97.5 = fsi_site_coefs[, "Q97.5", "SR.log2"]
)

# intercept
fsi_site_intercept <- data.table(
  Site = rownames(fsi_site_coefs[, , "Intercept"]),
  Intercept = fsi_site_coefs[, "Estimate", "Intercept"],
  Int_se = fsi_site_coefs[, "Est.Error", "Intercept"],
  Q2.5 = fsi_site_coefs[, "Q2.5", "Intercept"],
  Q97.5 = fsi_site_coefs[, "Q97.5", "Intercept"]
)

# posterior fixed effect
fsi_fix <- data.table(
  intercept = fixef(fit_fsi1)["Intercept", "Estimate"],
  slope = fixef(fit_fsi1)["SR.log2", "Estimate"]
)

# combine slopes with intercepts 
fsi_site_slopes <- merge.data.table(fsi_site_slopes, fsi_site_intercept[, .(Intercept, Site)], by = "Site")


## Add zone to dataset 
fsi_site_slopes <- setDT(merge.data.table(fsi_site_slopes, envr.preds[, .(Clim.class, Age, Elev, WA, FI, Nspace, Site)], 
                                          by = "Site"))
fsi_site_slopes$Clim.class <- factor(fsi_site_slopes$Clim.class, levels = c("Af", "Am", "Cfa", "Aw", "Bwh"))
fsi_site_slopes$Site <- factor(fsi_site_slopes$Site, 
                               level = c("Sabah","EBee","Sardinillia","AguaSalud","BEF-China",
                                         "CADE","MataDiv","UADY","BrazilDry","IDENT"))



## --- Interaction models with context covariates


## slope coefficients 
fsi_site_coefs_age <- coef(fit_fsi_age)$Site
fsi_site_coefs_wa <- coef(fit_fsi_wa)$Site
fsi_site_coefs_space <- coef(fit_fsi_space)$Site
fsi_site_coefs_fi <- coef(fit_fsi_fi)$Site

# interaction effect
fsi_envr <- data.table(
  Slope = c(fsi_site_coefs_age[, "Estimate", "SR.log2:Age"][1], 
            fsi_site_coefs_wa[, "Estimate", "SR.log2:WA"][1],
            fsi_site_coefs_space[, "Estimate", "SR.log2:Nspace"][1],
            fsi_site_coefs_fi[, "Estimate", "SR.log2:FI"][1]), 
  Slope_se = c(fsi_site_coefs_age[, "Est.Error", "SR.log2:Age"][1],
               fsi_site_coefs_wa[, "Est.Error", "SR.log2:WA"][1],
               fsi_site_coefs_space[, "Est.Error", "SR.log2:Nspace"][1],
               fsi_site_coefs_fi[, "Est.Error", "SR.log2:FI"][1]),
  Q2.5 = c(fsi_site_coefs_age[, "Q2.5", "SR.log2:Age"][1],
           fsi_site_coefs_wa[, "Q2.5", "SR.log2:WA"][1],
           fsi_site_coefs_space[, "Q2.5", "SR.log2:Nspace"][1],
           fsi_site_coefs_fi[, "Q2.5", "SR.log2:FI"][1]),
  Q97.5 = c(fsi_site_coefs_age[, "Q97.5", "SR.log2:Age"][1],
            fsi_site_coefs_wa[, "Q97.5", "SR.log2:WA"][1],
            fsi_site_coefs_space[, "Q97.5", "SR.log2:Nspace"][1],
            fsi_site_coefs_fi[, "Q97.5", "SR.log2:FI"][1]), 
  var = c("Age", "WA", "Nspace", "FI")
)

# fixed effect
fsi_envr_fix <- data.table(
  Slope = c(fsi_site_coefs_age[, "Estimate", "Age_sc"][1], 
            fsi_site_coefs_wa[, "Estimate", "WA_sc"][1],
            fsi_site_coefs_space[, "Estimate", "Nspace_sc"][1],
            fsi_site_coefs_fi[, "Estimate", "FI"][1]), 
  Slope_se = c(fsi_site_coefs_age[, "Est.Error", "Age_sc"][1],
               fsi_site_coefs_wa[, "Est.Error", "WA_sc"][1],
               fsi_site_coefs_space[, "Est.Error", "Nspace_sc"][1],
               fsi_site_coefs_fi[, "Est.Error", "FI"][1]),
  Q2.5 = c(fsi_site_coefs_age[, "Q2.5", "Age_sc"][1],
           fsi_site_coefs_wa[, "Q2.5", "WA_sc"][1],
           fsi_site_coefs_space[, "Q2.5", "Nspace_sc"][1],
           fsi_site_coefs_fi[, "Q2.5", "FI"][1]),
  Q97.5 = c(fsi_site_coefs_age[, "Q97.5", "Age_sc"][1],
            fsi_site_coefs_wa[, "Q97.5", "WA_sc"][1],
            fsi_site_coefs_space[, "Q97.5", "Nspace_sc"][1],
            fsi_site_coefs_fi[, "Q97.5", "FI"][1]), 
  var = c("Age_sc", "WA_sc", "Nspace_sc", "FI")
)


########################################################################
## plotting 
########################################################################


## set colors
# site
pal1 <- rev(met.brewer("VanGogh3", 12))
pal2 <- rev(met.brewer("OKeeffe2", 12))
pal.site <- c(pal1[1], pal1[3], pal1[5], pal1[7], pal1[9], 
              pal1[11], pal1[12], pal2[10], pal2[8], pal2[1])

# climate
pal1 <- rev(met.brewer("VanGogh3", 6))
pal2 <- rev(met.brewer("OKeeffe2", 5))
pal.clim <- c(pal1[1], pal1[3], pal1[6], pal2[3], pal2[1])


## --- Effect size base model 


## Generate confidence interval for the fixed slope 
newdat <- data.frame(SR.log2 = xseq, Site = NA)

ep <- posterior_epred(fit_fsi1, newdata = newdat, re_formula = NA)

# prediction 
pred_fix <- data.table(
  SR.log2 = xseq,
  fit = apply(ep, 2, mean),
  lwr = apply(ep, 2, quantile, probs = 0.025),
  upr = apply(ep, 2, quantile, probs = 0.975)
)

# plot standard coefficients across the sites


(fig4A <- ggplot(rediv, aes(x = SR.log2, y = FSI_sc, color = Clim.class)) +
    geom_abline(data = fsi_site_slopes, aes(intercept = Intercept, slope = Slope, color = Clim.class), linewidth = .8) +
    geom_ribbon(data = pred_fix, aes(x = SR.log2, ymin = lwr, ymax = upr), inherit.aes = FALSE, alpha = .2) +
    geom_line(data = pred_fix, aes(x = SR.log2, y = fit), inherit.aes = FALSE, colour = "black", linewidth = 1, linetype = 2) +
    scale_color_manual(values = pal.clim) +
    labs(y = expression(bold(paste("Stand Structural Integrity"))),
         x = expression(bold(paste("Species Richness (log"["2"], " scale)"))), 
         color = "Climate") +
    scale_x_continuous(breaks=c(0,1,2,3,4), labels=c("1","2","4","8","16")) + 
    theme_bw() + coord_cartesian(ylim = c(-.35,1), xlim = c(.1, 4)) + 
    theme(
      legend.position = c(.2,.7),
      legend.background = element_blank(),
      legend.key.spacing = unit(0.2, "cm"),
      legend.key.size = unit(0.2, "cm"), 
      legend.spacing = unit(0.1, "cm"),
      legend.text = element_text(size = 10, face = "bold"),
      legend.title = element_text(face = "bold", size = 12),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"), 
      panel.grid.minor = element_blank(), 
      axis.title = element_text(face = "bold", size = 14),
      axis.text = element_text(face = "bold", size = 12)))


(fig4B <- 
    ggplot(data = fsi_site_slopes, aes(x = Slope, y = Site)) +
    geom_errorbar(aes(xmin = Q2.5, xmax = Q97.5), width = 0) +  # Error bars
    geom_point(aes(color = Clim.class), size = 3) +  # Mean points
    geom_vline(xintercept = 0, linetype = 2) +
    scale_color_manual(values = pal.site) +
    labs( y = "", x = "Structural Integrity Slope",  color = "Climate") +
    theme_bw() +  
    theme(legend.position = c(.85, .2),
          legend.background = element_blank(),
          legend.key = element_blank(), 
          legend.key.height = unit(0.05, "cm"), 
          legend.key.width = unit(0.05, "cm"), 
          legend.key.size = unit(0.05, "cm"), 
          legend.spacing = unit(0.0005, "cm"),
          legend.key.spacing = unit(0.0005, "cm"),
          legend.spacing.y = unit(0.001, "cm"),
          legend.title = element_text(face = "bold", size = 12),
          legend.text = element_text(face = "bold", size = 11),
          axis.text.x = element_text(face = "bold", hjust = 1, size = 12), 
          axis.text.y = element_text(face = "bold", size = 12),
          axis.title = element_text(face = "bold", size = 12)))


(fig4.grid <- plot_grid(fig4A, fig4B, rel_heights = c(.45, .55)))

## Export panel as png
save_plot(fig4.grid, base_height = 6, base_width = 3,
          filename = "output/Panels/fig4.png", dpi = 600)


############################################################################################
## --- Effect size interaction effect


# extract interaction effect
ce_age_df <- ce_age[["SR.log2:Age_sc"]]
ce_wa_df <- ce_wa[["SR.log2:WA_sc"]]
ce_space_df <- ce_space[["SR.log2:Nspace_sc"]]
ce_fi_df <- ce_fi[["SR.log2:FI"]]



## estimate conditional effect 
ce_age_df |> group_by(effect2__) |> summarise(meanEst = mean(estimate__))
ce_wa_df |> group_by(effect2__) |> summarise(meanEst = mean(estimate__))
ce_space_df |> group_by(effect2__) |> summarise(meanEst = mean(estimate__))
ce_fi_df |> group_by(effect2__) |> summarise(meanEst = mean(estimate__))



# plot interaction effect size (i.e. Spp. richness x context)
(fig5D <- ggplot(data = fsi_envr_fix, aes(x = var, y = Slope)) + 
    geom_point() + geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width = .1) + 
    geom_hline(yintercept = 0, linetype = 2, color = "gray70") + 
    labs(x = "Context", y = expression(bold(paste("Effect Size (SS"["int"],")"))), 
         title = "Intercept") + 
    theme_bw() + 
    theme(axis.text = element_text(face = "bold", size = 12), 
          title = element_text(face = "bold", size = 10),
          axis.title = element_text(face = 'bold', size = 13)))

save_plot(fig5D, filename = "output/Panels/fig5D", dpi = 300, base_height = 3, base_width = 3)

## plot conditional effects

(fig5.age <- ggplot(ce_age_df, aes(x = SR.log2, y = estimate__)) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__, fill = factor(round(Age_sc,0))), alpha = 0.25, colour = NA) +
  geom_line(aes(color = factor(round(Age_sc,0))), linewidth = 1) +
  scale_color_manual(values = c(pal1[1],  pal1[6], pal2[1]),
                     labels = c("7", "13", "18")) +
  scale_fill_manual(values = c(pal1[1],  pal1[6], pal2[1]),
                    labels = c("7", "13", "18")) +
  labs( x = "",
        y = expression(bold(paste("Stand Structural Integrity"))), colour = "Age", fill = "Age") +
  theme_bw() + 
    theme(legend.position = c(.2, .8), 
          legend.background = element_blank(), 
          legend.title = element_text(face = "bold", size = 14),
          legend.text = element_text(face = "bold", size = 12),
          axis.text.x = element_text(face = "bold", hjust = 1, size = 14), 
          axis.text.y = element_text(face = "bold", size = 14),
          axis.title = element_text(face = "bold", size = 13)))

(fig5.wa <- ggplot(ce_wa_df, aes(x = SR.log2, y = estimate__)) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__, fill = factor(round(WA_sc,0))), alpha = 0.25, colour = NA) +
  geom_line(aes(color = factor(round(WA_sc,0))), linewidth = 1) +
  scale_color_manual(values = c(pal1[1],  pal1[6], pal2[1]),
                       labels = c("low", "mid", "high")) +
  scale_fill_manual(values = c(pal1[1],  pal1[6], pal2[1]),
                      labels = c("low", "mid", "high")) +
  labs( x = "", y = "", colour = "Water Availability", fill = "Water Availability") +
  theme_bw() + coord_cartesian(ylim = c(-.5,1.5)) +
  theme(legend.position = c(.33, .8), 
        legend.background = element_blank(), 
        legend.title = element_text(face = "bold", size = 14),
        legend.text = element_text(face = "bold", size = 12),
        axis.text.x = element_text(face = "bold", hjust = 1, size = 14), 
        axis.text.y = element_text(face = "bold", size = 14),
        axis.title = element_text(face = "bold", size = 13)))

(fig5.space <- ggplot(ce_space_df, aes(x = SR.log2, y = estimate__)) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__, fill = factor(round(Nspace_sc,0))), alpha = 0.25, colour = NA) +
  geom_line(aes(color = factor(round(Nspace_sc,0))), linewidth = 1) +
  scale_color_manual(values = c(pal1[1],  pal1[6], pal2[1]),
                       labels = c("low", "mid", "high")) +
  scale_fill_manual(values = c(pal1[1],  pal1[6], pal2[1]),
                      labels = c("low", "mid", "high")) +
  labs( x = expression(bold(paste("Species Richness (log"[2], " scale)"))),
        y = expression(bold(paste("Stand Structural Integrity"))), colour = "Growing Space", fill = "Growing Space") +
  theme_bw() + coord_cartesian(ylim = c(-.5,1.5)) +
    theme(legend.position = c(.33, .8), 
          legend.background = element_blank(), 
          legend.title = element_text(face = "bold", size = 14),
          legend.text = element_text(face = "bold", size = 12),
          axis.text.x = element_text(face = "bold", hjust = 1, size = 14), 
          axis.text.y = element_text(face = "bold", size = 14),
          axis.title = element_text(face = "bold", size = 13)))

(fig5.fi <- ggplot(ce_fi_df, aes(x = SR.log2, y = estimate__)) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__, fill = factor(round(FI,0))), alpha = 0.25, colour = NA) +
  geom_line(aes(color = factor(round(FI,0))), linewidth = 1) +
  scale_color_manual(values = c(pal1[1],  pal1[6], pal2[1]),
                       labels = c("low", "mid", "high")) +
  scale_fill_manual(values = c(pal1[1],  pal1[6], pal2[1]),
                      labels = c("low", "mid", "high")) +
  labs( x = expression(bold(paste("Species Richness (log"[2], " scale)"))),
        y = "", colour = "Soil Fertility", fill = "Soil Fertility") +
  theme_bw() + coord_cartesian(ylim = c(-.5,1.5)) +
    theme(legend.position = c(.2, .8), 
          legend.background = element_blank(), 
          legend.title = element_text(face = "bold", size = 14),
          legend.text = element_text(face = "bold", size = 12),
          axis.text.x = element_text(face = "bold", hjust = 1, size = 14), 
          axis.text.y = element_text(face = "bold", size = 14),
          axis.title = element_text(face = "bold", size = 13)))

## save and export
fig5.grid <- cowplot::plot_grid(fig5.age, fig5.wa, fig5.space, fig5.fi, nrow = 2)
save_plot(fig5.grid, filename = "output/Panels/fig5.png", dpi = 300, base_width = 8, base_height = 7)

## save RData to compute contrast
save.image("bayes_model.RData")
