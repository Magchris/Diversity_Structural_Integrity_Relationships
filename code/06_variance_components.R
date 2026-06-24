##---------------------------------------------------------------------------------------------------------------------
## Variance decomposition using Bayesian model
## Title: Code Script for Diversity/Structure Relationships across Tropical Biodiversity Experiments
## Author: Magnus Onyiriagwu Supervised by: Clara Zemp
##---------------------------------------------------------------------------------------------------------------------

## load required libraries
pkgs <- c("data.table", "arm",  "ggplot2", "cowplot",  "brms", "cmdstanr")
librarian::shelf(pkgs)

#######################################################################

## Import structural complexity and context datasets
rediv <- fread("data/datatable/rediv_fsi.csv")

# create site x composition variable 
rediv$Site.comp <- factor(paste(rediv$Site, rediv$comp, sep = "_"))

## Estimate Posterior Distribution 
fit_fsi_comp <- brm(
  formula = FSI_mean ~ SR.log2 + (1|Site.comp) + (1 + SR.log2 | Site) + (1|comp),
  data = rediv,
  family = gaussian(),
  prior = c(
    prior(normal(0, 2), class = "Intercept"),
    prior(normal(0, 1), class = "b", coef = "SR.log2"),
    prior(student_t(3, 0, 1), class = "sd", group = "Site"),
    prior(student_t(3, 0, 1), class = "sd", group = "comp"),
    prior(student_t(3, 0, 1), class = "sigma")
  ),
  chains = 4, iter = 6000, warmup = 3000, cores = 4, seed = 123, backend = "cmdstanr"
)



#######################################################################
# Extract posterior draws
draws <- as_draws_df(fit_fsi_comp)

## data vectors
dat <- rediv

sr_vec  <- dat$SR.log2
site_id <- dat$Site
comp_id <- dat$comp
site.comp_id <- dat$Site.comp

## levels
site_levels <- levels(dat$Site)
comp_levels <- levels(dat$comp)
site_comp_levels <- levels(dat$Site.comp)

## extract group-level effects from brms draws
get_re_matrix <- function(draws_df, prefix, levels_vec, coef_name = "Intercept") {
  out <- matrix(NA_real_, nrow = nrow(draws_df), ncol = length(levels_vec))
  colnames(out) <- levels_vec
  
  for (j in seq_along(levels_vec)) {
    nm <- paste0("r_", prefix, "[", levels_vec[j], ",", coef_name, "]")
    if (!nm %in% names(draws_df)) {
      stop("Could not find posterior draw column: ", nm)
    }
    out[, j] <- draws_df[[nm]]
  }
  out
}

# random intercepts for Site
re_site_int <- get_re_matrix(draws, "Site", site_levels, "Intercept")

# random slopes for SR.log2 within Site
re_site_slope <- get_re_matrix(draws, "Site", site_levels, "SR.log2")

# random intercepts for comp
re_comp_int <- get_re_matrix(draws, "comp", comp_levels, "Intercept")

# random intercepts for Site x comp
re_site_comp_int <- get_re_matrix(draws, "Site.comp", site_comp_levels, "Intercept")

# fixed richness coefficient
beta_sr <- draws$b_SR.log2

# residual SD
sigma_y <- draws$sigma

## compute posterior SD components draw by draw
comp_post <- lapply(seq_len(nrow(draws)), function(i) {
  
  # fixed effect contribution of species richness
  # place it on the response scale as variation across the observed richness gradient
  s_SR <- abs(beta_sr[i]) * sd(sr_vec, na.rm = TRUE)
  
  # site random intercept contribution
  s_SITE <- sd(re_site_int[i, ], na.rm = TRUE)
  
  # composition random intercept contribution
  s_COMP <- sd(re_comp_int[i, ], na.rm = TRUE)
  
  # expeirment x composition random intercept contribution
  s_SITECOMP <- sd(re_site_comp_int[i, ], na.rm = TRUE)
  
  
  # site x richness contribution:
  # site-specific slope deviations translated across observed plot richness values
  sr_site_dev <- re_site_slope[i, match(site_id, site_levels)] * sr_vec
  s_SITESR <- sd(sr_site_dev, na.rm = TRUE)
  
  # extra residual
  s_RES <- sigma_y[i]
  
  # total variance on SD-component scale
  var_total <- s_SR + s_SITE + s_COMP + s_SITECOMP + s_SITESR + s_RES
  
  c(
    s_SR     = s_SR,
    s_SITE   = s_SITE,
    s_COMP   = s_COMP,
    s_SITECOMP = s_COMP,
    s_SITESR = s_SITESR,
    s_RES    = s_RES,
    p_SR     = s_SR / var_total,
    p_SITE   = s_SITE / var_total,
    p_COMP   = s_COMP / var_total,
    p_SITECOMP   = s_SITECOMP / var_total,
    p_SITESR = s_SITESR / var_total,
    p_RES    = s_RES / var_total
  )
})

comp_post <- do.call(rbind, comp_post)
comp_post <- as.data.frame(comp_post)

## posterior summary function
summ_fun <- function(x) {
  c(
    mean  = round(mean(x, na.rm = TRUE), 3),
    sd    = round(sd(x, na.rm = TRUE), 3),
    q2.5  = unname(round(quantile(x, 0.025, na.rm = TRUE), 3)),
    q16   = unname(round(quantile(x, 0.16,  na.rm = TRUE), 3)),
    q50   = unname(round(quantile(x, 0.50,  na.rm = TRUE), 3)),
    q84   = unname(round(quantile(x, 0.84,  na.rm = TRUE), 3)),
    q97.5 = unname(round(quantile(x, 0.975, na.rm = TRUE), 3))
  )
}

varcomp_summary <- rbind(
  Species_richness    = summ_fun(comp_post$s_SR),
  Site                = summ_fun(comp_post$s_SITE),
  Species_composition = summ_fun(comp_post$s_COMP),
  SR_x_Site           = summ_fun(comp_post$s_SITESR),
  Comp_x_Site         = summ_fun(comp_post$s_SITECOMP),
  Residual            = summ_fun(comp_post$s_RES)
)

prop_summary <- rbind(
  Species_richness    = summ_fun(comp_post$p_SR),
  Site                = summ_fun(comp_post$p_SITE),
  Species_composition = summ_fun(comp_post$p_COMP),
  SR_x_Site           = summ_fun(comp_post$p_SITESR),
  Comp_x_Site         = summ_fun(comp_post$p_SITECOMP),
  Residual            = summ_fun(comp_post$p_RES)
)

round(varcomp_summary, 3)
round(prop_summary, 3)

## plot posterior SD components
plot_df <- as.data.frame(varcomp_summary)
plot_df$component <- rownames(plot_df)
plot_df <- plot_df[order(plot_df$mean), ]
plot_df$component <- factor(plot_df$component, levels = plot_df$component)
rownames(plot_df) <- NULL

## Plotting 
var.comp <- ggplot(plot_df, aes(x = mean, y = component)) +
  geom_errorbarh(aes(xmin = q2.5, xmax = q97.5), height = 0, linewidth = 0.6) +
  geom_errorbarh(aes(xmin = q16, xmax = q84), height = 0, linewidth = 2) +
  geom_point(size = 3) +
  geom_vline(xintercept = 0, colour = "grey50", linetype = 2) +
  labs(x = "Posterior SD FSI", y = NULL) +
  scale_y_discrete(
    limits = c("Site", "Residual", "Species_composition", "Comp_x_Site", "SR_x_Site", "Species_richness"),
    labels = c("Site" = "Expt", "Residual" = "Residual", "Species_composition" = "Comp", 
               "Comp_x_Site" = "Comp x Expt", "SR_x_Site" = "SR x Expt", "Species_richness" = "SR")) +
  theme_classic() +
  theme(
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(size = 11, face = "bold"),
    axis.title = element_text(face = "bold"))


## Export Image
save_plot(bayes.plot, filename = "output/Panels/bayesplot.png", dpi = 600)

