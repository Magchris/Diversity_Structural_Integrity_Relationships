##---------------------------------------------------------------------------------------------------------------------
## Structural Equation Modelling
## Title: Code Script for Diversity/Structure Relationships across Tropical Biodiversity Experiments
## Author: Magnus Onyiriagwu Supervised by: Clara Zemp
##---------------------------------------------------------------------------------------------------------------------

## load required libraries
pkgs <- c("data.table", "nls2", "dplyr", "lme4", "lmerTest", "cowplot", "piecewiseSEM")
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
rediv <- merge.data.table(rediv, envr.preds[, .(Clim.class, Site, Age, Elev, WA, temp, FI, Nspace)], by = "Site")
rediv$Clim.class <- factor(rediv$Clim.class, levels = c("Af", "Am", "Cfa", "Aw", "Bwh"))
rediv$Site <- factor(rediv$Site, levels = c("Sabah","EBee","Sardinillia","AguaSalud","BEF-China",
                                            "CADE","MataDiv","UADY","BrazilDry","IDENT"))

## log transform species richness and factorize species composition
rediv[, `:=` (SR.log2 = log(div.lev, 2))]

#######################################################################
## Compute composite structural integrity 


## fit SSCI composite model using experimental plots
start <- expand.grid(
  a = seq(0.1, 1, by = 0.1),
  b = seq(0.1, 3, by = 0.1)
)

FSI_mod <- nls2(FSI_mean ~ a * FMFrac_mean^b * FENL_mean, data = rediv, start = start)

summary(FSI_mod)
plot(rediv$FSI_mean ~ predict(FSI_mod))

beta.a <- summary(FSI_mod)$coefficients[1, 1]
beta.b <- summary(FSI_mod)$coefficients[2, 1]

## compute composite SSCI for experimental plots
rediv$FSI_comp <- beta.a * rediv$FMFrac_mean^beta.b * rediv$FENL_mean

rediv <- rediv[complete.cases(rediv)]


################################################################################
## structural equation modelling


# Stand age

mod1.1 <- psem(
  lmer(TopH ~ SR.log2:Age + (1|Site), data = rediv),
  lmer(FENL_mean ~ SR.log2:Age + TopH  + (1|Site), data = rediv),
  lmer(FMFrac_mean ~ SR.log2:Age + TopH  + (1|Site), data = rediv),
  lmer(FSI_comp ~ SR.log2:Age + FMFrac_mean + FENL_mean + TopH + (1|Site), data = rediv))
dSep(mod1.1)


# unmeasured correlated errors
summary(update(mod1, FENL_mean %~~% FMFrac_mean))


##############################

# Soil fertility 
mod2.1 <- psem(
  lmer(TopH ~ SR.log2:FI + (1|Site), data = rediv),
  lmer(FENL_mean ~ SR.log2:FI + TopH  + (1|Site), data = rediv),
  lmer(FMFrac_mean ~ SR.log2:FI + TopH  + (1|Site), data = rediv),
  lmer(FSI_comp ~ SR.log2:FI + FMFrac_mean + FENL_mean + TopH + (1|Site), data = rediv))
dSep(mod2.1)

# unmeasured correlated errors
summary(update(mod2.1, FENL_mean %~~% FMFrac_mean))


##############################

# Water availability 

mod3.1 <- psem(
  lmer(TopH ~ SR.log2:WA + (1|Site), data = rediv),
  lmer(FENL_mean ~ SR.log2:WA + TopH  + (1|Site), data = rediv),
  lmer(FMFrac_mean ~ SR.log2:WA + TopH  + (1|Site), data = rediv),
  lmer(FSI_comp ~ SR.log2:WA + FMFrac_mean + FENL_mean + TopH + (1|Site), data = rediv))
dSep(mod3.1)

# unmeasured correlated errors

summary(update(mod3.1, FENL_mean %~~% FMFrac_mean))


##############################

# Planting Density

mod4.1 <- psem(
  lmer(TopH ~ SR.log2:Nspace + (1|Site), data = rediv),
  lmer(FENL_mean ~ SR.log2:Nspace + TopH  + (1|Site), data = rediv),
  lmer(FMFrac_mean ~ SR.log2:Nspace + TopH  + (1|Site), data = rediv),
  lmer(FSI_comp ~ SR.log2:Nspace + FMFrac_mean + FENL_mean + TopH + (1|Site), data = rediv))
dSep(mod4.1)

# unmeasured correlated errors
summary(update(mod4.1, FENL_mean %~~% FMFrac_mean))


