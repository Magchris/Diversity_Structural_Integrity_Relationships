##---------------------------------------------------------------------------------------------------------------------
## Estimating Forest Structural Integrity via Bootstrapping Method
## Title: Diversity/Structure Relationships across Tropical Biodiversity Experiments
## Author: Magnus Onyiriagwu 
## Supervised by: Clara Zemp
##---------------------------------------------------------------------------------------------------------------------

## load required libraries
pkgs <- c("data.table", "dplyr", "MetBrewer", "cowplot", "DT", "multcompView")
librarian::shelf(pkgs)


#######################################################################
## Import experiment and reference datasets 
rediv <- fread("data/core_rediv.csv")
rediv.ref <- fread("data/RedivRef.csv")


## auxiliary information
envr.preds <- fread("data/auxdata/envr_preds.csv", drop = "V1") 

## Add zone to dataset 
rediv <- merge.data.table(rediv, envr.preds[, .(Clim.zone, Zone, Site)], by = "Site")

## log transform species richness
rediv[, `:=` (SR.log2 = log(div.lev, 2))]
rediv[, plotID := .I]


#######################################################################
## Functions 
#######################################################################

## compute FSI
Int_fun <- function(ref, exp){
  denom <- (ref + exp)
  ifelse(denom == 0, NA, 1 - ((ref - exp) / denom))}

## Bootstrap function
boot_fun <- function(exp_value, ref_values, B = 1000) {
 
  Bvals <- numeric(B)
  
  for (b in seq_len(B)) {
    ref_boot_mean <- mean(sample(ref_values, size = length(ref_values), replace = TRUE))
    Bvals[b] <- Int_fun(ref_boot_mean, exp_value)
  }
  
  list(
    mean = mean(Bvals, na.rm = TRUE), se  = sd(Bvals, na.rm = TRUE))
}

#######################################################################
## Bootstrap stand structural integrity 
df.fsi <- rediv[, {
  ref_vals <- rediv.ref[Site == .BY$Zone, SSCI_comp]
  out <- boot_fun(exp_value = SSCI_comp, ref_values = ref_vals, B = 1000)
  
  .(FSI_mean = out$mean, FSI_se = out$se)
}, by = .(plotID, Site, Zone)]


## Bootstrap effective number of layers (integrity)
df.fenl <- rediv[, {
  ref_vals <- rediv.ref[Site == .BY$Zone, ENL]
  out <- boot_fun(exp_value = ENL, ref_values = ref_vals, B = 1000)
  
  .(FENL_mean = out$mean, FENL_se = out$se)
}, by = .(plotID, Site, Zone)]


## Bootstrap the mean fractal dimension (integrity)
df.frac <- rediv[, {
  ref_vals <- rediv.ref[Site == .BY$Zone, meanFrac]
  out <- boot_fun(exp_value = meanFrac, ref_values = ref_vals, B = 1000)
  
  .(FMFrac_mean = out$mean, FMFrac_se = out$se)
}, by = .(plotID, Site, Zone)]


## Merge data 
rediv <- merge(rediv, df.fsi, by = c("plotID", "Site", "Zone"), all.x = TRUE)
rediv <- merge(rediv, df.fenl, by = c("plotID", "Site", "Zone"), all.x = TRUE)
rediv <- merge(rediv, df.frac, by = c("plotID", "Site", "Zone"), all.x = TRUE)


###################################################
## Standardize
rediv[, `:=` (SSCI_sc = SSCI, meanFrac_sc = meanFrac, ENL_sc = ENL, FSI_sc = FSI_mean, FMFrac_sc = FMFrac_mean,
              FENL_sc = FENL_mean, SSCIpot_sc = SSCIpot, ENLpot_sc = ENLpot, MFracpot_sc = MFracpot)]

expts <- unique(rediv$Site)
preds <- c("SSCI_sc", "ENL_sc", "meanFrac_sc", "FSI_sc", "FENL_sc",
           "FMFrac_sc", "SSCIpot_sc", "ENLpot_sc", "MFracpot_sc")

for (i in expts) {
  for (j in preds) {
    
    rediv[Site == i, (j) := scale(get(j))[,1]]
  }
}

## save the dataset
write.csv(rediv, "data/datatable/rediv_fsi.csv", row.names = FALSE)

################################################################################
## Comparing Integrity Attributes between sites with one-way ANOVA
################################################################################

## Add zone to dataset 
rediv <- merge.data.table(rediv, envr.preds[, .(Clim.class, Age, Elev, WA, FI, Nspace, Site)], by = "Site")
rediv$Clim.class <- factor(rediv$Clim.class, levels = c("Af", "Am", "Cfa", "Aw", "Bwh"))
rediv$Site <- factor(rediv$Site, levels = c("Sabah","EBee","Sardinillia","AguaSalud","BEF-China",
                                            "CADE","MataDiv","UADY","BrazilDry","IDENT"))


# remove hyphen in site names (else, fails in multcompletters) 
rediv[Site == "BEF-China", Site := "BEF_China"]

sign.ls <- list()
for (name in names(rediv[, .(FSI_mean)])){
  # remove row with NA or Inf
  site_data <- rediv[is.finite(get(name)) & !is.na(get(name))]
  
  # execute anova 
  formula <- as.formula(paste(name, " ~ Site"))
  test <- aov(formula, data = site_data)
  result <- summary(test)[[1]]$'Pr(>F)'[1]
  
  # PostHoc Test
  if (result < 0.05){
    print(paste0(name, " is significantly different among the sites"))
    result.tukey <- TukeyHSD(test)
    signif <- multcompLetters4(test, result.tukey)
    sig.letters <- as.data.frame.list(signif$Site)
    sig.letters$Site <- rownames(sig.letters)
    sig.letters <- sig.letters[,  c("Site", "Letters")]
    row.names(sig.letters) <- NULL
    sign.ls[[name]] <- sig.letters
  } else (print(paste0(name, " is not significantly different among the sites")))
} 

## return transformed site name to original
rediv$Site[rediv$Site == "BEF_China"] <- "BEF-China"

## extract letters for ssci
lt.fsi <- sign.ls[["FSI_mean"]]
lt.fsi$Site[lt.fsi$Site == "BEF_China"] <- "BEF-China"
lt.fsi$Site <- factor(lt.fsi$Site, levels = c("Sabah","EBee","Sardinillia","AguaSalud","BEF-China","CADE","MataDiv","UADY","BrazilDry","IDENT"))

## merge 
rediv <- merge.data.table(rediv, lt.fsi, by = "Site")


# create label dataframe
labels_fsi <- rediv |> group_by(Site) |>
  summarise(FSI = quantile(FSI_mean, probs = 1, na.rm = TRUE), 
            Letters = unique(Letters),
            .groups = "drop") |> setDT()

site_mean <- rediv |> group_by(Site) |>
  summarise(mean_FSI = mean(FSI_mean, na.rm = TRUE),
    .groups = "drop") |> mutate(x = seq_len(n()))

fsi_summary <- rediv |> group_by(Site) |>
  summarise(mean_FSI = mean(FSI_mean, na.rm = TRUE),
            se_FSI   = mean(FSI_se, na.rm = TRUE),  # or first(FSI_se)
            .groups = "drop")

################################################################################
## plotting
pal1 <- rev(met.brewer("VanGogh3", 6))
pal2 <- rev(met.brewer("OKeeffe2", 5))
pal <- c(pal1[1], pal1[3], pal1[6], pal2[3], pal2[1])

(fsi_boxplot <-  rediv |>
    ggplot(aes(Site, FSI_mean, fill = Clim.class))+
    geom_violin() + 
    geom_point(data = fsi_summary, aes(x = Site, y = mean_FSI),inherit.aes = FALSE,size = 2) +
    geom_hline(yintercept = 1, linetype = 2, color = "gray50") +
    labs(y = "Stand structural Integrity", x = "Experiments", fill = "Climate Zone") +
    geom_text(data = labels_fsi, aes(x = Site, y = FSI, label = Letters), 
              position = position_dodge(width = .5), inherit.aes = FALSE, fontface = "bold", size = 5) +
    scale_fill_manual(values = pal) +
    theme_bw() +
    theme(
      legend.position = "right", 
      legend.title = element_text(face = "bold", size = 13),
      legend.text = element_text(face = "bold", size = 12),
      axis.title = element_text(face = "bold", size = 12), 
      axis.text.x = element_text(face = "bold", size = 12, angle = 30, hjust = 1),
      strip.text = element_text(face = "bold", size = 12), 
      axis.text.y = element_text(face = "bold", size = 12)))

## ---- Export Image
cowplot::save_plot(fsi_boxplot, 
                   filename = "output/Panels/fsi_boxplot.png", 
                   base_height = 4, base_width = 6)
