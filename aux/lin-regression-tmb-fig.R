## header ---------------------------------------------------------------------

options(scipen = 999)
options(stringsAsFactors = F)
rm(list = ls())
library(ggplot2)
library(data.table)
library(this.path)
library(scriptName)

## settings -------------------------------------------------------------------

### fixed settings
dirBase <- dirname(this.dir())
dirdtMeta <- file.path(dirBase, "res")
dirPlots <- file.path(dirBase, "plots", "fig")
dir.create(dirPlots, recursive = T, showWarnings = F)

## clmnt ----------------------------------------------------------------------

### script name
myName <- current_filename()
cat("[", myName, "] ",
    "Antani, speriamo duri poco. ",
    "\n", sep = "")

### load the input table
pathInTmb <- list.files(dirdtMeta, pattern = "metadata-tmb.txt", full.names = T)

dtMeta <- fread(pathInTmb)
dtMeta[subtype == "metastatic melanoma", subtype := "MM"]

### calculate mean of xgboost, lightgbm, and tabnet
### normalised with a fixed value for the size of the target
### as in the original paper (41 Mbp)
dtMeta[, mean_TMB_comparison := rowMeans(.SD, na.rm = TRUE),
       .SDcols = tail(names(dtMeta), 3)]

## compare matched pipeline vs models by subtype ------------------------------

cat("[", myName, "] ",
    "Comparing matched pipeline TMB vs models by subtype. ",
    "\n", sep = "")

subtypes <- unique(dtMeta$subtype)
for (subtype_i in subtypes) {
  dtSubtype <- dtMeta[subtype == subtype_i]
  subtype_name <- gsub(" ", "_", subtype_i)
  
  ### xgboost
  lmFit <- lm(matched_TMB ~ xgm_TMB, data = dtSubtype)
  rSlope <- coef(lmFit)[2]
  rSquared <- summary(lmFit)$r.squared
  
  plotTolo <- ggplot(dtSubtype, aes(x = xgm_TMB, y = matched_TMB)) +
    geom_point(size = 4, alpha = 0.7) +
    geom_smooth(method = "lm", se = T, color = "black") +
    theme_minimal() +
    theme(text = element_text(size = 36),
          axis.text = element_text(size = 36),
          plot.margin = unit(c(0.5, 1.5, 0.5, 1.5), "cm"),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_text(color = "black",
                                     margin = margin(t = 10)),
          axis.text.y = element_text(color = "black",
                                     margin = margin(r = 10))) +
    ### coord_fixed(ratio = 1, xlim = c(0, 40), ylim = c(0, 40), expand = F) +
    labs(subtitle = bquote(beta[1] * " = " * .(formatC(rSlope, format = "f",
                                                       digits = 3)) * 
                             ", R"^2 * " = " * .(formatC(rSquared, format = "f",
                                                         digits = 3))),
         title = subtype_i)
  pdf(file = file.path(dirPlots, paste0("fit-regression-xgb-", 
                                        subtype_name, ".pdf")),
      width = 8, height = 8)
  print(plotTolo)
  dev.off()
  
}
