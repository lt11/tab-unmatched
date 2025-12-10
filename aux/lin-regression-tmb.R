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
dirPlots <- file.path(dirBase, "plots")
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

### linear model
lmFit <- lm(mean_TMB_comparison ~ TMB_internally_calculated, data = dtMeta)
rSlope <- coef(lmFit)[2]
rSquared <- summary(lmFit)$r.squared

### plot linear regression fit
plotTolo <- ggplot(dtMeta, aes(x = TMB_internally_calculated,
                               y = mean_TMB_comparison,
                               color = subtype)) +
  geom_point(size = 4, alpha = 0.7) +
  geom_smooth(method = "lm", se = T, color = "black") +
  theme_minimal() +
  theme(text = element_text(size = 36),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black")) +
  xlim(0, 350) +
  ylim(0, 350) +
  labs(x = "Reported Mean TMB [1/Mbp]",
       y = "Mean TMB [1/Mbp]",
       subtitle = bquote(beta[1] * " = " * .(round(rSlope, 3)) * 
                           ", R"^2 * " = " * .(round(rSquared, 3)))
  )

pdf(file = file.path(dirPlots, "mean-vs-reported.pdf"), width = 12, height = 8)
print(plotTolo)
dev.off()

### calculate IQR bounds (nostri)
Q1 <- quantile(dtMeta$mean_TMB_comparison, 0.25, na.rm = T)
Q3 <- quantile(dtMeta$mean_TMB_comparison, 0.75, na.rm = T)
valIQR <- Q3 - Q1
valLow <- Q1 - 1.5 * valIQR
valUp <- Q3 + 1.5 * valIQR
nNostri <- length(dtMeta[mean_TMB_comparison >= valLow
                         & mean_TMB_comparison <= valUp])

### calculate IQR bounds (loro)
Q1 <- quantile(dtMeta$TMB_internally_calculated, 0.25, na.rm = T)
Q3 <- quantile(dtMeta$TMB_internally_calculated, 0.75, na.rm = T)
valIQR <- Q3 - Q1
valLow <- Q1 - 1.5 * valIQR
valUp <- Q3 + 1.5 * valIQR
nLoro <- length(dtMeta[TMB_internally_calculated >= valLow
                       & TMB_internally_calculated <= valUp])

### subset the data to remove outliers (nostri)
dtMetaNoOutl <- dtMeta[mean_TMB_comparison >= valLow
                       & mean_TMB_comparison <= valUp]

### linear model
lmFit <- lm(mean_TMB_comparison ~ TMB_internally_calculated,
             data = dtMetaNoOutl)
rSlope <- coef(lmFit)[2]
rSquared <- summary(lmFit)$r.squared

plotTolo <- ggplot(dtMetaNoOutl, aes(x = TMB_internally_calculated,
                                     y = mean_TMB_comparison,
                                     color = subtype)) +
  geom_point(size = 4, alpha = 0.7) +
  geom_smooth(method = "lm", se = T, color = "black") +
  theme_minimal() +
  labs(x = "Reported Mean TMB [1/Mbp]",
       y = "Mean TMB [1/Mbp]",
       subtitle = bquote(beta[1] * " = " * .(round(rSlope, 3)) * 
                           ", R"^2 * " = " * .(round(rSquared, 3)))) +
  theme(text = element_text(size = 36),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black")) +
  xlim(0, 35) +
  ylim(0, 35)
  
pdf(file = file.path(dirPlots, "mean-vs-reported-no-outliers.pdf"),
    width = 12, height = 8)
print(plotTolo)
dev.off()

sd(dtMeta$mean_TMB_comparison)
sd(dtMeta$TMB_internally_calculated)
