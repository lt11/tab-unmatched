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

### cleaning/inspection
### str(dtMeta)

### renaming factor levels
dtMeta[race == "BLACK OR AFRICAN AMERICAN", race := "BLACK"]
dtMeta$race <- factor(dtMeta$race, levels = c("WHITE", "BLACK", "ASIAN"))
dtMeta <- dtMeta[!is.na(race)]

## free-axis plots -------------------------------------------------------------

### analysis persone normali
pathPdf <- file.path(dirPlots, "plot-tmb-matched.pdf")
plotPdf <- ggplot(dtMeta, aes(x = race, y = matched_TMB)) + 
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, color = "black"),
        legend.position = "none",
        text = element_text(size = 36),
        axis.text.y = element_text(color = "black")) +
  geom_boxplot(outlier.colour = "black", outlier.shape = 1,
               outlier.size = 4) +
  labs(x = "Ancestry", y = "TMB [N/Mbp]",
       title = "Matched analysis")
pdf(file = pathPdf, width = 8, height = 12)
print(plotPdf)
dev.off()

### xgboost
pathPdf <- file.path(dirPlots, "plot-tmb-xgboost.pdf")
plotPdf <- ggplot(dtMeta, aes(x = race, y = xgm_TMB)) + 
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, color = "black"),
        legend.position = "none",
        text = element_text(size = 36),
        axis.text.y = element_text(color = "black")) +
  geom_boxplot(outlier.colour = "black", outlier.shape = 1,
               outlier.size = 4) +
  labs(x = "Ancestry", y = "TMB [N/Mbp]",
       title = "XGBoost")
pdf(file = pathPdf, width = 8, height = 12)
print(plotPdf)
dev.off()

### lightgbm
pathPdf <- file.path(dirPlots, "plot-tmb-lightgbm.pdf")
plotPdf <- ggplot(dtMeta, aes(x = race, y = lgbm_TMB)) + 
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, color = "black"),
        legend.position = "none",
        text = element_text(size = 36),
        axis.text.y = element_text(color = "black")) +
  geom_boxplot(outlier.colour = "black", outlier.shape = 1,
               outlier.size = 4) +
  labs(x = "Ancestry", y = "TMB [N/Mbp]",
       title = "LightGBM")
pdf(file = pathPdf, width = 8, height = 12)
print(plotPdf)
dev.off()

### tabnet
pathPdf <- file.path(dirPlots, "plot-tmb-tabnet.pdf")
plotPdf <- ggplot(dtMeta, aes(x = race, y = tabnet_TMB)) + 
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, color = "black"),
        legend.position = "none",
        text = element_text(size = 36),
        axis.text.y = element_text(color = "black")) +
  geom_boxplot(outlier.colour = "black", outlier.shape = 1,
               outlier.size = 4) +
  labs(x = "Ancestry", y = "TMB [N/Mbp]",
       title = "TabNet")
pdf(file = pathPdf, width = 8, height = 12)
print(plotPdf)
dev.off()

## fixed-axis plots ------------------------------------------------------------

### analysis persone normali
pathPdf <- file.path(dirPlots, "plot-fixed-tmb-matched.pdf")
plotPdf <- ggplot(dtMeta, aes(x = race, y = matched_TMB)) + 
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, color = "black"),
        legend.position = "none",
        text = element_text(size = 36),
        axis.text.y = element_text(color = "black")) +
  geom_boxplot(outlier.colour = "black", outlier.shape = 1,
               outlier.size = 4) +
  labs(x = "Ancestry", y = "TMB [N/Mbp]",
       title = "Matched analysis") +
  coord_cartesian(ylim = c(0, 50))
pdf(file = pathPdf, width = 8, height = 12)
print(plotPdf)
dev.off()

### xgboost
pathPdf <- file.path(dirPlots, "plot-fixed-tmb-xgboost.pdf")
plotPdf <- ggplot(dtMeta, aes(x = race, y = xgm_TMB)) + 
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, color = "black"),
        legend.position = "none",
        text = element_text(size = 36),
        axis.text.y = element_text(color = "black")) +
  geom_boxplot(outlier.colour = "black", outlier.shape = 1,
               outlier.size = 4) +
  labs(x = "Ancestry", y = "TMB [N/Mbp]",
       title = "XGBoost") +
  coord_cartesian(ylim = c(0, 50))
pdf(file = pathPdf, width = 8, height = 12)
print(plotPdf)
dev.off()

### lightgbm
pathPdf <- file.path(dirPlots, "plot-fixed-tmb-lightgbm.pdf")
plotPdf <- ggplot(dtMeta, aes(x = race, y = lgbm_TMB)) + 
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, color = "black"),
        legend.position = "none",
        text = element_text(size = 36),
        axis.text.y = element_text(color = "black")) +
  geom_boxplot(outlier.colour = "black", outlier.shape = 1,
               outlier.size = 4) +
  labs(x = "Ancestry", y = "TMB [N/Mbp]",
       title = "LightGBM") +
  coord_cartesian(ylim = c(0, 50))
pdf(file = pathPdf, width = 8, height = 12)
print(plotPdf)
dev.off()

### tabnet
pathPdf <- file.path(dirPlots, "plot-fixed-tmb-tabnet.pdf")
plotPdf <- ggplot(dtMeta, aes(x = race, y = tabnet_TMB)) + 
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, color = "black"),
        legend.position = "none",
        text = element_text(size = 36),
        axis.text.y = element_text(color = "black")) +
  geom_boxplot(outlier.colour = "black", outlier.shape = 1,
               outlier.size = 4) +
  labs(x = "Ancestry", y = "TMB [N/Mbp]",
       title = "TabNet") +
  coord_cartesian(ylim = c(0, 50))
pdf(file = pathPdf, width = 8, height = 12)
print(plotPdf)
dev.off()
