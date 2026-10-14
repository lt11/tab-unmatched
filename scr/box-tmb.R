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

### test splits
testSplits <- c("test-mixtcga", "test-ffpe-dlbcl", "test-metmel", "test-mixanc")

### models: TMB column, file tag and plot title
dtModels <- data.table(
  column = c("matched_TMB", "xgm_TMB", "lgbm_TMB", "tabnet_TMB",
             "logreg_TMB", "hf_TMB"),
  tag = c("matched", "xgboost", "lightgbm", "tabnet",
          "logreg", "hardfilter"),
  title = c("Matched analysis", "XGBoost", "LightGBM", "TabNet",
            "Logistic regression", "Hard filters"))

## functions ------------------------------------------------------------------

### box plot of TMB by model
plotTmb <- function(dt, title) {
  ggplot(dt, aes(x = model, y = TMB)) +
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, color = "black"),
          legend.position = "none",
          text = element_text(size = 36),
          axis.text.y = element_text(color = "black"),
          plot.margin = margin(10, 10, 10, 80)) +
    geom_boxplot(outlier.colour = "black", outlier.shape = 1,
                 outlier.size = 4) +
    labs(x = "Model", y = "TMB [N/Mbp]",
         title = title)
}

## clmnt ----------------------------------------------------------------------

### script name
myName <- current_filename()
cat("[", myName, "] ",
    "Antani, speriamo duri poco. ",
    "\n", sep = "")

### load the input table
pathInTmb <- list.files(dirdtMeta, pattern = "metadata-tmb.txt", full.names = T)
dtMeta <- fread(pathInTmb)

### keep test splits only
dtMeta <- dtMeta[split %in% testSplits]

### long format: one row per sample and model
dtLong <- melt(dtMeta, id.vars = c("fastq_id", "split"),
               measure.vars = dtModels$column,
               variable.name = "model", value.name = "TMB")
dtLong$model <- factor(dtLong$model, levels = dtModels$column,
                       labels = dtModels$title)

## free-axis plots (log scale) -------------------------------------------------

### log(1 + x) keeps samples with TMB = 0
for (testSplit in testSplits) {
  pathPdf <- file.path(dirPlots, paste0("plot-tmb-", testSplit, ".pdf"))
  plotPdf <- plotTmb(dtLong[split == testSplit], testSplit) +
    scale_y_continuous(trans = "log1p",
                       breaks = c(0, 1, 3, 10, 30, 100, 300))
  pdf(file = pathPdf, width = 10, height = 12)
  print(plotPdf)
  dev.off()
}

## fixed-axis plots ------------------------------------------------------------

for (testSplit in testSplits) {
  pathPdf <- file.path(dirPlots, paste0("plot-fixed-tmb-", testSplit, ".pdf"))
  plotPdf <- plotTmb(dtLong[split == testSplit], testSplit) +
    coord_cartesian(ylim = c(0, 50))
  pdf(file = pathPdf, width = 10, height = 12)
  print(plotPdf)
  dev.off()
}
