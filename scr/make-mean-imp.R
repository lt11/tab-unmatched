## header ---------------------------------------------------------------------

options(scipen = 999)
options(stringsAsFactors = F)
rm(list = ls())
library(data.table)
library(this.path)
library(scriptName)

## settings -------------------------------------------------------------------

### fixed settings
args <- commandArgs(trailingOnly = T)
if (length(args) >= 1L) {
  dirBase <- normalizePath(args[1], mustWork = T)
} else {
  dirBase <- dirname(this.dir())
}
dirData <- file.path(dirBase, "res")
dirRuns <- list.files(dirBase, pattern = "^run-1[a-z]$",
                      full.names = T)
### processing order of the models; any other model goes last
strModels <- c("xgboost", "lgbm", "tabnet")
dtImp <- data.table()

## clmnt ----------------------------------------------------------------------

### script name
myName <- current_filename()
cat("[", myName, "] ",
    "Antani, speriamo duri poco. ",
    "\n", sep = "")

for (indL in dirRuns) {
  strRun <- basename(indL)
  cat("[", myName, "] ", "Processing run: ",
      strRun, "\n", sep = "")

  ### read data
  pathInImp <- file.path(indL, "res", "importance.txt")
  if (!file.exists(pathInImp)) {
    cat("[", myName, "] ",
        pathInImp, " not found.", "\n", sep = "")
    next
  }

  dtOne <- fread(pathInImp)
  dtOne[, Run := strRun]

  ## append importance --------------------------------------------------------

  dtImp <- rbindlist(list(dtImp, dtOne), use.names = T, fill = T)
}

if (nrow(dtImp) == 0) {
  stop("[", myName, "] no importance.txt found under ",
       dirBase, "/run-1?/res")
}

## mean and sd ----------------------------------------------------------------

dtMean <- dtImp[, .(N_runs = .N,
                    Mean = mean(Importance),
                    SD = sd(Importance)),
                by = .(Model, Feature)]
dtMean[, Model := factor(Model,
                         levels = union(strModels, sort(unique(Model))))]
setorder(dtMean, Model, -Mean)
dtMean[, Model := as.character(Model)]
dtMean[, Rank_mean := seq_len(.N), by = Model]

for (indM in unique(dtMean$Model)) {
  cat("[", myName, "] ", "Model: ", indM, "\n", sep = "")
  dtM <- dtMean[Model == indM]
  for (indR in seq_len(nrow(dtM))) {
    cat("[", myName, "] ",
        dtM$Rank_mean[indR], "\t",
        dtM$Feature[indR], "\t",
        "mean: ", signif(dtM$Mean[indR], 6), "\t",
        "sd: ", signif(dtM$SD[indR], 6), "\t",
        "n: ", dtM$N_runs[indR], "\n", sep = "")
  }
}

dir.create(dirData, showWarnings = F)
fwrite(dtMean, file = file.path(dirData, "mean-imp.txt"),
       append = F, sep = "\t", col.names = T)
cat("[", myName, "] ", "Written: ",
    file.path(dirData, "mean-imp.txt"), "\n", sep = "")
