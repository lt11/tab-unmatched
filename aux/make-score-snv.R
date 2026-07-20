## header ---------------------------------------------------------------------

options(scipen = 999)
options(stringsAsFactors = F)
rm(list = ls())
library(data.table)
library(this.path)
library(scriptName)

## settings -------------------------------------------------------------------

### fixed settings
dirBase <- dirname(this.dir())
dirData <- file.path(dirBase, "res")
strSets <- c("train", "validation", "test-melanoma", "test-mixtcga")
dtScores <- data.table()

## clmnt ----------------------------------------------------------------------

### script name
myName <- current_filename()
cat("[", myName, "] ",
    "Antani, speriamo duri poco. ",
    "\n", sep = "")

for (indL in strSets) {
  cat("[", myName, "] ", "Processing set: ",
      indL, "\n", sep = "")
  
  ### read data
  pathInPred <- list.files(dirData, pattern = paste0("preds-", indL),
                           full.names = T)
  if (length(pathInPred) == 0) {
    cat("[", myName, "] ",
        indL, " not found.", "\n", sep = "")
    next
  }
  
  dtPred <- fread(pathInPred)
  
  ## keep only tmb snvs -------------------------------------------------------
  
  dtPred <- dtPred[ontology_missense == 1 
                   | ontology_nonsense == 1]
  
  oneRow <- c(indL)
  nGermline <- length(which(dtPred$target == 0))
  nSomatic <- length(which(dtPred$target == 1))
  oneRow <- cbind(oneRow, nGermline, nSomatic)
  
  ### target = 1 -> somatic, target = 0 -> germline
  nG <- length(which(dtPred$xgm_preds == "germline"))
  nS <- length(which(dtPred$xgm_preds == "somatic"))
  nTP <- length(which(dtPred$target == 1
                      & dtPred$xgm_preds == "somatic"))
  nFN <- length(which(dtPred$target == 1
                      & dtPred$xgm_preds == "germline"))
  nTN <- length(which(dtPred$target == 0
                      & dtPred$xgm_preds == "germline"))
  nFP <- length(which(dtPred$target == 0
                      & dtPred$xgm_preds == "somatic"))
  valPrec <- nTP / (nTP + nFP)
  valReca <- nTP / (nTP + nFN)
  valAcc <- (nTP + nTN) / (nTP + nTN + nFP + nFN)
  valTNR <- nTN / (nTN + nFP)
  valNPV <- nTN / (nTN + nFN)
  oneRow <- c(oneRow, nG, nS, nTP, nTN, nFP, nFN,
              valPrec, valReca, valAcc, valTNR, valNPV)
  cat("[", myName, "] ",
      "XGBoost precision: ", valPrec, "\n", sep = "")
  cat("[", myName, "] ",
      "XGBoost recall: ", valReca, "\n", sep = "")
  cat("[", myName, "] ",
      "XGBoost accuracy: ", valAcc, "\n", sep = "")
  cat("[", myName, "] ",
      "XGBoost TNR: ", valTNR, "\n", sep = "")
  cat("[", myName, "] ",
      "XGBoost NPV: ", valNPV, "\n", sep = "")
  
  nG <- length(which(dtPred$lgbm_preds == "germline"))
  nS <- length(which(dtPred$lgbm_preds == "somatic"))
  nTP <- length(which(dtPred$target == 1
                      & dtPred$lgbm_preds == "somatic"))
  nFN <- length(which(dtPred$target == 1
                      & dtPred$lgbm_preds == "germline"))
  nTN <- length(which(dtPred$target == 0
                      & dtPred$lgbm_preds == "germline"))
  nFP <- length(which(dtPred$target == 0
                      & dtPred$lgbm_preds == "somatic"))
  valPrec <- nTP / (nTP + nFP)
  valReca <- nTP / (nTP + nFN)
  valAcc <- (nTP + nTN) / (nTP + nTN + nFP + nFN)
  valTNR <- nTN / (nTN + nFP)
  valNPV <- nTN / (nTN + nFN)
  oneRow <- c(oneRow, nG, nS, nTP, nTN, nFP, nFN,
              valPrec, valReca, valAcc, valTNR, valNPV)
  cat("[", myName, "] ",
      "LightGBM precision: ", valPrec, "\n", sep = "")
  cat("[", myName, "] ",
      "LightGBM recall: ", valReca, "\n", sep = "")
  cat("[", myName, "] ",
      "LightGBM accuracy: ", valAcc, "\n", sep = "")
  cat("[", myName, "] ",
      "LightGBM TNR: ", valTNR, "\n", sep = "")
  cat("[", myName, "] ",
      "LightGBM NPV: ", valNPV, "\n", sep = "")
  
  nG <- length(which(dtPred$tabnet_preds == "germline"))
  nS <- length(which(dtPred$tabnet_preds == "somatic"))
  nTP <- length(which(dtPred$target == 1
                      & dtPred$tabnet_preds == "somatic"))
  nFN <- length(which(dtPred$target == 1
                      & dtPred$tabnet_preds == "germline"))
  nTN <- length(which(dtPred$target == 0
                      & dtPred$tabnet_preds == "germline"))
  nFP <- length(which(dtPred$target == 0
                      & dtPred$tabnet_preds == "somatic"))
  valPrec <- nTP / (nTP + nFP)
  valReca <- nTP / (nTP + nFN)
  valAcc <- (nTP + nTN) / (nTP + nTN + nFP + nFN)
  valTNR <- nTN / (nTN + nFP)
  valNPV <- nTN / (nTN + nFN)
  oneRow <- c(oneRow, nG, nS, nTP, nTN, nFP, nFN,
              valPrec, valReca, valAcc, valTNR, valNPV)
  cat("[", myName, "] ",
      "TabNet precision: ", valPrec, "\n", sep = "")
  cat("[", myName, "] ",
      "TabNet recall: ", valReca, "\n", sep = "")
  cat("[", myName, "] ",
      "TabNet accuracy: ", valAcc, "\n", sep = "")
  cat("[", myName, "] ",
      "TabNet TNR: ", valTNR, "\n", sep = "")
  cat("[", myName, "] ",
      "TabNet NPV: ", valNPV, "\n", sep = "")

  if (!"logreg_preds" %in% colnames(dtPred)
      || all(is.na(dtPred$logreg_preds))) {
    nG <- NA
    nS <- NA
    nTP <- NA
    nFN <- NA
    nTN <- NA
    nFP <- NA
    valPrec <- NA
    valReca <- NA
    valAcc <- NA
    valTNR <- NA
    valNPV <- NA
  } else {
    nG <- length(which(dtPred$logreg_preds == "germline"))
    nS <- length(which(dtPred$logreg_preds == "somatic"))
    nTP <- length(which(dtPred$target == 1
                        & dtPred$logreg_preds == "somatic"))
    nFN <- length(which(dtPred$target == 1
                        & dtPred$logreg_preds == "germline"))
    nTN <- length(which(dtPred$target == 0
                        & dtPred$logreg_preds == "germline"))
    nFP <- length(which(dtPred$target == 0
                        & dtPred$logreg_preds == "somatic"))
    valPrec <- nTP / (nTP + nFP)
    valReca <- nTP / (nTP + nFN)
    valAcc <- (nTP + nTN) / (nTP + nTN + nFP + nFN)
    valTNR <- nTN / (nTN + nFP)
    valNPV <- nTN / (nTN + nFN)
  }
  oneRow <- c(oneRow, nG, nS, nTP, nTN, nFP, nFN,
              valPrec, valReca, valAcc, valTNR, valNPV)
  cat("[", myName, "] ",
      "LogiRegr precision: ", valPrec, "\n", sep = "")
  cat("[", myName, "] ",
      "LogiRegr recall: ", valReca, "\n", sep = "")
  cat("[", myName, "] ",
      "LogiRegr accuracy: ", valAcc, "\n", sep = "")
  cat("[", myName, "] ",
      "LogiRegr TNR: ", valTNR, "\n", sep = "")
  cat("[", myName, "] ",
      "LogiRegr NPV: ", valNPV, "\n", sep = "")

  if (!"hf_preds" %in% colnames(dtPred)
      || all(is.na(dtPred$hf_preds))) {
    nG <- NA
    nS <- NA
    nTP <- NA
    nFN <- NA
    nTN <- NA
    nFP <- NA
    valPrec <- NA
    valReca <- NA
    valAcc <- NA
    valTNR <- NA
    valNPV <- NA
  } else {
    nG <- length(which(dtPred$hf_preds == "germline"))
    nS <- length(which(dtPred$hf_preds == "somatic"))
    nTP <- length(which(dtPred$target == 1
                        & dtPred$hf_preds == "somatic"))
    nFN <- length(which(dtPred$target == 1
                        & dtPred$hf_preds == "germline"))
    nTN <- length(which(dtPred$target == 0
                        & dtPred$hf_preds == "germline"))
    nFP <- length(which(dtPred$target == 0
                        & dtPred$hf_preds == "somatic"))
    valPrec <- nTP / (nTP + nFP)
    valReca <- nTP / (nTP + nFN)
    valAcc <- (nTP + nTN) / (nTP + nTN + nFP + nFN)
    valTNR <- nTN / (nTN + nFP)
    valNPV <- nTN / (nTN + nFN)
  }
  oneRow <- c(oneRow, nG, nS, nTP, nTN, nFP, nFN,
              valPrec, valReca, valAcc, valTNR, valNPV)
  cat("[", myName, "] ",
      "HardFilt precision: ", valPrec, "\n", sep = "")
  cat("[", myName, "] ",
      "HardFilt recall: ", valReca, "\n", sep = "")
  cat("[", myName, "] ",
      "HardFilt accuracy: ", valAcc, "\n", sep = "")
  cat("[", myName, "] ",
      "HardFilt TNR: ", valTNR, "\n", sep = "")
  cat("[", myName, "] ",
      "HardFilt NPV: ", valNPV, "\n", sep = "")
  
  ## append metrics -----------------------------------------------------------
  
  dtScores <- rbindlist(list(dtScores, as.data.table(as.list(oneRow))),
                        use.names = F, fill = T)
}

colnames(dtScores) <- c(
  "Set",
  "Germline (gold standard)",
  "Somatic (gold standard)",
  "Germline (XGBoost)",
  "Somatic (XGBoost)",
  "TP (XGBoost)",
  "TN (XGBoost)",
  "FP (XGBoost)",
  "FN (XGBoost)",
  "Precision (XGBoost)",
  "Recall (XGBoost)",
  "Accuracy (XGBoost)",
  "TNR (XGBoost)",
  "NPV (XGBoost)",
  "Germline (LightGBM)",
  "Somatic (LightGBM)",
  "TP (LightGBM)",
  "TN (LightGBM)",
  "FP (LightGBM)",
  "FN (LightGBM)",
  "Precision (LightGBM)",
  "Recall (LightGBM)",
  "Accuracy (LightGBM)",
  "TNR (LightGBM)",
  "NPV (LightGBM)",
  "Germline (TabNet)",
  "Somatic (TabNet)",
  "TP (TabNet)",
  "TN (TabNet)",
  "FP (TabNet)",
  "FN (TabNet)",
  "Precision (TabNet)",
  "Recall (TabNet)",
  "Accuracy (TabNet)",
  "TNR (TabNet)",
  "NPV (TabNet)",
  "Germline (LogiRegr)",
  "Somatic (LogiRegr)",
  "TP (LogiRegr)",
  "TN (LogiRegr)",
  "FP (LogiRegr)",
  "FN (LogiRegr)",
  "Precision (LogiRegr)",
  "Recall (LogiRegr)",
  "Accuracy (LogiRegr)",
  "TNR (LogiRegr)",
  "NPV (LogiRegr)",
  "Germline (HardFilt)",
  "Somatic (HardFilt)",
  "TP (HardFilt)",
  "TN (HardFilt)",
  "FP (HardFilt)",
  "FN (HardFilt)",
  "Precision (HardFilt)",
  "Recall (HardFilt)",
  "Accuracy (HardFilt)",
  "TNR (HardFilt)",
  "NPV (HardFilt)")

fwrite(dtScores, file = file.path(dirData, "scores-snv.txt"),
       append = F, sep = "\t", col.names = T)
