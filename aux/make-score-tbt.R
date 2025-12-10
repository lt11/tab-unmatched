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
  allTumours <- unique(dtPred[, subtype])
  for (indT in allTumours) {
    oneRow <- c(indL)
    
    nGermline <- length(which(dtPred$target == 0 
                              & dtPred$subtype == indT))
    nSomatic <- length(which(dtPred$target == 1
                             & dtPred$subtype == indT))
    oneRow <- cbind(oneRow, indT, nGermline, nSomatic)
    
    ### target = 1 -> somatic, target = 0 -> germline
    nG <- length(which(dtPred$xgm_preds == "germline"
                       & dtPred$subtype == indT))
    nS <- length(which(dtPred$xgm_preds == "somatic"
                       & dtPred$subtype == indT))
    nTP <- length(which(dtPred$target == 1
                        & dtPred$xgm_preds == "somatic"
                        & dtPred$subtype == indT))
    nFN <- length(which(dtPred$target == 1
                        & dtPred$xgm_preds == "germline"
                        & dtPred$subtype == indT))
    nTN <- length(which(dtPred$target == 0
                        & dtPred$xgm_preds == "germline"
                        & dtPred$subtype == indT))
    nFP <- length(which(dtPred$target == 0
                        & dtPred$xgm_preds == "somatic"
                        & dtPred$subtype == indT))
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
    
    nG <- length(which(dtPred$lgbm_preds == "germline"
                       & dtPred$subtype == indT))
    nS <- length(which(dtPred$lgbm_preds == "somatic"
                       & dtPred$subtype == indT))
    nTP <- length(which(dtPred$target == 1
                        & dtPred$lgbm_preds == "somatic"
                        & dtPred$subtype == indT))
    nFN <- length(which(dtPred$target == 1
                        & dtPred$lgbm_preds == "germline"
                        & dtPred$subtype == indT))
    nTN <- length(which(dtPred$target == 0
                        & dtPred$lgbm_preds == "germline"
                        & dtPred$subtype == indT))
    nFP <- length(which(dtPred$target == 0
                        & dtPred$lgbm_preds == "somatic"
                        & dtPred$subtype == indT))
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
    
    nG <- length(which(dtPred$tabnet_preds == "germline"
                       & dtPred$subtype == indT))
    nS <- length(which(dtPred$tabnet_preds == "somatic"
                       & dtPred$subtype == indT))
    nTP <- length(which(dtPred$target == 1
                        & dtPred$tabnet_preds == "somatic"
                        & dtPred$subtype == indT))
    nFN <- length(which(dtPred$target == 1
                        & dtPred$tabnet_preds == "germline"
                        & dtPred$subtype == indT))
    nTN <- length(which(dtPred$target == 0
                        & dtPred$tabnet_preds == "germline"
                        & dtPred$subtype == indT))
    nFP <- length(which(dtPred$target == 0
                        & dtPred$tabnet_preds == "somatic"
                        & dtPred$subtype == indT))
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
    
    ## append metrics -----------------------------------------------------------
    
    dtScores <- rbindlist(list(dtScores, as.data.table(as.list(oneRow))),
                          use.names = F, fill = T)
    
  }
}
colnames(dtScores) <- c(
  "Set",
  "Data",
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
  "NPV (TabNet)")
dtScores[Data == "metastatic melanoma", Data := "MM"]

fwrite(dtScores, file = file.path(dirData, "scores-tbt.txt"),
       append = F, sep = "\t", col.names = T)
