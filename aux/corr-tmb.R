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
dirMeta <- file.path(dirBase, "res")

## clmnt ----------------------------------------------------------------------

### script name
myName <- current_filename()
cat("[", myName, "] ",
    "Antani, speriamo duri poco. ",
    "\n", sep = "")

### load the input table
pathInTmb <- list.files(dirMeta, pattern = "metadata-tmb.txt", full.names = T)
dtMeta <- fread(pathInTmb)

### clean
dtMeta <- dtMeta[matched_TMB != 0]

### renaming factor levels
dtMeta[race == "BLACK OR AFRICAN AMERICAN", race := "BLACK"]
dtMeta$race <- factor(dtMeta$race, levels = c("WHITE", "BLACK", "ASIAN"))

### normality tests
# shapiro.test(dtMeta[["matched_TMB"]])
# shapiro.test(dtMeta[["xgm_TMB"]])
# shapiro.test(dtMeta[["lgbm_TMB"]])
# shapiro.test(dtMeta[["tabnet_TMB"]])

## melanoma data ---------------------------------------------------------------

rValX <- cor(dtMeta[subtype == "metastatic melanoma", matched_TMB],
             dtMeta[subtype == "metastatic melanoma", xgm_TMB],
             method = "kendall")
rValL <- cor(dtMeta[subtype == "metastatic melanoma", matched_TMB],
             dtMeta[subtype == "metastatic melanoma", lgbm_TMB],
             method = "kendall")
rValT <- cor(dtMeta[subtype == "metastatic melanoma", matched_TMB],
             dtMeta[subtype == "metastatic melanoma", tabnet_TMB],
             method = "kendall")

## the tcga data ---------------------------------------------------------------

rValX <- cor(dtMeta[subtype == "BRCA" | subtype == "SARC" |  subtype == "UCEC",
                    matched_TMB],
             dtMeta[subtype == "BRCA" | subtype == "SARC" |  subtype == "UCEC",
                    xgm_TMB],
             method = "kendall")
rValL <- cor(dtMeta[subtype == "BRCA" | subtype == "SARC" |  subtype == "UCEC",
                    matched_TMB],
             dtMeta[subtype == "BRCA" | subtype == "SARC" |  subtype == "UCEC",
                    lgbm_TMB],
             method = "kendall")
rValT <- cor(dtMeta[subtype == "BRCA" | subtype == "SARC" |  subtype == "UCEC",
                    matched_TMB],
             dtMeta[subtype == "BRCA" | subtype == "SARC" |  subtype == "UCEC",
                    tabnet_TMB],
             method = "kendall")

## test matched tmb vs model tmb distributions ---------------------------------

ks.test(dtMeta[["matched_TMB"]],
        dtMeta[["xgm_TMB"]])

ks.test(dtMeta[subtype == "metastatic melanoma", matched_TMB],
        dtMeta[subtype == "metastatic melanoma", xgm_TMB])
ks.test(dtMeta[subtype == "metastatic melanoma", matched_TMB],
        dtMeta[subtype == "metastatic melanoma", lgbm_TMB])
ks.test(dtMeta[subtype == "metastatic melanoma", matched_TMB],
        dtMeta[subtype == "metastatic melanoma", tabnet_TMB])

ks.test(dtMeta[subtype == "BRCA" | subtype == "SARC" |  subtype == "UCEC",
               matched_TMB],
        dtMeta[subtype == "BRCA" | subtype == "SARC" |  subtype == "UCEC",
               xgm_TMB])
ks.test(dtMeta[subtype == "BRCA" | subtype == "SARC" |  subtype == "UCEC",
               matched_TMB],
        dtMeta[subtype == "BRCA" | subtype == "SARC" |  subtype == "UCEC",
               lgbm_TMB])
ks.test(dtMeta[subtype == "BRCA" | subtype == "SARC" |  subtype == "UCEC",
               matched_TMB],
        dtMeta[subtype == "BRCA" | subtype == "SARC" |  subtype == "UCEC",
               tabnet_TMB])

## test matched tmb vs model tmb by ethnic group -------------------------------

wilcox.test(dtMeta[race == "BLACK", matched_TMB],
            dtMeta[race == "BLACK", xgm_TMB])
wilcox.test(dtMeta[race == "WHITE", matched_TMB],
            dtMeta[race == "WHITE", xgm_TMB])
wilcox.test(dtMeta[race == "ASIAN", matched_TMB],
            dtMeta[race == "ASIAN", xgm_TMB])

wilcox.test(dtMeta[race == "BLACK", matched_TMB],
            dtMeta[race == "BLACK", lgbm_TMB])
wilcox.test(dtMeta[race == "WHITE", matched_TMB],
            dtMeta[race == "WHITE", lgbm_TMB])
wilcox.test(dtMeta[race == "ASIAN", matched_TMB],
            dtMeta[race == "ASIAN", lgbm_TMB])

wilcox.test(dtMeta[race == "BLACK", matched_TMB],
            dtMeta[race == "BLACK", tabnet_TMB])
wilcox.test(dtMeta[race == "WHITE", matched_TMB],
            dtMeta[race == "WHITE", tabnet_TMB])
wilcox.test(dtMeta[race == "ASIAN", matched_TMB],
            dtMeta[race == "ASIAN", tabnet_TMB])

### double check the last one
indA <- which(dtMeta$race == "ASIAN")
wilcox.test(dtMeta$matched_TMB[indA],
            dtMeta$tabnet_TMB[indA])

### test across different ethnicities: conceptually wrong
### since it does not take into account genetic background effects
