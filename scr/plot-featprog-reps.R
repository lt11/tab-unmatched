#!/usr/bin/env Rscript

## header ---------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(this.path)
})

## functions ------------------------------------------------------------------

MakeMetricPlot <- function(metricName) {
  d <- dtPlot[Metric_label == metricName]
  
  if (nrow(d) == 0L) {
    warning("No data for metric: ", metricName)
    return(NULL)
  }
  
  ggplot(d,
         aes(x = N_features,
             y = M_mean,
             colour = Str_model,
             shape  = Data_set)) +
    geom_errorbar(aes(ymin = M_mean - M_sd,
                      ymax = M_mean + M_sd),
                  width = 0.35,
                  alpha = 0.35) +
    geom_point(size = 4, alpha = 0.7) +
    geom_line(aes(group = interaction(Str_model, Data_set)), alpha = 0.7) +
    scale_color_manual(values = modelColours, drop = F) +
    scale_shape_manual(values = dataShapes, drop = F) +
    labs(x = "Number of features",
         y = metricName,
         colour = "Str_model",
         shape = "Data_set") +
    scale_x_continuous(breaks = seq(65, 90, 5), limits = c(65, 90)) +
    theme_bw() +
    theme(plot.margin = unit(c(0.5, 1., 0.5, 1.), "cm"),
          axis.title = element_text(size = 36),
          axis.text = element_text(size = 36),
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(color = "black"),
          axis.text.y = element_text(color = "black"),
          legend.position = "none")
}

## clmnt ----------------------------------------------------------------------

### input
baseDir <- this.dir()
fileIn <- file.path(baseDir, "stats-feat-by-feat-reps.tsv")
dtIn <- fread(fileIn)

## data preparation -----------------------------------------------------------

### build dtMetadata from column names
### we expect columns like:
###   auc_validation_xgboost
###   accuracy_test_melanoma_lightgbm
###   precision_test_mixtcga_tabnet
### pattern: <metric>_<dataset>_<model>
colIds <- c("N_run", "N_features")
colMeasures <- setdiff(names(dtIn), colIds)

dtMeta <- rbindlist(
  lapply(colMeasures, function(col) {
    parts <- strsplit(col, "_")[[1]]
    metric <- parts[1]
    model <- tail(parts, 1)
    strData <- paste(parts[2:(length(parts) - 1)], collapse = "_")
    data.table(Col_name = col,
               Metric_id = metric,
               Str_model = model,
               Dataset_original = strData)
  }),
  use.names = T,
  fill = T
)

### map nicer metric labels
mapLabelMetric <- c(auc = "AUC",
                    accuracy = "Accuracy",
                    precision = "Precision",
                    recall = "Recall",
                    tnr = "TNR",
                    npv = "NPV")
dtMeta[, Metric_label := mapLabelMetric[Metric_id]]

### remove underscores
labelMap <- c(validation = "validation",
              test_mm = "test melanoma",
              test_melanoma = "test melanoma",
              test_tcga = "test mix tcga",
              test_mixtcga = "test mix tcga")
dtMeta[, Data_set := labelMap[Dataset_original]]

### melt to long format
dtLong <- melt(dtIn,
               id.vars = colIds,
               measure.vars = dtMeta$Col_name,
               variable.name = "Col_name",
               value.name = "M_value")

### join metadata
dtLong <- merge(dtLong, dtMeta, by = "Col_name", all.x = T)

### aggregate replicates
### N_run encodes the replicate as <number><letter>, e.g. 70a, ..., 70e.
dtLong[, N_run_base := as.integer(sub("^([0-9]+).*", "\\1", N_run))]

dtPlot <- dtLong[, .(M_mean = mean(M_value, na.rm = T),
                     M_sd = sd(M_value, na.rm = T),
                     N_reps = sum(!is.na(M_value))),
                 by = .(N_run_base, N_features, Col_name, Metric_id,
                        Str_model, Dataset_original, Metric_label, Data_set)]

## set colour and shape encodings
modelColours <- c("xgboost" = "red",
                  "lightgbm" = "darkgreen",
                  "tabnet" = "blue")

### data sets shapes
dataShapes <- c("validation" = 16,     # circle
                "test melanoma" = 17,  # triangle
                "test mix tcga" = 15)  # square

### make factors for consistent legend order
dtPlot[, Str_model := factor(Str_model, levels = names(modelColours))]
dtPlot[, Data_set := factor(Data_set,  levels = names(dataShapes))]
dtPlot[, Metric_label := factor(Metric_label,
                                levels = c("AUC", "Accuracy", "Precision",
                                           "Recall", "TNR", "NPV"))]

## plots ----------------------------------------------------------------------

for (metricName in levels(dtPlot$Metric_label)) {
  p <- MakeMetricPlot(metricName)
  if (!is.null(p)) {
    fileOut <- paste0("plot-", tolower(metricName), ".pdf")
    ggsave(file.path(baseDir, fileOut), p, width = 8, height = 8)
  }
}
