#!/usr/bin/env Rscript

## header ---------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

## functions ------------------------------------------------------------------

### scales, labels and theme shared by the mean/SD and the raw plots
MetricPlotLayers <- function(metricName) {
  list(
    scale_color_manual(values = modelColours, drop = F),
    scale_shape_manual(values = dataShapes, drop = F),
    labs(x = "Number of features",
         y = metricName,
         colour = "Str_model",
         shape = "Data_set"),
    scale_x_continuous(breaks = xBreaks, limits = xLimits),
    theme_bw(),
    theme(plot.margin = unit(c(0.5, 1., 0.5, 1.), "cm"),
          axis.title = element_text(size = 36),
          axis.text = element_text(size = 36),
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(color = "black"),
          axis.text.y = element_text(color = "black"),
          legend.position = "none")
  )
}

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
    MetricPlotLayers(metricName)
}

### one point per replicate; a small horizontal jitter (fixed seed) keeps
### replicates at the same feature count from overlapping
MakeRawPlot <- function(metricName) {
  d <- dtLong[Metric_label == metricName & !is.na(M_value)]

  ggplot(d,
         aes(x = N_features,
             y = M_value,
             colour = Str_model,
             shape  = Data_set)) +
    geom_point(size = 3, alpha = 0.5,
               position = position_jitter(width = 0.2, height = 0, seed = 1)) +
    MetricPlotLayers(metricName)
}

### save p twice: free y axis, and y axis forced to [0, 1]; coord_cartesian
### zooms rather than dropping error bars that extend past the limits
SaveMetricPlots <- function(p, prefix, metricName) {
  fileOut <- paste0(prefix, tolower(metricName), "ylimfree.pdf")
  ggsave(file.path(dirOut, fileOut), p, width = 8, height = 8)
  fileOut <- paste0(prefix, tolower(metricName), "ylimfixed.pdf")
  pFixed <- p +
    scale_y_continuous(breaks = seq(0, 1, 0.2)) +
    coord_cartesian(ylim = c(0, 1), expand = T)
  ggsave(file.path(dirOut, fileOut), pFixed, width = 8, height = 8)
}

## clmnt ----------------------------------------------------------------------

### input and output paths
args <- commandArgs(trailingOnly = T)
if (length(args) != 2L) {
  stop("Usage: plot-featprog-reps.R <input.tsv> <output_dir>")
}
fileIn <- args[1]
if (!file.exists(fileIn)) {
  stop("Input file not found: ", fileIn)
}
dirOut <- args[2]
dir.create(dirOut, recursive = T, showWarnings = F)
dtIn <- fread(fileIn)

## data preparation -----------------------------------------------------------

### build dtMetadata from column names
### we expect columns like:
###   auc_test_melanoma_xgboost
###   accuracy_test_mixtcga_lightgbm
###   precision_test_ffpe_dlbcl_tabnet
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
labelMap <- c(test_mm = "test melanoma",
              test_melanoma = "test melanoma",
              test_tcga = "test mix tcga",
              test_mixtcga = "test mix tcga",
              test_mixanc = "test mix ancestry",
              test_ffpe_dlbcl = "test ffpe dlbcl")
dtMeta[, Data_set := labelMap[Dataset_original]]

### melt to long format
dtLong <- melt(dtIn,
               id.vars = colIds,
               measure.vars = dtMeta$Col_name,
               variable.name = "Col_name",
               value.name = "M_value")

### join metadata
dtLong <- merge(dtLong, dtMeta, by = "Col_name", all.x = T)

## set colour and shape encodings
modelColours <- c("xgboost" = "red",
                  "lightgbm" = "darkgreen",
                  "tabnet" = "blue")

### data sets shapes
dataShapes <- c("test melanoma" = 17,      # triangle
                "test mix tcga" = 15,      # square
                "test mix ancestry" = 18,  # diamond
                "test ffpe dlbcl" = 16)    # circle

### make factors for consistent legend order; done before aggregating so
### that both the raw and the aggregated tables carry them
dtLong[, Str_model := factor(Str_model, levels = names(modelColours))]
dtLong[, Data_set := factor(Data_set,  levels = names(dataShapes))]
dtLong[, Metric_label := factor(Metric_label,
                                levels = c("AUC", "Accuracy", "Precision",
                                           "Recall", "TNR", "NPV"))]

### aggregate replicates
### N_run encodes the replicate as <number><letter>, e.g. 70a, ..., 70e.
dtLong[, N_run_base := as.integer(sub("^([0-9]+).*", "\\1", N_run))]

dtPlot <- dtLong[, .(M_mean = mean(M_value, na.rm = T),
                     M_sd = sd(M_value, na.rm = T),
                     N_reps = sum(!is.na(M_value))),
                 by = .(N_run_base, N_features, Col_name, Metric_id,
                        Str_model, Dataset_original, Metric_label, Data_set)]

### the x axis spans the feature counts in the table, so that no run
### falls outside it
xLimits <- range(dtPlot$N_features, na.rm = T) + c(-1, 1)
xBreaks <- pretty(xLimits)

## plots ----------------------------------------------------------------------

for (metricName in levels(dtPlot$Metric_label)) {
  p <- MakeMetricPlot(metricName)
  if (!is.null(p)) {
    SaveMetricPlots(p, "plot-msd-", metricName)
    SaveMetricPlots(MakeRawPlot(metricName), "plot-raw-", metricName)
  }
}
