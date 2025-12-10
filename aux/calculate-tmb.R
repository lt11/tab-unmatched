## header ---------------------------------------------------------------------

options(scipen = 999)
options(stringsAsFactors = F)
rm(list = ls())
library(data.table)
library(this.path)

## settings -------------------------------------------------------------------

### folders
dirBase <- dirname(this.dir())
pathRes <- file.path(dirBase, "res")
pathTab <- file.path(dirBase, "tab")
pathOut <- file.path(pathRes, "metadata-tmb.txt")

### kit sizes
kitSize <- c(nimblegen_seqcap_ez_v3_kit = 151.7,
             hgsc_vcrome_kit =  45.1,
             agilent_custom_v2_kit =  56.8)

### the ANN values to count
tmbVariants <- c("missense_variant",
                 "rare_amino_acid_variant",
                 "initiator_codon_variant",
                 "stop_gained",
                 "stop_lost",
                 "start_lost",
                 "frameshift_variant",
                 "frameshift_variant&splice_acceptor_variant&splice_donor_variant&splice_region_variant&intron_variant",
                 "frameshift_variant&splice_acceptor_variant&splice_region_variant&intron_variant",
                 "frameshift_variant&splice_donor_variant&splice_region_variant&intron_variant",
                 "frameshift_variant&splice_region_variant",
                 "frameshift_variant&stop_gained",
                 "frameshift_variant&stop_gained&splice_region_variant",
                 "frameshift_variant&stop_lost",
                 "disruptive_inframe_deletion")

## clmnt ----------------------------------------------------------------------

### all preds* csv files in pathRes
pathCsv <- list.files(pathRes, pattern = "preds-.*\\.csv$", full.names = T)

### read all files and bind them
lsCsvFiles <- lapply(pathCsv, fread)
dtAllCsvFiles  <- rbindlist(lsCsvFiles, use.names = T, fill = T)

### rename the column for convenience
colnames(dtAllCsvFiles) <- gsub("\\s+", "_", colnames(dtAllCsvFiles))

### for each fastq_id count the somatic rows where ANN is in the tmbVariants
dtCounts <- dtAllCsvFiles[
  ANN %in% tmbVariants,
  .(
    matched_somatic = sum(target == 1, na.rm = T),
    xgm_somatic = sum(xgm_preds == "somatic", na.rm = T),
    lgbm_somatic = sum(lgbm_preds == "somatic", na.rm = T),
    tabnet_somatic = sum(tabnet_preds == "somatic", na.rm = T)
  ),
  by = fastq_id
]

### read the metadata
dtMetadata <- fread(file.path(pathTab, "metadata.csv"), stringsAsFactors = F)
colnames(dtMetadata) <- gsub("\\s+", "_", colnames(dtMetadata))

### inner‐join on fastq_id to keep only those present in both
dtMetadataCounts <- merge(dtMetadata, dtCounts, by = "fastq_id",
                          all = F) ### only fastq_id in both

### add kit_size
dtMetadataCounts[, kit_size := kitSize[kit]]

### add TMB columns for each model
dtMetadataCounts[, `:=`(
  matched_TMB = matched_somatic / kit_size,
  xgm_TMB = xgm_somatic / kit_size,
  lgbm_TMB = lgbm_somatic / kit_size,
  tabnet_TMB = tabnet_somatic / kit_size
)]

### add TMB columns to compare with the McLaughlin paper
### (wrong as it seems they used also synonymous variants) 
dtMetadataCounts[, `:=`(
  matched_TMB_comparison = matched_somatic / 41,
  xgm_TMB_comparison = xgm_somatic / 41,
  lgbm_TMB_comparison = lgbm_somatic / 41,
  tabnet_TMB_comparison = tabnet_somatic / 41
)]

### what follows is the calculation of the gTMB which
### we invented to compare against the original paper;
### gTMB is a wrong implementation of the TMB, and it turned out
### that it was not implemented in the original paper, despite
### gTMB is what they claim to have implemented in the methods;
### che coglioni del cazzo

# ### global counts without filtering by tmbVariants
# dtCountsGlob <- dtAllCsvFiles[
#   ,
#   .(
#     matched_somatic_glob = sum(target == 1, na.rm = T),
#     xgm_somatic_glob = sum(xgm_preds == "somatic", na.rm = T),
#     lgbm_somatic_glob = sum(lgbm_preds == "somatic", na.rm = T),
#     tabnet_somatic_glob = sum(tabnet_preds == "somatic", na.rm = T)
#   ),
#   by = fastq_id
# ]
# 
# ### inner‐join on fastq_id to keep only those present in both
# dtMetadataCounts <- merge(dtMetadataCounts, dtCountsGlob, by = "fastq_id",
#                           all = F) ### only fastq_id in both
# 
# ### add gTMB columns for each model
# dtMetadataCounts[, `:=`(
#   matched_TMB_glob = matched_somatic_glob / kit_size,
#   xgm_TMB_glob = xgm_somatic_glob / kit_size,
#   lgbm_TMB_glob = lgbm_somatic_glob / kit_size,
#   tabnet_TMB_glob = tabnet_somatic_glob / kit_size
# )]
# 
# ### add gTMB columns to compare with the McLaughlin paper
# dtMetadataCounts[, `:=`(
#   matched_TMB_glob_comparison = matched_somatic_glob / 41,
#   xgm_TMB_glob_comparison = xgm_somatic_glob / 41,
#   lgbm_TMB_glob_comparison = lgbm_somatic_glob / 41,
#   tabnet_TMB_glob_comparison = tabnet_somatic_glob / 41
# )]

### save
fwrite(x = dtMetadataCounts, file = pathOut, append = F, quote = F, sep = "\t")
