# scFlow Workshop - Part II - Merge Post-QC Samples
#  Combiz Khozoie <c.khozoie@imperial.ac.uk>

#   ____________________________________________________________________________
#   Initialization                                                          ####

n_cores <- future::availableCores(methods = "mc.cores")
options(mc.cores = n_cores)

##  ............................................................................
##  Load packages                                                           ####
library(scFlow)
library(SingleCellExperiment)

#   ____________________________________________________________________________
#   Specify Inputs                                                          ####

ensembl_path <- "~/Documents/workshopscflow/refs/ensembl_mappings.tsv"
sces_folder <- "~/Documents/workshopscflow/data/tidy/sces_post_qc"

#   ____________________________________________________________________________
#   Merge of Post-QC Samples                                                ####

# Retrieve list of available SCE folders
sce_l <- list.files(
  sces_folder,
  full.names   = TRUE,
  include.dirs = TRUE
  )

sce <- merge_sce(
  sce_l,
  ensembl_mapping_file = ensembl_path
)

#   ____________________________________________________________________________
#   Generate Post-Merge Report                                              ####

sce <- annotate_merged_sce(
  sce,
  plot_vars     = c(
    "total_features_by_counts","total_counts","pc_mito","pc_ribo"
    ),
  unique_id_var = "manifest",
  facet_vars    = c("seqdate", "diagnosis"),
  outlier_vars  = c("total_features_by_counts", "total_counts")
)

report_merged_sce(sce)

#   ____________________________________________________________________________
#   Save SCE                                                                ####

# Optional: filter samples
sce <- sce[, sce$manifest != "honiz"]

write_sce(
  sce,
  file.path(getwd(), "MS_Example_Merged")
)
