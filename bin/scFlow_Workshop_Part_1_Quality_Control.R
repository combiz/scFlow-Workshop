# scFlow Workshop - Part I - Quality Control
#  Combiz Khozoie <c.khozoie@imperial.ac.uk>

#   ____________________________________________________________________________
#   Initialization                                                          ####

n_cores <- future::availableCores(methods = "system")
options(mc.cores = n_cores)

##  ............................................................................
##  Load packages                                                           ####
library(scFlow)
library(parallel)
library(SingleCellExperiment)

#   ____________________________________________________________________________
#   Specify Inputs                                                          ####

ensembl_path <- "~/Documents/workshopscflow/refs/ensembl_mappings.tsv"

mat_path <- "~/Documents/workshopscflow/data/raw/dovim/outs/raw_feature_bc_matrix"
samplesheet_path <- "~/Documents/workshopscflow/conf/SampleSheet.tsv"

#   ____________________________________________________________________________
#   Import gene-cell matrix and metadata                                    ####

mat <- read_sparse_matrix(mat_path)

metadata <- read_metadata(
  unique_key       = "dovim",
  key_colname      = "manifest",
  samplesheet_path = samplesheet_path#,
  #col_classes     = list("aplevel" = "factor")
)

sce <- generate_sce(mat, metadata)

#   ____________________________________________________________________________
#   Ambient RNA profiling                                                   ####

# !! Time consuming; for the workshop, read the SCE below instead.

sce <- find_cells(
  sce,
  lower  = 100,
  retain = "auto",
  niters = 10000
)


sce <- read_sce(
  "~/Documents/workshopscflow/data/tidy/dovim_post_find_cells",
  read_metadata = TRUE
  )


#   ____________________________________________________________________________
#   Annotate SCE with Gene Info & QC Metrics                                ####

sce <- annotate_sce(
  sce                  = sce,
  min_library_size     = 100,
  max_library_size     = "adaptive",
  min_features         = 100,
  max_features         = "adaptive",
  max_mito             = "adaptive",
  min_ribo             = 0,
  max_ribo             = 1,
  min_counts           = 2,
  min_cells            = 2,
  drop_unmapped        = TRUE,
  drop_mito            = TRUE,
  drop_ribo            = FALSE,
  nmads                = 4.0,
  ensembl_mapping_file = ensembl_path,
  species              = "human"
)

sce <- filter_sce(
  sce,
  filter_genes = TRUE,
  filter_cells = TRUE
)

#   ____________________________________________________________________________
#   Identify and filter out doublets/multiplets                             ####

sce <- find_singlets(
  sce                 = sce,
  singlet_find_method = "doubletfinder",
  vars_to_regress_out = c("nCount_RNA", "pc_mito"),
  pca_dims            = 10,
  var_features        = 2000,
  doublet_rate        = 0,
  dpk                 = 8,
  pK                  = 0.02
)

sce <- filter_sce(
    sce,
    filter_genes = TRUE,
    filter_cells = TRUE
)

##  ............................................................................
##  Save Outputs                                                            ####

report_qc_sce(sce)

# Save SingleCellExperiment
write_sce(
  sce         = sce,
  folder_path = file.path(getwd(), paste0(unique(sce$manifest), "_sce"))
)
