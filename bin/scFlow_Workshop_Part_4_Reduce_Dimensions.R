# scFlow Workshop - Part IV - Reduce Dimensions
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

#sce_path <- "~/Documents/workshopscflow/MS_Example_Integrated"
sce_path <- file.path(getwd(), "MS_Example_Integrated")


sce <- read_sce(sce_path, read_metadata = TRUE)

#   ____________________________________________________________________________
#   Reduce Dimensions                                                       ####

##  ............................................................................
##  UMAP                                                                    ####
sce <- reduce_dims_sce(
  sce,
  input_reduced_dim    = c("PCA", "Liger"),
  reduction_methods    = c("UMAP","UMAP3D"),
  vars_to_regress_out  = c("nCount_RNA","pc_mito"),
  pca_dims             = 30,
  n_neighbors          = 35,
  n_components         = 2,
  init                 = "spectral",
  metric               = "euclidean",
  n_epochs             = 200,
  learning_rate        = 1,
  min_dist             = 0.4,
  spread               = 0.85,
  set_op_mix_ratio     = 1,
  local_connectivity   = 1,
  repulsion_strength   = 1,
  negative_sample_rate = 5,
  fast_sgd             = FALSE
)

plot_reduced_dim(
  sce,
  feature_dim = "manifest",
  reduced_dim = "UMAP_Liger",
  size = 0.7
)

plot_reduced_dim(
  sce,
  feature_dim = "manifest",
  reduced_dim = "UMAP_PCA",
  size = 2
)

##  ............................................................................
##  Alternative method: tSNE                                                ####
sce <- reduce_dims_sce(
  sce,
  input_reduced_dim        = c("PCA", "Liger"),
  reduction_methods        = "tSNE",
  vars_to_regress_out      = c("nCount_RNA","pc_mito"),
  dims                     = 2,
  initial_dims             = 30,
  perplexity               = 50,
  theta                    = 0.5,
  stop_lying_iter          = 250,
  mom_switch_iter          = 250,
  max_iter                 = 1000,
  pca_center               = TRUE,
  pca_scale                = FALSE,
  normalize                = TRUE,
  momentum                 = 0.5,
  final_momentum           = 0.8,
  eta                      = 1000,
  exaggeration_factor      = 12
)

plot_reduced_dim(
  sce,
  feature_dim = "manifest",
  reduced_dim = "tSNE_Liger",
  size        = 0.7
)

plot_reduced_dim(
  sce,
  feature_dim = "manifest",
  reduced_dim = "tSNE_PCA",
  size        = 0.7
)

#   ____________________________________________________________________________



##  ............................................................................
##  Save SCE                                                               ####

write_sce(
  sce,
  file.path(getwd(), "MS_Example_DimRed"),
  write_metadata = TRUE
)
