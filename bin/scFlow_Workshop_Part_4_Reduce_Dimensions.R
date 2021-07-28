# scFlow Workshop - Part IV - Reduce Dimensions
#  Combiz Khozoie <c.khozoie@imperial.ac.uk>

#   ____________________________________________________________________________
#   Initialization                                                          ####

n_cores <- future::availableCores(methods = "mc.cores")
options(mc.cores = n_cores)

##  ............................................................................
##  Load packages                                                           ####
library(scFlow)
library(parallel)
library(SingleCellExperiment)

#   ____________________________________________________________________________
#   Specify Inputs                                                          ####

sce_path <- "~/Documents/workshopscflow/MS_Example_Integrated"

sce <- read_sce(sce_path, read_metadata = TRUE)

#   ____________________________________________________________________________
#   Reduce Dimensions                                                       ####

##  ............................................................................
##  UMAP                                                                    ####
sce <- reduce_dims_sce(
  sce,
  input_reduced_dim = c("PCA", "Liger"),
  reduction_methods = c("UMAP","UMAP3D"),
  vars_to_regress_out = c("nCount_RNA","pc_mito"),
  umap_pca_dims = 30,
  umap_n_neighbors = 35,
  umap_n_components = 2,
  umap_init = "spectral",
  umap_metric = "euclidean",
  umap_n_epochs = 200,
  umap_learning_rate = 1,
  umap_min_dist = 0.4,
  umap_spread = 0.85,
  umap_set_op_mix_ratio = 1,
  umap_local_connectivity = 1,
  umap_repulsion_strength = 1,
  umap_negative_sample_rate = 5,
  umap_fast_sgd = FALSE
)

plot_reduced_dim(
  sce,
  feature_dim = "manifest",
  reduced_dim = "UMAP_Liger",
  size = 0.7
)

##  ............................................................................
##  Alternative method: tSNE                                                ####
sce <- reduce_dims_sce(
  sce,
  input_reduced_dim = c("PCA", "Liger"),
  reduction_methods = "tSNE",
  vars_to_regress_out = c("nCount_RNA","pc_mito"),
  tsne_dims = 2,
  tsne_initial_dims = 30,
  tsne_perplexity = 50,
  tsne_theta = 0.5,
  tsne_stop_lying_iter = 250,
  tsne_mom_switch_iter = 250,
  tsne_max_iter = 1000,
  tsne_pca_center = TRUE,
  tsne_pca_scale = FALSE,
  tsne_normalize = TRUE,
  tsne_momentum = 0.5,
  tsne_final_momentum = 0.8,
  tsne_eta = 1000,
  tsne_exaggeration_factor = 12
)

plot_reduced_dim(
  sce,
  feature_dim = "manifest",
  reduced_dim = "tSNE_Liger",
  size = 0.7
)


#   ____________________________________________________________________________



##  ............................................................................
##  Save SCE                                                               ####

write_sce(
  sce,
  "~/Documents/workshopscflow/MS_Example_DimRed",
  write_metadata = TRUE
)
