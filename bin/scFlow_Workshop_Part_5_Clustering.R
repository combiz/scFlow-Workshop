# scFlow Workshop - Part V - Clustering
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

options("scflow_reddimplot_pointsize" = 0.3)
options("scflow_reddimplot_alpha" = 0.5)

sce_path <- "~/Documents/workshopscflow/MS_Example_DimRed"

sce <- read_sce(sce_path, read_metadata = TRUE)

#   ____________________________________________________________________________
#   Clustering                                                              ####

sce <- cluster_sce(
  sce,
  cluster_method = "leiden",
  reduction_method = "UMAP_Liger",
  res = 0.01, # Lower is more coarse
  k = 100, # Higher is more coarse
  louvain_iter = 1
)

plot_reduced_dim(
  sce,
  feature_dim = "clusters",
  reduced_dim = "UMAP_Liger",
  label_clusters = TRUE,
  size = 0.3
)

plot_reduced_dim(
  sce,
  feature_dim = "clusters",
  reduced_dim = "UMAP_Liger",
  size = 0.3,
  highlight_feature = "4"
)

#   ____________________________________________________________________________
#   Generate Integration Report                                             ####

sce <- annotate_integrated_sce(
  sce,
  categorical_covariates = c("manifest", "diagnosis", "sex", "seqdate"),
  input_reduced_dim = "UMAP"
)

report_integrated_sce(sce)

plot_reduced_dim(
  sce,
  feature_dim = "seqdate",
  reduced_dim = "UMAP_PCA"
  #reduced_dim = "UMAP_Liger"
)


#   ____________________________________________________________________________
#   Cell-type Annotation                                                    ####


#   ____________________________________________________________________________
##  Save SCE                                                               ####

write_sce(
  sce,
  "~/Documents/workshopscflow/MS_Example_Clustered"
)
