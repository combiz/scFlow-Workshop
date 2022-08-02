# scFlow Workshop - Part III - Integration
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

sce_path <- file.path(getwd(), "MS_Example_Merged")
sce      <- read_sce(sce_path)

#   ____________________________________________________________________________
#   Integrate                                                              ####

sce <- integrate_sce(
  sce,
  method             = "Liger",
  unique_id_var      = "manifest",
  take_gene_union    = FALSE,
  remove_missing     = TRUE,
  num_genes          = 3000,
  combine            = "union",
  keep_unique        = FALSE,
  capitalize         = FALSE,
  use_cols           = TRUE,
  #k                 = 30,
  k                  = 15,
  lambda             = 5.0,
  thresh             = 0.0001,
  max_iters          = 100,
  nrep               = 1,
  rand_seed          = 1,
  knn_k              = 20,
  k2                 = 500,
  #prune_thresh      = 0.2,
  prune_thresh       = 0.1,
  ref_dataset        = NULL,
  min_cells          = 2,
  quantiles          = 50,
  nstart             = 10,
  resolution         = 1,
  dims_use           = NULL,
  dist_use           = "CR",
  center             = FALSE,
  small_clust_thresh = 0
)

##  ............................................................................
##  Save SCE                                                               ####

write_sce(
  sce,
  file.path(getwd(), "MS_Example_Integrated"),
  write_metadata = TRUE
)
