# scFlow Workshop - Part IX - Dirichlet Analysis
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

sce_path <- "~/Documents/workshopscflow/MS_Example_Final"

sce <- read_sce(sce_path)

#   ____________________________________________________________________________
#   Perform a Dirichlet model analysis of relative cell-type frequencies   ####

results <- model_celltype_freqs(
  sce,
  unique_id_var = "manifest",
  celltype_var = "cluster_celltype",
  dependent_var = "diagnosis",
  ref_class = "Control",
  var_order = c("Control", "MS")
)

## ............................................................................
## Save Outputs ####

report_celltype_model(results)
