# scFlow Workshop - Part VIII - Impacted Pathway Analysis
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
#gene_file <- "~/Documents/workshopscflow/data/tidy/DE_Results/Control_vs_diagnosisMS_DE.tsv"
gene_file <- file.path(getwd(), "Control_vs_diagnosisMS_DE.tsv")

#   ____________________________________________________________________________
#   Perform impacted pathway analysis                                       ####

enrichment_result <- find_impacted_pathways(
  gene_file           = gene_file,
  enrichment_tool     = 'WebGestaltR',
  enrichment_method   = 'ORA',
  enrichment_database = 'GO_Biological_Process'
)

#   ____________________________________________________________________________
##  Save Outputs                                                            ####

report_impacted_pathway(enrichment_result)
