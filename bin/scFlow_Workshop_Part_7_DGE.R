# scFlow Workshop - Part VII - Differential Gene Expression
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
options("scflow_species" = "human")

ensembl_path <- "~/Documents/junk/src/ensembl-ids/ensembl_mappings.tsv"

sce_path <- "~/Documents/workshopscflow/MS_Example_Final"

sce <- read_sce(sce_path)

#   ____________________________________________________________________________
#   Set up differential gene expression                                    ####

celltype <- "Oligo"

sce_subset <- sce[, sce$cluster_celltype == celltype]

#   ____________________________________________________________________________
#   Perform differential gene expression                                    ####

de_results <- perform_de(
  sce_subset,
  de_method = "MASTZLM",
  mast_method = "glmer",
  min_counts = 1,
  min_cells_pc = 0.1,
  rescale_numerics = TRUE,
  pseudobulk = FALSE,
  celltype_var = "cluster_celltype",
  sample_var = "manifest",
  dependent_var = "diagnosis",
  ref_class = "Control",
  confounding_vars = c("cngeneson", "pc_mito"),
  random_effects_var = "manifest",
  fc_threshold = 1.1,
  pval_cutoff = 0.05,
  force_run = FALSE,
  max_cores = NULL,
  ensembl_mapping_file = ensembl_path
)

#   ____________________________________________________________________________
##  Save Outputs                                                            ####

result <- names(de_results)[1]

report_de(
  de_results[[result]],
  report_file = paste0(result, "_scflow_de_report")
  )

write.table(
  de_results[[result]],
  file = file.path(getwd(), paste0(result, "_DE.tsv")),
  quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE
  )

#   ____________________________________________________________________________
#   Plots for DGE                                                           ####

p <- volcano_plot(
  de_results[[result]],
  fc_threshold = 1.189207,
  pval_cutoff = 0.05,
  n_label = 10
  )
p

p <- plot_violin(
  sce,
  group_var = "diagnosis",
  subset_var = "cluster_celltype",
  subset_group = "Oligo",
  gene = "GRIA4",
  var_order = c("Control", "MS"),
  alpha = 0.2,
  size = 0.06
)
p
