# scFlow Workshop - Part VI - Celltyping
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
options("scflow_reddimplot_pointsize" = 0.3)
options("scflow_reddimplot_alpha" = 0.5)

sce_path <- "~/Documents/workshopscflow/MS_Example_Clustered"
ctd_folder <- "~/Documents/junk/refs/ctd"

sce <- read_sce(sce_path)

#   ____________________________________________________________________________
#   Celltyping                                                              ####

sce <- map_celltypes_sce(
  sce,
  ctd_folder = ctd_folder,
  clusters_colname = "clusters",
  cells_to_sample = 10000,
  species = "human"
)

plot_reduced_dim(
  sce,
  feature_dim = "cluster_celltype",
  reduced_dim = "UMAP_Liger",
  label_clusters = TRUE
)

#   ____________________________________________________________________________
#   Report Celltype Metrics                                                 ####

sce <- annotate_celltype_metrics(
  sce,
  unique_id_var = "manifest",
  cluster_var = "clusters",
  celltype_var = "cluster_celltype",
  facet_vars = c("manifest","diagnosis","sex"),
  metric_vars = c("pc_mito","pc_ribo","total_counts","total_features_by_counts")
)

report_celltype_metrics(sce)

##  ............................................................................
##  OPTIONAL: Manually Revise Cell-types                                    ####

write_celltype_mappings(sce, folder_path = getwd())

# To revise, edit celltype_mappings.tsv and reload below
celltype_mappings <- read_celltype_mappings(
  file.path(getwd(), "celltype_mappings_revised.tsv")
  )

sce <- map_custom_celltypes(
  sce,
  mappings = celltype_mappings,
  clusters_colname = "clusters"
)

#   ____________________________________________________________________________
#   Inspect Marker Genes                                                    ####

sce@metadata$markers$cluster_celltype$marker_plot

plot_reduced_dim_gene(
  sce,
  gene = "PLP1",
  reduced_dim = "UMAP_Liger",
  size = 1
)

sce@metadata$markers$clusters$marker_plot

plot_reduced_dim_gene(
  sce,
  gene = "PLXDC2",
  reduced_dim = "UMAP_Liger",
  size = 1
)


#   ____________________________________________________________________________
##  Save SCE                                                               ####

write_sce(
  sce,
  "~/Documents/workshopscflow/MS_Example_Final"
)
