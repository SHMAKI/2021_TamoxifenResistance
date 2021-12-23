# before staring R, putting a .Renviron file containing the line RETICULATE_PYTHON="/path/to/python3" into the user's home directory
# https://stackoverflow.com/questions/50145643/unable-to-change-python-path-in-reticulate
library(reticulate) # for UMAP clustering
py_config() 


library(tidyverse)
library(Seurat)
#monocle 3 alpha version
#https://cole-trapnell-lab.github.io/monocle3/monocle3_docs/#constructing-single-cell-trajectories
#devtools::install_github('cole-trapnell-lab/monocle3')
library(monocle3)
library(ggsci)
# library("reshape2")
# library("Matrix")
# library("ggplot2")
# library(clusterProfiler)
# library(org.Hs.eg.db)
# library(data.table)
library(scales)
library(ggpubr)
source("src/theme.R")

#DIR6 <- "results6_monocle3_n1500_remove_mito_RPratio"
#dir.create(DIR6)
Seu_allreg <- readRDS("processed_data/Seu_allreg.rds")
df_filteredgenes <- readRDS("processed_data/df_filteredgenes.rds")

#Generate a cell_data_set from Seurat obj
umi_matrix <- as.matrix(Seurat::GetAssayData(object = Seu_allreg, slot = "counts"))
sample_info <- data.frame(Seu_allreg@meta.data,
                          week = rownames(Seu_allreg@meta.data), 
                          stringsAsFactors = F)

sample_info$week[grep("Control", rownames(sample_info))] <- 0
sample_info$week[grep("TAM3weeks", rownames(sample_info))] <- 3
sample_info$week[grep("TAM6weeks", rownames(sample_info))] <- 6
sample_info$week[grep("TAM9weeks", rownames(sample_info))] <- 9

#--------Fig.2b correlation in each time point--------
dat5_CTRL <- umi_matrix[,grep("Control", colnames(umi_matrix))]
dat5_W3 <- umi_matrix[,grep("TAM3weeks", colnames(umi_matrix))]
dat5_W6 <- umi_matrix[,grep("TAM6weeks", colnames(umi_matrix))]
dat5_W9 <- umi_matrix[,grep("TAM9weeks", colnames(umi_matrix))]
#calc cor
dat5_CTRL_cor <- cor(dat5_CTRL, method = "p")
dat5_W3_cor <- cor(dat5_W3, method = "p")
dat5_W6_cor <- cor(dat5_W6, method = "p")
dat5_W9_cor <- cor(dat5_W9, method = "p")
#X[upper.tri(X)]で上三角成分
dat5_CTRL_cor_v <- dat5_CTRL_cor[upper.tri(dat5_CTRL_cor)]
dat5_W3_cor_v <- dat5_W3_cor[upper.tri(dat5_W3_cor)]
dat5_W6_cor_v <- dat5_W6_cor[upper.tri(dat5_W6_cor)]
dat5_W9_cor_v <- dat5_W9_cor[upper.tri(dat5_W9_cor)]

DAT_cor <- data.frame(value=c(dat5_CTRL_cor_v, dat5_W3_cor_v, dat5_W6_cor_v, dat5_W9_cor_v),
                      week=c(rep("0", length(dat5_CTRL_cor_v)), 
                             rep("3", length(dat5_W3_cor_v)),
                             rep("6", length(dat5_W6_cor_v)),
                             rep("9", length(dat5_W9_cor_v))
                      )
)

n_fun <- function(x){
  return(data.frame(y = median(x)+0.01, label = paste0(round(median(x),4))))
}
g <- ggplot(DAT_cor, aes(x=week, y=value)) +
  geom_boxplot(width=.5, fill="lightgray",outlier.fill = NA, outlier.color = NA, color="black") +
  stat_summary(fun.data = n_fun, geom="text", color="red") +
  theme_bw() + ng1 + ylab("Correlation coefficient\nbetween cells") + xlab("week")
my_comparisons <- list( c("0", "3"), c("3", "6"), c("6", "9") ) # Add pairwise comparisons p-value
g <- g + stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", bracket.size = 0.5, step.increase = 0.02, tip.length = 0.01)
g <- g + coord_cartesian(ylim=c(0.8, 1.1))
ggsave(file =paste0("Fig/Fig_2b.pdf"), plot=g, width=4, height=4)

#--------create monocle obj-----
gene_annotation <- df_filteredgenes
colnames(gene_annotation)[colnames(gene_annotation) == "gene_name"] <- "gene_short_name"
rownames(gene_annotation) <- rownames(umi_matrix)

cds <- new_cell_data_set(as(umi_matrix, "sparseMatrix"),
                         cell_metadata = sample_info,
                         gene_metadata = gene_annotation)

# Pre-process the data
cds <- preprocess_cds(cds, num_dim = 100, residual_model_formula_str = "~ S.Score + G2M.Score")
pdf(paste0("intermediate_figs/plot_pc_variance_explained.pdf"))
plot_pc_variance_explained(cds)
dev.off()

#--------Fig.2c summarize cell cycle phase--------
cell_group_df_for_rate = tibble::tibble(cell=row.names(colData(cds)), clusters=cds@clusters$UMAP$clusters, week=colData(cds)$week, Phase = colData(cds)$Phase)
POP_phase <- cell_group_df_for_rate %>% group_by(week, Phase, .drop = F) %>% summarise (n = n()) %>% mutate(freq = n / sum(n))
g <- ggplot(POP_phase,  aes(x = week, y = freq, fill = Phase))
g <- g + geom_bar(stat = "identity", position = "fill") + theme_bw(base_size = 16)
g <- g + geom_text(aes(label=ifelse(freq >= 0.05, paste0(round(freq, 2) * 100, "%"), "")),
                   position=position_stack(vjust=0.5), colour="black") + ylab("Frequency (%)")
g <- g + scale_y_continuous(labels = scales::percent) + ng1 
g <- g + scale_color_manual(values = c("skyblue", "#CC6677", "darkgreen")) + scale_fill_manual(values = c("skyblue", "#CC6677", "darkgreen"))
g_horiz <- g + coord_flip() + guides(fill = guide_legend(title = "cell cycle", override.aes = list(size = 4)))
ggsave(g_horiz, filename = paste0("Fig/Fig2c.pdf"), w=8, h=4)


#cds <- reduce_dimension(cds, max_components = 3, cores = 10, reduction_method = "UMAP")
cds <- reduce_dimension(cds, max_components = 3, reduction_method = "UMAP")
# cluster cells
cds <- cluster_cells(cds, partition_qval = 5e-2, reduction_method = "UMAP", k = 25, resolution=5e-2, verbose = T, random_seed = 0)


#check how many partitions there are (I had two partitions):
levels(cds@clusters$UMAP$partitions)
#set them all to "1"
cds@clusters$UMAP$partitions[cds@clusters$UMAP$partitions == "2"] <- "1"
#checking again still says there are two types of partitions, "1" and "2"
levels(cds@clusters$UMAP$partitions)
#but there are no cells in partition "2":
length(cds@clusters$UMAP$partitions[cds@clusters$UMAP$partitions == "2"])
#they are all in partition 1:
length(cds@clusters$UMAP$partitions[cds@clusters$UMAP$partitions == "1"])
#check plot:
plot_cells(cds, color_cells_by = "partition")

# 3d UMAP plot in Fig. S5
# plot_cells_3d(cds, color_cells_by="week")#, color_cells_by="partition", group_cells_by="partition")#, show_trajectory_graph = F)# + scale_color_nejm()

fig2_d <- plot_cells(cds, color_cells_by = "Phase", show_trajectory_graph = F, group_label_size = 0) + theme_bw() + ng1 + scale_color_manual(values = c("skyblue", "#CC6677", "darkgreen")) + guides(color = guide_legend(title = "cell cycle", override.aes = list(size = 4)))
ggsave(fig2_d, filename = paste0("Fig/Fig_2d.pdf"),w=5,h=4)

fig2_e <- plot_cells(cds, color_cells_by = "week", show_trajectory_graph = F, group_label_size = 0) + theme_bw() + ng1 + scale_color_nejm()
ggsave(fig2_e, filename = paste0("Fig/Fig_2e.pdf"),w=5,h=4)

fig2_f <- plot_cells(cds, color_cells_by = "cluster", show_trajectory_graph = F, group_label_size = 8) + theme_bw() + ng1 + scale_color_npg() + guides(color = guide_legend(title = "subgroup", override.aes = list(size = 4)))
ggsave(fig2_f, filename = paste0("Fig/Fig_2f.pdf"),w=5.5,h=4)


# learn graph
cds = learn_graph(cds, use_partition = F, close_loop = T, verbose = T, list(minimal_branch_len = 12,
                                                                            euclidean_distance_ratio = 2.5,
                                                                            geodesic_distance_ratio = 1/2,
                                                                            orthogonal_proj_tip = F))
cds = order_cells(cds, reduction_method = "UMAP", verbose = T)#, root_pr_nodes=get_earliest_principal_node(cds_1125))
plot_cells(cds)

fig3_a <- plot_cells(cds_1125, color_cells_by = "cluster", show_trajectory_graph = T, group_label_size = 0, label_leaves=FALSE) + theme_bw() + ng1 + scale_color_npg() + guides(color = guide_legend(title = "subgroup", override.aes = list(size = 4)))
ggsave(fig3_a, filename = paste0(DIR, "Fig3_a.pdf"),w=5.5,h=4)

fig3_b <- plot_cells(cds_1125, color_cells_by = "week", show_trajectory_graph = T, group_label_size = 0, label_leaves=FALSE) + theme_bw() + ng1 + scale_color_nejm()
ggsave(fig3_b, filename = paste0(DIR, "Fig3_b.pdf"),w=5,h=4)

fig3_c <- plot_cells(cds_1125, color_cells_by = "pseudotime", label_branch_points=TRUE, label_leaves=FALSE) + theme_bw() + ng1
ggsave(fig3_c, filename = paste0(DIR, "Fig3_c.pdf"),w=5.5,h=4)



pdf(paste0("intermediate_figs/UMAP_raw.pdf"), w=7, h=6)
plot_cells(cds, show_trajectory_graph = F) + ng1
dev.off()
pdf(paste0(DIR6, "/UMAP_week2.pdf"), w=7, h=6)
plot_cells(cds, color_cells_by="week", label_cell_groups=FALSE, show_trajectory_graph = F) + scale_color_nejm() + ng1 #direct cell_metadata
dev.off()
pdf(paste0(DIR6, "/UMAP_CCPhase_fromSeurat.pdf"), w=7, h=6)
plot_cells(cds, color_cells_by="Phase", label_cell_groups=FALSE, show_trajectory_graph = F)  + scale_color_aaas() + ng1 #direct cell_metadata
dev.off()
saveRDS(cds, file = paste0(DIR6, "/cds.rds"))

#Group cells into clusters
#regires louvain python pachages
cds <- cluster_cells(cds, partition_qval = 1) #, partition_qval = 1, k = 30)#, partition_qval = 0.01)
pdf(paste0(DIR6, "/UMAP_Monoclepartition.pdf"))
plot_cells(cds, color_cells_by="partition", group_cells_by="partition", show_trajectory_graph = F)
dev.off()
pdf(paste0(DIR6, "/UMAP_MonocleCluster.pdf"), w=7, h=6)
plot_cells(cds, color_cells_by="cluster", group_label_size = 10,  show_trajectory_graph = F,
           group_cells_by="cluster") + ng1
dev.off()

pdf(paste0(DIR6, "/UMAP_MonocleCluster2.pdf"), w=7, h=6)
plot_cells(cds, color_cells_by="cluster", 
           group_cells_by="cluster", group_label_size = 10, show_trajectory_graph = F) + ng1
dev.off()

pdf(paste0(DIR6, "/UMAP_MonocleCluster_F.pdf"), w=7, h=6)
plot_cells(cds, color_cells_by="cluster", 
           group_cells_by="cluster", 
           label_cell_groups=FALSE, show_trajectory_graph = F)+ ng1
dev.off()

#Find marker genes expressed by each cluster
# marker_test_res = top_markers(cds, group_cells_by="cluster", reference_cells=1000, cores=4)
# top_specific_markers = marker_test_res %>% 
#   #dplyr::filter(fraction_expressing >= 0.10) %>% 
#   dplyr::filter(marker_test_q_value < 0.05) %>% 
#   dplyr::filter(!(str_detect(gene_short_name, "^MT-") | str_detect(gene_short_name, "^RPS") | str_detect(gene_short_name, "^RPL") | (gene_short_name == "") )) %>%
#   group_by(cell_group) %>% top_n(1, -marker_test_q_value)
# top_specific_marker_ids = unique(top_specific_markers %>% pull(gene_id))
# pdf(paste0(DIR6, "/markergenes_fromMonocle.pdf"))
# plot_genes_by_group(cds,
#                     top_specific_marker_ids,
#                     group_cells_by="cluster",
#                     ordering_type="maximal_on_diag",
#                     max.size=3)
# dev.off()

#Constructing single-cell trajectories
cds <- learn_graph(cds, #verbose = T,
                   learn_graph_control = list(minimal_branch_len = 10, 
                                              euclidean_distance_ratio = 1,
                                              geodesic_distance_ratio = 1/3,
                                              orthogonal_proj_tip = T
                                              )
                   )
pdf(paste0(DIR6, "/branch_raw.pdf"))
plot_cells(cds,
           color_cells_by = "week",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE, label_cell_groups=FALSE
            ) + scale_color_aaas()
dev.off()
##Order the cells in pseudotime
pdf(paste0(DIR6, "/branch_raw_branch.pdf"), w=7, h=6)
plot_cells(cds,
           color_cells_by = "week",
           label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE,　
           graph_label_size=1.5) + scale_color_nejm() + ng1
dev.off()
##Order the cells in pseudotime
pdf(paste0(DIR6, "/branch_raw_cluter.pdf"), w=7, h=6)
plot_cells(cds,
           color_cells_by = "cluster",
           label_cell_groups=T,
           label_leaves=F,
           label_branch_points=F, group_label_size = 12　
           ) + ng1 
dev.off()




#save.image()

# #choose the root
cds = order_cells(cds) #manual click w0 -> pseudo time = 0 
 
# # a helper function to identify the root principal points:
get_earliest_principal_node <- function(cds, week="0"){
  cell_ids <- which(colData(cds)[, "week"] == week)

  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                              (which.max(table(closest_vertex[cell_ids,]))))]
  root_pr_nodes
}
#cds = order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))

pdf(paste0(DIR6, "/branch_pseudotime_after.pdf"), w=8.5, h=6)
plot_cells(cds,
           color_cells_by = "pseudotime",
#           label_cell_groups=FALSE,
           label_leaves=FALSE,
#           label_branch_points=FALSE,
           group_label_size = 12,  graph_label_size=3
           ) + ng1
dev.off()

#Working with 3D trajectories
cds_3d = reduce_dimension(cds, max_components = 3)
cds_3d = cluster_cells(cds_3d,partition_qval = 1)#, k = 30, partition_qval = 0.01)
cds_3d = learn_graph(cds_3d, use_partition = F)#, learn_graph_control = list(minimal_branch_len = 3))
cds_3d = order_cells(cds_3d)#, root_pr_nodes=get_earliest_principal_node(cds))
cds_3d_plot_obj = plot_cells_3d(cds_3d, color_cells_by="week", show_trajectory_graph = F, cell_size = 100)# + scale_color_nejm()
cds_3d_plot_obj
cds_3d_plot_obj2 = plot_cells_3d(cds_3d, color_cells_by="cluster")#, show_trajectory_graph = F)# + scale_color_nejm()
cds_3d_plot_obj2
#animation https://www.medicalmed.press/2018/03/05/rstudio-ggplot2-animation/#i
devtools::install_github("dgrtwo/gganimate")
install.packages("rgl")
library(plotly)
scene = list(camera = list(eye = list(x = -1.25, y = 1.25, z = 1.25)))
scene2 = list(camera = list(eye = list(x = -1, y = 1, z = 1)))
ggplotly(cds_3d_plot_obj) %>% 
  layout(scene = scene)
play3d(cds_3d_plot_obj)


##Order the cells in pseudotime
pdf(paste0(DIR6, "/branch_raw3dto2d.pdf"))
plot_cells(cds_3d,
           color_cells_by = "week",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE, label_cell_groups=FALSE
) + scale_color_aaas()
dev.off()
##Order the cells in pseudotime
pdf(paste0(DIR6, "/branch_raw_branch3dto2d.pdf"), w=7, h=6)
plot_cells(cds_3d,
           color_cells_by = "week",
           label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE,　
           graph_label_size=1.5) + scale_color_nejm() + ng1
dev.off()
##Order the cells in pseudotime
pdf(paste0(DIR6, "/branch_raw_cluter3dto2d.pdf"), w=7, h=6)
plot_cells(cds_3d,
           color_cells_by = "cluster",
           label_cell_groups=T,
           label_leaves=F,
           label_branch_points=F, group_label_size = 12　
) + ng1 
dev.off()
pdf(paste0(DIR6, "/branch_pseudotime_after3dto2d.pdf"), w=8.5, h=6)
plot_cells(cds_3d,
           color_cells_by = "pseudotime",
           #           label_cell_groups=FALSE,
           label_leaves=FALSE,
           #           label_branch_points=FALSE,
           group_label_size = 12,  graph_label_size=3
) + ng1
dev.off()


saveRDS(cds, "cds_2d.rds")
saveRDS(cds_3d, "cds_3d.rds")

cds <- cds_3d
DIR6 <- "results6_monocle3_n1500_remove_mito_RPratio/3d"
dir.create(DIR6)
#Analyzing branches in single-cell trajectories
ciliated_cds_pr_test_res = graph_test(cds, neighbor_graph="principal_graph")#, cores=4)
pr_deg_ids = row.names(subset(ciliated_cds_pr_test_res, q_value < 0.05))

# plot_cells(cds, genes=c("ANXA2", "GAPDH", "LC3B", "CCND1"),
#            show_trajectory_graph=FALSE,
#            label_cell_groups=FALSE,
#            label_leaves=FALSE)

#As before, we can collect the trajectory-variable genes into modules:
gene_module_df = monocle3:::find_gene_modules(cds[pr_deg_ids,], resolution=c(0,10^seq(-6,-1)))
gene_module_df2 <- gene_module_df
colnames(gene_module_df2)[1] <- "ensembl_gene_id"
results.merge <- merge(gene_module_df2, res, by = "ensembl_gene_id", incomparables=NA)
fwrite(results.merge, paste0(DIR6, "/gene_module_df.txt"), sep = "\t")

newDF <- data.frame(colData(cds),clusters = clusters(cds))
g <- ggplot(newDF, aes(x=clusters, y=percent.mt)) + geom_boxplot() + theme_bw()+ ng1
ggsave(g, filename = paste0(DIR6, "/mito.rate.pdf"))
#Here we plot the aggregate module scores within each group of cell types as annotated by Packer & Zhu et al:
# cell_group_df = tibble::tibble(cell=row.names(colData(cds)), cell_group=colData(cds)$week)
# agg_mat = aggregate_gene_expression(cds, gene_module_df, cell_group_df)
# row.names(agg_mat) = stringr::str_c("Module ", row.names(agg_mat))
# pheatmap::pheatmap(agg_mat,
#                    scale="column", clustering_method="ward.D2",
#                    cluster_cols = F, filename = paste0(DIR6, "/pheatmap_week.pdf"))
# 
# cell_group_df = tibble::tibble(cell=row.names(colData(cds)), cell_group=colData(cds)$Phase)
# agg_mat = aggregate_gene_expression(cds, gene_module_df, cell_group_df)
# row.names(agg_mat) = stringr::str_c("Module ", row.names(agg_mat))
# pheatmap::pheatmap(agg_mat,
#                    scale="column", clustering_method="ward.D2", 
#                    filename = paste0(DIR6, "/pheatmap_Phase.pdf"))

cell_group_df = tibble::tibble(cell=row.names(colData(cds)), cell_group=cds@clusters$UMAP$clusters)
agg_mat = aggregate_gene_expression(cds, gene_module_df, cell_group_df)
row.names(agg_mat) = stringr::str_c("Module ", row.names(agg_mat))
pheatmap::pheatmap(agg_mat,
                   scale="column", clustering_method="ward.D2",
                   filename = paste0(DIR6, "/pheatmap_UMAP_clusters.pdf"))






#use ComplexHeatmap
#library(devtools)
#install_github("jokergoo/ComplexHeatmap")
library(ComplexHeatmap)
library(RColorBrewer)
COL <- colorRampPalette(c("#0068b7","white","magenta"))(n=31)
COLP <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)
ht1 <- Heatmap(t(scale(as.matrix(agg_mat))),
               clustering_distance_rows = "pearson",
               clustering_distance_columns = "pearson",
               clustering_method_rows = "ward.D2", 
               clustering_method_columns = "ward.D2", 
               col = COLP,
               heatmap_legend_param = list(title="Scaled\nexpression"),
               column_title = "co-regulated gene modules",
               column_title_side = "bottom",
               row_dend_reorder = TRUE,
               column_dend_reorder = TRUE
)
ht3 <- Heatmap(t(as.matrix(agg_mat)),
               clustering_distance_rows = "pearson",
               clustering_distance_columns = "pearson",
               clustering_method_rows = "ward.D2", 
               clustering_method_columns = "ward.D2", 
               col = COLP,
               heatmap_legend_param = list(title="Scaled\nexpression"),
               column_title = "co-regulated gene modules",
               column_title_side = "bottom",
               row_dend_reorder = TRUE,
               column_dend_reorder = TRUE
)
COLo <- colorRampPalette(c("white", "green4"))(n=15)#(c("#0068b7","white","magenta"))(n=31)
clu.cells <- as.matrix(POP3[,-1]) * 100
rownames(clu.cells) <- 1:nrow(clu.cells)
ht2 <- Heatmap(clu.cells,
               col =COLo,
               cluster_columns = F,
               heatmap_legend_param = list(title="Frequency (%)"), 
               row_title = "Cluster", 
               column_title = "week",
               column_title_side = "bottom",
               row_title_side = "right"#, border = "gray"
)

pdf(paste0(DIR6,"/ComlexHeatmap.pdf"), w=16, h=4)
ht1 + ht2
dev.off()

pdf(paste0(DIR6,"/ComlexHeatmap_notscaled.pdf"), w=16, h=4)
ht3 + ht2
dev.off()



Module_Matrix <- t(scale(as.matrix(agg_mat)))
cor(Module_Matrix) %>% hist()
MAT <- cor(Module_Matrix)
Module_Matrix_temp <- Module_Matrix
MAT_temp <- cor(Module_Matrix_temp)
matmax <- MAT[upper.tri(MAT)] %>% max
while (matmax > 0.8) {
  toRemove <- which(MAT_temp==matmax, arr.ind = TRUE)
  w <- which.max(c(sum(gene_module_df2$module == toRemove[1,1]),
              sum(gene_module_df2$module == toRemove[1,2]))
  )
  MAT_temp <- MAT_temp[-toRemove[1,w], -toRemove[1,w]]
  matmax <- MAT_temp[upper.tri(MAT_temp)] %>% max
}
pat <- rownames(MAT_temp) %>% str_split(., pattern = "\ ", simplify = T)
pat <- pat[,2] %>% as.numeric

ht4 <- Heatmap(Module_Matrix[,pat],
               clustering_distance_rows = "pearson",
               clustering_distance_columns = "pearson",
               clustering_method_rows = "ward.D2", 
               clustering_method_columns = "ward.D2", 
               col = COLP,
               heatmap_legend_param = list(title="Scaled\nexpression"),
               column_title = "co-regulated gene modules",
               column_title_side = "bottom",
               row_dend_reorder = TRUE,
               column_dend_reorder = TRUE
)


pdf(paste0(DIR6,"/ComlexHeatmap_reduced.pdf"), w=16, h=4)
ht4 + ht2
dev.off()

#pca ランダムフォレストによる次元削減
# install.packages("caret")
# library(caret)
# pca <- preProcess(data_train, method = "pca", pcaComp = 10)



qval <- 0.05 #0.0001
#Cluster3と2のDEG_monocle
res <<- res
ng1 <<- ng1
func_calc_DEG <- function(cds, x, qval) {
  cds_subset = cds[, clusters(cds) %in% x]
  subset_pr_test_res = graph_test(cds_subset, cores=4, neighbor_graph	= "principal_graph")
  pr_deg_ids = row.names(subset(subset_pr_test_res, q_value < qval))
  gene_module_df = monocle3:::find_gene_modules(cds_subset[pr_deg_ids,], resolution=c(0,10^seq(-6,-1)))#, partition_qval = qval)
  cell_group_df = tibble::tibble(cell=row.names(colData(cds_subset)), 
                                 cell_group=clusters(cds_subset)[colnames(cds_subset)])
  agg_mat = aggregate_gene_expression(cds_subset, gene_module_df, cell_group_df)
  module_dendro = hclust(dist(agg_mat))
  gene_module_df$module <- factor(gene_module_df$module, levels = row.names(agg_mat)[module_dendro$order])
  gene_module_df2 <- gene_module_df
  colnames(gene_module_df2)[1] <- "ensembl_gene_id"
  results.merge <- merge(gene_module_df2, res, by = "ensembl_gene_id", incomparables=NA)
  return(list(results.merge, cds_subset, gene_module_df))
}

#x <- c(1,4)
#x <- c(3,4)
#x <- c(3,8)
#x <- c(2,9)
#x <- c(2,4,9)
#x <- c(2,8)
x <- c(7,8,10)
x <- c(3:6)
x <- c(2, 6)
test1 <- func_calc_DEG(cds, c(7,8,10), qval)
test2 <- func_calc_DEG(cds, c(3:6), qval)
fwrite(test1[[1]], file = paste0(DIR6, "/subset",  paste(c(7,8,10), collapse="_"), ".txt"), sep = "\t")
fwrite(test2[[1]], file = paste0(DIR6, "/subset",  paste(c(3:6), collapse="_"), ".txt"), sep = "\t")

M <- max(as.numeric(test1[[3]]$module))

pdf(paste0(DIR6, "/subset", paste(c(7,8,10), collapse="_"), ".pdf"), w=10, h=8)
plot_cells(test1[[2]],
           genes=test1[[3]],
           label_branch_points = F,
           label_leaves = F,
           label_roots = F,
           cell_size = 2
) + ng1 + theme(axis.text.x = element_blank(), axis.text.y =  element_blank())
dev.off()

# TEST <- plot_cells_3d(test[[2]],
#            genes=test[[3]],
# #           label_branch_points = F,
# #           label_leaves = F,
# #           label_roots = F,
#            cell_size = 2
# )


pdf(paste0(DIR6, "/subset_all.pdf"), w=60, h=60)
plot_cells(cds,
           genes=gene_module_df,
           label_branch_points = F,
           label_leaves = F,
           label_roots = F,
           cell_size = 4
) + ng1
dev.off()



#Cluster 7と10のDEG_DEsingle
library(DEsingle)
rawcount <- cds[, clusters(cds) %in% c(7, 10)]@assays$data$counts
group <- factor(as.character(clusters(cds[, clusters(cds) %in% c(7, 10)])))
results <- DEsingle(counts = rawcount, group = group, parallel = T) #takes very long time about 20 sec/gene
results.classified <- DEtype(results = results, threshold = 0.05)
results.classified2 <- data.frame(results.classified, ensembl_gene_id=rownames(results.classified))
results.merge <- merge(results.classified2, res, by = "ensembl_gene_id", incomparables=NA)
NAME <-  "10_vs_7"
write.table(results.merge, file = paste0(DIR6, "/", NAME, "_classified_DEsingle.txt"), sep = "\t", col.names = NA, quote = F)
saveRDS(results.merge, file = paste0(DIR6, "/", NAME, "_DEsingle.rds"))

#Cluster 3と2のDEG_DEsingle
rawcount <- cds_subset32@assays$data$counts
group <- factor(as.character(clusters(cds_subset32)))
results <- DEsingle(counts = rawcount, group = group, parallel = T) #takes very long time about 20 sec/gene
results.classified <- DEtype(results = results, threshold = 0.05)
results.classified2 <- data.frame(results.classified, ensembl_gene_id=rownames(results.classified))
results.merge <- merge(results.classified2, res, by = "ensembl_gene_id", incomparables=NA)
NAME <-  "3_vs_2"
write.table(results.merge, file = paste0(DIR6, "/", NAME, "_classified_DEsingle.txt"), sep = "\t", col.names = NA, quote = F)
saveRDS(results.merge, file = paste0(DIR6, "/", NAME, "_DEsingle.rds"))


#Cluster 5と10のDEG_DEsingle
library(DEsingle)
rawcount <- cds[, clusters(cds) %in% c(5, 10)]@assays$data$counts
group <- factor(as.character(clusters(cds[, clusters(cds) %in% c(5, 10)])))
results <- DEsingle(counts = rawcount, group = group, parallel = T) #takes very long time about 20 sec/gene
results.classified <- DEtype(results = results, threshold = 0.05)
results.classified2 <- data.frame(results.classified, ensembl_gene_id=rownames(results.classified))
results.merge <- merge(results.classified2, res, by = "ensembl_gene_id", incomparables=NA)
NAME <-  "10_vs_5"
write.table(results.merge, file = paste0(DIR6, "/", NAME, "_classified_DEsingle.txt"), sep = "\t", col.names = NA, quote = F)
saveRDS(results.merge, file = paste0(DIR6, "/", NAME, "_DEsingle.rds"))

#############190801ここから


# 
# 
# install_github("wjawaid/enrichR")
# library(enrichR)
# 
# dbs <- enrichR::listEnrichrDbs()
# if (is.null(dbs)) websiteLive <- FALSE
# if (websiteLive) head(dbs)
# 
# dbs <- c("GO_Molecular_Function_2015", "GO_Cellular_Component_2015", "GO_Biological_Process_2015")
# if (websiteLive) enriched <- enrichr(c("Runx1", "Gfi1", "Gfi1b", "Spi1", "Gata1", "Kdr"), databases = "GO_Biological_Process_2018")
# if (websiteLive) enriched[["GO_Biological_Process_2015"]]
# enrichr(results.merge$ensembl_gene_id)
# 
