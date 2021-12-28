# before staring R, putting a .Renviron file containing the line RETICULATE_PYTHON="/path/to/python3" into the user's home directory
# https://stackoverflow.com/questions/50145643/unable-to-change-python-path-in-reticulate
library(reticulate) # for UMAP clustering
# py_config() 
library(tidyverse)
library(Seurat)
# monocle 3 alpha version
# https://cole-trapnell-lab.github.io/monocle3/monocle3_docs/#constructing-single-cell-trajectories
# devtools::install_github('cole-trapnell-lab/monocle3')
library(monocle3)
library(ggsci)
library(scales)
library(ggpubr)
library("ComplexHeatmap")

source("src/theme.R")

# import Seurat obj
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
g <- g + stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", bracket.size = 0.5, step.increase = 0.03, tip.length = 0.01)
g <- g + coord_cartesian(ylim=c(0.8, 1.1))

# ggsave(file =paste0("Fig/Fig2b.pdf"), plot=g, width=4, height=4)

#--------create monocle obj-----
gene_annotation <- df_filteredgenes
colnames(gene_annotation)[colnames(gene_annotation) == "gene_name"] <- "gene_short_name"
rownames(gene_annotation) <- rownames(umi_matrix)

cds <- new_cell_data_set(as(umi_matrix, "sparseMatrix"),
                         cell_metadata = sample_info,
                         gene_metadata = gene_annotation)

# Pre-process the data
cds <- preprocess_cds(cds, num_dim = 100, residual_model_formula_str = "~ S.Score + G2M.Score")

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

# ggsave(g_horiz, filename = paste0("Fig/Fig2c.pdf"), w=8, h=4)


#--------Fig.2def UMAP and clustering--------
# output is slightly different from the paper because we didn't use 'umap.fast_sgd = FALSE' and 'cores = 1'
cds <- reduce_dimension(cds, max_components = 3, cores = 18, reduction_method = "UMAP")

# cluster cells
cds <- cluster_cells(cds, partition_qval = 5e-2, reduction_method = "UMAP", k = 25, resolution=5e-2, verbose = T, random_seed = 0)

#check how many partitions there are (I had two partitions):
levels(cds@clusters$UMAP$partitions)
# set them all to "1"
cds@clusters$UMAP$partitions[cds@clusters$UMAP$partitions == "2"] <- "1"

# 3d UMAP plot in Fig. S5
# plot_cells_3d(cds, color_cells_by="week")#, color_cells_by="partition", group_cells_by="partition")#, show_trajectory_graph = F)# + scale_color_nejm()

fig2_d <- plot_cells(cds, color_cells_by = "Phase", show_trajectory_graph = F, group_label_size = 0) + theme_bw() + ng1 + scale_color_manual(values = c("skyblue", "#CC6677", "darkgreen")) + guides(color = guide_legend(title = "cell cycle", override.aes = list(size = 4)))
# ggsave(fig2_d, filename = paste0("Fig/Fig2d.pdf"),w=5,h=4)

fig2_e <- plot_cells(cds, color_cells_by = "week", show_trajectory_graph = F, group_label_size = 0) + theme_bw() + ng1 + scale_color_nejm()
# ggsave(fig2_e, filename = paste0("Fig/Fig2e.pdf"),w=5,h=4)

fig2_f <- plot_cells(cds, color_cells_by = "cluster", show_trajectory_graph = F, group_label_size = 8) + theme_bw() + ng1 + scale_color_npg() + guides(color = guide_legend(title = "subgroup", override.aes = list(size = 4)))
# ggsave(fig2_f, filename = paste0("Fig/Fig2f.pdf"),w=5.5,h=4)

#-----------Fig.2g module analysis-------
pr_graph_test_res <- graph_test(cds,  neighbor_graph="knn", expression_family = "quasipoisson", verbose = T, cores = 18, k=50)
pr_deg_ids <- row.names(subset(pr_graph_test_res, q_value < 0.05 & morans_I > 0.05))

# We will learn more about graph_test() in the differential expression analysis section later.
# We can take all the genes that vary across this set of cells and group those that have similar patterns of expression into modules:
gene_module_df <- find_gene_modules(cds[pr_deg_ids,], max_components = 3, umap.fast_sgd	= F, umap.min_dist = .1,
                                    verbose = T, cores = 18, k=75, random_seed = 0, louvain_iter=5)
gene_module_df %>% group_by(module) %>% summarize(n=n()) %>% pull(n) %>% min() 
# colnames(gene_module_df2)[1] <- "ensembl_gene_id"
# fwrite(gene_module_df, paste0("intermediate_figs/gene_module_df.txt"), sep = "\t")

cell_group_df = tibble::tibble(cell=row.names(colData(cds)), cell_group=cds@clusters$UMAP$clusters)
agg_mat = aggregate_gene_expression(cds, gene_module_df, cell_group_df)

# create heatmap
ht1 <- Heatmap(t(scale(as.matrix(agg_mat))), row_names_gp = (gpar(fontsize = 18)), column_names_gp = (gpar(fontsize = 18)),
               clustering_distance_columns =  "euclidean",
               clustering_method_columns = "ward.D2",
               col = COL,
               heatmap_legend_param = list(title="Scaled\nexpression"),
               column_title = "co-regulated gene modules", 
               column_title_side = "bottom",
               column_names_rot = 0, border = "gray60",
)

cell_group_df_for_rate = tibble::tibble(cell=row.names(colData(cds)), clusters=cds@clusters$UMAP$clusters, week=colData(cds)$week, Phase = colData(cds)$Phase)
POP <- cell_group_df_for_rate %>% group_by(week, clusters, .drop = F) %>% summarise (n = n()) %>% mutate(freq = n / sum(n))
POP3 <- POP %>% dplyr::select(clusters, week, freq) %>% tidyr::spread(week, freq)
COLo <- colorRampPalette(c("white", "green4"))(n=15)
clu.cells <- as.matrix(POP3[,-1]) * 100
rownames(clu.cells) <- 1:nrow(clu.cells)
ht2 <- Heatmap(clu.cells, row_names_gp = (gpar(fontsize = 18)), column_names_gp = (gpar(fontsize = 18)),
               col =COLo,
               cluster_columns = F,
               cluster_rows = F,
               heatmap_legend_param = list(title="Frequency (%)"), 
               column_title = "week", column_names_rot = 0,
               column_title_side = "bottom", border = "gray60",
               row_title = "cell subgroups", row_title_side = "left", row_names_side = "left", row_order = c(1, 6, 2, 3, 4, 5),
)

# pdf(paste0("Fig/Fig2g_ComlexHeatmap.pdf"), w=7, h=3)
# ht2 + ht1
# dev.off()

#--------Fig.3abcde pseudotime--------
# learn graph
cds = learn_graph(cds, use_partition = F, close_loop = T, verbose = T, list(minimal_branch_len = 12,
                                                                            euclidean_distance_ratio = 2.5,
                                                                            geodesic_distance_ratio = 1/2,
                                                                            orthogonal_proj_tip = F))
# #choose the root and calc pseudotime
cds = order_cells(cds, reduction_method = "UMAP", verbose = T) #manual click w0 -> pseudo time = 0
plot_cells(cds)

fig3_a <- plot_cells(cds, color_cells_by = "cluster", show_trajectory_graph = T, group_label_size = 0, label_leaves=FALSE) + theme_bw() + ng1 + scale_color_npg() + guides(color = guide_legend(title = "subgroup", override.aes = list(size = 4)))
# ggsave(fig3_a, filename = paste0("Fig/Fig3a.pdf"),w=5.5,h=4)

fig3_b <- plot_cells(cds, color_cells_by = "week", show_trajectory_graph = T, group_label_size = 0, label_leaves=FALSE) + theme_bw() + ng1 + scale_color_nejm()
# ggsave(fig3_b, filename = paste0("Fig/Fig3b.pdf"),w=5,h=4)

fig3_c <- plot_cells(cds, color_cells_by = "pseudotime", label_branch_points=TRUE, label_leaves=FALSE) + theme_bw() + ng1
# ggsave(fig3_c, filename = paste0("Fig/Fig3c.pdf"),w=5.5,h=4)

# df pseudotime - cluster - week
Pseudotime.dat <- data.frame(cds@colData, pseudotime = cds@principal_graph_aux$UMAP$pseudotime, Cluster = cds@clusters$UMAP$clusters)

g <- ggplot(Pseudotime.dat, aes(x=week, y= pseudotime)) +  geom_violin(fill="white", scale = "width") +
  geom_boxplot(width=.2, fill="lightgray",outlier.fill = NA, outlier.color = NA, color="black") + 
  theme_bw() + scale_fill_aaas() + ng1 + ylab("pseudotime")+ xlab("week")+ylim(c(0,16))
my_comparisons <- list( c("0", "3"), c("3", "6"), c("6", "9") ) # Add pairwise comparisons p-value
g <- g + stat_compare_means(comparisons = my_comparisons)
ggsave(file =paste0("Fig/Fig3d.pdf"), plot=g, width=4, height=4)

g <- ggplot(Pseudotime.dat, aes(x=Cluster, y= pseudotime)) + geom_violin(fill="white", scale = "width") +
  geom_boxplot(width=.2, fill="lightgray",outlier.fill = NA, outlier.color = NA, color="black") + 
  theme_bw() + ng1 + ylab("pseudotime")+ xlab("subgroup")
ggsave(file =paste0("Fig/Fig3e.pdf"), plot=g, width=4, height=4)


#--------Fig.3g pseudotime-2345-tree-------
cds_cluster2345 <- cds[, cds@clusters$UMAP$clusters %in% c(2, 3, 4, 5)]

subset_pr_test_res <- graph_test(cds_cluster2345, neighbor_graph="principal_graph", cores=18, verbose = T)
pr_deg_ids <- row.names(subset(subset_pr_test_res, q_value < 0.05 & morans_I > 0.05)) 
gene_module_df_tree <- find_gene_modules(cds_cluster2345[pr_deg_ids,],  verbose = T, cores = 18, k=25, max_components = 3, resolution = 0.03)
agg_mat_tree <- aggregate_gene_expression(cds_cluster2345, gene_module_df_tree)
module_dendro_tree <- hclust(dist(agg_mat_tree))
p <- plot_cells(cds_cluster2345,
                genes=gene_module_df_tree, cell_size = 1,
                label_cell_groups=TRUE,
                show_trajectory_graph=FALSE, group_label_size = 8) + theme_bw() + ng1
ggsave(p, filename = paste0("Fig/Fig3g.pdf"), width = 12, height = 7)
