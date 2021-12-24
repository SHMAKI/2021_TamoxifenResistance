library(tidyverse)
library(Seurat)
library(rtracklayer)

# import count data
datT <- t(read.csv("data/DrMagi_F0551_GA06071-4.genematrix.csv", header=T, stringsAsFactors=F, row.names = 1))
# import annotation data
annot <- read.csv("data/DrMagi_F0551_GA06071-4.report.csv", header=T, stringsAsFactors=F, skip = 8)
# delete final row "Union"
annot <- annot[annot$Barcode != "UNION",]

# cell labeling to datT
# sum(colnames(datT) ==  annot$Barcode)
Type <- paste0(str_split(annot$Type, "_",simplify = T)[,1], str_split(annot$Type, "_",simplify = T)[,2])
COL <- paste0(Type[match(colnames(datT),annot$Barcode)], "_", seq(1, nrow(annot),1))
colnames(datT) <- COL

# convert gene id: ENSG -> Entrez Gene ID
gtf <- rtracklayer::import('/Users/s_magi/ReferenceGenome/Homo_sapiens.GRCh38.81.gtf') #use release81
match(rownames(datT), as.data.frame(gtf)$gene_id) %>% is.na() %>% sum()
gtf_Ens_gene <- subset(as.data.frame(gtf), as.data.frame(gtf)$type == "gene")
gtf_Ens_gene <- gtf_Ens_gene[,c("gene_id","gene_name","gene_biotype")]

# Set up object CPM > x in more than y % of total cells
x <- 5
y <- 5
nCount_all <- apply(datT, 2, sum)
merged_CPM <- datT
for(i in 1:length(nCount_all)) {
  merged_CPM[,i] <- datT[,i] * 1E6/ nCount_all[i]
}
count_merged_CPM <- apply(merged_CPM, 1, function(t) {sum(t>x)/length(t)*100 > y})
rm(merged_CPM)
# sum(count_merged_CPM)/ length(count_merged_CPM) of total cells were kept for analysis.
datT_filtered <- datT[count_merged_CPM,]
df_filteredgenes <- gtf_Ens_gene[match(rownames(datT_filtered), gtf_Ens_gene$gene_id),]
# saveRDS(df_filteredgenes, "processed_data/df_filteredgenes.rds")
# if it is true, gene ids in the data are completely included in gtf_Ens_gene.
# sort(match(gtf_Ens_gene$gene_id, rownames(datT_filtered))) %>% unique() %>% length == nrow(datT_filtered)

# cell cycle regression https://satijalab.org/seurat/cell_cycle_vignette.html
# a list of cell cycle markers, from Tirosh et al, 2015
cc.genes <- readLines(con = "data/regev_lab_cell_cycle_genes.txt")
# We can segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes[1:43]
s.genes[s.genes == "MLF1IP"] <- "CENPU"
g2m.genes <- cc.genes[44:97]
s.genes_ensembl <- gtf_Ens_gene[match(s.genes, gtf_Ens_gene$gene_name),]$gene_id %>% na.omit()
g2m.genes_ensembl <- gtf_Ens_gene[match(g2m.genes, gtf_Ens_gene$gene_name),]$gene_id %>% na.omit()

# create Seurat object
Seu <- CreateSeuratObject(counts = datT_filtered, project = "TamR", min.cells = 5)
# Visualize QC metrics as a violin plot
mitogenes.id <- df_filteredgenes[grep("^MT-", df_filteredgenes$gene_name),]$gene_id
gtf_Ens_gene$gene_name
mitogenes.id.pattern <- paste(mitogenes.id,  collapse = "|")
rps.id <- df_filteredgenes[grep("^RPS[1-9]", df_filteredgenes$gene_name),]$gene_id
rpl.id <- df_filteredgenes[grep("^RPL[1-9]", df_filteredgenes$gene_name),]$gene_id
rps.id.pattern <- paste(rps.id,  collapse = "|")
rpl.id.pattern <- paste(rpl.id,  collapse = "|")
rp.id.pattern <- paste(rps.id.pattern, rpl.id.pattern,  sep = "|")
Seu[["percent.mt"]] <- PercentageFeatureSet(Seu, pattern = mitogenes.id.pattern)
Seu[["percent.rp"]] <- PercentageFeatureSet(Seu, pattern = rp.id.pattern)
VlnPlot(Seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rp"), ncol = 4);ggsave(paste0("intermediate_figs/Vlnplot.pdf"), width = 12, h=6)
plot1 <- FeatureScatter(Seu, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Seu, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(Seu, feature1 = "percent.mt", feature2 = "nFeature_RNA")
plot4 <- FeatureScatter(Seu, feature1 = "percent.mt", feature2 = "percent.rp")
plot5 <- FeatureScatter(Seu, feature1 = "percent.rp", feature2 = "nFeature_RNA")
plot1;ggsave(paste0("intermediate_figs/QCplot1.pdf"))
plot2;ggsave(paste0("intermediate_figs/QCplot2.pdf"))
plot3;ggsave(paste0("intermediate_figs/QCplot3.pdf"))
plot4;ggsave(paste0("intermediate_figs/QCplot4.pdf"))
plot5;ggsave(paste0("intermediate_figs/QCplot5.pdf"))

# data filtering
Feature.min <- 1500
Feature.max <- Inf
Mito.max <- 25
Seu <- subset(Seu, subset = nFeature_RNA >= Feature.min & nFeature_RNA < Feature.max & percent.mt < Mito.max)
Seu <- NormalizeData(Seu, verbose = FALSE)
Seu <- FindVariableFeatures(Seu, selection.method = "vst")
Seu <- ScaleData(Seu, features = rownames(Seu))
Seu <- RunPCA(Seu, features = VariableFeatures(Seu), ndims.print = 1:4, nfeatures.print = 10)

# saveRDS(Seu, paste0("processed_data/Seu.rds"))

# cell cycle regression
Seu_allreg <- CellCycleScoring(Seu, s.features = s.genes_ensembl, g2m.features = g2m.genes_ensembl, set.ident = TRUE)
# saveRDS(Seu_allreg, paste0("processed_data/Seu_allreg.rds"))
