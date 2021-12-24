library(tidyverse)
library(DESeq2)
library(ggpubr) #ggmaplot
library(pheatmap) #pheatmap
library(rtracklayer) #readGFF
library(amap) # Dist
library(ggdendro) #ggdendrogram

source("src/theme.R")

# count files
tbl_ensembl_pe <- read.delim("data/counts_ensembl_from_pe.tsv", skip=1, row.names=1)
tbl_ensembl_sr <- read.delim("data/counts_ensembl_from_sr.tsv", skip=1, row.names=1)
# check all gene IDs are the same
sum(rownames(tbl_ensembl_pe) == rownames(tbl_ensembl_sr)) ==  nrow(tbl_ensembl_sr)

tbl_en <- cbind(tbl_ensembl_pe, tbl_ensembl_sr[, grep("Sample", colnames(tbl_ensembl_sr))])
tbl_en_annot <- tbl_en[, -grep("Sample", colnames(tbl_en))]
tbl_en_dat <- tbl_en[, grep("Sample", colnames(tbl_en))]

# import gtf file which can be downloaded from ensembl.org
gtf_Ens <- readGFF("~/ReferenceGenome/Homo_sapiens.GRCh38.96.gtf")
gtf_Ens_gene <- subset(gtf_Ens, gtf_Ens$type == "gene")
gtf_Ens_gene_1 <- gtf_Ens_gene[,c("gene_id","gene_name","gene_biotype")]

# match samplename with experimental condition
sampleCD <- read.csv("data/samplecondition.csv")
samplename <- colnames(tbl_en_dat)

MATCH <- rep(NA, nrow(sampleCD))
MATCHname <- rep(NA, nrow(sampleCD))
for(i in 1:nrow(sampleCD)) {
  if(length(grep(as.character(sampleCD$Sample.ID)[i], samplename) != 0)) {
    MATCH[i] <- grep(as.character(sampleCD$Sample.ID)[i], samplename)
    MATCHname[i] <- grep(as.character(sampleCD$Sample.ID)[i], samplename, value = T)
  }
}
# check duplicates or losses
sum(!is.na(sampleCD2$MATCH)) == length(samplename)

sampleCD2 <- cbind(sampleCD, MATCH, MATCHname)
sampleCD3 <- na.omit(sampleCD2)
sampleCD3 <- sampleCD3[order(sampleCD3$MATCH),]

# thresholds
TH <- 0.001
LFCTH <- 0.5
FC <- LFCTH

FCtable <- matrix(data = 0, nrow = nrow(tbl_en_dat), ncol = 13)
padjtable <- matrix(data = 1, nrow = nrow(tbl_en_dat), ncol = 12)
baseMeantable <- matrix(data = 0, nrow = nrow(tbl_en_dat), ncol = 12)
FCSEtable <- matrix(data = 0, nrow = nrow(tbl_en_dat), ncol = 12)
res_list <- vector("list", 12)

#Deseq2 in each time points
for(i in 1:12) {
  SUB <- subset(sampleCD3, week == i & treatment != "H")
  featureX <- tbl_en_dat[, SUB$MATCH]
  groupC <- data.frame(con = factor(SUB$treatment))
  dds <- DESeqDataSetFromMatrix(countData = featureX, colData = groupC, design = ~ con)
  dds$con <- relevel(dds$con, "C")
  dds <- DESeq(dds)
  res <- results(dds)
  res_list[i] <- res
    
  FCtable[,i+1] <- res$log2FoldChange
  padjtable[,i] <- res$padj
  baseMeantable[,i] <- res$baseMean
  FCSEtable[,i] <- res$lfcSE
}

rownames(FCtable) <- rownames(tbl_en_dat)
rownames(padjtable) <- rownames(tbl_en_dat)
rownames(baseMeantable) <- rownames(tbl_en_dat)
rownames(FCSEtable) <- rownames(tbl_en_dat)
colnames(FCtable) <- paste0("W", 0:12)
colnames(padjtable) <- paste0("W", 1:12)
colnames(baseMeantable) <- paste0("W", 1:12)
colnames(FCSEtable) <- paste0("W", 1:12)

FCtableX <- data.frame(FCtable, gtf_Ens_gene_1[match(rownames(FCtable), gtf_Ens_gene_1$gene_id),])
padjtableX <- data.frame(padjtable, gtf_Ens_gene_1[match(rownames(padjtable), gtf_Ens_gene_1$gene_id),])
baseMeantableX <- data.frame(baseMeantable, gtf_Ens_gene_1[match(rownames(baseMeantable), gtf_Ens_gene_1$gene_id),])
FCSEtableX <- data.frame(FCSEtable, gtf_Ens_gene_1[match(rownames(FCSEtable), gtf_Ens_gene_1$gene_id),])

# fig 1f---------
# thresholds
FDRtest <- c(1E-3)
log2FCtest <- c(0.5)
nClu <- c(6)
AtLeast <- 2

COL <- colorRampPalette(c("#0068b7","white","magenta"))(n=31)

for (x in 1:length(log2FCtest)) {
  log2FC <- log2FCtest[x]
  #DIRx <- paste0("FC",log2FC)
  #dir.create(DIRx)
  maxFCweekTF <- abs(FCtable) > log2FC

  for (y in 1:length(FDRtest)) {
    FDR <- FDRtest[y]
    # DIRy <- paste0("FDR<",FDR)
    # DIRxy <- paste0(DIRx, "/", DIRy)
    # dir.create(DIRxy)
    maxFCpadjTF <- padjtable < FDR
    TFtable <- na.omit(as.data.frame(maxFCpadjTF * maxFCweekTF[,2:13]))
    SELECT <- apply(TFtable, 1, sum) > AtLeast
    names(SELECT)[SELECT]
    FCtable3 <- FCtable[rownames(FCtable) %in% names(SELECT)[SELECT],]
    padjtable3 <- FCtable[rownames(padjtable) %in% names(SELECT)[SELECT],]
    MATCH <- match(rownames(FCtable3), gtf_Ens_gene_1$gene_id)
    FCtable4 <- cbind(as.data.frame(FCtable3), gtf_Ens_gene_1[MATCH,])
    # write.csv(FCtable4, file=paste0(DIRxy,"/allDEGs.csv"))
    #extract protein-coding gemes
    FCtable5 <- na.omit(FCtable3[which(FCtable4$gene_biotype == "protein_coding"),])
    
    for (z in 1:length(nClu)) {
      Clu <- nClu[z]
      # DIRz <- paste0("clu", Clu)
      DIRxyz <- "Fig/Fig1f" #paste0(DIRx, "/", DIRy, "/", DIRz)
      dir.create(DIRxyz)
      #clustering
      scaleFCtable5 <- t(apply(FCtable5, 1, scale))
      OBJ <- pheatmap(FCtable5,
                      scale = "row",
                      color = COL,
                      cluster_cols = F,
                      cutree_rows = Clu,
                      clustering_distance_rows = "correlation", # "euclidean", 
                      clustering_method = "ward.D2", #"ward.D2"
                      cellwidth = 10, cellheight = 0.1, fontsize = 8, 
                      filename = paste(DIRxyz, "/heatmap_", Clu, ".pdf", sep="")
      )
      
      LargeE <- cbind(scaleFCtable5, cluster = cutree(OBJ$tree_row, k = Clu))
      MATCH <- match(rownames(scaleFCtable5), gtf_Ens_gene_1$gene_id)
      LargeF <- cbind(as.data.frame(LargeE), gtf_Ens_gene_1[MATCH,])
      write.csv(LargeF, file = paste(DIRxyz, "/scaledDATA_",Clu, ".csv", sep=""))
      
      for(w in 1:Clu) {
        smallE <- LargeF[which(LargeF$cluster == w),]
        smallE2 <- smallE[,1:13]
        DIR2 <- paste(DIRxyz, "/cluster_", w, sep="")
        dir.create(DIR2)
        pdf(paste(DIR2,"/cluster", w, "_linegraph.pdf",sep=""), height = 6, width = 6)
        matplot(t(smallE2), type="l", col = alpha("gray", 0.1), lty="dashed",axes=F,xlab="Week", ylab="gene expression (TAM/Ctrl, z-score)", cex=2)
        par(new=T)
        matlines(apply(smallE2, 2, median), type="l", col = "red",axes=F, xlab="", ylab="")
        axis(2) #add y-axis
        axis(side=1,at=1:ncol(smallE2),labels=c(0:12)) #add x-axis
        dev.off()
        
        tp <- 0:12
        df_tidy <- t(smallE2) %>% as.data.frame() %>% mutate(Time = tp) %>% 
          tidyr::gather(Gene, Z_score, -Time)
        g <- ggplot(df_tidy, aes(x=Time, y=Z_score)) + geom_line(color = "orange", alpha = 0.2, aes(group=Gene))
        g <- g + stat_summary(geom="line", fun.y = "median", color="red")
        g <- g + scale_x_continuous(breaks = tp) + xlab("Time (weeks)")  + theme_bw(base_size = 15) + ng1
        ggsave(g, filename = paste(DIR2,"/cluster", w, "_linegraph.pdf",sep=""), height = 4, width = 6)
        
      }
    }
  }
}

# fig1e and sfig3----------
DIST_p <- Dist(scale(t(FCtable5)), method = "correlation")
h <- hclust(DIST_p, method="ward.D2")
g <- ggdendrogram(h, segments = TRUE, labels = TRUE, 
                  leaf_labels = TRUE, rotate = 45, theme_dendro = TRUE) + theme(panel.border = element_blank())
# ggsave(filename = "Fig/Fig1e.pdf", g, width = 4, height=6)

Dist_vec <- c(as.matrix(DIST_p)[2,1],
              as.matrix(DIST_p)[3,2],
              as.matrix(DIST_p)[4,3],
              as.matrix(DIST_p)[5,4],
              as.matrix(DIST_p)[6,5],
              as.matrix(DIST_p)[7,6],
              as.matrix(DIST_p)[8,7],
              as.matrix(DIST_p)[9,8],
              as.matrix(DIST_p)[10,9],
              as.matrix(DIST_p)[11,10],
              as.matrix(DIST_p)[12,11],
              as.matrix(DIST_p)[13,12]
              )

# ### growth_rate from Fig.1B
df_growth <- read.csv("../../Figure_1_bcd/data/growthdata.csv")
df_growth_TAM_mean <- df_growth %>% dplyr::filter(Condition=="TAM") %>% group_by(Week) %>% summarize(mean=mean(Value))
abs_derivative_grate <- abs(df_growth_TAM_mean$mean[1:10] - df_growth_TAM_mean$mean[2:11])

df <- data.frame(week = 1:12,
                 growth_derivative = c(abs_derivative_grate, NA, NA),
                 distance = Dist_vec)

p1 <- ggplot(df, aes(x=week)) +
  geom_line(aes(y=growth_derivative,color="growth rate"), size=1) + geom_point(aes(y=growth_derivative,color="growth rate"), size=2) +
  geom_line(aes(y=distance, color="gene expression"), size=1) + geom_point(aes(y=distance, color="gene expression"), size=2) +
  theme_bw() + scale_x_continuous(breaks = 1:12) + ng1
p1 <- p1 + scale_y_continuous(name = "Absolute differences of\ngrowth rate from a week before",
                              sec.axis = sec_axis(~.*second_rate, name = "Distance from a week before"))
p1 <- p1 + theme(legend.position = c(.75,.8), legend.direction = "vertical")
# ggsave(p1, filename = "Fig/FigS3.pdf", w=8,h=4)
