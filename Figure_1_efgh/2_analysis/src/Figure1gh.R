library(tidyverse)
library(tidyheatmap) #devtools::install_github("jbengler/tidyheatmap")

DIRtgmine <- "data/EAfromTargetimine"
files <- list.files(DIRtgmine, pattern = ".tsv$", recursive = T)

list_kegg <- list()
list_reactome <- list()
v <- 1
w <- 1
for(i in 1:length(files)) {
  module <- str_split(string = files[i], pattern = " ", simplify = T)[1]
  cluster <- str_split(string = module, pattern = "/")[[1]] %>% tail(n=1)
  tmp_file <- read_tsv(paste0(DIRtgmine, "/", files[i]), col_names = F) %>% mutate(cluster = cluster)
  type <- str_split(string = files[i], pattern = " ")[[1]] %>% tail(n=1)
  
  if(type == "KEGG.tsv") {
    list_kegg[[v]] <- tmp_file
    v <- v + 1
  } else if(type == "Reactome.tsv") {
    list_reactome[[w]] <- tmp_file
    w <- w + 1
  }
}

df_kegg <- do.call("rbind", list_kegg)
df_reactome <- do.call("rbind", list_reactome)

# define graph func
DrawEnrichmentHeatmap <- function(temp_df, top_th=5, slim_th=.7, filename="") {
  temp_df2 <- temp_df %>% mutate(terms_id =paste0(X1, " (",X4,")"))
  temp_df2_topsig <- temp_df2 %>% arrange(cluster) %>% group_by(cluster) %>% dplyr::filter(X2 < 1)

  # remove redundant terms
  list_slim <- list()
  x <- 1
  for(tmp_module in unique(temp_df2_topsig$cluster)) {
    tmp_df <- temp_df2_topsig %>% dplyr::filter(cluster == tmp_module)
    tmp_vec1 <- str_split(tmp_df %>% pull(X3), pattern = ",", simplify = T)
    
    rem_id_term <- c()
    rem_num <- c()
    if (nrow(tmp_vec1)<2){
      list_slim[[x]] <- tmp_df
      x <- x + 1
      next
    }
    for (i in 1:(nrow(tmp_vec1)-1)) {
      if (i %in% rem_num) {next}
      for (j in (i+1):nrow(tmp_vec1)) {
        if (j %in% rem_num) {next}
        vi <- tmp_vec1[i,][tmp_vec1[i,] != ""]
        vj <- tmp_vec1[j,][tmp_vec1[j,] != ""]
        if (sum(vj %in% vi) / length(vj) > slim_th | sum(vi %in% vj) / length(vi) > slim_th) { #どちらかがどちらかに一定割合以上含まれる場合
          wmax = c(i,j)
          rem_num_tmp <- wmax[which.max(tmp_df$X2[wmax])]
          rem_id_term <- c(rem_id_term, tmp_df$terms_id[rem_num_tmp])
          rem_num <- c(rem_num, rem_num_tmp)
          if (rem_num_tmp == i) {break}
        }
      }
    }
    tmp_df_filter <- tmp_df %>% dplyr::filter(!(terms_id %in% rem_id_term)) %>% dplyr::top_n(n = top_th, -X2)
    list_slim[[x]] <- tmp_df_filter
    x <- x + 1
  }
  
  temp_df2_slim <- do.call("rbind", list_slim)
  use_id_term <- temp_df2_slim %>% dplyr::pull(terms_id) %>% unique()
  temp_df3 <- temp_df2_topsig %>% dplyr::select(-X3) %>% spread(key= cluster,value = X2) %>% dplyr::filter(terms_id %in% use_id_term)
  
  temp_df3$terms_id <- factor(temp_df3$terms_id, levels = use_id_term)
  temp_df3 %>% arrange(terms_id)
  
  temp_df4 <- temp_df3 %>% replace(is.na(.),1) %>%
    tidyr::gather(key = cluster, value = qvalue, colnames(temp_df3)[!(colnames(temp_df3) %in% c("X1","X4","terms_id"))])
  colnames(temp_df4)[1:2] <- c("terms", "id")
  temp_df4 <- temp_df4 %>% mutate(log10qval=-log10(as.numeric(qvalue)))
  temp_df4$cluster <- factor(temp_df4$cluster)
  temp_df4$terms_id <- factor(temp_df4$terms_id, levels=use_id_term)
  #draw_heatmap
  temp_df4 %>% 
    arrange(terms_id, cluster) %>% tidy_heatmap(., border_color = "grey60",
                                                rows = terms_id,
                                                columns = cluster,
                                                values = log10qval,
                                                scale = "none",
                                                colors = c("#ffffff","#390655"),
                                                color_legend_min = 0,
                                                color_legend_max = 4,
                                                cellheight = 16, cellwidth = 16, fontsize = 14, angle_col = 0,
                                                #filename = paste0(DIRtgmine, "/tgmine_df_reactome", top_th, "_slim", slim_th, ".pdf")
                                                filename = filename
    )
}

DrawEnrichmentHeatmap(df_reactome, top_th=5, slim_th=.7, filename="Fig/Fig1g.pdf")
DrawEnrichmentHeatmap(df_kegg, top_th=5, slim_th=.7, filename="Fig/Fig1h.pdf")
