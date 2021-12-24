library(tidyverse)
library(tidyheatmap) #devtools::install_github("jbengler/tidyheatmap")

DIRtgmine <- "data/EAfromTargetimine/"
module_level<- c(5,9,11,8,10,3,7,6,4,2,12,1)
files <- list.files("data/EAfromTargetimine", pattern = ".tsv$")

list_kegg <- list()
list_reactome <- list()
v <- 1
w <- 1
for(i in 1:length(files)) {
  module <- str_split(string = files[i], pattern = " ", simplify = T)[1]
  moduleNo <-  str_extract_all(module, "[0-9.]+") %>% as.numeric
  tmp_file <- read_tsv(paste0(DIRtgmine, "/", files[i]), col_names = F) %>% mutate(module = moduleNo)
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

top_th <- 10
slim_th <- .7

# define graph func
DrawEnrichmentHeatmap_module <- function(temp_df, top_th=5, slim_th=.7, module_level,filename="") {
  temp_df2 <- temp_df %>% mutate(terms_id =paste0(X1, " (",X4,")"))
  temp_df2$module <- factor(temp_df2$module, levels=module_level)
  temp_df2_topsig <- temp_df2 %>% arrange(module) %>% group_by(module) %>% dplyr::filter(X2 < 1)
  
  # remove redundant terms
  list_slim <- list()
  x <- 1
  for(tmp_module in unique(temp_df2_topsig$module)) {
    tmp_df <- temp_df2_topsig %>% dplyr::filter(module == tmp_module)
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
        if (sum(vj %in% vi) / length(vj) > slim_th | sum(vi %in% vj) / length(vi) > slim_th) { 
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
  temp_df3 <- temp_df2_topsig %>% dplyr::select(-X3) %>% spread(key= module,value = X2) %>% dplyr::filter(terms_id %in% use_id_term)
  
  temp_df3$terms_id <- factor(temp_df3$terms_id, levels = use_id_term)
  temp_df3 %>% arrange(terms_id)
  
  temp_df4 <- temp_df3 %>% replace(is.na(.),1) %>%
    tidyr::gather(key = module, value = qvalue, colnames(temp_df3)[!(colnames(temp_df3) %in% c("X1","X4","terms_id"))])
  colnames(temp_df4)[1:2] <- c("terms", "id")
  temp_df4 <- temp_df4 %>% mutate(log10qval=-log10(as.numeric(qvalue)))
  temp_df4$module <- factor(temp_df4$module, levels=module_level)
  temp_df4$terms_id <- factor(temp_df4$terms_id, levels=use_id_term)
  #draw_heatmap
  temp_df4 %>% 
    arrange(terms_id, module) %>% tidy_heatmap(., border_color = "grey60",
                                                rows = terms_id,
                                                columns = module,
                                                values = log10qval,
                                                scale = "none",
                                                colors = c("#ffffff","#390655"),
                                                color_legend_min = 0,
                                                color_legend_max = 4,
                                                cellheight = 16, cellwidth = 16, fontsize = 14, angle_col = 0,
                                                filename = filename
    )
}

# DrawEnrichmentHeatmap_module(df_reactome, top_th=top_th, 
#                              slim_th=slim_th, module_level=module_level,
#                              filename="Fig/Fig2g.pdf")
# DrawEnrichmentHeatmap_module(df_kegg, top_th=top_th,
#                              slim_th=slim_th,
#                              module_level=module_level,
#                              filename="Fig/FigS6b.pdf")









