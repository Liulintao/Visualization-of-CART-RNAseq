library(tidyverse)
library(pheatmap)
library(gplots)
#library(fastcluster)

dir2 = "~/Desktop/p05_heatmap_tliu/v1"
dir = "~/Desktop/p05_heatmap_tliu/v2"
gene_list = read_csv(file.path(dir, "Gene_List2.csv"))
DGE = read_tsv(file.path(dir2, "Th1-VS-Th9.DEseq2_Method.GeneDiffExpFilter.xls"))
fpkm = read_tsv(file.path(dir2, "AllSamples.GeneExpression.FPKM.xls"))

draw_heatmap2 = function(fpkm_table, output, showRow, width, height)
{
  htable = fpkm_table %>%
    dplyr::select(Th1_1_FPKM, Th1_2_FPKM, Th1_3_FPKM, Th9_1_FPKM, Th9_2_FPKM, Th9_3_FPKM) %>%
    mutate(Th1_1_FPKM=log2(Th1_1_FPKM + 1),
           Th1_2_FPKM=log2(Th1_2_FPKM + 1),
           Th1_3_FPKM=log2(Th1_3_FPKM + 1),
           Th9_1_FPKM=log2(Th9_1_FPKM + 1),
           Th9_2_FPKM=log2(Th9_2_FPKM + 1),
           Th9_3_FPKM=log2(Th9_3_FPKM + 1)) %>%
    as.matrix() %>%
    t() %>% scale() %>% t()
  rownames(htable) = fpkm_table$SymbolID
  breaksList = seq(-2, 2, by = 0.1)
  
  if(showRow & nrow(htable)<=60)
  {
    labRow = rownames(htable)
  }
  else
  {
    labRow = FALSE
  }
  
  pdf(output, width, height)
  heatmap.2(htable, col="bluered", trace="none", 
            cexRow = 0.7, 
            cexCol = 1, 
            labRow = labRow,
            #lhei = c(1,6),
            lhei = c(1,1.9),
            #lhei = c(1,2),
            #lhei = c(1,3.6),
            #lhei = c(1,5),
            #lhei = c(1,3),
            scale="row", Rowv = TRUE, Colv="Rowv", 
            notecol="black",      # change font color of cell labels to black
            density.info="none",  # turns off density plot inside color legend
            margins =c(9,9),     # widens margins around plot
            distfun = function(x) dist(x, method="euclidean"),
            hclustfun = function(x) hclust(x, method="complete"))
  dev.off()
}

draw_heatmap = function(symbol_list, name, width, height)
{
  symbol_list = toupper(symbol_list)
  symbol_list = symbol_list[!is.na(symbol_list)]
  temp_fpkm = fpkm %>%
    filter(SymbolID %in% symbol_list) 
  temp_fpkm_dge = temp_fpkm %>%
    inner_join(DGE, by=c("SymbolID"="Symbol"))
  
  if(length(setdiff(symbol_list, temp_fpkm[['SymbolID']]))!=0)
  {
    print(paste(setdiff(symbol_list, temp_fpkm[['SymbolID']]), "in list", name, "was not expressed."))
  }
  
  output = file.path(dir, paste("htmap_",name, sep=""))
  draw_heatmap2(temp_fpkm_dge, paste(output, "_padj.pdf", sep=""), TRUE, width, height)
}

draw_heatmap(gene_list[[1]], colnames(gene_list[1]), 5, 10)
draw_heatmap(gene_list[[2]], colnames(gene_list[2]), 5, 4.2)
draw_heatmap(gene_list[[3]], colnames(gene_list[3]), 5, 4.5)
draw_heatmap(gene_list[[4]], colnames(gene_list[4]), 5, 7)
draw_heatmap(gene_list[[5]], colnames(gene_list[5]), 5, 9)
draw_heatmap(gene_list[[6]], colnames(gene_list[6]), 5, 6)
