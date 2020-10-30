
# QC_Plots

Cell_QC_Plots = function( so, plotfile = NULL){
## plots distributions of single cell counts, number of genes and percent mitochondria
library(Seurat)
library(ggplot2)
pdf(plotfile)
par(mfrow=c(2,2))

hist(so$nCount_RNA, breaks = 100, col = "dodgerblue", main = "Count Depth", xlab = "") 
abline(v = median(so$nCount_RNA), col = "red", lwd = 2)

hist(so$nFeature_RNA,breaks = 100, col = "dodgerblue", main = "# of Genes", xlab = "") 
abline(v = median(so$nFeature_RNA), col = "red", lwd = 2)

hist(so$percent.mt, breaks = 100 , col = "dodgerblue", main = "% Mitochondria", xlab = "") 
abline( v = median(so$percent.mt), col = "red", lwd = 2 )

## Plot Mitochondira by counts
df = data.frame("Count" = so$nCount_RNA,
                "Genes" = so$nFeature_RNA,
                "Mitochondria" = so$percent.mt)

g1 = ggplot(df, aes(x = Count, y = Genes) ) +
  geom_point(size = 1, alpha = 0.6, aes(col = Mitochondria )) +
  scale_color_gradient(low = "black", high = "red") + 
  theme_bw()


print(g1)
dev.off()
}




# library(dplyr)
# library(ggplot2)

Joined_QC = function(mega, QC_dir, plot_name) {
  ## Plots QC metrics across samples to find problematic samples
  library(ggplot2)
  library(dplyr)
  df = data.frame("Sample" = mega$Sample, "nCount_RNA" = log2(mega$nCount_RNA), "percent.mt" = mega$percent.mt )
  ords = df %>% 
    group_by(Sample) %>% 
    summarise("Med" = median(nCount_RNA)) %>% 
    arrange(desc(Med)) %>% 
    .$Sample
  df$Sample = factor(df$Sample, ordered = TRUE, levels = ords )
  
  g1 = ggplot(df, aes(x = Sample, y = nCount_RNA, group = Sample )) +
    geom_violin( aes(fill = Sample ) ) +
    geom_hline(yintercept =  median(df$nCount_RNA), linetype = "dashed", col = "black") + 
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none" ) + 
    ylab("Log2_nCount_RNA") 
  
  g2 = ggplot(df, aes(x = Sample, y = percent.mt, group = Sample )) +
    geom_violin( aes(fill = Sample ) ) +
    geom_hline(yintercept =  median(df$percent.mt), linetype = "dashed", col = "black") + 
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none" ) + 
    ylab("percent.mt") 
  
  cowplot::plot_grid(g1,g2, ncol= 1)
  ggsave(last_plot(), filename = file.path(QC_dir, plot_name), device = "pdf")
}




Prob_Clusts = function(so_big, QC_dir, percent_threshold = 0.7 ) {
  ## Find and plot problematic clusters driven primarily by a single sample
  library(dplyr)
  library(ggplot2)
  samp_clusts = data.frame( table(so_big$seurat_clusters, so_big$Sample))
  rela = samp_clusts %>% group_by(Var1) %>% mutate("Rel" = Freq/sum(Freq))
  prob_clusts = rela$Var1[rela$Rel > percent_threshold ]
  
  samp_clusts = samp_clusts[samp_clusts$Var1 %in% prob_clusts, ]
  
  g1 = ggplot(samp_clusts, aes(x = Var1 , y = Freq, fill = Var2 )) +
    geom_bar(stat = "identity", position = "fill" ) + 
    theme_bw() +
    ggtitle("Problematic Clusters ( >70% of cluster from single sample)")
  ggsave(last_plot(),device = "pdf", filename = file.path(QC_dir, "Problematic_Clusters.pdf"))
}







AutoCellType = function(mega, GO, QC_dir ){
  ## GO is a list of lists (celltypes and their marker genes)
  cluster.averages <- AverageExpression(object = mega, verbose = FALSE, use.scale = TRUE)
  clustmat =  cluster.averages$RNA
  
  resSS = GSVA::gsva(as.matrix(clustmat), GO, method="ssgsea", ssgsea.norm = FALSE, min.sz= 3, verbose = FALSE )
  df = as.data.frame(matrix(nrow= ncol(clustmat), ncol = 3 ) )
  colnames(df) = c("Cluster", "SSGSEA","Score" )
  rownames(df) = colnames(clustmat)
  
  vals = c()
  for (i in 1:ncol(resSS)){
    vals =c(vals, rownames(resSS)[which.max(resSS[ ,i])] )
  }
  
  df$SSGSEA = vals
  df$Score = apply(resSS, 2, max)
  df$Cluster = colnames(clustmat)
  mega$auto_celltype = df$SSGSEA[match(mega$seurat_clusters, df$Cluster)]
  g1 = DimPlot(mega, group.by= "auto_celltype", label = TRUE)
  ggplot2::ggsave(g1, device = "png", filename=file.path(QC_dir, "Auto_celltype.png"), width = 10, height = 7, dpi = 100)
  
  return(mega)
}
