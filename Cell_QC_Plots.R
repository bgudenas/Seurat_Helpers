
# QC_Plots

Cell_QC_Plots = function( so, plotfile = NULL){
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




library(dplyr)
library(ggplot2)

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
