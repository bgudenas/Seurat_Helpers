
# QC_Plots

Cell_QC_Plots = function( seurat_object, plotfile = NULL){
library(Seurat)
  library(ggplot2)
pdf(plotfile)
par(mfrow=c(2,2))

hist(mega$nCount_RNA, breaks = 100, col = "dodgerblue", main = "Count Depth", xlab = "") 
abline(v = median(mega$nCount_RNA), col = "red", lwd = 2)

hist(mega$nFeature_RNA,breaks = 100, col = "dodgerblue", main = "# of Genes", xlab = "") 
abline(v = median(mega$nFeature_RNA), col = "red", lwd = 2)

hist(mega$percent.mt, breaks = 100 , col = "dodgerblue", main = "% Mitochondria", xlab = "") 
abline( v = median(mega$percent.mt), col = "red", lwd = 2 )

## Plot Mitochondira by counts
df = data.frame("Count" = mega$nCount_RNA,
                "Genes" = mega$nFeature_RNA,
                "Mitochondria" = mega$percent.mt)

g1 = ggplot(df, aes(x = Count, y = Genes) ) +
  geom_point(size = 1, alpha = 0.6, aes(col = Mitochondria )) +
  scale_color_gradient(low = "black", high = "red") + 
  theme_bw()


print(g1)
dev.off()
}