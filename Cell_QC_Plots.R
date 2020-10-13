
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