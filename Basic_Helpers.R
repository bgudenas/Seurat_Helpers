PlotSeurat = function(plot, name, path ){

p1 = p1 + ggtitle(name)
#p1 = AugmentPlot(plot = p1)
ggsave(p1, device = "pdf", filename = paste0(path, name, ".pdf"), height = 12, width = 10)

}



AutoCellType = function(clustmat, GO ){

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

return(df)
}

 
ProcSeurat = function(mega, vars.to.regress=NULL ) {

mega <- NormalizeData(mega, normalization.method = "LogNormalize", scale.factor = 10000)

if (!is.null(vars.to.regress) ) {
mega <- ScaleData(mega, vars.to.regress = vars.to.regress )
} else {
mega <- ScaleData( object = mega, verbose = FALSE )
}
mega <- FindVariableFeatures(mega, nfeatures = 3000 )
mega <- RunPCA(mega, npcs = 41, ndims.print = 1:5, nfeatures.print = 5 )

mega <- FindNeighbors(mega, reduction = "pca", dims = 1:40 )
mega <- FindClusters(mega,  n.start = 100, random.seed=54321 )
mega <- RunUMAP(object = mega, reduction = "pca", dims = 1:40)

return(mega)
}
