

CellType_Transfer = function(ref_so, query_so) {
  ## transfer cell type labels from one seurat dataset to another

  ref_so <- FindVariableFeatures(object = ref_so, nfeatures = 2000 )
  g.anchors <- FindTransferAnchors(reference = ref_so, query = query_so , 
                                   dims = 1:30)
  predictions <- TransferData(anchorset = g.anchors, refdata = ref_so$cell_type, 
                              dims = 1:30)
  ref_so <- AddMetaData(object = query_so, metadata = predictions)
  
}


PlotSeurat = function(plot, name, path ){
  p1 = p1 + ggtitle(name)
  p1 = AugmentPlot(plot = p1)
  ggsave(p1, device = "pdf", filename = paste0(path, name, ".pdf"), height = 12, width = 10, useDingbats = FALSE)
  
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
  ggplot2::ggsave(g1, device = "png", filename=file.path(QC_dir, "Auto_celltype.png"), width = 10, height = 7)
  
  return(mega)
}


ProcSeurat = function(mega, vars.to.regress=NULL, nDim = 20, nGenes = 2000 ) {
  
  mega <- NormalizeData(mega, normalization.method = "LogNormalize", scale.factor = 10000)
  
  if (!is.null(vars.to.regress) ) {
    mega <- ScaleData(mega, vars.to.regress = vars.to.regress )
  } else {
    mega <- ScaleData( object = mega, verbose = FALSE )
  }
  mega <- FindVariableFeatures(mega, nfeatures = nGenes )
  mega <- RunPCA(mega, npcs = 75, ndims.print = 1:5, nfeatures.print = 5 )
  
  mega <- FindNeighbors(mega, reduction = "pca", dims = 1:nDim  )
  mega <- FindClusters(mega,  n.start = 100, random.seed=54321, resolution = 0.5  )
  mega <- RunUMAP(object = mega, reduction = "pca", dims = 1:nDim )
  
  return(mega)
}
