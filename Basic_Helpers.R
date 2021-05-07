

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
