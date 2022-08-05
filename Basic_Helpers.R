

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

clean_HVG = function(HVG, map){
  message(paste0("HVG length = ", length(HVG)))
  map = map[map$gene_biotype == "protein_coding", ]
  sex_genes = c("Xist","Eif2s3y","Tsix")
  HVG = HVG[!(HVG %in% sex_genes)]
  HVG = HVG[!grepl("^Rpl", HVG)] ## remove ribosomal genes
  HVG = HVG[!grepl("^Mt-", HVG)] ## remove mitochondrial genes
  HVG = HVG[HVG %in% map$external_gene_name ] ## only protein coding
  
  message(paste0("HVG length post-filter = ", length(HVG)))
  return(HVG)
}


ProcSeurat = function(mega, vars.to.regress=NULL, nDim = 20, nGenes = 2000, nResolution = 0.8, map) {
  
  mega <- NormalizeData(mega, normalization.method = "LogNormalize", scale.factor = 10000)
  
  if (!is.null(vars.to.regress) ) {
    mega <- ScaleData(mega, vars.to.regress = vars.to.regress, block.size=2000)
  } else {
    mega <- ScaleData( object = mega, verbose = FALSE, block.size=2000)
  }
  mega <- FindVariableFeatures(mega, nfeatures = nGenes+1000 ) ## add buffer to account for genes lost
  HVG = VariableFeatures(mega)
  HVG = clean_HVG(HVG, map)
  HVG = HVG[1:nGenes] # trim to nGenes
  
  mega <- RunPCA(mega, npcs = 100, features = HVG, ndims.print = 1:5, nfeatures.print = 5 )
  
  mega <- FindNeighbors(mega, reduction = "pca", dims = 1:nDim  )
  mega <- FindClusters(mega,  n.start = 100, random.seed=54321, resolution = nResolution  )
  mega <- RunUMAP(object = mega, reduction = "pca", dims = 1:nDim, n.neighbors=45, n.epochs=500)
  
  return(mega)
}
