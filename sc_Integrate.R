## Description: loop through samples and remove cells with low counts/genes and doublets within each sample. 
## Then integrate all samples and do standard seurat pre-processing. Plot standard QC's along the way.


sc_Integrate = function( samps, ## sample names equal in length to count_paths
                         count_paths , ## paths to cellranger dirs
                         out_data_path, ## out name for RDS file
                         nDims = 30, ## PCA dims
                         mt_filt=20, ## percent mitochondria filter
                         min_genes = 400,
                         CC = TRUE,
                         rem_Xist = FALSE, ## remove Xist from HVG
                         drop_cells = NULL, ## cell names to remove before processing
                         normData = NULL, ## seurat object which bypasses integration IE just want to test drop_cells, PCA/Cluster params
                         GO=NULL, ## GO  to do automated cell type annotation
                         markers=NULL, ## gene markers to plot
                         project = "Proj",
                         QC_dir ="." ## where to put QC plots
                         ) {
  ## TODO make package to avoid this
  source("/home/bgudenas/src/Seurat_Helpers/Cell_QC_Plots.R")
  source("/home/bgudenas/src/Seurat_Helpers/FindDoublets.R")
  library(Seurat)
  library(ggplot2)
  library(stringr)
  
  stopifnot(length(samps) == length(count_paths))
  if ( is.null(normData)){

# Loop through each sample and apply adaptive thresholds ------------------
  print(paste("Samples = ", length(samps) ))
  so_list = vector( length = length(samps), mode = "list")
  for ( i in 1:length(samps)) {
    
    counts <- Read10X(data.dir = count_paths[i] )
    so =  CreateSeuratObject(counts = counts, project = samps[i], min.cells = 20, min.features = min_genes)
    so$Sample = rep(samps[i], ncol(so) )
    
    if (sum(grepl("MT-", rownames(so))) == 0 ){ ## check if Mt genes are mouse or human
      so[["percent.mt"]] <- PercentageFeatureSet(object = so, pattern = "^Mt")
    } else {
      so[["percent.mt"]] <- PercentageFeatureSet(object = so, pattern = "^MT-")
    }
    dir.create("Sample_QC")
    Cell_QC_Plots(so, plotfile = file.path("Sample_QC", paste0(samps[i],".pdf") ))
    
    # Adaptive QC thresholds --------------------------------------------------
    doublets = quiet(Find_Doublets( count_paths[i] ))
    singlet_names = names(doublets)[doublets == "Singlet" ]
    print(paste0("Doublet_Filter=", sum(doublets == "Doublet")))
    
    qc.low_lib = scater::isOutlier(so$nCount_RNA, log = TRUE,  type="lower")
    print(paste0("Low_Count_Filter=", sum(qc.low_lib)))
    
    qc.low_genes = scater::isOutlier(so$nFeature_RNA, log = TRUE,  type="lower")
    print(paste0("Low_Gene_Filter=", sum(qc.low_genes)))
    
    cell_filt = !( qc.low_lib & qc.low_genes )
    keep_cells = unique(colnames(so)[cell_filt ] )
    keep_cells = keep_cells[keep_cells %in% singlet_names ]
    so = subset(so, cells = keep_cells  )
    
    so_list[[i]] = so
    print(paste(samps[i],"=", ncol(so_list[[i]])))
    
  }
  
  so_big <- merge(so_list[[1]], y = so_list[-1], add.cell.ids = samps, project = project )
  rm(so_list, counts)
  table(so_big$orig.ident)
  Joined_QC(so_big, QC_dir, "Joined_QC.pdf")
  
  so_big <- subset(x = so_big, subset = percent.mt < mt_filt )
  ## Joint low count filter
  keep_cells = names(so_big$orig.ident)[! scater::isOutlier(so_big$nCount_RNA, log = TRUE,  type="lower") ]
  ## Joint High count filter
  keep_cells = c(keep_cells,
                 names(so_big$orig.ident)[! scater::isOutlier(so_big$nCount_RNA,  type="higher") ]
  )
  keep_cells = keep_cells[duplicated(keep_cells)] ## cells passing both filters will be duplicated in vector
  so_big = subset(so_big, cells = keep_cells  ) 
  
  } else { so_big = readRDS( normData )}
  
  if ( !is.null( drop_cells )){
    keep_cells = colnames(so_big)[!(colnames(so_big) %in% drop_cells )]
    num_drops = ncol(so_big) - length(keep_cells)
    print(paste0("Removing drop cells", num_drops ))
    so_big = subset(so_big, cells = keep_cells  ) 
  }
  
  Joined_QC(so_big, QC_dir, "Joined_Post_filt_QC.pdf")
  nCells = ncol(so_big)
  print(paste("Total cells =", nCells ))
  print(paste("Total genes =", nrow(so_big)))
  print(as.data.frame(table(so_big$Sample)))
  
  so_big <- NormalizeData(so_big, normalization.method = "LogNormalize", scale.factor = 10000)
  if (sum(cc.genes$s.genes %in% rownames(so_big)) == 0 ){
    cc.genes.updated.2019$s.genes = stringr::str_to_title(cc.genes.updated.2019$s.genes)
    cc.genes.updated.2019$g2m.genes = stringr::str_to_title(cc.genes.updated.2019$g2m.genes)
  }
  if ( CC == TRUE ){
  so_big <- CellCycleScoring(so_big, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = TRUE)
  so_big$CC.Difference <- so_big$S.Score - so_big$G2M.Score
  so_big <- quiet( ScaleData(so_big, vars.to.regress = "CC.Difference" ) )
  } else { so_big <- quiet( ScaleData(so_big  ) )  }
  
  so_big <- FindVariableFeatures(object = so_big, nfeatures = 2000, selection.method = "vst")
  HVG = VariableFeatures(object = so_big)
  if (rem_Xist == TRUE) { HVG = HVG[ HVG != "Xist" ]}
  so_big <- RunPCA(object = so_big,  verbose = FALSE, features = HVG)
  ggsave(ElbowPlot(so_big, ndims = 50), device = "pdf", filename = file.path(QC_dir, "Integrated_Elbow_Plot.pdf"))
  
  so_big <- RunUMAP(object = so_big, reduction = "pca", dims = 1:nDims, n.epochs = 500 )
  so_big <- FindNeighbors(object = so_big, reduction = "pca", dims = 1:nDims)
  so_big <- FindClusters(so_big, n.start =  100, resolution = 0.6, random.seed = 54321, group.singletons = FALSE) ## decrease resolution for broader clusters
  
  saveRDS(so_big, out_data_path )
  
  Prob_Clusts(so_big, QC_dir) ## plots clusters driven by single samples
  
  g1 = DimPlot(so_big, group.by = "Sample" ) + ggtitle(paste0("Samples= ",length(unique(so_big$Sample)) ))
  ggsave(g1, device = "pdf", filename = file.path(QC_dir, "UMAP_Samples.pdf"))
  so_big$Log_nCount = log( so_big$nCount_RNA )
  g1 = FeaturePlot(so_big, features = c("Log_nCount", "percent.mt") ) + ggtitle(paste0("Cells = ", nCells ))
  ggsave(g1, device = "pdf", filename = file.path(QC_dir, "UMAP_QC_metrics.pdf"))
  
  g1 = DimPlot(so_big, group.by = "Phase" ) + ggtitle(paste0("Cells = ", nCells ))
  ggsave(g1, device = "pdf", filename = file.path(QC_dir, "UMAP_Phase.pdf"))
  
  g1 = DimPlot(so_big, label = TRUE ) + ggtitle(paste0("Cells = ", nCells ))
  ggsave(g1, device = "pdf", filename = file.path(QC_dir, "UMAP_Clusters.pdf"))
  
  if (!is.null(markers)){
    g1 = FeaturePlot(so_big, min.cutoff = "q10", features = markers )
    ggsave(g1, device = "pdf", filename = file.path(QC_dir, "UMAP_Markers.pdf"))
  }
}


quiet <- function(x) { 
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(x)) 
} 
