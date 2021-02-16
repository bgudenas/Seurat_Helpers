# DoubletFinder wrapper-----------------------------------------------------------
## description: standard seurat pre-processing followed by doublet calling using DoubletFinder
## Returns a character vector of (Singlet | Doublet ) named by cell name


Find_Doublets = function(count_path,
                    mt_filt = 10,
                    nDims = 15,
                    min_count = NULL,
                    min_gene = NULL){
library(Seurat)
library(DoubletFinder)
  sink(tempfile()) 
gene_counts = Read10X(count_path )
so = CreateSeuratObject(gene_counts, min.cells = 5, min.features = 200 )

## check if mouse cells
if (sum(grepl("MT-", rownames(so))) == 0 ){
  so[["percent.mt"]] = PercentageFeatureSet(object = so, pattern = "^Mt")
} else {
  so[["percent.mt"]] = PercentageFeatureSet(object = so, pattern = "^MT-")
}
## hard mitochondria filter
so = subset(so, subset = percent.mt < mt_filt  )

# Adaptive QC thresholds --------------------------------------------------
if (is.null(min_count) | is.null(min_gene)) {
qc.low_lib = scater::isOutlier(so$nCount_RNA, log = TRUE,  type="lower")
qc.low_genes = scater::isOutlier(so$nFeature_RNA, log = TRUE,  type="lower")

print(table(qc.low_lib, qc.low_genes))
keep_cells = colnames(so)[!( qc.low_lib | qc.low_genes)]
so = subset(so, cells = keep_cells  )
} else {
  keep_cells = colnames(so)[so$nCount_RNA > min_count ]
  print(paste0("LOW Counts removing = ", ncol(so) - length(keep_cells) ))
  so_big = subset(so, cells = keep_cells  ) 
  
  keep_cells = colnames(so)[so$nFeature_RNA > min_gene ]
  print(paste0("LOW Genes removing = ", ncol(so) - length(keep_cells) ))
  so_big = subset(so, cells = keep_cells  ) 
  
}

## standard pre-processing
so = NormalizeData(so)
so = FindVariableFeatures(so, selection.method = "vst", nfeatures = 2000)
so = ScaleData(so)
so = RunPCA(so)
so = RunUMAP(so, dims = 1:nDims)

## DoubletFinder
sweep.res.list_so = paramSweep_v3(so, PCs = 1:nDims, sct = FALSE)
sweep.stats_so = summarizeSweep(sweep.res.list_so, GT = FALSE)
bcmvn_so = find.pK(sweep.stats_so)

pK =  as.numeric(as.character(bcmvn_so$pK[which.max(bcmvn_so$BCmetric)]))
print(paste0("pK =", pK ))

doublet_rate = (nrow(so@meta.data)/500 )* 0.004 ## (0.4%  doublet rate per 500 cells )
nExp_poi = round(doublet_rate * nrow(so@meta.data))  ## 4% doublet rate

so = doubletFinder_v3(so, PCs = 1:nDims, pN = 0.25, pK = pK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
DF_col =  which(grepl("DF.classifications", colnames(so@meta.data) ))

dub_calls = so@meta.data[ ,DF_col]
names(dub_calls) = colnames(so)

sink()
return(dub_calls)


}
