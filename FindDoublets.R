# DoubletFinder wrapper-----------------------------------------------------------
## description: standard seurat pre-processing followed by doublet calling using DoubletFinder
## Returns a character vector of (Singlet | Doublet ) named by cell name


Find_Doublets = function(count_path,
                    mt_filt = 15,
                    nDims = 15 ){
library(Seurat)
library(DoubletFinder)
  sink(tempfile()) 
gene_counts = Read10X(count_path )
so = CreateSeuratObject(gene_counts, min.cells = 50)

## check if mouse cells
if (sum(grepl("MT-", rownames(so))) == 0 ){
  so[["percent.mt"]] = PercentageFeatureSet(object = so, pattern = "^Mt")
} else {
  so[["percent.mt"]] = PercentageFeatureSet(object = so, pattern = "^MT-")
}

# Adaptive QC thresholds --------------------------------------------------
qc.low_lib = scater::isOutlier(so$nCount_RNA, log = TRUE,  type="lower")
qc.low_genes = scater::isOutlier(so$nFeature_RNA, log = TRUE,  type="lower")

print(table(qc.low_lib, qc.low_genes))

keep_cells = colnames(so)[!( qc.low_lib | qc.low_genes)]
so = subset(so, cells = keep_cells  )
## hard mitochondria filter
so = subset(so, subset = percent.mt < mt_filt  )
## standard pre-processing
so = NormalizeData(so)
so = FindVariableFeatures(so, selection.method = "vst", nfeatures = 2000)
so = ScaleData(so)
so = RunPCA(so)
so = RunUMAP(so, dims = 1:nDims)
so = FindNeighbors(object = so, reduction = "pca", dims = 1:nDims)
so = FindClusters(so, n.start =  100 , resolution = 0.5)

## DoubletFinder
sweep.res.list_so = paramSweep_v3(so, PCs = 1:nDims, sct = FALSE)
sweep.stats_so = summarizeSweep(sweep.res.list_so, GT = FALSE)
bcmvn_so = find.pK(sweep.stats_so)

pK =  as.numeric(as.character(bcmvn_so$pK[which.max(bcmvn_so$BCmetric)]))

homotypic.prop = modelHomotypic(so$seurat_clusters)          
nExp_poi = round(0.04 * nrow(so@meta.data))  ## 5% doublet rate
nExp_poi.adj = round(nExp_poi*(1-homotypic.prop))

so = doubletFinder_v3(so, PCs = 1:nDims, pN = 0.25, pK = pK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
DF_col =  which(grepl("DF.classifications", colnames(so@meta.data) ))

dub_calls = so@meta.data[ ,DF_col]
names(dub_calls) = colnames(so)

sink()
return(dub_calls)


}
