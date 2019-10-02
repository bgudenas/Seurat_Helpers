
source("./Seurat_Helpers/Basic_Helpers.R")

Celltype_Assign = function(mega ) {
  
cluster.averages <- AverageExpression(object = mega, verbose = FALSE, use.scale = TRUE)
clustmat =  cluster.averages$RNA
GO  = readRDS("./Data/scRNA_celltypes.rds")
cellAssignments = AutoCellType(clustmat, GO)
clusters  = mega@active.ident

mega@active.ident = as.factor(cellAssignments$SSGSEA [match(mega@active.ident, cellAssignments$Cluster) ] )
names(mega@active.ident) = names(clusters)

p1 = DimPlot(object = mega, label = TRUE,  label.size = 5, pt.size = 0.2 ) + ggtitle("Automated Cell types")
p2 = DimPlot(object = mega, label = TRUE, group.by = "Timepoint",  label.size = 5, pt.size = 0.2 ) + ggtitle( "Timepoints")
p3 = CombinePlots(plots = list(p1, p2), legend = "none")
pcount = sum(grepl("UMAP_cellTypes", list.files("~/Proj/scPB/Figures/"))) + 1
ggsave(p3, device = "pdf", filename = paste0("~/Proj/scPB/Figures/UMAP_cellTypes_", pcount, ".pdf"), width = 14, height = 10)

drops = c("Reticulocytes", "Nuocytes", "Platelets","Erythroid−like and erythroid precursor cells","Basophils","Macrophages", 
          "Oligodendrocyte progenitor cells","Ependymal cells", "Immature nerons", 
          "Glutaminergic neurons", "GABAergic neurons", "Cajal-Retzius cells", "Choroid plexus cells", "Myeloid−derived suppressor cells")

dropcells = names(mega@active.ident)[ (mega@active.ident %in% drops) ]
dropcells = unlist(lapply(stringr::str_split(dropcells, "_"), "[[", 1))
names(dropcells) =  mega$Timepoint[(mega@active.ident %in% drops) ]
mega@active.ident = clusters

df = data.frame("Cells" = dropcells, "Timepoint" = names(dropcells) )
write.table(df, "./Data/DropCells.csv", append = TRUE, sep = ",")

}

Reform_Seurat = function(samps ) {
  sampIDs = c()
  dropcells = read.csv("./Data/DropCells.csv", row.names=NULL)
  
  for (i in 1:length(samps )) {
    ID = unlist(lapply(stringr::str_split(dirname(samps[i]), "/"), "[[", 3 ))
    print(ID)
    ## Break on "-"; TP is which element has 3 chars
    TP = unlist(stringr::str_split(ID, "-" ))[unlist(lapply(stringr::str_split(ID, "-" ), nchar)) == 3]
    sampIDs = c(sampIDs, TP)
    counts <- Read10X(data.dir = samps[i] )
    
    drops = dropcells$Cells[dropcells$Timepoint %in% TP ]
    print(paste("Removed", sum( colnames(counts) %in% drops), "in", TP ))
    
    counts =counts[ , !(colnames(counts) %in% drops) ]
    s1 = CreateSeuratObject(counts  =  counts, min.cells = 40, min.features = 500, project = TP )
    rm(counts)
    s1[["percent.mt"]] <- PercentageFeatureSet(object = s1, pattern = "^mt-")    
    s1$Timepoint = TP
    
    if ( !exists("mega") ){
      mega = s1
    } else {
      mega = merge( mega, y=s1, project = "Pineal_Atlas" )
    }
  }
  
  mega = subset(x = mega, subset = nCount_RNA > 5000  &  nFeature_RNA > 1500 )
  mega = subset(x = mega, subset = nCount_RNA < 70000 )
  mega = subset(x = mega, subset = percent.mt < 12 )
  mega = ProcSeurat(mega)
  
  ##Change to simple cap for mouse
  s.genes = tolower(cc.genes$s.genes)
  s.genes = unlist(lapply(s.genes, simpleCap ))
  g2m.genes = tolower(cc.genes$g2m.genes)
  g2m.genes = unlist(lapply(g2m.genes, simpleCap ))
  mega <- CellCycleScoring(mega, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  mega$CC.Difference <- mega$S.Score - mega$G2M.Score
  mega = ProcSeurat(mega,  "CC.Difference")
  
  return(mega)
}
