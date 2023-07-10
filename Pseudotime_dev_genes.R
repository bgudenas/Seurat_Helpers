

Add_pseuodtimes_seurat = function(so,
                                  slingshot_sce){
  library(ggplot2)
  library(RColorBrewer)
  
  slings = colnames(ps@colData)[grepl("sling", colnames(ps@colData))]
  sling_mat = matrix(nrow = ncol(so), ncol = length(slings), data = NA)
  colnames(sling_mat) = slings
  rownames(sling_mat) = colnames(so)
  
  for (i in 1:length(slings)){
    
    pseudo_vec = ps@colData[ ,slings[i]]
    names(pseudo_vec) = rownames(ps@colData)
    sling_mat[ ,i] = pseudo_vec[match(colnames(so), names(pseudo_vec))]
  }
  
  # add sling_mat to seurat object ------------------------------------------
  so@meta.data = cbind(so@meta.data, sling_mat)
  colos <- colorRampPalette(brewer.pal(11,'Spectral'))(100)
  
  for (i in slings){
    df = slingCurves(ps)[[1]]$s[slingCurves(ps)[[1]]$ord, ]
    FeaturePlot(so, features = i,raster = TRUE, raster.dpi=c(1012,1012)) +
      scale_colour_gradientn(colors = colos)
    ggsave(last_plot(), device = "pdf", filename=paste0("../Figures/Slingshot/", i, ".pdf"), width=8, height=7)
  }
  
  so$mean_slingshot = rowMeans(sling_mat, na.rm=TRUE)
  FeaturePlot(so, features = "mean_slingshot", raster = TRUE, raster.dpi=c(1012,1012)) +
    scale_colour_gradientn(colors = colos)
  ggsave(last_plot(), device = "pdf", filename=paste0("../Figures/Slingshot/", "mean_slingshot", ".pdf"), width=8, height=7)
return(so)
  }
  
  



find_dev_genes = function(so,
                          sling_name = "slingPseudotime_1",
                          map_path="~/Annots/Annotables/Mm10.rds",
                          outout_rf_model="../Data/Slingshot/RF_model_ps1.rds"){
 library(Seurat)
 library(randomForest)
 
 map = readRDS(map_path)
 message(paste("Pseudotime vector = ", sling_name))
 
 map = map[map$gene_biotype == "protein_coding", ]
 map = map[!duplicated(map$hsapiens_homolog_associated_gene_name) & !duplicated(map$hsapiens_homolog_associated_gene_name, fromLast = TRUE),  ]
 drops = map$external_gene_name[grepl("^Rpl|^Rps|^mt-", map$external_gene_name)] ## remove ribosomal/ mitocondrial genes
 map = map[!(map$external_gene_name %in% drops), ]
 pc_genes = c(map$external_gene_name,  map$hsapiens_homolog_associated_gene_name)
 
 pseudo_vec = so@meta.data[ ,colnames(so@meta.data) == sling_name]
 names(pseudo_vec) = colnames(so)
 pseudo_vec = pseudo_vec[!is.na(pseudo_vec)]
 
 mega = FindVariableFeatures(mega, nfeatures = 12000) ## limit search space to protein coding genes in the top 12k most variable just for comp speed
 HVG = VariableFeatures(mega)
 pc_genes = pc_genes[pc_genes %in% HVG]
 
# make a normalized expression matrix of cells in pseudotime vecto --------
# we round the matrix digits to 2 decimals to reduce overall size
 gene_mat = as.matrix(round(so@assays$RNA@data[rownames(so) %in% pc_genes, names(pseudo_vec)], 2)) 

message("Final gene matrix size")
cell_vec = sample(1:ncol(gene_mat), 0.5*ncol(gene_mat))
print(dim(gene_mat[ ,cell_vec]))

# make random forest to derive feature importance scores ------------------
 start.time <- Sys.time()

 rf_mod = randomForest::randomForest(x = t(gene_mat[ ,cell_vec]), 
                                     y = pseudo_vec[cell_vec], 
                                    #importance = TRUE,
                                     ntree=1000)
 end.time <- Sys.time()
 time.taken <- end.time - start.time
 time.taken
 saveRDS(rf_mod, outout_rf_model)
}




# Plotting functions ------------------------------------------------------


plot_pseuodtime_genes = function(so,
                                 sling_name = "slingPseudotime_1",
                                 ymax = 2,
                                 genes = c("Otx2","Pax6", "Crx", "Bsx"),
                                 outname = "../Figures/Slingshot/Slingshot1_genes.pdf"){
  library(ggplot2)
 # library(tidyquant)
  library(RColorBrewer)
  
  gene_mat = matrix(nrow = ncol(so), ncol = length(genes),data = 0)
  colnames(gene_mat) = genes
  
  df = data.frame("Barcode" = colnames(so) )
  df = cbind(df, gene_mat)
  for ( i in genes){
    df[ ,i] = as.numeric(so@assays$RNA@data[i, ])
  }
  
  ldf = reshape2::melt(df, "Barcode")
  pseudo_vec = so@meta.data[ ,colnames(so@meta.data) == sling_name]
  
  ldf$Pseudotime = as.numeric(pseudo_vec)
  ldf = ldf[!is.na(ldf$Pseudotime), ]
  
  ldf$scaled_Pseudotime = (ldf$Pseudotime - min(ldf$Pseudotime))/(max(ldf$Pseudotime) - min(ldf$Pseudotime))
  ldf$Celltype = so$Celltype[match(ldf$Barcode, colnames(so))]
  ldf$Timepoint = so$Timepoint[match(ldf$Barcode, colnames(so))]
  
  g1 = ggplot(ldf, aes(x = scaled_Pseudotime, y = value, group.by=variable)) +
  #  geom_smooth(aes(col = variable, fill = variable), span = 0.9, se = TRUE, n=2000) +
    #coord_cartesian(ylim = c(0, 2.5)) +
   # ylim(0, ymax) +
    geom_smooth(aes(col = variable, fill = variable), method = "loess", span = 0.2, se = FALSE) +
   # geom_smooth(aes(col = variable, fill = variable), method ="lm", formula = y ~ splines::bs(x, 4))
    theme_bw() +
    th +
    ylab("Expression")
  
  g2 = ggplot(ldf, aes(x = scaled_Pseudotime, y = Celltype, color = Celltype)) +
    geom_boxplot(outlier.shape = NA, aes(fill = Celltype), alpha = 0.8) +
    # scale_color_manual(values = colos ) +
    theme_bw() +
    th +
    xlab("Pseudotime") +
    theme(legend.position = "none")
  
  g3 = ggplot(ldf, aes(x = scaled_Pseudotime, y = Timepoint, color = Timepoint)) +
    geom_boxplot(outlier.shape = NA, aes(fill = Timepoint, alpha = 0.8)) +
    # scale_color_manual(values = colos ) +
    theme_bw() +
    th +
    xlab("Pseudotime") +
    theme(legend.position = "none")
  
  g4 = cowplot::plot_grid(plotlist = list(g1, g2,g3),
                          nrow = 3,
                          align = "v",
                          axis ="lr",
                          rel_heights = c(3,1.5,1.5))
  
  ggsave(g4, device = "pdf", filename = outname, width=9, height=6)
}
