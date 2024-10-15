
#Input_path="../Data/Retina_Atlas_Cleaned.rds"
#GS = readRDS("~/Proj/Adult_CB/MCD_PR_UBC_Gene_sets.rds")

Plot_permuted_gene_sets = function(seurat_path, 
                                   gene_set_path,
                                   plot=TRUE,
                                   plot_dir_path,
                                   downsample = FALSE,
                                   out_name="Gene_set_permuted_cutoffs",
				                           nPerm=1000,
				                           pos_thresh = 0.25){
  library(ggplot2)
  library(dplyr)
  library(Seurat)
  th <- theme(text = element_text(size=10, family = "Helvetica" ),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              legend.position = "none",
              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
map = readRDS("~/Annots/Annotables/Mm10.rds")
if (is.character(seurat_path)){ mega = readRDS(seurat_path)
} else { mega = seurat_path }

if (is.character(gene_set_path)){ GS = readRDS(gene_set_path)
} else { GS = gene_set_path }
message(paste0("Gene sets= "), length(GS))
message(paste("Permutations = ", nPerm))

if (downsample == TRUE){
  library(scuttle)
  library(SingleCellExperiment)
  message("DownSampling ------------")
  sce = as.SingleCellExperiment(mega)
  assay(sce, "logcounts") <- scuttle::normalizeCounts(sce, log = TRUE, downsample = TRUE)
  mega@assays$RNA@data = assay(sce, "logcounts")
  
}
if (sum(unlist(GS[[1]]) %in% rownames(mega))/length(unlist(GS[[1]])) < 0.3 ){ ## if less than 30% of genes present, convert to mouse
for ( i in names(GS)){
  genes = map$external_gene_name[ match(unlist(GS[[i]]), map$hsapiens_homolog_associated_gene_name) ]
  genes = genes[!is.na(genes)]
  GS[[i]] = genes
}
}

GS_filt = remove_short_lists(GS, rownames(mega))
print(paste0("Final GS length = ", length(GS_filt)))
GS = GS_filt

mega = suppressWarnings(AddModuleScore(mega, features = GS, name = names(GS)))
GS_names = paste0(names(GS), 1:length(GS))

#FeaturePlot(mega, features = GS_names, max.cutoff = "q95", min.cutoff = "q05")
#ggsave(last_plot(), device = "pdf", filename = file.path(plot_dir_path, paste0("UMAP_GeneSets_raw_", out_name, ".pdf")), width = 14, height = 14)

cutoffs = c()
dummy_df = data.frame(
  g_length = 0,
  cutoff = 0)
write.table(dummy_df,
              file = "Cutoff_vals_tmp.csv",
              sep = ",",
              row.names = FALSE)

for (i in names(GS)){
  print(i)
  g_length = length(GS[[i]])
  check_vals = read.csv("Cutoff_vals_tmp.csv")
  if ( g_length %in% check_vals$g_length){
    message("GS length found -- using stored value ")
    cutoff_val = check_vals$cutoff[check_vals$g_length == g_length]
  } else {
  cutoff_val = median(permute_AddModule(mega, 
                           gene_set = GS[[i]],
                           nPerm=nPerm))
  # Append the new rows to the existing CSV file
  new_rows = data.frame("g_length" = g_length,
               "cutoff" = cutoff_val)
  write.table(new_rows, file = "Cutoff_vals_tmp.csv",
              sep = ",", 
              row.names = FALSE,
              col.names = FALSE,
              append = TRUE)
  }
  cutoffs = c(cutoffs, cutoff_val)
}
suppressWarnings(file.remove("Cutoff_vals_tmp.csv", showWarnings=FALSE))
names(cutoffs) = names(GS)

for ( i in 1:length(cutoffs)){
  feat = paste0(names(cutoffs)[i], i)
  
  nes = mega@meta.data[ ,feat]
  nes = nes - cutoffs[[i]]
  mega@meta.data[ ,feat] = nes
}
if (plot == TRUE){
<<<<<<< HEAD
FeaturePlot(mega, features = GS_names, min.cutoff = "q01", max.cutoff = "q90", raster = TRUE, raster.dpi=c(1012,1012), reduction="umap")
ggsave(last_plot(), device = "pdf", filename = file.path(plot_dir_path, paste0("UMAP_GeneSets_permuted_",out_name, ".pdf")), width = 14, height = 14)
=======
FeaturePlot(mega,
            features = GS_names,
            pt.size = 1.5,
            min.cutoff = "q01", max.cutoff = "q90",
            raster = TRUE, raster.dpi=c(500,500),
            reduction="umap")
ggsave(last_plot(), device = "pdf", filename = file.path(plot_dir_path, paste0("UMAP_GeneSets_permuted_",out_name, ".pdf")), width = 13, height = 6)
>>>>>>> fd6ac7c (Updating funcs with new seurat)
}

df = as.data.frame( mega@meta.data[ ,GS_names])
colnames(df) = GS_names
df$Celltype = mega$Celltype
rownames(df) = colnames(mega)
saveRDS(df, file.path(plot_dir_path, paste0(out_name,".rds")))

if (plot == TRUE ){
ldf = reshape2::melt(df, "Celltype")

plot_list = list()
for ( i in unique(ldf$variable)){
  
  pdf = ldf[ldf$variable == i, ]
  
  sigs = pdf %>% 
    group_by(Celltype, variable) %>% 
    summarise("pos_perc" = sum(value > 0)/n()) %>% 
    filter(pos_perc > pos_thresh )
  
  pdf$Positive = "NO"
  pdf$Positive[pdf$Celltype %in% sigs$Celltype] = "YES"
  
  colos = setNames( c("grey","firebrick"),c( "NO", "YES") )
  
  ords = pdf %>% 
    group_by(Celltype) %>% 
    summarise("Avg" = mean(value )) %>% 
    arrange(desc(Avg)) %>% 
    .$Celltype %>% 
    unique
  pdf$Celltype = factor(pdf$Celltype, levels = ords)
    
g1 = ggplot(pdf, aes(x= Celltype, y = value )) +
  geom_hline(yintercept = 0, col="black", lwd = 1, linetype = "dashed") +
  geom_violin(aes(fill = Positive), adjust=0.75) +
  scale_fill_manual(values = colos ) +
  ylab("Enrichment score") +
  theme_bw() +
  th  +
  ggtitle( paste0(i, " (Positive >=", pos_thresh*100, "%)"))

plot_list[[i]] = g1
  
}
cowplot::plot_grid(plotlist = plot_list)
ggsave(last_plot(), device = "pdf", filename = file.path(plot_dir_path, paste0("Violins_GeneSets_permuted_", out_name, ".pdf")), width = 12, height = 12)
}
}

permute_AddModule = function(gene_set, so, nPerm=1000){
  set.seed(54321)
  #perm_mat = matrix(nrow = nPerm, ncol = ncol(so), data = 0)
  perm_vec = vector(mode = "list", length = nPerm)
  g_size = sum(gene_set %in% rownames(so))
  elig_genes = rownames(so)
  for (i in 1:nPerm){
  #  print(i)
    rand_genes = sample(rownames(so), size = g_size, replace = FALSE)
    so = suppressWarnings(AddModuleScore(so, features = list(rand_genes) ))
    
    perm_vec[[i]] = quantile( as.numeric(so$Cluster1), 0.95)
    #perm_mat[i, ] = as.numeric(so$Cluster1)
  }
  return(unlist(perm_vec))
}



#Plot_permuted_gene_sets(seurat_path="../Data/Retina_Atlas_Cleaned.rds", 
#                                   gene_set_path="~/Proj/Adult_CB/MCD_PR_UBC_Gene_sets.rds", 
#                                   plot_dir_path="../QC/",
#                                   out_name="Gene_set_permuted_cutoffs.rds")
Plot_ssgsea_gene_sets = function(seurat_path, 
                                 gene_set_path, 
                                 plot_dir_path,
                                 out_name="Gene_set_permuted_cutoffs"){
  library(ggplot2)
  library(dplyr)
  library(Seurat)
  th <- theme(text = element_text(size=10, family = "Helvetica" ),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              legend.position = "none",
              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  map = readRDS("~/Annots/Annotables/Mm10.rds")
  mega = readRDS(seurat_path)
  GS = readRDS(gene_set_path)
  
  ## convert genes to mouse if needed
  if (sum(unlist(GS[[1]]) %in% rownames(mega))/length(unlist(GS[[1]])) < 0.3 ){ 
    for ( i in names(GS)){
      genes = map$external_gene_name[ match(unlist(GS[[i]]), map$hsapiens_homolog_associated_gene_name) ]
      genes = genes[!is.na(genes)]
      GS[[i]] = genes
    }
  }
  if ( ncol(mega) > 50000 ){
    message("Chunking cells for SSGSEA")
    elements_per_chunk = 50000
    chunks = split(1:ncol(mega), ceiling(seq_along(1:ncol(mega))/elements_per_chunk))
    out_mat = matrix(nrow = ncol(mega), ncol = length(GS), data = 0)
    for (i in chunks ){
      tmp_mat = as.matrix(mega@assays$RNA@data[ ,i])
      ssgMat = GSVA::gsva(expr = tmp_mat,
                          gset.idx.list = GS,
                          method = "ssgsea",
                          ssgsea.norm = FALSE,
                          verbose=FALSE)
      out_mat[i[[1]], ] = t(ssgMat)
      
    }
  } else {
    ssgMat = GSVA::gsva(expr = tmp_mat,
                        gset.idx.list = GS,
                        method = "ssgsea",
                        ssgsea.norm = FALSE,
                        verbose=FALSE)
  }
  colnames(out_mat) =  rownames(ssgMat)
  saveRDS(out_mat, paste0("../Data/ssgsea_", out_name, ".rds"))
  
  # Plotting time -----------------------------------------------------------
  df = as.data.frame(out_mat)
  df$Celltype = mega$Celltype
  ldf = reshape2::melt(df)
  plot_list = list()
  for ( i in unique(ldf$variable)){
    
    pdf = ldf[ldf$variable == i, ]
    
    ords = pdf %>% 
      group_by(Celltype) %>% 
      summarise("Avg" = mean(value )) %>% 
      arrange(desc(Avg)) %>% 
      .$Celltype %>% 
      unique
    pdf$Celltype = factor(pdf$Celltype, levels = ords)
    
    g1 = ggplot(pdf, aes(x= Celltype, y = value )) +
      geom_hline(yintercept = 0, col="black", lwd = 1, linetype = "dashed") +
      geom_violin( adjust=0.75) +
      ylab("Enrichment score") +
      theme_bw() +
      th  +
      ggtitle( paste0(i, "-ssgsea"))
    
    plot_list[[i]] = g1
    
  }
  cowplot::plot_grid(plotlist = plot_list)
  ggsave(last_plot(), device = "pdf", filename = file.path(plot_dir_path, paste0("Violins_ssgsea_", out_name, ".pdf")), width = 12, height = 12)
}



Plot_gs_loess = function(df, y_index = 1, x_index = 2, out_plot_path = "./GSets_loess.pdf" ){
  library(ggplot2)
  library(ggside)
  
  df$y_var = df[ ,y_index]
  df$x_var = df[ ,x_index]
  
  g1 = ggplot(df, aes(x = x_var,
                      y = y_var, colour = Celltype)) +
    geom_smooth( span=0.25, method = "loess", se = FALSE) +
    geom_xsideboxplot(aes(y = Celltype, fill= Celltype), color = "black", orientation = "y", outlier.shape = NA, width= 1, varwidth = FALSE) +
    geom_ysideboxplot(aes(x = Celltype, fill= Celltype),color = "black", orientation = "x", outlier.shape = NA, width= 1, varwidth = FALSE) +
    scale_xsidey_discrete() +
    scale_ysidex_discrete() +
    ylab(label = colnames(df)[1]) +
    xlab(label = colnames(df)[2]) +
    theme_bw() +
    theme(legend.position = "none")
  
  ggsave(g1, device = "pdf",filename=out_plot_path, width = 10, height = 8)
  return(g1)
}




remove_short_lists = function(gs, so_genes){
  gs_filt = list()
  for ( i in names(gs)){
    genes = gs[[i]]
    genes = genes[genes %in% so_genes ]
    if ( length(genes) > 12){
      gs_filt[[i]] =  gs[[i]]
    }
  }
  return(gs_filt)
}
