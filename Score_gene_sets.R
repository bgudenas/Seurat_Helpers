
#Input_path="../Data/Retina_Atlas_Cleaned.rds"
#GS = readRDS("~/Proj/Adult_CB/MCD_PR_UBC_Gene_sets.rds")

Plot_permuted_gene_sets = function(seurat_path, 
                                   gene_set_path, 
                                   plot_dir_path,
                                   out_name="Gene_set_permuted_cutoffs.rds",
				   nPerm=100){
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
message(paste0("Gene sets= "), length(GS))
message(paste("Permutations = ", nPerm))

for ( i in names(GS)){
  genes = map$external_gene_name[ match(unlist(GS[[i]]), map$hsapiens_homolog_associated_gene_name) ]
  genes = genes[!is.na(genes)]
  GS[[i]] = genes
}

mega = AddModuleScore(mega, features = GS, name = names(GS))
GS_names = paste0(names(GS), 1:length(GS))

FeaturePlot(mega, features = GS_names, max.cutoff = "q95", min.cutoff = "q05")
ggsave(last_plot(), device = "pdf", filename = file.path(plot_dir_path, "UMAP_GeneSets.pdf"), width = 14, height = 14)


cutoffs = c()
for (i in names(GS)){
  print(i)
  cutoffs = c(cutoffs, median(permute_AddModule(mega, 
					gene_set = GS[[i]],
					nPerm=nPerm)
					))
}
names(cutoffs) = names(GS)
saveRDS(cutoffs, file.path("../Data/", out_name))

for ( i in 1:length(cutoffs)){
  feat = paste0(names(cutoffs)[i], i)
  
  nes = mega@meta.data[ ,feat]
  nes = nes - cutoffs[[i]]
  mega@meta.data[ ,feat] = nes
}

FeaturePlot(mega, features = GS_names, min.cutoff = "q01", max.cutoff = "q90")
ggsave(last_plot(), device = "pdf", filename = file.path(plot_dir_path, "UMAP_GeneSets_permuted.pdf"), width = 14, height = 14)

df = as.data.frame( mega@meta.data[ ,GS_names])
df$Celltype = mega$Celltype
ldf = reshape2::melt(df, "Celltype")

sigs = ldf %>% 
  group_by(Celltype, variable) %>% 
  summarise("pos_perc" = sum(value > 0)/n()) %>% 
  filter(pos_perc > .25 )

plot_list = list()
for ( i in unique(ldf$variable)){
  
  colos = setNames( rep("grey",length(unique(ldf$Celltype))), unique(ldf$Celltype) )
  colos[names(colos) %in% sigs$Celltype[sigs$variable == i ]] = "firebrick"
  
  pdf = ldf[ldf$variable == i, ]
g1 = ggplot(pdf, aes(x= Celltype, y = value )) +
  geom_hline(yintercept = 0, col="black", lwd = 1, linetype = "dashed") +
  geom_violin(aes(fill = Celltype), adjust=0.75) +
  scale_fill_manual(values = colos ) +
  ylab("Enrichment score") +
  theme_bw() +
  th  +
  ggtitle(i)

plot_list[[i]] = g1
  
}
cowplot::plot_grid(plotlist = plot_list)
ggsave(last_plot(), device = "pdf", filename = file.path(plot_dir_path, "Violins_GeneSets_permuted.pdf"), width = 12, height = 12)
}



permute_AddModule = function(gene_set, so, nPerm=100){
  set.seed(54321)
  #perm_mat = matrix(nrow = nPerm, ncol = ncol(so), data = 0)
  perm_vec = vector(mode = "list", length = nPerm)
  g_size = sum(gene_set %in% rownames(so))
  elig_genes = rownames(so)
  for (i in 1:nPerm){
    print(i)
    rand_genes = sample(rownames(so), size = g_size, replace = FALSE)
    so = AddModuleScore(so, features = list(rand_genes))
    
    perm_vec[[i]] = quantile( as.numeric(so$Cluster1), 0.95)
    #perm_mat[i, ] = as.numeric(so$Cluster1)
  }
  return(unlist(perm_vec))
}



#Plot_permuted_gene_sets(seurat_path="../Data/Retina_Atlas_Cleaned.rds", 
#                                   gene_set_path="~/Proj/Adult_CB/MCD_PR_UBC_Gene_sets.rds", 
#                                   plot_dir_path="../QC/",
#                                   out_name="Gene_set_permuted_cutoffs.rds")
