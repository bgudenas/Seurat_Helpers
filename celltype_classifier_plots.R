
# model_path="../../Atlas/Data/Xgboost_celltype_model.rds"
# atlas_path="../../Atlas/Data/Pineal_Atlas_Final.rds"
# tumor_path= "../Data/snPB_integrated_all_intron_MTfilt.rds"
# plot_name ="Pineal_celltypes"

classify_tumors = function(model_path,
                           atlas_path,
                           tumors,
                           plot_name,
                           plot_dir="../Figures/ML",
                           data_dir="../Data/ML/",
                           threshold=0.5,
                           subgroup_order=NULL,
                           nPerm=10) {
set.seed(54321)
library(Seurat)
library(stringr)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
th <- theme_bw() +
  theme(text = element_text(size=16, family = "Helvetica" ), panel.grid.major = element_blank(), panel.grid.minor = element_blank() )
library(xgboost)
source("~/src/Seurat_Helpers/celltype_classifier.R")
map = readRDS("~/Annots/Annotables/Mm10.rds")

## read data
xg_mod = readRDS(model_path)
atlas = readRDS(atlas_path)

mm_feats = xg_mod$feature_names
## check if human gene names
if (sum(mm_feats %in% map$hsapiens_homolog_associated_gene_name) > round(0.75*length(mm_feats), 0)){
  message("Human gene symbols detected -------")
  hs_feats = mm_feats
} else{
hs_feats = map$hsapiens_homolog_associated_gene_name[match(mm_feats, map$external_gene_name)]
}

hs_feats = hs_feats[hs_feats %in% rownames(tumors)]

bin_tumors = make_binary(tumors, bin_genes = hs_feats)

genes_present = map$external_gene_name[match(colnames(bin_tumors), map$hsapiens_homolog_associated_gene_name)]
genes_missing = mm_feats[!(mm_feats %in% genes_present)]

miss_mat = matrix(nrow = nrow(bin_tumors), ncol = length(genes_missing), data = 0) ## manually add in missing genes all as 0 (not expressed)
colnames(miss_mat) = genes_missing

colnames(bin_tumors) = genes_present

bin_tumors = cbind(bin_tumors, miss_mat)
print(quantile(rowSums(bin_tumors)))

bin_tumors = bin_tumors[ ,mm_feats]  ## arrange in same order

preds <- predict(xg_mod, bin_tumors, reshape = TRUE)

# derive null predictions -------------------------------------------------
message("Calibrating null predictions")
rand_mat = preds
rand_mat[rand_mat > 0 ] = 0 ## set matrix to zero to initialize loop
# TODO testing null input instead of random -- doesnt change much
# null_preds <- predict(xg_mod, rand_mat, reshape = TRUE)
# rand_mat = null_preds

for ( i in 1:nPerm){
tmp_tumors = bin_tumors
old_labels = colnames(tmp_tumors)
tmp_tumors = tmp_tumors[ ,sample(1:ncol(tmp_tumors))]
colnames(tmp_tumors) = old_labels
null_preds <- predict(xg_mod, tmp_tumors, reshape = TRUE)

rand_mat = rand_mat + null_preds
}
rand_mat = rand_mat/nPerm

annots = data.frame("Cell" = colnames(atlas),
                    "Celltype" = atlas$Celltype)

label_2_num = setNames(as.numeric(as.factor(unique(annots$Celltype)))-1,  
                       unique(annots$Celltype))

colnames(preds) = 1:ncol(preds)

for (i in 1:ncol(preds)){
  colnames(preds)[i] =  names(label_2_num)[label_2_num == (i-1)]
}
rownames(preds) = rownames(bin_tumors)

calibrated_preds = preds - rand_mat
calibrated_preds[calibrated_preds < 0 ] = 0

output = list("pred_mat" = preds, "ID" = tumors$ID, "Subgroup" = tumors$Subgroup, "calibrated_mat" = calibrated_preds)
saveRDS(output, paste0( data_dir, plot_name, "_cell_predictions.rds"))


plot_pred_heatmap(pred_list=output,
                             out_plot_dir=plot_dir,
                             out_name=paste0(plot_name, "_heatmap"),
                             threshold=0.6,
                             subgroup_order=subgroup_order)

plot_pred_heatmap(pred_list=output,
                  calibrated = TRUE,
                  out_plot_dir=plot_dir,
                  out_name=paste0(plot_name, "_heatmap_calibrated"),
                  threshold=0.4,
                  subgroup_order=subgroup_order)
  
max_preds = colnames(preds)[apply(preds, 1, which.max)]

output = list("pred_mat" = preds,
              "final_prediction" = max_preds,
              "ID" = tumors$ID,
              "Subgroup" = tumors$Subgroup,
              "calibrated_mat" = calibrated_preds)

saveRDS(output, paste0( data_dir, plot_name, "_cell_predictions.rds"))


plot_pred_heatmap(pred_list=output,
                             out_plot_dir=plot_dir,
                             out_name=paste0(plot_name, "_heatmap"),
                             threshold=0.5,
                             subgroup_order=subgroup_order)

plot_pred_heatmap(pred_list=output,
                  calibrated = TRUE,
                  out_plot_dir=plot_dir,
                  out_name=paste0(plot_name, "_heatmap_calibrated"),
                  threshold=0.4,
                  subgroup_order=subgroup_order)
  
plot_group_pred_heatmap (pred_list=output,
                  calibrated = FALSE,
                  out_plot_dir=plot_dir,
                  out_name=paste0(plot_name, "_Grouped_heatmap"),
                  threshold=0.5,
                  subgroup_order= NULL)

plot_group_pred_heatmap (pred_list=output,
                         calibrated = TRUE,
                         out_plot_dir=plot_dir,
                         out_name=paste0(plot_name, "_Grouped_heatmap_calibrated"),
                         threshold=0.5,
                         subgroup_order= NULL)


celltype_order = c(sort(unique(max_preds)), "Unknown_low_conf")

max_score = apply(preds, 1, max)
max_preds[max_score <= threshold] = "Unknown_low_conf"

ldf = as.data.frame(table(max_preds, tumors$Subgroup)) %>% 
  group_by(Var2) %>% 
  summarise("Percent" = Freq/sum(Freq), "Celltype" = max_preds)
ldf$Celltype = factor(ldf$Celltype, levels = celltype_order)

n <- length(unique(ldf$Celltype)) -1
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

colos = setNames(c(sample(col_vector, n), "grey"),
                  unique(ldf$Celltype))


ggplot(ldf, aes(x=Celltype, y = Percent)) +
  geom_bar( aes(fill = Celltype), stat="identity") +
  scale_fill_manual(values = colos ) +
  facet_wrap( ~ Var2) +
  theme_bw() +
  th +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#ggsave(last_plot(), 
#       filename = file.path(plot_dir, paste0(plot_name, "_predictions_subgroup.pdf")),
#       device="pdf", width = 12, height = 8)
}



# Plot a heatmap per sample of average classifier confidence calls --------

plot_pred_heatmap = function(pred_list,
                             calibrated=FALSE,
                             out_plot_dir="../Figures/ML",
                             out_name,
                             threshold=0.6,
                             subgroup_order=NULL){
  library(pheatmap)
  dir.create(out_plot_dir, showWarnings = FALSE)
  
  if (calibrated == TRUE){ 
    pmat = pred_list$calibrated_mat 
  } else { 
    pmat = pred_list$pred_mat
  }
  
  IDs = pred_list$ID
  subgroup = pred_list$Subgroup
  
  stopifnot(length(IDs) == nrow(pmat))
  
  
  samp_mat = matrix(nrow = length(unique(IDs)),
                    ncol = ncol(pmat),
                    data = 0)
  colnames(samp_mat) = colnames(pmat)
  rownames(samp_mat) = unique(IDs)
  
  for (i in unique(IDs)){
    samp_mat[i, ] = colMeans(pmat[IDs == i, ])
  }
  
  anno_cols = data.frame(row.names =  unique(IDs),
                         "Subgroup" = subgroup[match( unique(IDs), IDs)])
  if (!is.null(subgroup_order)){
    anno_cols$Subgroup = factor(anno_cols$Subgroup, levels = subgroup_order)
  }
  
 # ords = hclust(as.dist(1-cor(t(samp_mat))), method = "ward.D2")
  #ords = hclust(as.dist(t(samp_mat)), method = "ward.D2")
 # samp_mat = samp_mat[order(anno_cols$Subgroup, ords$order), ]
  samp_mat = samp_mat[order(anno_cols$Subgroup), ]
  cols = colorRampPalette(colors = c("blue4","blue","white","red","red4") )(100)
  
  samp_mat[samp_mat > threshold] = threshold
  p1 = pheatmap(samp_mat,
                annotation_row = anno_cols,
                cluster_rows = FALSE,
                cluster_cols = FALSE,
                cellwidth = 10,
                cellheight = 10,
                color = cols )
  out_file = file.path(out_plot_dir, paste0(out_name, ".pdf"))
  
  pdf(out_file, width = 12, height = 10)
  print(p1)
  dev.off()
  message("Heatmap created --------")
  return(p1)
}
# plot_pred_heatmap(pred_data_path = "../Data/ML/Pineal_celltypes_cell_predictions.rds",
#                   out_name="PB_pineal",
#                   subgroup_order=c("PB-miRNA1","PB-miRNA2", "PB-MYC/FOXR2","PB-RB", "PPTID",  
#                                    "PC","PAT", "PTPR")
# )



plot_group_pred_heatmap = function(pred_list,
                                   calibrated=FALSE,
                                   out_plot_dir="../Figures/ML",
                                   out_name,
                                   threshold=0.6,
                                   subgroup_order=NULL){
  library(pheatmap)
  dir.create(out_plot_dir, showWarnings = FALSE)
  
  if (calibrated == TRUE){ 
    pmat = pred_list$calibrated_mat 
  } else { 
    pmat = pred_list$pred_mat
  }
  
  IDs = pred_list$Celltype
  
  stopifnot(length(IDs) == nrow(pmat))
  
  
  samp_mat = matrix(nrow = length(unique(IDs)),
                    ncol = ncol(pmat),
                    data = 0)
  colnames(samp_mat) = colnames(pmat)
  rownames(samp_mat) = unique(IDs)
  
  for (i in unique(IDs)){
    samp_mat[i, ] = colMeans(pmat[IDs == i, ])
  }
  
  anno_cols = data.frame(row.names =  unique(IDs),
                         "Celltype" = unique(IDs)
  )
  # if (!is.null(subgroup_order)){
  #   anno_cols$Subgroup = factor(anno_cols$Subgroup, levels = subgroup_order)
  # }
  # 
  # ords = hclust(as.dist(1-cor(t(samp_mat))), method = "ward.D2")
  #ords = hclust(as.dist(t(samp_mat)), method = "ward.D2")
  # samp_mat = samp_mat[order(anno_cols$Subgroup, ords$order), ]
  # samp_mat = samp_mat[order(anno_cols$Subgroup), ]
  cols = colorRampPalette(colors = c("blue4","blue","white","red","red4") )(100)
  
  samp_mat[samp_mat > threshold] = threshold
  p1 = pheatmap(samp_mat,
                annotation_row = anno_cols,
                cluster_rows = FALSE,
                cluster_cols = FALSE,
                cellwidth = 10,
                cellheight = 10,
                color = cols )
  out_file = file.path(out_plot_dir, paste0(out_name, ".pdf"))
  
  pdf(out_file, width = 12, height = 10)
  print(p1)
  dev.off()
  message("Heatmap created --------")
  return(p1)
}



Pick_Threshold = function(pred_matrix, threshold=0.33){
  pred_matrix[pred_matrix < threshold] = NA
  pick_celltype = apply(pred_matrix, 1, which.max)
  pick_celltype[ lapply(pick_celltype, length) == 0 ] = NA
  celltype_col = unlist(pick_celltype)
  final  = colnames(pred_matrix)[celltype_col]
  final[is.na(final)] = "Unknown"
  return(final)
ggsave(last_plot(), 
       filename = file.path(plot_dir, paste0(plot_name, "_predictions_subgroup.pdf")),
       device="pdf", width = 12, height = 8)
}



# Plot a heatmap per sample of average classifier confidence calls --------

plot_pred_heatmap = function(pred_list,
                             calibrated=FALSE,
                             out_plot_dir="../Figures/ML",
                             out_name,
                             threshold=0.6,
                             subgroup_order=NULL){
  library(pheatmap)
  dir.create(out_plot_dir, showWarnings = FALSE)
  
  if (calibrated == TRUE){ 
    pmat = pred_list$calibrated_mat 
  } else { 
    pmat = pred_list$pred_mat
  }
  
  IDs = pred_list$ID
  subgroup = pred_list$Subgroup
  
  stopifnot(length(IDs) == nrow(pmat))
  
  
  samp_mat = matrix(nrow = length(unique(IDs)),
                    ncol = ncol(pmat),
                    data = 0)
  colnames(samp_mat) = colnames(pmat)
  rownames(samp_mat) = unique(IDs)
  
  for (i in unique(IDs)){
    samp_mat[i, ] = colMeans(pmat[IDs == i, ])
  }
  
  anno_cols = data.frame(row.names =  unique(IDs),
                         "Subgroup" = subgroup[match( unique(IDs), IDs)])
  if (!is.null(subgroup_order)){
    anno_cols$Subgroup = factor(anno_cols$Subgroup, levels = subgroup_order)
  }
  
 # ords = hclust(as.dist(1-cor(t(samp_mat))), method = "ward.D2")
  #ords = hclust(as.dist(t(samp_mat)), method = "ward.D2")
 # samp_mat = samp_mat[order(anno_cols$Subgroup, ords$order), ]
  samp_mat = samp_mat[order(anno_cols$Subgroup), ]
  cols = colorRampPalette(colors = c("blue4","blue","white","red","red4") )(100)
  
  samp_mat[samp_mat > threshold] = threshold
  p1 = pheatmap(samp_mat,
                annotation_row = anno_cols,
                cluster_rows = FALSE,
                cluster_cols = FALSE,
                cellwidth = 10,
                cellheight = 10,
                color = cols )
  out_file = file.path(out_plot_dir, paste0(out_name, ".pdf"))
  
  pdf(out_file, width = 12, height = 10)
  print(p1)
  dev.off()
  message("Heatmap created --------")
  return(p1)
}
# plot_pred_heatmap(pred_data_path = "../Data/ML/Pineal_celltypes_cell_predictions.rds",
#                   out_name="PB_pineal",
#                   subgroup_order=c("PB-miRNA1","PB-miRNA2", "PB-MYC/FOXR2","PB-RB", "PPTID",  
#                                    "PC","PAT", "PTPR")
# )
