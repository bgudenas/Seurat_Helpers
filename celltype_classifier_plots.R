



# model_path="../../Atlas/Data/Xgboost_celltype_model.rds"
# atlas_path="../../Atlas/Data/Pineal_Atlas_Final.rds"
# tumor_path= "../Data/snPB_integrated_all_intron_MTfilt.rds"
# plot_name ="Pineal_celltypes"

classify_tumors = function(model_path,
                           atlas_path,
                           tumors,
                           plot_name,
                           plot_dir="../Figures/",
                           threshold=0.5) {

library(Seurat)
library(stringr)
library(ggplot2)
library(dplyr)
th <- theme_bw() +
  theme(text = element_text(size=16, family = "Helvetica" ), panel.grid.major = element_blank(), panel.grid.minor = element_blank() )
library(xgboost)
source("~/src/Seurat_Helpers/celltype_classifier.R")
map = readRDS("~/Annots/Annotables/Mm10.rds")

## read data
xg_mod = readRDS(model_path)
atlas = readRDS(atlas_path)

mm_feats = xg_mod$feature_names
hs_feats = map$hsapiens_homolog_associated_gene_name[match(mm_feats, map$external_gene_name)]

hs_feats = hs_feats[hs_feats %in% rownames(tumors)]

bin_tumors = make_binary(tumors, pc_genes = hs_feats)

genes_present = map$external_gene_name[match(colnames(bin_tumors), map$hsapiens_homolog_associated_gene_name)]
genes_missing = mm_feats[!(mm_feats %in% genes_present)]

miss_mat = matrix(nrow = nrow(bin_tumors), ncol = length(genes_missing), data = 0) ## manually add in missing genes all as 0 (not expressed)
colnames(miss_mat) = genes_missing

colnames(bin_tumors) = genes_present

bin_tumors = cbind(bin_tumors, miss_mat)
print(quantile(rowSums(bin_tumors)))
bin_tumors = bin_tumors[ ,mm_feats]

preds <- predict(xg_mod, bin_tumors, reshape = TRUE)

annots = data.frame("Cell" = colnames(atlas),
                    "Celltype" = atlas$Celltype)
label_2_num = setNames(as.numeric(as.factor(unique(annots$Celltype)))-1,
                       unique(annots$Celltype))

colnames(preds) = 1:ncol(preds)

for (i in 1:ncol(preds)){
  colnames(preds)[i] =  names(label_2_num)[label_2_num == (i-1)]
}

max_preds = colnames(preds)[apply(preds, 1, which.max)]

celltype_order = c(sort(unique(max_preds)), "Unknown_low_conf")
# too_low = rowSums(bin_tumors) < 30
# max_preds[too_low] = "Unknown_low_conf"
max_score = apply(preds, 1, max)
max_preds[max_score <= threshold] = "Unknown_low_conf"

ldf = as.data.frame(table(max_preds, tumors$Subgroup)) %>% 
  group_by(Var2) %>% 
  summarise("Percent" = Freq/sum(Freq), "Celltype" = max_preds)
ldf$Celltype = factor(ldf$Celltype, levels = celltype_order)


library(RColorBrewer)
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

ggsave(last_plot(), 
       filename = file.path(plot_dir, paste0(plot_name, "_predictions_subgroup.pdf")),
       device="pdf", width = 12, height = 8)

# 
# samp_mat = c()
# for ( i in unique(tumors$ID)){
#  cells = tumors$ID==i
#  samp_preds = colMeans(preds[cells, ])
#  samp_mat = cbind(samp_mat, samp_preds)
# }
# colnames(samp_mat) =  unique(tumors$ID)
# ldf = reshape2::melt(samp_mat)
# ldf$Subgroup = tumors$Subgroup[match(ldf$Var2, tumors$ID)]
# 
# ldf = ldf[order(ldf$Subgroup, ldf$Var2), ]
# 
# ldf$Var2 = factor(ldf$Var2, levels =unique(ldf$Var2))
# ggplot(ldf, aes(x=Var1, y = value)) +
#   geom_bar( aes(fill = Subgroup), stat="identity") +
#   facet_wrap( ~ Var2) +
#   theme_bw() +
#   th +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# 
# ggsave(last_plot(),
#        filename = file.path(plot_dir, paste0(plot_name, "_predictions_sample.pdf")),
#        device="pdf", width = 14, height = 10)

}
