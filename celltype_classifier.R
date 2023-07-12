#so =  readRDS("../Data/Pineal_Atlas_Final.rds")

Make_sc_Classifier = function(so,
                              map_path = "~/Annots/Annotables/Mm10.rds",
                              feats_path = "../Data/ML/ML_features.rds",
                              ML_dir ="../Data/ML/",
                              best_param=NULL,
			                        nfeats=100,
			                        down_prop = 0.01,
			                        downsample = FALSE,
                              gene_list_path="/home/bgudenas/Proj/PB_origins/CB_PN_RT_Atlas/Data/genes_expressed/Genes_present_5percent_2datasets.rds"){
                              feats_path = "../Data/ML_features.rds",
                              out_feats = "../Data/ML_features.rds",
                              best_param=NULL){
  
shhh <- suppressPackageStartupMessages
shhh(library(Seurat))
shhh(library(ggplot2))
th <- theme(text = element_text(size=14, family = "Helvetica" ), panel.grid.major = element_blank(), panel.grid.minor = element_blank() )
shhh(library(dplyr))
shhh(library(xgboost))
dir.create(ML_dir, showWarnings = FALSE)

map = readRDS(map_path)
map = map[map$gene_biotype == "protein_coding", ]
#GO = readRDS("~/Annots/GSEA/MisgDB_c2_c5_c6_c8_h1.rds")
# HS = GO[["HSIAO_HOUSEKEEPING_GENES"]]
# mm_HS = map$external_gene_name[match(HS, map$hsapiens_homolog_associated_gene_name)]
# mm_HS = mm_HS[!is.na(mm_HS)]

map = map[!duplicated(map$hsapiens_homolog_associated_gene_name) & !duplicated(map$hsapiens_homolog_associated_gene_name, fromLast = TRUE),  ]
drops = map$external_gene_name[grepl("^Rpl|^Rps|^mt-", map$external_gene_name)] ## remove ribosomal/ mitocondrial genes
map = map[!(map$external_gene_name %in% drops), ]

if (!is.null(gene_list_path)){
 ## filter genes
  gene_list = readRDS(gene_list_path)
  map = map[ map$external_gene_name %in% gene_list | map$hsapiens_homolog_associated_gene_name %in% gene_list, ]
}

pc_genes = c(map$external_gene_name, map$hsapiens_homolog_associated_gene_name)
pc_genes = pc_genes[pc_genes != ""]
#pc_genes = map$external_gene_name[map$hsapiens_homolog_associated_gene_name != ""]
pc_genes = pc_genes[pc_genes %in% rownames(so)]
message(paste0("PC genes loaded = ",length(pc_genes)))

if (downsample == TRUE){
message("Downsampling -----------")
library(scuttle)
library(SingleCellExperiment)
sce = as.SingleCellExperiment(so)
assay(sce, "logcounts") <- normalizeCounts(sce, log = TRUE, downsample = TRUE, down.prop = down_prop)
mega@assays$RNA@data = assay(sce, "logcounts")
}


if (!file.exists(feats_path)){
message("Finding marker genes")
  so@active.ident = factor(so$Celltype)
  feats = FindAllMarkers(so, test.use="LR", only.pos=TRUE, min.diff.pct = 0.1 )
  saveRDS(feats, feats_path) 
} else {
 feats = readRDS(feats_path)
}

feats$pct.diff = feats$pct.1 -  feats$pct.2
feats = feats[feats$p_val_adj <= 0.01 , ]

ct_feats = feats %>%
  filter(gene %in% pc_genes ) %>% 
  group_by(cluster) %>% 
  top_n(wt = pct.diff, nfeats) %>% 
  top_n(wt = avg_log2FC, 100) %>% 
  .$gene %>% 
  unique

message("Number of Features = ", length(ct_feats))

bin_mat = make_binary(so, bin_genes = ct_feats)

#bin_mat = bin_mat[ ,colnames(bin_mat) %in% ct_feats]

annots = data.frame("Cell" = rownames(bin_mat),
                    "Celltype" = so$Celltype)
annots$Celltype[is.na(annots$Celltype)] = "unknown"

label_2_num = setNames(as.numeric(as.factor(unique(annots$Celltype)))-1,
                       unique(annots$Celltype))
annots$Label = label_2_num[match(annots$Celltype, names(label_2_num))]

## select 80% of each celltype for training
set.seed(54321)

train_cells = annots %>% 
  group_by(Celltype) %>% 
  sample_frac( 0.80 ) %>% 
  .$Cell

train_data = bin_mat[train_cells, ]
train_labels = annots$Label[match(train_cells, annots$Cell)]
train_celltypes = annots$Celltype[match(train_cells, annots$Cell)]

test_data = bin_mat[!(rownames(bin_mat) %in% train_cells), ]
test_labels = annots$Label[match(rownames(test_data), annots$Cell)]
test_celltypes = annots$Celltype[match(rownames(test_data), annots$Cell)]

#stopifnot(  all(table(c(test_labels, train_labels)) == table(annots$Celltype)) ) ## make sure celltypes of train + test == original annots
stopifnot( sum(rownames(test_data) %in% rownames(train_data)) == 0 ) ## check there is no overlap

train = list("data" = train_data, "label" = train_labels, "celltype" = train_celltypes)
test = list("data" = test_data, "label" = test_labels, "celltype" = test_celltypes)

model_output = paste0(ML_dir, "/Xgboost_celltype_model_", nfeats, ".rds")

if (!file.exists(model_output)){
  message("Starting XGBoost")
  pred_model = param_sweep_xgboost(train, label_2_num, weight = FALSE, best_param=best_param)
  pred_model = param_sweep_xgboost(train, label_2_num, best_param=best_param)
  saveRDS(pred_model, model_output)
} else {
  pred_model = readRDS(model_output)
}
message("Model finished")

preds <- predict(pred_model, test$data, reshape = TRUE)
colnames(preds) = 1:ncol(preds)

for (i in 1:ncol(preds)){
 colnames(preds)[i] =  names(label_2_num)[label_2_num == (i-1)]
}

max_preds = colnames(preds)[apply(preds, 1, which.max)]

conf_mat = caret::confusionMatrix(data = factor(max_preds, levels=unique(test$celltype)), 
                reference = factor(test$celltype,  levels=unique(test$celltype)),
)
print(conf_mat)

#cut_val = find_cutoff(pred_matrix = preds, reference_labels = test$celltype)

output = list("xgboost_model" = pred_model,
              "Downsample" = downsample,
              "nFeatures_per_celltype" = nfeats,
              "Cutoff"= 0.5,
              "conf_matrix" = conf_mat)
saveRDS(output, paste0(ML_dir, "/XGBoost_testset_output_", nfeats,"_DS_", downsample, ".rds"))

}


make_binary = function(so, q_threshold=0.8, bin_genes){
  set.seed(54321)
  ## derive q_threshold from a subset of genes (pc genes & orthologues for speed)
  map = readRDS("~/Annots/Annotables/Mm10.rds")
  map = map[map$gene_biotype == "protein_coding", ]
  map = map[!duplicated(map$hsapiens_homolog_associated_gene_name) & !duplicated(map$hsapiens_homolog_associated_gene_name, fromLast = TRUE),  ]
  genes = c(map$external_gene_name, map$hsapiens_homolog_associated_gene_name)
  genes = genes[genes != ""]
  genes = genes[genes %in% rownames(so)]
  ## if number of cells is greater than 50000, than randomly select 50k for threshold identification
  if (ncol(so) > 80000 ){
    cell_vec = sample(1:ncol(so), 80000)
  ## if number of cells is greater than 50000, than randomly select 50k for threshold identificaiton
  if (ncol(so) > 50000 ){
    cell_vec = sample(1:ncol(so), 50000)
  } else {cell_vec = 1:ncol(so)}
  
  if (grepl("originalexp", names(so@assays))){
    datExpr = so@assays$originalexp@data[genes ,cell_vec]
    ## find 80th percentile in normalized data (unless it is already saved in )
    threshold = find_threshold(datExpr,
                               q_threshold,
                               saved_threshold="../Data/q80_normData_bin_threshold.rds")
    ncounts = as.matrix(so@assays$originalexp@data[bin_genes, ])
    
  } else {
    datExpr = so@assays$RNA@data[genes ,cell_vec]
    threshold = find_threshold(datExpr,
                               q_threshold,
                               saved_threshold="../Data/q80_normData_bin_threshold.rds")
          ncounts = as.matrix(so@assays$RNA@data[bin_genes, ])
  }
   print(threshold)
   ncounts[ncounts > threshold] = 1
   ncounts[ncounts < 1 ] = 0
   ncounts = t(ncounts)
   return(ncounts)
}

param_sweep_xgboost = function(train, label_2_num, best_param = NULL, weight=TRUE, nround=500){
  
  ## create label weights for imbalanced classes min(class_sizes)/class_sizes
  weight_map = min(table(train$label))/table(train$label)
  weightsData = as.numeric(weight_map[ match(train$label, names(weight_map)) ])
  
param_sweep_xgboost = function(train, label_2_num, best_param = NULL, nround=400){
  
  if (is.null(best_param)){
  message("Starting XGboost grid-search")
  best_param = list()
  best_seednumber = 54321
  best_logloss = Inf
  best_logloss_index = 0
  
  for (iter in 1:50) {
    param <- list(objective = "multi:softprob",
                  eval_metric = "mlogloss",
                  num_class = length(label_2_num),
                  max_depth = sample(3:10, 1), # Typical values: 3-10
                  eta = sample(seq(0.01, 0.2, 0.01), 1),  #Typical values:0.01-0.2
                  gamma = sample(seq(0.0, 0.2, 0.1), 1),
                  scale_pos_weight = 1, #sample(seq(0, 2, 0.5), 1),
                  subsample = sample(seq(0.5, 1, 0.1), 1), #Typical values: 0.5-1
                  colsample_bytree = sample(seq(0.5, 1, 0.1), 1), #Typical values: 0.5-1
                  min_child_weight = sample(1:10, 1),
                  max_delta_step = sample(1:10, 1),
                  nthread=10
    )
    cv.nround = 500
    cv.nfold = 5
    seed.number = sample.int(10000, 1)[[1]]
    set.seed(seed.number)
    
    if (weight == TRUE){
    mdcv <- xgb.cv(data = train$data, 
                   label = train$label,
                   weight = weightsData,
                   params = param,
                   nfold = cv.nfold, 
                   nrounds = cv.nround,
                   verbose = FALSE, 
                   early_stopping_rounds=5, 
                   maximize=FALSE )
    } else if (weight == FALSE) {
      print("Not using class weights ----------------")
      mdcv <- xgb.cv(data = train$data, 
                     label = train$label,
                     params = param,
                     nfold = cv.nfold, 
                     nrounds = cv.nround,
                     verbose = FALSE, 
                     early_stopping_rounds=5, 
                     maximize=FALSE )
    }
    
    min_logloss =  min(mdcv$evaluation_log$test_mlogloss_mean)
    min_logloss_index = which.min(mdcv$evaluation_log$test_mlogloss_mean)
    
    if (min_logloss < best_logloss) {
      best_logloss = min_logloss
      best_logloss_index = min_logloss_index
      best_seednumber = seed.number
      best_param = param
    }
  }
  
  nround = best_logloss_index
  set.seed(best_seednumber)
  print(best_param)
  print("best_round")
  print(nround)
  output = list("best_param" = best_param, "best_round" = nround)
  }
  message("Making Final XGBoost model ---------")
  if (weight == TRUE){
  } 
  message("Making Final XGBoost model ---------")
  md <- xgboost(data=train$data,
                label = train$label,
                weight = weightsData,
		            params=best_param,
 	            	nrounds=nround )
  } else if (weight == FALSE){
    print("Not using class weights ----------------")
    md <- xgboost(data=train$data,
                  label = train$label,
                  params=best_param,
                  nrounds=nround )
    }
return(md)
}


find_threshold = function(datExpr, 
                          q_threshold, 
                          saved_threshold="../Data/q80_normData_bin_threshold.rds") {
  
  if (!file.exists(saved_threshold)){
    message("FInding threshold ----")
    vec = vector(mode = "numeric", length = ncol(datExpr))
    for (i in 1:ncol(datExpr)){
      vec[i] = quantile(as.numeric(datExpr[ ,i]), q_threshold)
    }
    threshold = median(vec)
    saveRDS(threshold, saved_threshold)
  } else {
    message("Saved threshold found ----")
    threshold = readRDS(saved_threshold)
  }
  return(threshold)
}


#### function to find classifier confidence threshold for low-confidence calls
find_cutoff = function(pred_matrix, reference_labels){
  stopifnot(length(reference_labels) == nrow(pred_matrix)) ## check pred matrix matches ref labels
  
  cuts = seq(0,1, 0.01)
  accs = vector(mode = "numeric", length = 100)
  
  for ( i in 1:length(cuts)){
    cut_val = cuts[i]
    pred_matrix[pred_matrix <= cut_val] = NA
    under_conf = rowSums(is.na(pred_matrix)) == ncol(pred_matrix)
    
    max_preds = colnames(pred_matrix)[ unlist(apply(pred_matrix, 1, which.max))]
    max_preds[under_conf] = "Unknown_low_conf"
    
    f_levs = c(unique(reference_labels), "Unknown_low_conf")
    conf_mat = caret::confusionMatrix(data = factor(max_preds, levels=f_levs), 
                                      reference = factor(reference_labels,  levels=f_levs),
    )
    accs[i] = conf_mat$overall[names(conf_mat$overall) == "Accuracy"]
  }
  highest_cut = cuts[which.max((1:length(accs))[accs == max(accs)])]
  
  message(paste0("Highest cut before accuracy dropoff = ", highest_cut))
  return(highest_cut)
}





# 
# feats$pos = 0
# 
# for (i in 1:nrow(feats)){
#   
#   geneName = feats$gene[i]
#   celltype = feats$cluster[i]
#   otx2 = ifelse(atlas@assays$RNA@data[geneName, ] > 0, 1, 0)
#   ratio = sum(otx2[atlas$Celltype == celltype ])/sum(otx2[atlas$Celltype != celltype ])
#   feats$pos[i] = ratio
# }
