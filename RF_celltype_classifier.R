Make_RF_Classifier = function(so,
                              out_dir,
                              features,
                              down_prop = NULL,
                              out_name){
  set.seed(54321)
  shhh <- suppressPackageStartupMessages
  shhh(library(Seurat))
  shhh(library(doParallel))
  shhh(library(dplyr))
  shhh(library(scuttle))
  shhh(library(SingleCellExperiment))
  shhh(library(caret))
  shhh(library(parallel))
  shhh(library(randomForest))
  
  annots = data.frame("Cell" = colnames(so),
                      "Celltype" = so$Celltype)
  
  keeps = annots %>% 
    group_by(Celltype) %>% 
    summarise(count = n()) %>% 
    filter(count > 400 ) %>% 
    .$Celltype
  
  annots = annots %>% 
    filter(Celltype %in% keeps)
  annots$Celltype = factor(annots$Celltype)
  ## select 70/30 split of each celltype for training/testing
  set.seed(54321)
  
  train_annots = annots %>% 
    group_by(Celltype) %>% 
    sample_frac( 0.7 ) %>% 
    mutate(bins = ntile(row_number(), 5)) %>% # ntile() splits the data into 5 equal-sized bins
    ungroup()
  
  train_cells = train_annots$Cell
  test_cells = annots$Cell[!(annots$Cell %in% train_cells)]
  
  # Calculate class weights
  class_counts <- table(train_annots$Celltype)
  class_weights <- 1 / (class_counts) ## TRY removing sqrt
  class_weights <- class_weights / sum(class_weights)
  names(class_weights) <- levels(train_annots$Celltype)
  cell_weights = class_weights[match(train_annots$Celltype, names(class_weights))]
    
  if (!is.null(down_prop)){
  print("Downsampling ------------")
  sce = as.SingleCellExperiment(so)
  assay(sce, "counts") <- normalizeCounts(sce, log = FALSE, downsample = TRUE, down.prop = down_prop)
  assay(sce, "logcounts") <- normalizeCounts(sce, log = TRUE, downsample = TRUE, down.prop = down_prop)
  so[["RNA"]]$counts = assay(sce, "counts")
  so[["RNA"]]$data = assay(sce, "logcounts")
  }
  ## make training/testing data
  train_data = round(t(as.matrix(so[["RNA"]]$data[features, train_cells])), 2)
  test_data = round(t(as.matrix(so[["RNA"]]$data[features, test_cells])), 2)
  stopifnot(length(intersect(rownames(train_data), rownames(test_data))) == 0) ## confirm no overlap of cells
  rm(so)
  
  train_celltypes = annots$Celltype[match(train_cells, annots$Cell)]
  test_celltypes = annots$Celltype[match(test_cells, annots$Cell)]
  
  cl <- makeCluster(6)     # Create cluster with number of cores
  registerDoParallel(cl)           # Register the parallel backend
  
  # Set up cross-validation control
  train_control <- trainControl(method = "cv",
                                number = 5,
                                index = createMultiFolds(train_annots$bins, k = 5, times = 1),
                                summaryFunction = multiClassSummary,
                                allowParallel = TRUE)
  
  mtry_top = floor((1/3)*ncol(train_data))*2
  mtry_bottom = floor((1/3)*ncol(train_data))*0.5
  mtry_values = round(seq(from = mtry_bottom, to = mtry_top, length.out = 4), 0)
  tune_grid = data.frame(.mtry = as.integer(mtry_values))
  
  # Train a randomForest model
  print("TRAINING ----------------------")
  model <- train(               
    x = train_data,                   # Dataset
    y = train_celltypes,
    weights = cell_weights,
    ntree= 1500,
    method = "rf",                 # Random Forest method
    trControl = train_control,     # Cross-validation settings
    tuneGrid = tune_grid              # Tune over 6 possible values of mtry
    )
  
  # Stop the cluster after training
  stopCluster(cl)
  
  # Predict on the new dataset
  predictions <- predict(model, newdata = test_data)
  
  conf_mat = caret::confusionMatrix(data = predictions, 
                                    reference = test_celltypes,
  )
  if (is.null(down_prop)) down_prop = "no_DS"
  print(conf_mat)
  print("Writing out data")
  #out_log = file.path(out_dir, paste0("Confusion_matrix", "_dp", down_prop,"_", Sys.Date(),".csv"))
  #write.csv(conf_mat$table, out_log)
  
  output = list("RF_model" = model,
                "Downsample" = down_prop,
                "conf_matrix" = conf_mat)
  out_data = file.path(out_dir, out_name)
  saveRDS(output, out_data)
}
