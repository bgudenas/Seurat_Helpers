
#Write_count_matrix.R ./Data/Pineal_Atlas_Final.rds ./Data/Pineal_Atlas_Final_counts.txt
args = commandArgs( trailingOnly=TRUE )

if (length(args)==0) {
  stop("Two command line arguments needed -- 1.seurat_object_path 2.output_name_txt", call.=FALSE)
}
library(Seurat)
so = readRDS(args[1])
message("Seurat object loaded")
message(paste0("dimensions=", dim(so)))

message(paste("Outname=", args[2]))

<<<<<<< HEAD

if (args[3] == "data"){
print("normalized data")
counts_df <- so@assays$RNA@data
} else if (grepl("originalexp", names(so@assays))){
=======
if (grepl("originalexp", names(so@assays))){
>>>>>>> 25952928a71987dd6bbce47066433df8306250cd
  counts_df = so@assays$originalexp@counts
} else {
counts_df <- so@assays$RNA@counts
}
<<<<<<< HEAD

=======
>>>>>>> 25952928a71987dd6bbce47066433df8306250cd
fname=args[2]
# make a column vector and write to a file
# mat_frame <- 1:counts_df@Dim[2]
mat_frame = c("Cell", rownames(so))
#mat_frame <- R.utils::insert(mat_frame, 1, c("Cell"))
mat_frame <- matrix(mat_frame, 1)
# colnames(mat_frame) = c("Cell", colnames(so))
write.table(mat_frame, file=fname, row.names = FALSE, col.names=FALSE, sep="\t", quote = FALSE)

# rows per chunk, this could be determined based on cells but just fixed below
# e.g. for 40000 cells I used 1000
elements_per_chunk = 1000
l = split(1:ncol(counts_df), ceiling(seq_along(1:ncol(counts_df))/elements_per_chunk))

# Write large data frame to csv in chunks
count = 1
for(i in l){
  print(count)
  write.table(round(t(as.matrix(counts_df[ ,i ])),2), file=fname, sep = "\t", row.names = TRUE,
              append = TRUE, quote = FALSE, col.names=FALSE)
  count = count + 1
}
print("DONE WRITING")


