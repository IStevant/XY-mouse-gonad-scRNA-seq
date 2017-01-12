###########################################
#                                         #
#          Load needed libraries          #
#                                         #
###########################################

# GitHub version of Monocle2 Jan2017
# install.packages("devtools")
# devtools::install_github("cole-trapnell-lab/monocle-release@monocle2")
# biocLite(c("DDRTree", "pheatmap"))

library(monocle)

###########################################
#                                         #
#                 Load data               #
#                                         #
###########################################

load(file="../Data/XY_single-cell_RPKM_table.Robj")


###########################################
#                                         #
#            Monocle Analysis             #
#                                         #
###########################################

# Prepare tables for monocle object
conds<- substr(colnames(males), 14,18)
maleCells_expr_matrix <- as.data.frame(males)
maleCells_sample_sheet <- data.frame(cells=names(males), stages=conds, cellType=as.vector(all_males_obj_j2_tree@ident))
rownames(maleCells_sample_sheet)<- names(males)
maleCells_gene_annotation <- as.data.frame(rownames(males))
rownames(maleCells_gene_annotation)<- rownames(males)
colnames(maleCells_gene_annotation)<- "genes"
pd <- new("AnnotatedDataFrame", data = maleCells_sample_sheet)
fd <- new("AnnotatedDataFrame", data = maleCells_gene_annotation)

# First create a CellDataSet from the relative expression levels
maleCells <- newCellDataSet(as.matrix(maleCells_expr_matrix),
phenoData = pd,
featureData = fd)
# Next, use it to estimate RNA counts from RPKM
rpc_matrix <- relative2abs(maleCells)
# Now, make a new CellDataSet using the RNA counts
maleCells <- newCellDataSet(as(as.matrix(rpc_matrix), "sparseMatrix"),
phenoData = pd,
featureData = fd,
lowerDetectionLimit=1,
expressionFamily=negbinomial())

maleCells <- estimateSizeFactors(maleCells)
maleCells <- estimateDispersions(maleCells)
disp_table <- dispersionTable(maleCells)

maleCells <- detectGenes(maleCells, min_expr = 0.1)

#######################################
# Unsupervised Ordering
#######################################

# Use genes from SEURAT to performe cell ordering
ordering_genes <- all_males_obj_j1_opt.sig.genes

maleCells <- setOrderingFilter(maleCells, ordering_genes)
plot_ordering_genes(maleCells)
maleCells <- reduceDimension(maleCells, max_components=4)
maleCells <- reduceDimension(maleCells, reduction_method="ICA")

maleCells <- orderCells(maleCells, reverse=FALSE)

pdf("monocle2_trajectory_pseudotime_ICA_Monocle_genes.pdf")
plot_cell_trajectory(maleCells, color_by="stages")
plot_cell_trajectory(maleCells, color_by="cellType")
dev.off()
