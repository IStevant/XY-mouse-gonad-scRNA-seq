###########################################
#                                         #
#          Load needed libraries          #
#                                         #
###########################################

# GitHub version of Monocle2 Jan2017
# install.packages("devtools")
# devtools::install_github("cole-trapnell-lab/monocle-release@monocle2")
# biocLite(c("DDRTree", "pheatmap"))

library("monocle")

###########################################
#                                         #
#                 Load data               #
#                                         #
###########################################

dataset <- read.csv(file="../Data/sample_test_monocle2.csv", row.names=1,check.names=FALSE)

###########################################
#                                         #
#            Monocle Analysis             #
#                                         #
###########################################

# Prepare tables for monocle object
conds<- substr(colnames(dataset), 14,18)
sampleCells_expr_matrix <- as.data.frame(dataset)
sampleCells_sample_sheet <- data.frame(cells=names(dataset), stages=conds))
rownames(sampleCells_sample_sheet)<- names(dataset)
sampleCells_gene_annotation <- as.data.frame(rownames(dataset))
rownames(sampleCells_gene_annotation)<- rownames(dataset)
colnames(sampleCells_gene_annotation)<- "genes"
pd <- new("AnnotatedDataFrame", data = sampleCells_sample_sheet)
fd <- new("AnnotatedDataFrame", data = sampleCells_gene_annotation)

# First create a CellDataSet from the relative expression levels
sampleCells <- newCellDataSet(as.matrix(sampleCells_expr_matrix),
phenoData = pd,
featureData = fd)
# Next, use it to estimate RNA counts from RPKM
rpc_matrix <- relative2abs(sampleCells)
# Now, make a new CellDataSet using the RNA counts
sampleCells <- newCellDataSet(as(as.matrix(rpc_matrix), "sparseMatrix"),
phenoData = pd,
featureData = fd,
lowerDetectionLimit=1,
expressionFamily=negbinomial())

sampleCells <- estimateSizeFactors(sampleCells)
sampleCells <- estimateDispersions(sampleCells)
disp_table <- dispersionTable(sampleCells)

sampleCells <- detectGenes(sampleCells, min_expr = 0.1)

#######################################
# Unsupervised Ordering
#######################################


disp_table <- dispersionTable(sampleCells)
ordering_genes <- subset(disp_table, mean_expression >0 & dispersion_empirical > 2 * dispersion_fit)$gene_id

sampleCells <- setOrderingFilter(sampleCells, ordering_genes)
sampleCells <- reduceDimension(sampleCells)

sampleCells <- orderCells(sampleCells, reverse=FALSE)

plot_cell_trajectory(sampleCells, color_by="stages")

BEAM_res <- BEAM(sampleCells, branch_point=1, cores = 1)
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("genes", "pval", "qval")]

plot_genes_branched_heatmap(sampleCells[row.names(subset(BEAM_res, qval < 1e-4)),],
branch_point = 1,
num_clusters = 4,
cores = 1,
use_gene_short_name = T,
show_rownames = T)

