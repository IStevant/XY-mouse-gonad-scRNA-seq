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

# Use genes from SEURAT to performe cell ordering (cf script 1.SEURAT_clustering all_cells.R)
ordering_genes <- all_males_obj_j1_opt.sig.genes

maleCells <- setOrderingFilter(maleCells, ordering_genes)

maleCells <- reduceDimension(maleCells, max_components=4)
maleCells <- orderCells(maleCells, reverse=FALSE)
plot_cell_trajectory(maleCells, color_by="stages")
plot_cell_trajectory(maleCells, color_by="cellType")


plot_cell_trajectory(maleCells, color_by="State")
maleCells <- orderCells(maleCells, root_state=6)
plot_cell_trajectory(maleCells, color_by="Pseudotime")


maleCells_expressed_genes <-  row.names(subset(fData(maleCells), num_cells_expressed >= 3))
maleCells_filtered <- maleCells[maleCells_expressed_genes,]
my_genes <- row.names(subset(fData(maleCells_filtered), genes %in% c("Wt1", "Pdgfra", "Trim71")))
cds_subset <- maleCells_filtered[my_genes,]
plot_genes_in_pseudotime(cds_subset, color_by="stages")


testis_genes <- row.names(subset(fData(maleCells), genes %in% c("Pdgfra", "Sox9", "Wnt5a", "Wt1", "Mro", "Dmrt1", "Foxp1", "Gadd45g")))

plot_genes_branched_pseudotime(maleCells[testis_genes,],
branch_point=2,
color_by="cellType",
ncol=2)


BEAM_res <- BEAM(maleCells, branch_point=2, cores = 3)
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]

diff_test_res <- differentialGeneTest(maleCells[expressed_genes,],fullModelFormulaStr="~stages")
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))


maleCells <- setOrderingFilter(maleCells, ordering_genes)
estimateDispersions(maleCells)
plot_ordering_genes(maleCells)

maleCells <- reduceDimension(maleCells, max_components=2)

maleCells <- orderCells(maleCells, reverse=FALSE)

plot_cell_trajectory(maleCells, color_by="stages")

