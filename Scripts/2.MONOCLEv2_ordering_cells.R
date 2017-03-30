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
library("RColorBrewer")

###########################################
#                                         #
#                 Load data               #
#                                         #
###########################################

load(file="../Data/XY_single-cell_RPKM_table.Robj")
load(file="../Data/clustered_male_tsne_100000_iter_1-4_600_1-6.Robj")
load(file="../Data/SEURAT_all_cells_sig_genes.Robj")


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
featureData = fd,
lowerDetectionLimit=1,
expressionFamily=negbinomial.size())


# Next, use it to estimate RNA counts from RPKM
rpc_matrix <- relative2abs(maleCells)
# Now, make a new CellDataSet using the RNA counts
maleCells <- newCellDataSet(as(as.matrix(rpc_matrix), "sparseMatrix"),
phenoData = pd,
featureData = fd,
lowerDetectionLimit=1,
expressionFamily=negbinomial.size())

maleCells <- estimateSizeFactors(maleCells)
maleCells <- estimateDispersions(maleCells)
disp_table <- dispersionTable(maleCells)

maleCells <- detectGenes(maleCells, min_expr = 0.1)

#######################################
# Unsupervised Ordering
#######################################


# disp_table <- dispersionTable(maleCells)
# ordering_genes <- subset(disp_table, mean_expression > 0 & dispersion_empirical > 2 * dispersion_fit)$gene_id

# Use genes from SEURAT to performe cell ordering (cf script 1.SEURAT_clustering all_cells.R)
ordering_genes <- all_males_obj_j1_opt.sig.genes
maleCells <- setOrderingFilter(maleCells, ordering_genes)
maleCells <- reduceDimension(maleCells, max_components=3, norm_method ="none")
maleCells <- orderCells(maleCells, reverse=FALSE)
plot_cell_trajectory(maleCells, color_by="cellType")

pdf("lineage.pdf")
plot_cell_trajectory(maleCells, color_by="cellType")
plot_cell_trajectory(maleCells, color_by="stages")
dev.off()


plot_cell_trajectory(maleCells, color_by="State")

maleCells <- orderCells(maleCells, root_state=7)
plot_cell_trajectory(maleCells, color_by="Pseudotime")

progenitor_cells <- colnames(males[,pData(maleCells)$State==7])
save(progenitor_cells, file="../Data/monocle_progenitor_cells_before_branching.Robj")

maleCells_filtered <- detectGenes(maleCells, min_expr = 0.5)
maleCells_filtered <- maleCells_filtered[fData(maleCells_filtered)$num_cells_expressed > 5, ]

BEAM_res <- BEAM(maleCells_filtered, branch_point=1, cores = 3)
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("genes", "pval", "qval")]


save(maleCells_filtered, file="../Data/Monocle_filtered_genes_for_BEAM.Robj")
save(BEAM_res, file="../Data/Monocle_BEAM_results.Robj")

load(file="../Data/Monocle_filtered_genes_for_BEAM.Robj")
load(file="../Data/Monocle_BEAM_results.Robj")

# Color palette for the heatmap
cool = rainbow(50, start=rgb2hsv(col2rgb('cyan'))[1], end=rgb2hsv(col2rgb('blue'))[1])
warm = rainbow(50, start=rgb2hsv(col2rgb('red'))[1], end=rgb2hsv(col2rgb('yellow'))[1])
cols = c(rev(cool), rev(warm))
mypalette <- colorRampPalette(cols)(62)


# Heatmap of the gene expression kinetics
plot_genes_branched_heatmap(maleCells_filtered[row.names(subset(BEAM_res, qval < 0.001)),],
	branch_point = 1,
	num_clusters = 5,
	hclust_method="ward.D",
	cores = 3,
	use_gene_short_name = T,
	show_rownames = T,
	hmcols = mypalette
)
