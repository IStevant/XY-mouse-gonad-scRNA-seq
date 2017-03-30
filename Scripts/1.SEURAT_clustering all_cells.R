#############################################################################################################################
# IMPORTANT : the SEURAT in the paper was produced with a version which is not available anymore...
# The script below was written for this version.
# I propose to install the closest version available in github (install_github("satijalab/seurat@da6cd08"))).
# With this version, the results don't change, the clustering remain the same, but the t-SNE does not look exactly the same than in the paper...
# You will have to adapt the "guse" parametter for the clustering (5 instead of 11) to obtain a similar result.
# I saved the Robj from the original clustering to reproduce the paper figure (../Data/clustered_male_tsne_100000_iter_1-4_600_1-6.Robj)
#############################################################################################################################


###########################################
#                                         #
#          Load needed libraries          #
#                                         #
###########################################

# library(devtools)
# install_github("satijalab/seurat@da6cd08")
# library(Seurat) 


# load library to une rep_along
library("taRifx")
# load Seurat
library("Seurat")
# load viridis color palette
library("viridis")

# modified SEURAT functions to fix bugs and add features
source("1.Heatmap_function_SEURAT_for_clustering all_cells.R")

###########################################
#                                         #
#                 Load data               #
#                                         #
###########################################

load(file="../Data/XY_single-cell_RPKM_table.Robj")

###########################################
#                                         #
#             Seurat Analysis             #
#                                         #
###########################################

# Log normal transformation
all_males_data=log(males+1)

# Create SEURAT object
all_males_obj=new("seurat",raw.data=all_males_data)

# Genes considered as expressed if detected in at least 3 cells
# Condition name in the 3rd field of the sample name
# Cells with at least 1000 expressed genes
# Expression threshold to 1 RPKM
all_males_obj=setup(all_males_obj,project="all_males",min.cells = 3,names.field = 3,names.delim = "_",min.genes = 1000,is.expr=1)

# Select genes with variance > 2 and mean expression >2
all_males_obj=mean.var.plot(all_males_obj,y.cutoff=2,x.low.cutoff=2,fxn.x = expMean,fxn.y = logVarDivMean)
length(all_males_obj@var.genes)

# Run pca on this subset of genes
all_males_obj=pca(all_males_obj,do.print=FALSE)

# Evaluate which PCs contain informative genes
all_males_obj_j1=jackStraw(all_males_obj,num.replicate = 1000, do.print = FALSE)

# Plot tests and p-values for PCs 1 to 9
# For more details, see https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4325543/
jackStrawPlot(all_males_obj_j1, PCs = 1:9)

# Project the all the genes in the previous PCA space
all_males_obj_j1_opt=project.pca(all_males_obj_j1,do.print=FALSE)

# Select the 600 first genes from the PCs selected as significant by Jackstraw
all_males_obj_j1_opt.sig.genes=pca.sig.genes(all_males_obj_j1_opt,c(1:4),pval.cut = 1e-5, max.per.pc=600)

# 1379 genes selected
length(all_males_obj_j1_opt.sig.genes)
save(all_males_obj_j1_opt.sig.genes, file="../Data/SEURAT_all_cells_sig_genes.Robj")
save(all_males_obj_j1_opt.sig.genes, file="../Data/SEURAT_all_cells_sig_genes.Robj")

# Run PCA on these 1379 genes
all_males_obj_j1_opt=pca(all_males_obj_j1_opt,pc.genes=all_males_obj_j1_opt.sig.genes, do.print = FALSE)

# Evaluate which PCs contain informative genes
all_males_obj_j1_opt=jackStraw(all_males_obj_j1_opt,num.replicate = 400,do.print = FALSE)

# Plot tests and p-values for PCs 1 to 19
jackStrawPlot(all_males_obj_j1_opt, PCs = 1:19)

# t-SNE on the 6 first PCs, with 100K iteration (takes hours to compute)
all_males_obj_j2_tsne=run_tsne(all_males_obj_j1_opt,dims.use = c(1:6), max_iter=100000)

# faster t-SNE, will produce more quickly a graph that is slightly different in shape but clusters should be the same
# all_males_obj_j2_tsne=run_tsne(all_males_obj_j1_opt,dims.use = c(1:6), max_iter=10000, do.fast=TRUE)

# Plot t-SNE -> Raw graph used to generate Figure 2a of the paper
tsne.plot(all_males_obj_j2_tsne,pt.size = 4)

# Save/Load object to analyse later
# save(all_males_obj_j2_tsne,file="../Data/male_tsne_100000_iter_1-4_600_1-6.Robj")
# Original clustering (cf comment on top of the script)
load(file="../Data/male_tsne_100000_iter_1-4_600_1-6.Robj")

# Define cell clusters with DBscan using 11 as epsilon neighborhood
# G.use is 11 for the for the original clustering, but will be 5 with the available SEURAT version I propose to install (cf comment on top of the script)
all_males_obj_j2_tsne_cl=DBclust_dimension(all_males_obj_j2_tsne,1,2,reduction.use = "tsne",G.use = 11,set.ident = TRUE)
tsne.plot(all_males_obj_j2_tsne_cl,pt.size = 3)

# Ordering clusters by similarities
all_males_obj_j2_tree=buildClusterTree(all_males_obj_j2_tsne_cl,do.reorder = TRUE,reorder.numeric = TRUE,pcs.use = c(1:3))

# t-sne plot colored by final clusters -> Used for Figure 2a of the paper
tsne.plot(all_males_obj_j2_tree,do.label = TRUE,label.pt.size = 1.5)

# Save/Load object to analyse later
save(all_males_obj_j2_tree,file="../Data/clustered_male_tsne_100000_iter_1-4_600_1-6.Robj")
load(file="../Data/clustered_male_tsne_100000_iter_1-4_600_1-6.Robj")

# Calculate differentially expressed genes between clusters using Tobit
markers.all=find_all_markers(all_males_obj_j2_tree, do.print = TRUE, test.use="tobit")

# Correct p-values for False Discovery Rate
FDR.adjusted.pvals <- p.adjust(markers.all$p_val, method = "BH")
markers.all.adj_pval <- data.frame(markers.all, FDR=FDR.adjusted.pvals)

# Keep only over-expressed genes -> Dataset 1 of the paper
markers.use=subset(markers.all.adj_pval,FDR<0.05 & avg_diff>0)$gene

# Save/Load object to analyse later
save(markers.use,file="../Data/DE_genes_391_cells_4_clusters.Robj")
load(file="../Data/DE_genes_391_cells_4_clusters.Robj")

# Heatmap of the over-expressed genes (scaled expression)
males_data <- doHeatMap(all_males_obj_j2_tree,genes.use = markers.use,slim.col.label = TRUE,remove.key = TRUE,cexRow=0.5, do.return=TRUE)

# Heatmap with no-scaled data and embryonic stages on top of the heatmap
males_data <- doHeatMap_2(all_males_obj_j2_tree,genes.use = markers.use,slim.col.label = TRUE,remove.key = TRUE,cexRow=0.5, do.return=TRUE)

###########################################
#                                         #
#           Additional graphs             #
#                                         #
###########################################

# heatmap of the marker genes -> Figure 2d

marker_genes <- c(
	"Trim71",
	"Tbx18",
	"Hoxb6",
	"Hoxa5",
	"Nr2f1",
	"Lhx9",
	"Tcf21",
	"Pdgfra",
	"Wnt5a",
	"Arx",
	"Gli1",
	"Sry",
	"Gadd45g",
	"Mro",
	"Aard",
	"Amh",
	"Dhh",
	"Hsd17b3",
	"Cyp11a1",
	"Cyp17a1",
	"Star",
	"Insl3",
	"Jack3"
)

males_data <- doHeatMap_2(all_males_obj_j2_tree,genes.use = marker_genes,slim.col.label = TRUE,remove.key = TRUE,cexRow=0.5, do.return=TRUE)


# t-SNE colored by gene expression level (Sup. Fig. 6)

tsne_gene_exp <- function(tsne_result,gene){
	p <- ggplot(tsne_result, aes(tSNE_1,tSNE_2)) +
	geom_point(shape = 21, stroke=0.5, aes(fill=as.numeric(log(males[gene,]+1))), size = 3) +
	# scale_fill_gradient2(high="darkred", low="yellow")+
	scale_fill_gradient(low = "white", high = "darkred")+
	theme_bw() +
	ggtitle(gene) +
	theme(
		plot.title = element_text(size=18, face="bold.italic"),
		axis.text=element_text(size=12),
		axis.title=element_text(size=16),
		legend.text = element_text(size =12),
		legend.title=element_blank()
	)
	print(p)

}


seurat_tsne <- all_males_obj_j2_tsne@tsne.rot

tsne_gene_exp(seurat_tsne,"Trim71")
tsne_gene_exp(seurat_tsne,"Wnt5a")
tsne_gene_exp(seurat_tsne,"Aard")
tsne_gene_exp(seurat_tsne,"Insl3")



# Plot individual genes (violin)

males <- log(males+1)

plot_gene_violin <- function(data, title) {
	mean_exp <- data
	clusters <- substr(names(data),14,18)
	plot_mean <- data.frame(
			t(mean_exp),
			clusters=clusters
		)
	colnames(plot_mean) <- c("expr", "clusters")
	g <- ggplot(plot_mean, aes(clusters, expr)) +
	geom_violin(aes(fill=factor(clusters)), adjust = .6, scale = "width") + 
	geom_jitter(height = 0, size=0.5) +
	stat_smooth(aes(x = clusters, y = as.numeric(as.character(expr)), group=1), method="loess", colour="black", se=FALSE) + 
	theme_bw() +
	xlab("Cells") +
	ylab("log(RPKM+1)") +
	ggtitle(title) +
	scale_colour_manual(
	values="black"
	) +
	scale_fill_manual(
	values=stagePalette,
	name=""
	) +
	theme(
	axis.text=element_text(size=18),
	axis.title=element_text(size=18),
	legend.text = element_text(size =18),
	legend.title = element_text(size =18 ,face="bold"),
	plot.title = element_text(size=18, face="bold.italic", hjust = 0.5),
	axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank()
	)
	print(g)
}


males_stages <- as.vector(substr(all_males_obj_j2_tsne@ident,1,5))
stagePalette <- c("#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c")
# stagePalette <- c("#abd9e9", "#ffffbf", "#fdae61", "#d7191c")


graph_violin <- function(genes){
	for (gene in genes){
		print(gene)
		subset_genes <- males[gene,]
		plot_gene_violin(subset_genes, gene)
	}
}

load(file="../Data/XY_single-cell_RPKM_table.Robj")
males <- log(males+1)

graph_violin(c(
	"Arx"
	)
)