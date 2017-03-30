#############################################################################################################################
# IMPORTANT : The results were generated with the available github SEURAT version ("satijalab/seurat@da6cd08").
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
library("amap")
library("miscTools")
library("reshape2")
library("ggplot2")


# modified SEURAT functions to fix bugs and add features
source("1.Heatmap_function_SEURAT_for_clustering all_cells.R")

###########################################
#                                         #
#                 Load data               #
#                                         #
###########################################

load(file="../Data/XY_single-cell_RPKM_table.Robj")
load(file="../Data/clustered_male_tsne_100000_iter_1-4_600_1-6.Robj")
# progenitor cells before branching (monocle analysis)
load(file="../Data/monocle_progenitor_cells_before_branching.Robj")

###########################################
#                                         #
#             Seurat Analysis             #
#                                         #
###########################################

# Select sertoli cells from cluster 3 (seurat analysis with all the cells)
sertoli_cells <- names(all_males_obj_j2_tree@ident[all_males_obj_j2_tree@ident==3])
# Regroup progenitor cells from monocle and Sertoli cells from Seurat analysis
sertoli_cell_lineage <- unique(c(progenitor_cells, sertoli_cells))
# Expression matrix of the selected cells
Sertoli_cells <- males[,sertoli_cell_lineage]
# log-transform the datta
sertoli_data=log(Sertoli_cells+1)

# Create SEURAT object
sertoli_obj=new("seurat",raw.data=sertoli_data)

# Genes considered as expressed if detected in at least 3 cells
# Condition name in the 3rd field of the sample name
# Cells with at least 1000 expressed genes
# Expression threshold to 1 RPKM
sertoli_obj=setup(sertoli_obj,project="sertoli",min.cells = 3,names.field = 3,names.delim = "_",min.genes = 1000,is.expr=1,) 

# Select genes with variance > 1.5 and mean expression >1
sertoli_obj=mean.var.plot(sertoli_obj,y.cutoff = 1.5,x.low.cutoff = 1,fxn.x = expMean,fxn.y = logVarDivMean)
length(sertoli_obj@var.genes)

# Run pca on this subset of genes
sertoli_obj=pca(sertoli_obj, do.print = FALSE)

# Evaluate which PCs contain informative genes
sertoli_obj_j1=jackStraw(sertoli_obj,num.replicate = 1000, do.print = FALSE)

# Plot tests and p-values for PCs 1 to 9
# For more details, see https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4325543/
jackStrawPlot(sertoli_obj_j1,PCs = 1:12)

# Project the all the genes in the previous PCA space
sertoli_obj_j1_opt=project.pca(sertoli_obj_j1,do.print=FALSE)

# Select the 600 first genes from the PCs selected as significant by Jackstraw
sertoli_obj_j1_opt.sig.genes=pca.sig.genes(sertoli_obj_j1_opt,c(1:3),pval.cut = 1e-5, max.per.pc=600)
length(sertoli_obj_j1_opt.sig.genes)

# Run PCA on these genes
sertoli_obj_j1_opt=pca(sertoli_obj_j1_opt,pc.genes=sertoli_obj_j1_opt.sig.genes, do.print = FALSE)

# Evaluate which PCs contain informative genes
sertoli_obj_j2=jackStraw(sertoli_obj_j1_opt,num.replicate = 400,do.print = FALSE)
jackStrawPlot(sertoli_obj_j2,PCs = 1:10)

# t-SNE on the 4 first PCs, with 5K iteration (takes hours to compute)
sertoli_obj_j1_tsne=run_tsne(sertoli_obj_j1_opt,dims.use = c(1:4), max_iter=5000, perplexity=25)

# Plot t-SNE -> Raw graph used to generate Figure 4a of the paper
tsne.plot(sertoli_obj_j1_tsne,pt.size = 4)

# Save/Load object to analyse later
save(sertoli_obj_j1_tsne,file="../Data/sertoli_e11_tsne_5000_iter_1-3_600_1-4_p25.Robj")
load("../Data/sertoli_e11_tsne_5000_iter_1-3_600_1-4_p25.Robj")

# Define cell clusters with DBscan using 2 as epsilon neighborhood
sertoli_obj_j2_tsne_cl=DBclust_dimension(sertoli_obj_j1_tsne,1,2,reduction.use = "tsne",G.use =2,set.ident = TRUE)
tsne.plot(sertoli_obj_j2_tsne_cl,pt.size = 4)

# Ordering groups by similarities
sertoli_obj_j2_tree=buildClusterTree(sertoli_obj_j2_tsne_cl,do.reorder = TRUE,reorder.numeric = TRUE,pcs.use = c(1:5))

# switch branches according to embryonic stages
ident <- sertoli_obj_j2_tree@ident
levels(ident) <- c(3,4,5,1,2)
ident = factor(ident,levels(ident)[c(4, 5, 1, 2, 3)])
sertoli_obj_j2_tree@ident <- ident

# t-sne plot colored by final clusters -> Used for Figure 4a of the paper
tsne.plot(sertoli_obj_j2_tree,do.label = TRUE,label.pt.size = 1.5)


# Calculate differentially expressed genes between clusters using Tobit
markers.all=find_all_markers(sertoli_obj_j2_tree, do.print = TRUE, test.use="tobit")

# Correct p-values for False Discovery Rate
markers.all <- 
	cbind(
		markers.all,
		FDR=p.adjust(markers.all$p_val, method = "BH")
	)


# Save/Load object to analyse later
save(markers.all, file="../Data/DE_genes_Sertoli_cells.Robj")
load(file="../Data/DE_genes_Sertoli_cells.Robj")

# Keep only over-expressed genes
markers.exp=subset(markers.all,FDR<0.05)
markers.exp=subset(markers.exp,avg_diff>1.7)

# Select the top specific genes for each group
genes.viz_1 <- head(markers.exp[markers.exp$cluster==1,"gene"],50)
genes.viz_2 <- head(markers.exp[markers.exp$cluster==2,"gene"],50)
genes.viz_3 <- head(markers.exp[markers.exp$cluster==3,"gene"],50)
genes.viz_4 <- head(markers.exp[markers.exp$cluster==4,"gene"],50)
genes.viz_5 <- head(markers.exp[markers.exp$cluster==5,"gene"],50)


markers.use <- c(
	as.character(genes.viz_1),
	as.character(genes.viz_2),
	as.character(genes.viz_3),
	as.character(genes.viz_4),
	as.character(genes.viz_5)
	)

# Heatmap with no-scaled data and embryonic stages on top of the heatmap (Figure 4b)
sertoli_heatmap <- doHeatMap_2(sertoli_obj_j2_tree,genes.use = markers.use,slim.col.label = TRUE,remove.key = TRUE,cexRow=0.5, do.return=TRUE)


###########################################
#                                         #
#       DE genes pattern clustering       #
#                                         #
###########################################

markers.exp=subset(markers.all,FDR<0.01)

cell_cluster <- sertoli_obj_j2_tree@ident

# Get cells from each clusters
cell_cluster1 <- names(cell_cluster[cell_cluster==1])
cell_cluster2 <- names(cell_cluster[cell_cluster==2])
cell_cluster3 <- names(cell_cluster[cell_cluster==3])
cell_cluster4 <- names(cell_cluster[cell_cluster==4])
cell_cluster5 <- names(cell_cluster[cell_cluster==5])

# Mean gene expression per cluster
# exp_cluster1 <- rowMeans(log(males[,colnames(males) %in% cell_cluster1]+1))
# exp_cluster2 <- rowMeans(log(males[,colnames(males) %in% cell_cluster2]+1))
# exp_cluster3 <- rowMeans(log(males[,colnames(males) %in% cell_cluster3]+1))
# exp_cluster4 <- rowMeans(log(males[,colnames(males) %in% cell_cluster4]+1))
# exp_cluster5 <- rowMeans(log(males[,colnames(males) %in% cell_cluster5]+1))

exp_cluster1 <- rowMeans(males[,colnames(males) %in% cell_cluster1])
exp_cluster2 <- rowMeans(males[,colnames(males) %in% cell_cluster2])
exp_cluster3 <- rowMeans(males[,colnames(males) %in% cell_cluster3])
exp_cluster4 <- rowMeans(males[,colnames(males) %in% cell_cluster4])
exp_cluster5 <- rowMeans(males[,colnames(males) %in% cell_cluster5])

# Generate a matrix of gene mean expression per cluster
mean_by_cluster <- 	cbind(
	cluster1=exp_cluster1,
	cluster2=exp_cluster2,
	cluster3=exp_cluster3,
	cluster4=exp_cluster4,
	cluster5=exp_cluster5
	)

mean_by_cluster_DE_genes <- log(mean_by_cluster[rownames(mean_by_cluster) %in% rownames(markers.exp),]+1)
# mean_by_cluster_DE_genes <- mean_by_cluster[rownames(mean_by_cluster) %in% rownames(markers.all),]


write.table(mean_by_cluster_DE_genes, "/home/zazooo/mean_DE_genes_by_cluster.txt", sep="\t")
mean_by_cluster_DE_genes <- read.table("/home/zazooo/ownCloud/PhD/scRNAseq/Analysis/seq_1+2/Male/Sertoli_clustering_170206/mean_DE_genes_by_cluster.txt", sep="\t")

clPalette <- c("#303633", "#8be8cb", "#7ea2aa", "#888da7", "#9c7a97")

gene_clusters <- function(data, clusterNb){
	# Run clustering
	d <- Kmeans(data, centers=clusterNb, iter.max=100000, method="pearson")

	# Open a pdf to record the graphs
	pdf("graph_file.pdf", width=3, height=3)

	# For each cluster, generate the graphs
	for (cluster in c(1:clusterNb)) {
		print(paste("Compute cluster ", cluster))
		group <- names(d$cluster[d$cluster==cluster])
		exp_group <- data[rownames(data) %in% group,]

		# format the K-means output for ggplot. Melt call the columns Var1 (genes), Var2 (conditions), value (expression value per gene)
		melt_dataframe <- melt(exp_group)
		# rename columns for clarity
		colnames(melt_dataframe) <- c("genes", "conditions", "value")
		# in my data, conditions are "cluster1", "cluster2"..., so I pick "1", "2" and make it as numeric to get a continous x axis (easier to deal with)
		melt_dataframe$conditions <- as.numeric(substr(melt_dataframe$conditions,8,9))

		# I calculate the median of expression of the genes per condition to draw it on the graph
		exp_group_mean <- melt(colMedians(exp_group))
		# transform conditions ("cluster1", "cluster2"...) as numeric values
		exp_group_mean$Var1 <- as.numeric(substr(exp_group_mean$Var1,8,9))
		# rename columns for clarity
		colnames(exp_group_mean) <- c("conditions", "value")

		# plot (jitter plot) gene expression per condition with the median to highlight the trajectory
		g <-ggplot() +
		geom_jitter(data = melt_dataframe, aes(x = conditions, y = value, group = conditions, color=factor(conditions)), size=2.5)+
		theme_bw()+
		geom_line(data = exp_group_mean, aes(x=conditions, y=value), color="black", size=1)+
		xlab(" ")+
		ylab(" ")+
		scale_x_continuous(breaks=c(1:5), labels=c("Prog.", "Sc I", "Sc II", "Sc III", "Sc IV")) +
		scale_color_manual(
		values=alpha(clPalette, 0.7),
		name=""
		) +
		ggtitle(paste(length(rownames(exp_group)), " genes", sep="")) +
		theme(
			axis.text=element_text(size=18),
			axis.title=element_text(size=18),
			legend.text = element_text(size =18),
			legend.title = element_text(size =18 ,face="bold"),
			plot.title = element_text(size=18, face="bold", hjust = 0.5),
			legend.position="none",
			axis.text.x = element_text(angle = 45, hjust = 1)
		)

		# write file containing gene name contained in the current k-mean cluster
		write.table(rownames(exp_group), paste("genes_cluster_", cluster,".txt", sep=""))

		# If needed, export graph one by one as png
		# png(paste("path/to/your/folder/genes_cluster_", cluster,".png", sep=""), width = 300, height = 280)
		print(g)
		# dev.off()

		}

		# Close pdf with all the graphs
		dev.off()
}


gene_clusters(as.matrix(mean_by_cluster_DE_genes), 15)





















clPalette <- c("#dbdbdb", "#c269ff", "#9541ff", "#5300db", "#260070")

stagePalette <- c("#303633", "#8be8cb", "#7ea2aa", "#888da7", "#9c7a97")






























plot_genes <- function(data) {
	mean_exp <- data
	clusters <- substr(colnames(data),23,24)
	# print(colnames(data))
	# clusters <- early_clusters

	# print(seq(names(data)))
	plot_mean <- melt(data)
	plot_mean <- cbind(
		plot_mean,
		cluster=substr(plot_mean$Var2,25,25)
	)

	# print(plot_mean)

	# plot_mean$cells <- factor(plot_mean$cells, levels = plot_mean$stages)

	# print(colnames(data))

	g <- ggplot(plot_mean, aes(cluster, as.numeric(as.character(value)))) +
	# ggplot(plot_mean, aes(as.numeric(cells), as.numeric(as.character(mean)))) +
	# geom_point(shape=21, size=3.5, stroke=0.5, aes(fill=factor(stages))) + 
	geom_violin(aes(fill=factor(cluster)), adjust = 1, scale = "width") + 
	# geom_jitter(height = 0, size=1) +
	geom_jitter(height = 0, size=0.5) +
	stat_smooth(aes(x = cluster, y = as.numeric(as.character(value)), group=1), method="loess", colour="black", se=FALSE) + 
	theme_bw() +
	xlab(" ") +
	ylab("log(RPKM+1)") +
	ggtitle(unique(plot_mean$Var1)) +
	scale_colour_manual(
	values="black"
	) +
	scale_fill_manual(
	values=stagePalette,
	name=""
	) +
	scale_x_discrete(labels=c("Prog.", "Sc I", "Sc II", "Sc III", "Sc IV")) +
	theme(
	axis.text=element_text(size=18),
	axis.title=element_text(size=18),
	legend.text = element_text(size =18),
	legend.title = element_text(size =18 ,face="bold"),
	plot.title = element_text(size=18, face="bold.italic", hjust = 0.5),
	legend.position="none",
	axis.text.x = element_text(angle = 45, hjust = 1)
)

	# stat_smooth() +
	# geom_point(aes(colour = factor(clusters))) +
	# ylim(0,4.5)

	# ggplot(plot_mean, aes(as.numeric(as.character(cells)), as.numeric(as.character(mean)))) +
	# stat_smooth() +
	# geom_point(aes(colour = factor(clusters)))

		# png(paste("/home/zazooo/",unique(plot_mean$Var1),".png", sep=""), width = 280, height = 280)
		# print(g)
  # 		dev.off()

		pdf(paste("/home/zazooo/",unique(plot_mean$Var1),".pdf", sep=""), width = 3.5, height = 3.5)
		print(g)
  		dev.off()
}

sertoli_data <- sertoli_obj_j2_tree@raw.data
subset_genes <- sertoli_data["Tex15",]
colnames(subset_genes) <- paste(colnames(subset_genes), sertoli_obj_j2_tree@ident, sep="_")
males_subset <- as.matrix(na.omit(subset_genes))
plot_genes(males_subset)



subset_genes <- sertoli_data["Wnt4",]
colnames(subset_genes) <- paste(colnames(subset_genes), sertoli_obj_j2_tree@ident, sep="_")
males_subset <- as.matrix(na.omit(subset_genes))
plot_genes(males_subset)

subset_genes <- sertoli_data["Sry",]
colnames(subset_genes) <- paste(colnames(subset_genes), sertoli_obj_j2_tree@ident, sep="_")
males_subset <- as.matrix(na.omit(subset_genes))
plot_genes(males_subset)


subset_genes <- sertoli_data["Sox9",]
colnames(subset_genes) <- paste(colnames(subset_genes), sertoli_obj_j2_tree@ident, sep="_")
males_subset <- as.matrix(na.omit(subset_genes))
plot_genes(males_subset)


subset_genes <- sertoli_data["Wt1",]
colnames(subset_genes) <- paste(colnames(subset_genes), sertoli_obj_j2_tree@ident, sep="_")
males_subset <- as.matrix(na.omit(subset_genes))
plot_genes(males_subset)


subset_genes <- sertoli_data["Dmrt1",]
colnames(subset_genes) <- paste(colnames(subset_genes), sertoli_obj_j2_tree@ident, sep="_")
males_subset <- as.matrix(na.omit(subset_genes))
plot_genes(males_subset)

subset_genes <- sertoli_data["Nr5a1",]
colnames(subset_genes) <- paste(colnames(subset_genes), sertoli_obj_j2_tree@ident, sep="_")
males_subset <- as.matrix(na.omit(subset_genes))
plot_genes(males_subset)


subset_genes <- sertoli_data["Amh",]
colnames(subset_genes) <- paste(colnames(subset_genes), sertoli_obj_j2_tree@ident, sep="_")
males_subset <- as.matrix(na.omit(subset_genes))
plot_genes(males_subset)


subset_genes <- sertoli_data["Hsd17b1",]
colnames(subset_genes) <- paste(colnames(subset_genes), sertoli_obj_j2_tree@ident, sep="_")
males_subset <- as.matrix(na.omit(subset_genes))
plot_genes(males_subset)


subset_genes <- sertoli_data["Hsd17b3",]
colnames(subset_genes) <- paste(colnames(subset_genes), sertoli_obj_j2_tree@ident, sep="_")
males_subset <- as.matrix(na.omit(subset_genes))
plot_genes(males_subset)


subset_genes <- sertoli_data["Fshr",]
colnames(subset_genes) <- paste(colnames(subset_genes), sertoli_obj_j2_tree@ident, sep="_")
males_subset <- as.matrix(na.omit(subset_genes))
plot_genes(males_subset)


subset_genes <- sertoli_data["Dhh",]
colnames(subset_genes) <- paste(colnames(subset_genes), sertoli_obj_j2_tree@ident, sep="_")
males_subset <- as.matrix(na.omit(subset_genes))
plot_genes(males_subset)


subset_genes <- sertoli_data["Gadd45g",]
colnames(subset_genes) <- paste(colnames(subset_genes), sertoli_obj_j2_tree@ident, sep="_")
males_subset <- as.matrix(na.omit(subset_genes))
plot_genes(males_subset)


subset_genes <- sertoli_data["Mro",]
colnames(subset_genes) <- paste(colnames(subset_genes), sertoli_obj_j2_tree@ident, sep="_")
males_subset <- as.matrix(na.omit(subset_genes))
plot_genes(males_subset)





subset_genes <- sertoli_data["Tbx18",]
colnames(subset_genes) <- paste(colnames(subset_genes), sertoli_obj_j2_tree@ident, sep="_")
males_subset <- as.matrix(na.omit(subset_genes))
plot_genes(males_subset)



subset_genes <- sertoli_data["Lin28a",]
colnames(subset_genes) <- paste(colnames(subset_genes), sertoli_obj_j2_tree@ident, sep="_")
males_subset <- as.matrix(na.omit(subset_genes))
plot_genes(males_subset)



subset_genes <- sertoli_data["Nr2f1",]
colnames(subset_genes) <- paste(colnames(subset_genes), sertoli_obj_j2_tree@ident, sep="_")
males_subset <- as.matrix(na.omit(subset_genes))
plot_genes(males_subset)








subset_genes <- sertoli_data["Sry",]
colnames(subset_genes) <- paste(colnames(subset_genes), sertoli_obj_j2_tree@ident, sep="_")
males_subset <- as.matrix(na.omit(subset_genes))
plot_genes(males_subset)


subset_genes <- sertoli_data["Runx1",]
colnames(subset_genes) <- paste(colnames(subset_genes), sertoli_obj_j2_tree@ident, sep="_")
males_subset <- as.matrix(na.omit(subset_genes))
plot_genes(males_subset)


subset_genes <- sertoli_data["Nr0b1",]
colnames(subset_genes) <- paste(colnames(subset_genes), sertoli_obj_j2_tree@ident, sep="_")
males_subset <- as.matrix(na.omit(subset_genes))
plot_genes(males_subset)



subset_genes <- sertoli_data["Nr2f2",]
colnames(subset_genes) <- paste(colnames(subset_genes), sertoli_obj_j2_tree@ident, sep="_")
males_subset <- as.matrix(na.omit(subset_genes))
plot_genes(males_subset)



subset_genes <- sertoli_data["Fgf9",]
colnames(subset_genes) <- paste(colnames(subset_genes), sertoli_obj_j2_tree@ident, sep="_")
males_subset <- as.matrix(na.omit(subset_genes))
plot_genes(males_subset)


subset_genes <- sertoli_data["Mall",]
colnames(subset_genes) <- paste(colnames(subset_genes), sertoli_obj_j2_tree@ident, sep="_")
males_subset <- as.matrix(na.omit(subset_genes))
plot_genes(males_subset)


subset_genes <- sertoli_data["Lgr5",]
colnames(subset_genes) <- paste(colnames(subset_genes), sertoli_obj_j2_tree@ident, sep="_")
males_subset <- as.matrix(na.omit(subset_genes))
plot_genes(males_subset)




subset_genes <- sertoli_data["Igfbp7",]
colnames(subset_genes) <- paste(colnames(subset_genes), sertoli_obj_j2_tree@ident, sep="_")
males_subset <- as.matrix(na.omit(subset_genes))
plot_genes(males_subset)


subset_genes <- sertoli_data["Cst9",]
colnames(subset_genes) <- paste(colnames(subset_genes), sertoli_obj_j2_tree@ident, sep="_")
males_subset <- as.matrix(na.omit(subset_genes))
plot_genes(males_subset)


subset_genes <- sertoli_data["Ptgds",]
colnames(subset_genes) <- paste(colnames(subset_genes), sertoli_obj_j2_tree@ident, sep="_")
males_subset <- as.matrix(na.omit(subset_genes))
plot_genes(males_subset)

subset_genes <- sertoli_data["Dcn",]
colnames(subset_genes) <- paste(colnames(subset_genes), sertoli_obj_j2_tree@ident, sep="_")
males_subset <- as.matrix(na.omit(subset_genes))
plot_genes(males_subset)


subset_genes <- sertoli_data["Sry",]
colnames(subset_genes) <- paste(colnames(subset_genes), sertoli_obj_j2_tree@ident, sep="_")
males_subset <- as.matrix(na.omit(subset_genes))
plot_genes(males_subset)


subset_genes <- sertoli_data["Wt1",]
colnames(subset_genes) <- paste(colnames(subset_genes), sertoli_obj_j2_tree@ident, sep="_")
males_subset <- as.matrix(na.omit(subset_genes))
plot_genes(males_subset)


subset_genes <- sertoli_data["Fgf9",]
colnames(subset_genes) <- paste(colnames(subset_genes), sertoli_obj_j2_tree@ident, sep="_")
males_subset <- as.matrix(na.omit(subset_genes))
plot_genes(males_subset)


subset_genes <- sertoli_data["Sfrp1",]
colnames(subset_genes) <- paste(colnames(subset_genes), sertoli_obj_j2_tree@ident, sep="_")
males_subset <- as.matrix(na.omit(subset_genes))
plot_genes(males_subset)


subset_genes <- sertoli_data["Sfrp2",]
colnames(subset_genes) <- paste(colnames(subset_genes), sertoli_obj_j2_tree@ident, sep="_")
males_subset <- as.matrix(na.omit(subset_genes))
plot_genes(males_subset)


subset_genes <- sertoli_data["Nr0b1",]
colnames(subset_genes) <- paste(colnames(subset_genes), sertoli_obj_j2_tree@ident, sep="_")
males_subset <- as.matrix(na.omit(subset_genes))
plot_genes(males_subset)


subset_genes <- sertoli_data["Arx",]
colnames(subset_genes) <- paste(colnames(subset_genes), sertoli_obj_j2_tree@ident, sep="_")
males_subset <- as.matrix(na.omit(subset_genes))
plot_genes(males_subset)
