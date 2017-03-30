###########################################
#                                         #
#          Load needed libraries          #
#                                         #
###########################################

library("data.table") 
library("ggplot2")
library("gplots")
library("viridis")

###########################################
#                                         #
#                 Load data               #
#                                         #
###########################################

load(file="../Data/XY_single-cell_RPKM_table.Robj")
load(file="../Data/male_tsne_100000_iter_1-4_600_1-6.Robj")
load(file="../Data/clustered_male_tsne_100000_iter_1-4_600_1-6.Robj")

# Use cell clustering from SEURAT to exctract Cluster2 cells (cf script 1.SEURAT_clustering all_cells.R)
cluster2_cells <- names(all_males_obj_j2_tree@ident[all_males_obj_j2_tree@ident==2])
cluster1_cells <- names(all_males_obj_j2_tree@ident[all_males_obj_j2_tree@ident==1])

cluster2_expr <- as.matrix(all_males_obj_j2_tree@raw.data[,c(cluster1_cells,cluster2_cells)])
cluster2_expr <- cluster2_expr[rowSums(cluster2_expr)>0,]

stage <- as.numeric(substr(c(cluster1_cells,cluster2_cells),15,18))

###########################################
#                                         #
#         Perform linear regression       #
#                                         #
###########################################

lmPvals <- data.frame(
	pval=numeric(),
	coef=numeric()
	)

for (i in 1:dim(cluster2_expr)[1]) {
    apelm.sum <- summary(lm(cluster2_expr[i, ] ~ stage))
    pval <- pf(apelm.sum$fstatistic[1], apelm.sum$fstatistic[2], apelm.sum$fstatistic[3], 
        lower.tail = FALSE)
    mypval <- data.frame(
    	pval=pval,  
    	coef=coef(apelm.sum)[2]
    	)
    rownames(mypval) <- rownames(cluster2_expr)[i]
    lmPvals <- rbind(
    	lmPvals,
    	mypval
    )
}

hist(lmPvals$pval, breaks=50)

# p.adjust
FDR.adjusted.pvals <- p.adjust(lmPvals$pval, method = "BH")
FDR.adjusted.pvals <- cbind(
	lmPvals,
	adj.pval=FDR.adjusted.pvals
	)

hist(FDR.adjusted.pvals$adj.pval, breaks=50)
sig.FDR.adjusted.pvals <- as.data.frame(FDR.adjusted.pvals[FDR.adjusted.pvals$adj.pval<0.01,])

write.csv(sig.FDR.adjusted.pvals, "lm_prog.csv")
sig.FDR.adjusted.pvals <- read.csv("lm_prog.csv", header=T, row.names = 1)

# sig.FDR.adjusted.pvals <- as.data.frame(sig.FDR.adjusted.pvals[sig.FDR.adjusted.pvals$adj.pval<0.01,])

up_genes <- rownames(sig.FDR.adjusted.pvals[sig.FDR.adjusted.pvals$coef>0,])
down_genes <- rownames(sig.FDR.adjusted.pvals[sig.FDR.adjusted.pvals$coef<0,])
all_genes <- rownames(sig.FDR.adjusted.pvals)

###########################################
#                                         #
#        Plot up-down genes heatmap       #
#                                         #
###########################################

up <- cluster2_expr[up_genes,]
down <- cluster2_expr[down_genes,]

stagePalette <- c("#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c")


col <- factor(stage)
levels(col) <- stagePalette
col <- as.vector(col)

heatmap.2(up, trace="none", col=c(viridis(50),rep("#FDE725FF", 50)), Colv="FALSE", ColSideColors=col, dendrogram='none', labRow = FALSE, labCol = FALSE)

heatmap.2(down, trace="none", col=c(viridis(50),rep("#FDE725FF", 50)), Colv="FALSE", ColSideColors=col, dendrogram='none', labRow = FALSE, labCol = FALSE)


###########################################
#                                         #
#      Plot individual genes (violin)     #
#                                         #
###########################################


plot_gene_violin <- function(data, title) {
	mean_exp <- data
	stages <- substr(names(data),14,18)
	plot_mean <- as.data.frame(
		cbind(
			cells=seq(names(data)),
			mean=mean_exp,
			stages=stages
			)
		)
	g <- ggplot(plot_mean, aes(stages, as.numeric(as.character(mean)))) +
	geom_violin(aes(fill=factor(stages)), adjust = .6, scale = "width") + 
	geom_jitter(height = 0, size=0.5) +
	stat_smooth(aes(x = stages, y = as.numeric(as.character(mean)), group=1), method="loess", colour="black", se=FALSE) + 
	theme_bw() +
	xlab("Cells") +
	ylab("Expression (log(RPKM+1))") +
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
		# axis.title=element_text(size=18),
		axis.title=element_blank(),
		legend.text = element_text(size =18),
		legend.title = element_text(size =18 ,face="bold"),
		plot.title = element_text(size=24, face="bold.italic", hjust = 0.5),
		axis.text.x = element_text(angle = 45, hjust = 1),
		legend.position="none"
	)
	print(g)
}




males_stages <- as.vector(substr(all_males_obj_j2_tsne@ident,1,5))
stagePalette <- c("#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c")


graph_violin <- function(genes){
	for (gene in genes){
		print(gene)
		subset_genes <- cluster2_expr[gene,]
		plot_gene_violin(subset_genes, gene)
	}
}

graph_violin(c(
	"Hmga2",
	"Foxp1",
	"Gata6",
	"Wt1",
	"Pdgfra",
	"Arx",
	"Tcf21",
	"Wnt5a",
	"Sfrp1",
	"Sfrp2"
	)
)
