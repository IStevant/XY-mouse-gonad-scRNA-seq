###########################################
#                                         #
#          Load needed libraries          #
#                                         #
###########################################

library("data.table") 
library("ggplot2") 

###########################################
#                                         #
#                 Load data               #
#                                         #
###########################################

# Use cell clustering from SEURAT to exctract Cluster2 cells (cf script 1.SEURAT_clustering all_cells.R)
cluster2_cells <- names(all_males_obj_j2_tree@ident[all_males_obj_j2_tree@ident==2])

cluster2_expr <- as.matrix(all_males_obj_j2_tree@raw.data[,cluster2_cells])
cluster2_expr <- cluster2_expr[rowSums(cluster2_expr)>0,]

stage <- as.numeric(substr(cluster2_cells,15,18))

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
sig.FDR.adjusted.pvals <- FDR.adjusted.pvals[FDR.adjusted.pvals$adj.pval<0.05,]


###########################################
#                                         #
#       Plot mean exp genes per stage     #
#                                         #
###########################################

A <- setDT(df[df$factor == "a", ])[, mean(value)]


mean_exp <- as.numeric(as.character(rowMeans(cluster2_expr)))
stages <- substr(colnames(cluster2_expr),14,18)

plot_mean <- as.data.frame(
	cbind(
		mean=mean_exp,
		stages=as.factor(stages)
		)
	)

ggplot(plot_mean, aes(stages, as.numeric(as.character(mean)))) +
		# ggplot(plot_mean, aes(as.numeric(cells), as.numeric(as.character(mean)))) +
		# geom_point(shape=21, size=3.5, stroke=0.5, aes(fill=factor(stages))) + 
		# geom_violin(aes(fill=factor(stages)), adjust = .6, scale = "width") + 
		geom_jitter(aes(fill=factor(stages)),height = 0) +
		stat_smooth(aes(x = stages, y = as.numeric(as.character(mean)), group=1), method="loess", colour="black") + 
		theme_bw() +
		xlab("Cells") +
		ylab("mean expression (log(RPKM+1))") +
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
			plot.title = element_text(size=18, face="bold")
		)
)



plot_genes <- function(data, title) {
	mean_exp <- as.numeric(as.character(rowMeans(data)))
	stages <- substr(colnames(data),14,18)

	plot_mean <- as.data.frame(
		cbind(
			mean=mean_exp,
			stages=as.factor(stages)
			)
		)




	g <- ggplot(plot_mean, aes(stages, as.numeric(as.character(mean)))) +
		# ggplot(plot_mean, aes(as.numeric(cells), as.numeric(as.character(mean)))) +
		# geom_point(shape=21, size=3.5, stroke=0.5, aes(fill=factor(stages))) + 
		geom_violin(aes(fill=factor(stages)), adjust = .6, scale = "width") + 
		geom_jitter(height = 0, size=1) +
		stat_smooth(aes(x = stages, y = as.numeric(as.character(mean)), group=1), method="loess", colour="black") + 
		theme_bw() +
		xlab("Cells") +
		ylab("mean expression (log(RPKM+1))") +
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
			plot.title = element_text(size=18, face="bold")
		)
	)
	print(g)
}


males_subset <- cluster2_expr[rownames(up_genes),]

plot_genes(males_subset, "Significantly up-regulated genes (228 genes)")


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
	axis.title=element_text(size=18),
	legend.text = element_text(size =18),
	legend.title = element_text(size =18 ,face="bold"),
	plot.title = element_text(size=18, face="bold.italic")
)
	print(g)
}




males_stages <- as.vector(substr(all_males_obj_j2_tsne@ident,1,5))
stagePalette <- c("#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c")
stagePalette <- c("#abd9e9", "#ffffbf", "#fdae61", "#d7191c")


graph_violin <- function(genes){
	for (gene in genes){
		print(gene)
		subset_genes <- cluster2_expr[gene,]
		plot_gene_violin(subset_genes, gene)
	}
}


graph_violin(c(
	"Trim71",
	"Lin28a",
	"Epha7",
	"Snurf",
	"Igf2",
	"Wt1"
	)
)
