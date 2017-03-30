set.ifnull=function(x,y) {
if(is.null(x)) x=y
return(x)
}

setGeneric("doHeatMap_2", function(object,cells.use=NULL,genes.use=NULL,disp.min=-2.5,disp.max=2.5,draw.line=TRUE,do.return=FALSE,order.by.ident=TRUE,col.use=viridis(100),slim.col.label=FALSE,group.by=NULL,remove.key=FALSE,...) standardGeneric("doHeatMap_2"))
#' @export
setMethod("doHeatMap_2","seurat",
function(object,cells.use=NULL,genes.use=NULL,disp.min=-2.5,disp.max=2.5,draw.line=TRUE,do.return=FALSE,order.by.ident=TRUE,col.use=viridis(100),slim.col.label=FALSE,group.by=NULL,remove.key=FALSE,...) {
cells.use=set.ifnull(cells.use,object@cell.names)
genes.use=ainb(genes.use,rownames(object@scale.data))
cells.use=ainb(cells.use,object@cell.names)
cells.ident=object@ident[cells.use]
if (!is.null(group.by)) cells.ident=factor(fetch.data(object,group.by)[,1])
cells.ident=factor(cells.ident,labels = ainb(levels(cells.ident),cells.ident))
if (order.by.ident) {
cells.use=cells.use[order(cells.ident)]
}
data.use=as.matrix(object@data[genes.use,cells.use])
data.use=minmax(data.use,min=0,max=6)
vline.use=NULL;
colsep.use=NULL
hmFunction=heatmap.2
if (remove.key) hmFunction=heatmap2NoKey
if (draw.line) {
colsep.use=cumsum(table(cells.ident))
}
if(slim.col.label && order.by.ident) {
col.lab=rep("",length(cells.use))
col.lab[round(cumsum(table(cells.ident))-table(cells.ident)/2)+1]=levels(cells.ident)
# hmFunction(data.use,Rowv=NA,Colv=NA,trace = "none",col=col.use,colsep = colsep.use,labCol=col.lab,cexCol=0.2+1/log10(length(unique(cells.ident))),...)
conditions <- substr(colnames(data.use),14,18)
# print(conditions)
colours <- c("#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c")
# colours <- c("#abd9e9", "#ffffbf", "#fdae61", "#d7191c")
cond <- factor(conditions)
cols <- factor(conditions)
levels(cols) <- colours
cols <- as.vector(cols)
# print(cols)
# hmFunction(data.use,trace = "none",col=col.use,colsep = colsep.use,labCol=col.lab,cexCol=0.2+1/log10(length(unique(cells.ident))), ColSideColors=cols,...)
heatmap.2(data.use, scale='none', trace="none", keysize = 1, Colv=FALSE, Rowv=FALSE, dendrogram="none", ColSideColors=cols, col=viridis(100),colsep = colsep.use,labCol=col.lab,cexCol=0.2+1/log10(length(unique(cells.ident))))
par(xpd=TRUE, mar=c (0, 0, 0, 0))
legend("left", legend=levels(cond),fill=colours)
}
else {
conditions <- substr(colnames(data.use),14,18)
print(conditions)
colours <- c("#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c")
cond <- factor(conditions)
col <- factor(conditions)
levels(col) <- colours
col <- as.vector(col)
# hmFunction(data.use,Rowv=NA,Colv=NA,trace = "none",col=col.use,colsep = colsep.use,...)
hmFunction(data.use, Colv=NA,trace = "none",col=col.use,colsep = colsep.use, ColSideColors=col, ...)
par(xpd=TRUE, mar=c (3, 4, 4, 3))
legend("topleft", inset=c(-0.1,0.1),legend=levels(cond),fill=colours)

heatmap.2(data.use, scale='none', trace="none", keysize = 1, Colv=FALSE, Rowv=FALSE, dendrogram="none", ColSideColors=col, col=viridis(100))
par(xpd=TRUE, mar=c (3, 4, 4, 3))
legend("topleft", inset=c(-0.1,0.1),legend=levels(cond),fill=colours)
}
if (do.return) {
return(data.use)
}
}
)