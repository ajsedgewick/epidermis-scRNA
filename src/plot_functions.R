library(NMF)
library(ggplot2)
library(reshape)

#gg_color_hue <- function(n) {
#  hues = seq(15, 375, length = n + 1)
#  hcl(h = hues, l = 65, c = 100)[1:n]
#}
#gg.colors <- gg_color_hue(8)
#cluster.colors <- c("green4", gg.colors[3], "greenyellow", "gold",
#                    gg.colors[5], "purple2", "peru", gg.colors[1])

# k=11 cluster colors
#cluster.colors <- c("green4", gg.colors[3], "greenyellow", "gold",
#                   gg.colors[5], "purple2", gg.colors[2], gg.colors[4], gg.colors[6], "peru", gg.colors[1])
cluster.colors <- c('#008B00FF', '#7CAE00FF', '#FFA500FF', '#FFD700FF', '#00BFC4FF',
                    '#308ED2FF', '#7771B5FF', '#912CEEFF',  '#8B5A2BFF', '#CD853FFF','#FF5347FF')

make.FC.mat <- function(de.list){
  cluster.FC <- lapply(de.list, function(x) data.frame(Gene=rownames(x), logFC=x$logFC, stringsAsFactors=F))
  cluster.FC <- Reduce(function(df1, df2) merge(df1, df2, by='Gene', all=T), cluster.FC)
  rownames(cluster.FC) <- cluster.FC$Gene
  cluster.FC <- cluster.FC[,-1]
  names(cluster.FC) <- names(de.list)
  return(cluster.FC)
}

# Take merged FC matrix, filter by colinds and threshold FC for plotting (optional abs(FC))
format.FC <- function(fc.mat, colinds, FC.thresh=1.5, abs.flag=F, cluster.order=1:length(colinds)){
  out.FC <- fc.mat[,colinds]
  out.FC <- out.FC[,cluster.order]
  colnames(out.FC) <- paste0("cl", 1:length(colinds))
  if(abs.flag){
    out.FC <- out.FC[apply(abs(out.FC), 1, max, na.rm=T) > FC.thresh,]
  } else {
    out.FC <- out.FC[apply(out.FC, 1, max, na.rm=T) > FC.thresh,]
  }
  out.FC[is.na(out.FC)] <- 0
  return(out.FC)
}

cluster.plot.fc <- function(cluster.fc, outfile){
fc.melt <- melt(cluster.fc, id.vars=c("Gene","Index"), measure.vars=paste("cl", 1:8, sep="."),
                variable_name = "Cluster")

p <- ggplot(fc.melt, aes(Index, FC, col=Cluster)) + geom_point(size=.5) +
    scale_color_manual(values = rev(cluster.colors)) + facet_grid(Cluster~., scales = "free") + theme_bw()
  ggsave(outfile, plot=p)
}


tissue.plot.fc <- function(data.fc, outfile, color=colorRampPalette(c("#2C6493", "#FFFFFF", "#E41A1C"))(50)){
  aheatmap(data.fc, color = color, breaks=0, hclustfun="average", distfun = "pearson",
           scale="none", Colv = NA, annCol=data.frame(cluster=1:ncol(data.fc)), 
           annColors=list(cluster=cluster.colors),
           filename= outfile)
}

pso.plot.fc <- function(data.fc, outfile, color=colorRampPalette(c("#2C6493", "#FFFFFF", "#E41A1C"))(50)){
  aheatmap(data.fc, color = color, breaks=0, hclustfun="average", distfun = "pearson",
           scale="none", Colv = NA, annCol=data.frame(cluster=1:ncol(data.fc)), annColors=list(cluster=cluster.colors),
           width = 8, height = 10
           ,filename = outfile)
}


