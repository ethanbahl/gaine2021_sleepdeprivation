########## Supplemental Figures: Normalization.
library(EDASeq) # EDASeq_2.20.0
library(scatterplot3d) # scatterplot3d_0.3-41
library(RColorBrewer) # RColorBrewer_1.1-2
library(svglite) # svglite_1.2.3

### grab the 3 sets.
set = data.gene$set
dataWithin = data.gene$dataWithin
dataNorm = data.gene$dataNorm

colors = RColorBrewer::brewer.pal(12, "Paired")[c(3,4,9,10)]
legend.cols = colors
group_color = pData(set)$group
plot.cols = colors[group_color]

#svglite("figures/suppfig_edaseq.svg", height=16, width=18, system_fonts = list(sans = "Arial"))

cex.lab = 1.75
cex.main = 2
cex.axis=1.3

m = matrix(
    c(
        1,2,3,
        4,5,6,
        7,8,9,
        10,11,12
    ),
    byrow=TRUE,
    nrow=4
)
layout(mat = m, heights = c(0.04, rep(0.24, 3)))
par(mar=c(0,0,0,0))#, family="Arial")

plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
text(x = 0.525, y = 0.5, paste("raw"), 
     cex = 3, font=2, col = "black")
segments(x0=0.1, x1=1, y0=0.2, y1=0.2, lwd=2.5, col="black")

plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
text(x = 0.525, y = 0.5, paste("GC normalization"), 
     cex = 3, font=2, col = "black")
segments(x0=0.1, x1=1, y0=0.2, y1=0.2, lwd=2.5, col="black")

plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
text(x = 0.525, y = 0.5, paste("GC / depth normalization"), 
     cex = 3, font=2, col = "black")  
segments(x0=0.1, x1=1, y0=0.2, y1=0.2, lwd=2.5, col="black")

par(
  mar=c(5,5,4,2)
)

# GC content parameters.
ylab = "log(gene counts)"
xlab = "GC content"
ylim=c(2,7)
xlim=c(min(fData(set)$gc),max(fData(set)$gc))
gc.col = "#333333"
c1 = alpha(gc.col, 0.8)
c2 = alpha(gc.col, 0.2)
c3 = alpha(gc.col, 0.8)

########## GC bias plots.
biasPlot(set, "gc", log=T, col=plot.cols, lwd=1.4, xlim=xlim, ylim=ylim, main="GC content distribution", xlab=xlab, ylab=ylab, cex.main=cex.main, cex.lab=cex.lab, cex.axis=cex.axis)
#mtext(side=1, line=2, xlab, col="black", font=2, cex=1.2)
#mtext(side=2, line=3, "Y-axis label, bold, bigger", col="orange", font=2, cex=1.2)
legend("topleft", legend=levels(group_color), fill=legend.cols, cex=1.1, bty="n")
boxplot(fData(set)$gc, add=TRUE, horizontal=TRUE, col=c2, medcol=c3, whiskcol=c1, staplecol=c3, boxcol=c3, outcol=c3, pch=23, cex=0, xaxt="n", at=3)
text(x=median(fData(data.gene$set)$gc), y=2.5, labels=c("GC distribution of all genes tested"), col=c1, cex=cex.axis)

biasPlot(dataWithin, "gc", log=T, col=plot.cols, lwd=1.4, xlim=xlim, ylim=ylim, main="GC content distribution", xlab=xlab, ylab=ylab, cex.main=cex.main, cex.lab=cex.lab, cex.axis=cex.axis)
legend("topleft", legend=levels(group_color), fill=legend.cols, bty="n", cex=1.1)
boxplot(fData(set)$gc, add=TRUE, horizontal=TRUE, col=c2, medcol=c3, whiskcol=c1, staplecol=c3, boxcol=c3, outcol=c3, pch=23, cex=0, xaxt="n", at=3)
text(x=median(fData(set)$gc), y=2.5, labels=c("GC distribution of all genes tested"), col=c1, cex=cex.axis)

biasPlot(dataNorm, "gc", log=T, col=plot.cols, lwd=1.4, xlim=xlim, ylim=ylim, main="GC content distribution", xlab=xlab, ylab=ylab, cex.main=cex.main, cex.lab=cex.lab, cex.axis=cex.axis)
legend("topleft", legend=levels(group_color), fill=legend.cols, cex=1.1, bty="n")
boxplot(fData(set)$gc, add=TRUE, horizontal=TRUE, col=c2, medcol=c3, whiskcol=c1, staplecol=c3, boxcol=c3, outcol=c3, pch=23, cex=0, xaxt="n", at=3)
text(x=median(fData(set)$gc), y=2.5, labels=c("GC distribution of all genes tested"), col=c1, cex=cex.axis)

########## RLE
ylim=c(-1, 1)
plotRLE(set, col=plot.cols, ylim=ylim, main=paste0("RLE plot"), cex=0, xaxt="n", cex.main=cex.main, cex.lab=cex.lab, cex.axis=cex.axis, ylab="relative log expression")
legend("topleft", legend=levels(group_color), fill=legend.cols, cex=1.1, bty="n")

plotRLE(dataWithin, col=plot.cols, ylim=ylim, main=paste0("RLE plot"), cex=0, xaxt="n", cex.main=cex.main, cex.lab=cex.lab, cex.axis=cex.axis, ylab="relative log expression")
legend("topleft", legend=levels(group_color), fill=legend.cols, cex=1.1, bty="n")

plotRLE(dataNorm, col=plot.cols, ylim=ylim, main=paste0("RLE plot"), cex=0, xaxt="n", cex.main=cex.main, cex.lab=cex.lab, cex.axis=cex.axis, ylab="relative log expression")
legend("topleft", legend=levels(group_color), fill=legend.cols, cex=1.1, bty="n")


########## PCA
Y <- apply(log(counts(set)+1), 1, function(y) scale(y, center=TRUE, scale=FALSE))
s1 <- svd(Y)

Y <- apply(log(normCounts(dataWithin)+1), 1, function(y) scale(y, center=TRUE, scale=FALSE))
s2 <- svd(Y)

Y <- apply(log(normCounts(dataNorm)+1), 1, function(y) scale(y, center=TRUE, scale=FALSE))
s3 <- svd(Y)


xlim=ylim=zlim = c(-1,1)

s = s1
s$u = s$u[,-2]
scatterplot3d(x=s$u[,1], y=s$u[,2], z=s$u[,3], pch=16, color=plot.cols, main="3D PCA plot", cex.symbols=2, angle=55, type="h", xlab="PC1", ylab="", zlab="PC4", xlim=xlim, ylim=ylim, zlim=zlim, cex.main=cex.main, cex.lab=cex.lab-0.65, cex.axis=cex.axis-0.4)
text(x = 3.5, y = -0.85, labels="PC3", srt = 55, cex=cex.lab-0.05)
mtext(text="PC2 not shown, correlation with batch > 0.9", side=1, line=3, at = c(0), cex=0.8)

s = s2
s$u = s$u[,-2]
scatterplot3d(x=s$u[,1], y=s$u[,2], z=s$u[,3], pch=16, color=plot.cols, main="3D PCA plot", cex.symbols=2, angle=55, type="h", xlab="PC1", ylab="", zlab="PC4", xlim=xlim, ylim=ylim, zlim=zlim, cex.main=cex.main, cex.lab=cex.lab-0.65, cex.axis=cex.axis-0.4)
text(x = 3.5, y = -0.85, labels="PC3", srt = 55, cex=cex.lab-0.05)
mtext(text="PC2 not shown, correlation with batch > 0.9", side=1, line=3, at = c(0), cex=0.8)

s = s3
s$u = s$u[,-1]
scatterplot3d(x=s$u[,1], y=s$u[,2], z=s$u[,3], pch=16, color=plot.cols, main="3D PCA plot", cex.symbols=2, angle=55, type="h", xlab="PC2", ylab="", zlab="PC4", xlim=xlim, ylim=ylim, zlim=zlim, cex.main=cex.main, cex.lab=cex.lab-0.65, cex.axis=cex.axis-0.4)
text(x = 3.5, y = -0.85, labels="PC3", srt = 55, cex=cex.lab-0.05)
mtext(text="PC1 not shown, correlation with batch > 0.9", side=1, line=3, at = c(0), cex=0.8)

#dev.off()
