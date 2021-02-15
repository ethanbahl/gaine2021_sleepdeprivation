########## Volcano plot.
library(RColorBrewer) # RColorBrewer_1.1-2
library(svglite) # svglite_1.2.3

### plot setup.
tt = results.gene$top
ylim = c(0, max(-log10(tt[,"FDR"])))
xlim=c(-1.3,1.3)
cex.lab = 1.75
cex.main = 2
cex.axis=1.3
sig = tt[,"FDR"] <= 0.1

### save plot to SVG.
# svglite("figures/volcano_rdbu.svg", height=8, width=8, system_fonts = list(sans = "Arial"))
par(mar=c(5,5,4,2))

### plot insignificant genes.
plot(tt[-which(sig),"logFC"], -log10(tt[-which(sig), "FDR"]),
    pch=16, cex=0.7,
    xlab="logFC (relative to control)",
    ylab=expression('-log'[10]*'(FDR)'),
    main="volcano plot",
    xlim=xlim,
    ylim=ylim,
    col=scales::alpha("#333333", 0.3),
    cex.main=cex.main,
    cex.lab=cex.lab,
    cex.axis=cex.axis
)

### add colored points for significant genes.
pal = brewer.pal(9, "RdBu")[c(1,9)]


points(tt[sig, "logFC"],-log10(tt[sig, "FDR"]),
    pch=16, cex=1,
    col=scales::alpha(ifelse(tt[sig,"logFC"] > 0, pal[1], pal[2]), 0.7)
)

### add gridlines.
abline(v=0, lty=2)
abline(h= -log10(0.1), lty=2)

### add legend.
legend("topright", legend=c("upregulated in SD", "downregulated in SD", "not significant"), fill=c(pal, "#333333"), bty="n", cex=1.1)

# dev.off()