########## Heatmap of differentially expressed genes.
library(EDASeq) # EDASeq_2.20.0
library(edgeR) # edgeR_3.28.0
library(RUVSeq) # RUVSeq_1.20.0
library(RColorBrewer) # RColorBrewer_1.1-2
library(ComplexHeatmap) # ComplexHeatmap_2.2.0
library(circlize) # circlize_0.4.8
library(svglite) # svglite_1.2.3

tt = results.gene$top
set = results.gene$set

### get the normalized counts.
de.heat = normCounts(set)
### convert to logCPM.
de.heat = cpm(de.heat, prior.count=2, log=TRUE)
### pick the most significant genes (FDR<=1e-2) with the strongest effect sizes (>0.5) for visualization.
de.heat = de.heat[rownames(tt[tt$FDR<=0.01 & abs(tt$effect.size)>0.5 & !is.na(tt$effect.size),]),]
### assign gene symbols to rownames.
rownames(de.heat) = fData(set)[rownames(de.heat), "gene_name"]

### split by batch for scaling.
de.heat.batch1 = de.heat[, pData(set)$batch == "batch1"]
de.heat.batch2 = de.heat[, pData(set)$batch == "batch2"]
de.heat.batch1 = t(scale(t(de.heat.batch1)))
de.heat.batch2 = t(scale(t(de.heat.batch2)))
### fuse back together.
de.heat = cbind(de.heat.batch1, de.heat.batch2)[,colnames(de.heat)]

### for splitting the color-coded columns by group.
col_split = pData(set)[,"group", drop=FALSE]
rownames(col_split) = rownames(pData(set))

### for splitting the rows by direction of effect.
row_split = data.frame(direction=ifelse(tt[tt$FDR<=0.01 & abs(tt$effect.size)>0.5 & !is.na(tt$effect.size), "logFC"] < 0, "downregulated by\nsleep deprivation", "upregulated by\nsleep deprivation" ))

### color scale.
color.cells = colorRamp2(seq(c(-1)* max(abs(c(de.heat))), max(abs(c(de.heat))), length = 101), colorRampPalette(rev(brewer.pal(11, "RdBu")))(101))
color.samples.palette = brewer.pal(12, "Paired")[c(3,4,9,10)]
color.samples = color.samples.palette[pData(set)$group]

set.seed(777)
h = Heatmap(
    de.heat,
    col=color.cells, 
    
    # columns.
    cluster_columns=FALSE,
    column_split = col_split,
    column_gap = unit(c(2,4,2), "mm"),
    column_title_gp = gpar(col="black", fill = color.samples.palette[c(1,2,3,4)], font = 2, fontsize=12),
    column_names_rot = 45,
    column_names_gp = gpar(col="#333333", fontsize=10),
    column_names_centered = TRUE,

    # rows.
    row_title_gp = gpar(col="#333333", font =2, fontsize=14),
    row_names_gp = gpar(col="#333333", fontsize=9, font=2, hjust="center", vjust="center"),
    row_split = row_split,
    row_gap = unit(2, "mm"),
    row_names_side = "right",
    row_dend_side = "right",
    row_title_side = "left",
    cluster_rows=FALSE,
    show_row_names=TRUE,

    width= unit(5, "in"),
    height=unit(12, "in"),
    show_heatmap_legend=TRUE,
    border=TRUE,
    rect_gp = gpar(col= "#FFFFFF", lwd=1),

    heatmap_legend_param = list(
        title = "scaled logCPM\n(normalized)", at = seq(-2,2,1), 
        labels = as.character(seq(-2,2,1))
    )
)

#svglite("figures/heatmap_rdbu.svg", height=13, width=8, system_fonts = list(sans = "Arial"))
h
#dev.off()