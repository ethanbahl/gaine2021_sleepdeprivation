########## Differential Expression Analysis.

library(EDASeq) # EDASeq_2.20.0
library(edgeR) # edgeR_3.28.0

########## load data.
### read count data.
count.data = as.matrix(read.csv("data/counts_all.csv", row.names=1, header=TRUE))
### read sample information, set up factors/levels.
sample.info = read.csv("data/sample_info.csv", row.names=1, header=TRUE, stringsAsFactors=TRUE)
### read gene information.
gene.info = read.csv("data/gene_annotation.csv", row.names=1, header=TRUE, stringsAsFactors=FALSE)

########## create unfiltered SeqExpressionSet object.
# create seqExpressionSet.
set.all = newSeqExpressionSet(counts=count.data, phenoData=sample.info, featureData=gene.info)

set = set.all

########## gene filtering.
### filtering with 'filterByExpr'.
y <- DGEList(counts(set), samples=pData(set), genes=fData(set))
design_filter = model.matrix(~ batch + exposure, data=pData(set))
keep = filterByExpr(y, design=design_filter)
y <- y[keep, , keep.lib.sizes=FALSE]
set = set[keep,]; rm(keep)

########## EDASeq.
### EDASeq normalization for GC content, and then sequencing depth.
dataWithin <- withinLaneNormalization(set, "gc", which="full", offset=TRUE)
dataNorm <- betweenLaneNormalization(dataWithin, which="upper", offset=TRUE)

data.gene = list(set=set, dataWithin=dataWithin, dataNorm=dataNorm)

########## edgeR (QL).
### construct the DGEList and add offsets from EDASeq.
y = DGEList(counts=counts(data.gene$dataNorm), samples=pData(data.gene$dataNorm), genes=fData(data.gene$dataNorm))
y$offset = -offst(data.gene$dataNorm)

### create design matrix.
design = model.matrix(~ batch + exposure, data=pData(data.gene$dataNorm))

### estimate dispersion and fit model.
y = estimateDisp(y, design=design, robust=TRUE)
fit = glmQLFit(y, design, robust=TRUE)

### conduct statistical test. 
test = glmQLFTest(fit)

### get differential expression test results.
top = topTags(test, Inf)$table
top$effect.size = sign(top$logFC) * (abs(top$logFC) / (sqrt(1/(top$logCPM + y[rownames(top),]$trended.dispersion + abs(min(top$logCPM)))))   )

### add first pass to 'out' variable.
results.gene = list(set=data.gene$dataNorm, y=y, fit=fit, test=test, top=top)

########## enrichment of positive controls.
### Fisher test.
with(subset(top, vecsey.tested==TRUE), fisher.test(table(FDR<=0.1, vecsey.significant==TRUE)))