########## getting positive control information from publicly available data.
### PMID: 22930738, using GEO2R (default settings) to get statistics on all genes.

# read data.
pos = read.table("data/pmid_22930738_geo2r_results.txt", sep="\t", header=TRUE, stringsAsFactors=F, quote="")

# some IDs map to multiple genes. They're separated by "///" in the 'Gene.symbol' column. String-split these.
split.genes = strsplit(pos$Gene.symbol, split="///")

# pull the genes mapped to a signficant ID at least once.
pos.sig = unique(unlist(split.genes[pos$adj.P.Val<=0.05]))

# add a boolean column to our gene annotation dataframe indicating whether a gene was significant in the previous study.
gene.info$vecsey.significant = gene.info$gene_name %in% pos.sig
# and add an indication of whether or not the gene was present in the file from the previous study.
gene.info$vecsey.tested = gene.info$gene_name %in% unlist(split.genes)

# multiple IDs (rows) can map to the same gene, giving multiple lines of evidence for that gene.
# but a gene can map to both significant and non-significant IDs.
# here, for each gene in our data, we get the number of IDs (from the previous study results) which were significant, and how many were non-significant.
# the column 'vecsey.significant' column we added above will be TRUE if the gene had at least one significant ID.
# FYI, there's probably a more elegant coding solution for this. The code snippet can take a while to run.
system.time({
    pos.nper = do.call('rbind', lapply(gene.info$gene_name, function(x) tapply(split.genes, pos$adj.P.Val <= 0.05, function(x2) sum(unlist(x2) == x))))
})

# add the resulting ID counts to our annotation.
gene.info = cbind(gene.info, vecsey.nsig = pos.nper[,2], vecsey.ninsig =pos.nper[,1])