# TODO: Move this into multi-sample comparisons and finish https://bioconductor.org/packages/release/workflows/vignettes/RNAseq123/inst/doc/limmaWorkflow.html

library(edgeR)
batch <- sce$batch
timepoint <- factor(sce$timepoint)
group_1 <- factor(sce$group)
group <- factor(paste0(sce$group, ".", sce$timepoint))
x <- DGEList(as.matrix(counts(sce)), samples = colData(sce), group = group)
cpm <- cpm(x)
lcpm <- cpm(x, log = TRUE)

L <- mean(x$samples$lib.size) * 1e-6
M <- median(x$samples$lib.size) * 1e-6
c(L, M)

table(rowSums(x$counts==0)==ncol(x))

keep.exprs <- filterByExpr(x, group=x$group)
x <- x[keep.exprs,, keep.lib.sizes=FALSE]
dim(x)

x <- calcNormFactors(x, method = "TMM")
x$samples$norm.factors


lcpm <- cpm(x, log=TRUE)

par(mfrow=c(2,2))

col.timepoint <- timepoint
levels(col.timepoint) <-  brewer.pal(nlevels(col.timepoint), "Set1")
col.timepoint <- as.character(col.timepoint)
plotMDS(lcpm, labels=group, col=col.timepoint)
title(main="timepoint")

col.group_1 <- group_1
levels(col.group_1) <- brewer.pal(nlevels(col.group_1), "Set2")
col.group_1 <- as.character(col.group_1)
plotMDS(lcpm, labels=group, col=col.group_1)
title(main="group_1")

col.group <- group
levels(col.group) <-  brewer.pal(nlevels(col.group), "Set3")
col.group <- as.character(col.group)
plotMDS(lcpm, labels=group, col=col.group)
title(main="group")

col.batch <- batch
levels(col.batch) <-  brewer.pal(nlevels(col.batch), "Dark2")
col.batch <- as.character(col.batch)
plotMDS(lcpm, labels=group, col=col.batch)
title(main="batch")
