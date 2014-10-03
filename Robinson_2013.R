
## ----barseq_library, results="hide", cache=FALSE-------------------------
library(methods)
library(plyr)

setClass("barcode.counts",
    representation(
        counts = "matrix",
        design = "data.frame",
        comparisons = "list",
        pvalue = "data.frame",
        qvalue = "data.frame",
        logFC = "data.frame",
        conditions = "character"
    )
)

setMethod("initialize", "barcode.counts",
    function(.Object, counts, design) {
        .Object@counts = counts

        if (NROW(design) != NCOL(counts)) {
            stop("Design must have as many rows as there are count columns")
        }

        # check for legal replicate names
        legal = c("tag", "biological", "technical", "condition", "time")
        for (n in colnames(design)) {
            if (!(n %in% legal)) {
                stop(paste(c("Illegal column name in design matrix",
                             "legal names are", legal), collapse=", "))
            }
        }
        
        .Object@design = design
        .Object@conditions = ""
        .Object@comparisons = list()
        .Object@pvalue = data.frame()

        .Object
    }
)

### combine.counts ###
## the combine.counts methods adds up and down tags for each sample, which
## is recommended before performing any analysis

setGeneric("combine.counts",
            function(.Object, ...) standardGeneric("combine.counts"))

setMethod("combine.counts", "barcode.counts",
    function(.Object, combine.by="tag") {
        if (length(.Object@comparisons) > 0) {
            stop(paste("Don't perform combining after differential abundance",
                       "comparisons have already been made"))
        }

        new.design = .Object@design[!(names(.Object@design) %in% combine.by)]
        interacter = do.call(interaction, c(new.design, list(drop=TRUE)))

        # sum each group
        .Object@counts = t(apply(.Object@counts, 1, function(row) {
            tapply(row, interacter, sum) 
        }))
    
        .Object@design = as.data.frame(apply(new.design, 2, function(x) tapply(x, interacter, function(x) x[1])))
        .Object
    })


add.data.frame = function(d, name, col) {
    # add a column to a data frame if the data frame isn't empty, otherwise
    # recreate it and add it as the first column
    if (NROW(d) == 0) {
        d = data.frame(a=col)
        colnames(d) = name
    }
    else {
        d[[name]] = col
    }
    d
}

setGeneric("differential.abundance",
            function(.Object, ...) standardGeneric("differential.abundance"))

setMethod("differential.abundance", "barcode.counts",
    function(.Object, method, condition1, condition2, ...) {
        if (length(.Object@conditions) == 2 &&
            (!all(c(condition1, condition2) == .Object@conditions))) {
            stop(paste("Already performed comparison between conditions",
                        .Object@conditions[1], "and", .Object@conditions[2]))
        }
        for (cond in list(condition1, condition2)) {
            if (!(cond %in% .Object@design$condition)) {
                stop(paste("Condition", cond, "not found in design matrix"))
            }
        }
        library(qvalue)
        treatment = .Object@design$condition
        if (method == "DESeq") {
            # perform a negative binomial test using the DESeq package
            require(DESeq)
            cds = newCountDataSet(.Object@counts, treatment)
	        cds = estimateSizeFactors(cds)
	        cds = estimateDispersions(cds, fitType="local", ...)
	        result = nbinomTest(cds, condition1, condition2)
	        .Object@pvalue = add.data.frame(.Object@pvalue, method,
	                                        result$pval)
	        .Object@qvalue = add.data.frame(.Object@qvalue, method,
	                                            qvalue(result$pval)$qvalue)
	        .Object@logFC = add.data.frame(.Object@logFC, method,
	                                            log2(result$foldChange))
        } else if (method == "edgeR") {
            require(edgeR)
            d = DGEList(counts=.Object@counts, group=treatment)
            d = estimateCommonDisp(d)
            d = estimateTagwiseDisp(d)
            result = exactTest(d, pair=c(condition1, condition2))
	        .Object@pvalue = add.data.frame(.Object@pvalue, method,
	                                        result$table$PValue)
	        .Object@qvalue = add.data.frame(.Object@qvalue, method,
	                                        qvalue(result$table$PValue)$qvalue)
	        .Object@logFC = add.data.frame(.Object@logFC, method,
	                                            result$table$logFC)
        } else if (method == "baySeq") {
            require(baySeq)
            groups = list(NDE=rep(1, length(treatment),
                          DE=treatment))
            CD = new("countData", data=.Object@counts,
                        replicates=treatment, groups=groups)
            CD@libsizes = getLibsizes(CD)
            CD@annotation = data.frame(name=rownames(.Object@counts))
            CD = getPriors.NB(CD, samplesize = 1000, estimation = "QL",
                                cl=NULL)
            result = getLikelihoods.NB(CD, pET='BIC', cl=NULL)
        } else {
            stop(paste("Only methods for differential abundance are DESeq,",
                       "edgeR and baySeq"))
        }
        .Object@comparisons[[method]] = result
        .Object@conditions = c(condition1, condition2)
        .Object
    }
)


setGeneric("filter.counts",
            function(.Object, ...) standardGeneric("filter.counts"))

setMethod("filter.counts", "barcode.counts",
    function(.Object, min.reads=NULL, max.reads=NULL) {
        if (length(.Object@comparisons) > 0) {
            stop(paste("Don't perform filtering after differential abundance",
                       "comparisons have already been made"))
        }
    
        # for now, we'll keep it simple
        if (!is.null(min.reads)) {
            .Object@counts = .Object@counts[rowSums(.Object@counts)
                                                    >= min.reads, ]
        }
    
        if (!is.null(max.reads)) {
            .Object@counts = .Object@counts[rowSums(.Object@counts)
                                                    <= max.reads, ]
        }
        .Object
    }
)

setMethod("[", c("barcode.counts"),
    function(x, i, j, ..., drop=TRUE) {
    x@counts = x@counts[, i]
    x@design = x@design[i, ]
    if (drop == TRUE) {
        x@design = droplevels(x@design)
    }
    x
})


## ----methods, cache=FALSE------------------------------------------------
require(methods)


## ----import_data---------------------------------------------------------
# import Bar-Seq data and the class methods
load("input/replicates.RData")

original = new("barcode.counts", counts=experiment.counts, design=experiment.design)
original.combined = combine.counts(original, "tag")

num.strains = sum(rowSums(original@counts) > 0)


## ----process_data, dependson="import_data"-------------------------------
counts = filter.counts(original.combined, min.reads=100, max.reads=1000000)


## ----counts_obj, dependson="process_data"--------------------------------
# "counts" means something different later: save this
counts.obj = counts


## ----levenshtein_matrices------------------------------------------------
# matrices using various levenshtein distances
counts.0 = read.table("input/counts/BarNone_0m_counts.txt", header=TRUE, row.names=1)
counts.1 = read.table("input/counts/BarNone_1m_counts.txt", header=TRUE, row.names=1)
counts.2 = read.table("input/counts/BarNone_2m_counts.txt", header=TRUE, row.names=1)
counts.3 = read.table("input/counts/BarNone_3m_counts.txt", header=TRUE, row.names=1)


## ----pool_counts, dependson="import_data"--------------------------------
pooled = combine.counts(counts, "technical")
pooled = pooled[pooled@design$condition != "t0"]


## ----subsample_proportions, dependson="pool_counts"----------------------
#subsample.proportions = c(seq(0.002, .1, .001), seq(.103, 1, .003))
subsample.proportions = c(seq(0.0025, .1, .0025), seq(.105, 1, .005))


## ----edgeR, dependson="pool_counts"--------------------------------------
pooled.counts = differential.abundance(pooled, "YPD", "YPGAL", method="edgeR")


## ----DESeq, dependson="edgeR"--------------------------------------------
pooled.counts = differential.abundance(pooled.counts, "YPD", "YPGAL", method="DESeq")


## ----significant.each, dependson="DESeq"---------------------------------
significant.edgeR = pooled.counts@qvalue$edgeR < .05
significant.DESeq = pooled.counts@qvalue$DESeq < .05


## ----GSEA_prep-----------------------------------------------------------
library(GSEABase)
library(data.table)
library(org.Sc.sgd.db)

frame = toTable(org.Sc.sgdGO)
frame = frame[frame$Ontology == "BP", ]
goframeData = data.frame(frame$go_id, frame$Evidence, frame$systematic_name)

goFrame=GOFrame(goframeData,organism="Saccharomyces cerevisiae")
goAllFrame=GOAllFrame(goFrame)
gsc <- GeneSetCollection(goAllFrame, setType = GOCollection())


## ----setup_wilcox, dependson=c("GSEA_prep", "DESeq")---------------------
library(GO.db)

GO = as.list(GOTERM)

# filter
all.genes = rownames(pooled.counts@counts)
combined.names = sapply(gsc, function(g) paste(sort(intersect(g@geneIds, all.genes)), collapse="|"))
set.sizes = sapply(gsc, function(g) length(g@geneIds))
gsc = gsc[!duplicated(combined.names) & set.sizes > 3]

gsc.terms = GO[names(gsc)]
ontologies = sapply(gsc.terms, Ontology)

gene.contained = sapply(gsc, function(g) all.genes %in% g@geneIds)

# filter- must be BP and have at least 4 genes in the sample
cond = ontologies == "BP" & colSums(gene.contained) >= 4
gsc.filtered = gsc[cond]
gene.contained = gene.contained[, cond]
gsc.terms = gsc.terms[cond]


## ----save_gene_contained, dependson="setup_wilcox"-----------------------
save(gene.contained, file="gene_contained.Rdata")


## ----wilcoxon_test, dependson="setup_wilcox"-----------------------------
wilcoxon.pvalues = apply(gene.contained, 2, function(g) wilcox.test(pooled.counts@logFC$edgeR[g], pooled.counts@logFC$edgeR[!g])$p.value)


## ----load_barseq---------------------------------------------------------
load("input/BarSeq_experiment_080513.Rdata")

# this defines depth as over the entire analysis; the paper
# uses it within condition (half as much)
summaries$depth = summaries$depth / 2


## ----barseq_processed_graph, dependson="load_barseq", warning=FALSE------
library(data.table)
library(ggplot2)

x = summaries[, general.design:=substr(design, 1, 6)]
# combine the tech.1 and tech.2 levels
summaries$general.design[summaries$general.design == "tech.2"] = "tech.1"
summaries$general.design = factor(summaries$general.design)

summaries$bio.reps = as.factor(c(2, 2, 3, 3, 4, 4, 4)[as.numeric(summaries$general.design)])
summaries$tech.reps = as.factor(c(2, 1, 2, 1, 2, 1, 1)[as.numeric(summaries$general.design)])
summaries$tech.reps = relevel(summaries$tech.reps, 2)

summaries$FDR[is.na(summaries$FDR)] = 0

# load the GSEA data from wilcox_significant.Rdata, which has the GSEA results
# of all subsamples
load("input/wilcox_significant.Rdata")

summaries$GSEA = wilcox.significant

predict.spline = function(x, y, df=20) {
    main.spl = smooth.spline(x[order(x)], y[order(x)], df=df, all.knots=TRUE)
    predict(main.spl)$y
}

summaries = summaries[order(summaries$general.design, summaries$depth), ]
x = summaries[, fit:=predict.spline(depth, significant), by=general.design]
x = summaries[, MSE.fit:=predict.spline(depth, MSE), by=general.design]
x = summaries[, fdr.fit:=predict.spline(depth, FDR), by=general.design]
x = summaries[, GSEA.fit:=predict.spline(depth, GSEA), by=general.design]

# problematic artifact of the spline fitting:
# some spline fits get confused right at the end of the line and drop precipitously
# this occurs when one subdesign (set of replicates) goes on further than another
# we trim precipitous drops right at the end
trim.bad.end = function(x, threshold) {
    # replace it with NA
    diffs = diff(x)
    change.indices = which(abs(x) > threshold)
    # only look at last 5%
    change.indices = change.indices[change.indices > length(x) * .9]
    if (length(change.indices) > 0) {
        x[change.indices[1]:length(x)] = NA
    }
    x
}

summaries = summaries[, fit:=trim.bad.end(fit, 100), by=c("bio.reps", "tech.reps")]
summaries = summaries[, MSE.fit:=trim.bad.end(MSE.fit, .1), by=c("bio.reps", "tech.reps")]
summaries = summaries[, fdr.fit:=trim.bad.end(fdr.fit, .02), by=c("bio.reps", "tech.reps")]
summaries = summaries[, GSEA.fit:=trim.bad.end(GSEA.fit, 50), by=c("bio.reps", "tech.reps")]

bio.rep.colors = scale_colour_manual(name="Biological Replicates", values=c("blue", "purple", "red"))
tech.rep.lty = scale_linetype_manual(name="Technical Replicates", values=1:2)

base.plot = ggplot(summaries, aes(x=depth, y=fit, col=bio.reps, lty=tech.reps)) + xlab("Read Depth per Condition") + ylab("Number of Significant Strains") + tech.rep.lty + bio.rep.colors + theme(axis.text.x = element_text(angle = 30, hjust = 1))
p1 = base.plot + geom_line() + xlim(0, 25e6)
p3 = base.plot + geom_line(aes(y=GSEA.fit)) + xlim(0, 25e6) + ylab("Number of Significant Gene Sets")


## ----by_logFC, dependson="barseq_processed_graph"------------------------
library(reshape)
percents.logFC = data.table(melt(summaries[, list(depth, bio.reps, tech.reps,
                                       percent, percent.logFC1.5, percent.logFC2)],
                                       id=c("depth", "bio.reps", "tech.reps")))
levels(percents.logFC$variable) = c("All", "1.5", "2")
x = percents.logFC[, percent.fit:=predict.spline(depth, value),
                by=c("bio.reps", "tech.reps", "variable")]


## ----power_data, dependson="by_logFC"------------------------------------
# judge the power of each method by where the spline curve is at the recommended 6 million
full.dt = summaries[general.design == "Full", ]
bio3.dt = summaries[general.design == "Bio.3.", ]
bio2.dt = summaries[general.design == "Bio.2.", ]

full.sig = approx(full.dt$depth, full.dt$fit, 6e6)$y
bio3.sig = approx(bio3.dt$depth, bio3.dt$fit, 6e6)$y
bio2.sig = approx(bio2.dt$depth, bio2.dt$fit, 6e6)$y

percents.logFC.report = percents.logFC[bio.reps == 4 & tech.reps == 1, ]

power.FC.depth = function(FC, depth) {
	library(data.table)  # need to reimport because this is used in Sexpr
	ss = percents.logFC.report[variable == FC, ]
	approx(ss$depth, ss$percent.fit, depth)$y
}


## ----misc_data-----------------------------------------------------------
total.reads = 185.2e6

total.mapped = sum(original@counts)
percent.mapped = total.mapped / total.reads
HO.index = which.max(rowSums(original@counts))
percent.HO = sum(original@counts[HO.index, ]) / sum(original@counts)

minimum.reads = 100

num.strains = sum(rowSums(original@counts) > 0)
num.outliers = sum(rowSums(original@counts) < minimum.reads & rowSums(original@counts) > 0)
num.strains.filtered = NROW(counts@counts)

num.non.HO = sum(original@counts[-HO.index, ])


## ----uptag_downtag_num_fold_change---------------------------------------
# use the filtered for this one- don't count under 100-
# but don't combine it by tag (of course)
original.filtered = filter.counts(original, min.reads=100, max.reads=1000000)

num_up = rowSums(original.filtered@counts[, original.filtered@design$tag == 1]) + 1
num_dn = rowSums(original.filtered@counts[, original.filtered@design$tag == 2]) + 1
num_up = num_up / sum(num_up)
num_dn = num_dn / sum(num_dn)

num.fold.change = function(b) sum(abs(log((num_up) / (num_dn), b)) > 1)

num.not.t0 = sum(counts@counts[, counts@design$condition != "t0"])


## ----percentvar, dependson="count_obj"-----------------------------------
library(xtable)
library(edgeR)

# note that eigenR2 is not currently in CRAN and has to be installed
library(eigenR2)

row.sums = rowSums(counts.obj@counts)

not.t0 = counts.obj@design$condition != "t0"

Treatment = counts.obj@design$condition[not.t0, drop=TRUE]
Biological = counts.obj@design$biological[not.t0, drop=TRUE]

library(DESeq)

normalize = function(v) v / sum(v)

counts.m = counts.obj@counts[, not.t0]

TMM.factors = colSums(counts.m) * calcNormFactors(counts.m, method="TMM")
DESeq.factors = estimateSizeFactorsForMatrix(counts.m)

normalized.matrix.sizeFactor = log(t(t(counts.m + 1) / DESeq.factors))
normalized.matrix.TMM = log(t(t(counts.m + 1) / TMM.factors))
normalized.matrix.voom = voom(counts.m)$E

percents.explained = function(input.m, min.count=0) {
	m = input.m[row.sums >= min.count, ]
	
	eigenR2.treatment = eigenR2(dat=m, model=model.matrix(~ 1 + Treatment), adjust=TRUE)
	eigenR2.bio = eigenR2(dat=m, model=model.matrix(~ 1 + Biological), adjust=TRUE)
	
	100 * diff(c(0, eigenR2.treatment$eigenR2, eigenR2.bio$eigenR2, 1))
}

quantiles = seq(0, .75, .05)
bins = quantile(rowSums(counts.obj@counts), quantiles)

matrices = list(normalized.matrix.TMM, normalized.matrix.voom, normalized.matrix.sizeFactor)
var.data = do.call(rbind, lapply(matrices, function(m) {
	data.frame(t(sapply(1:length(bins), function(i) {
    		percents.explained(m, bins[i])
	})))
}))

library(reshape)

colnames(var.data) = c("Treatment", "Biological", "Technical")
var.data$Method = rep(c("TMM (edgeR)", "Voom (limma)", "estimateSizeFactors (DESeq)"), each=length(bins))

var.data$quantile = quantiles
var.melt = melt(var.data, id=c("quantile", "Method"))

reported.percents = percents.explained(normalized.matrix.TMM, quantile(row.sums, .1))


## ----remove_scipen-------------------------------------------------------
options(scipen = 10000)


## ----figure_1b, dependson="import_data", error=FALSE, include=FALSE, dev=c('tiff', 'jpeg'), dev.args=list(bg="white"), fig.width=5,fig.height=5, dpi=500----
library(gplots)
library(colorRamps)

filtered = filter.counts(original, min.reads=100, max.reads=1000000)

cormatrix = cor(counts.obj@counts, method="spearman")
sample.names = paste(counts.obj@design$condition, c("A", "B", "C", "D", "A", "B", "C", "D", "A", "B")[counts.obj@design$biological], 1:2[rep(1:2, 10)], sep="_")
rownames(cormatrix) = sample.names
colnames(cormatrix) = sample.names

par(mar=c(0, 0, 0, 0))
heatmap.2(cormatrix, dendrogram="none", Rowv=NA, Colv=NA, trace="n", col = blue2yellow(1000), cexRow=.65, cexCol=.65, key=TRUE, density.info="none")


## ----volcano_plot, dependson="DESeq"-------------------------------------
library(scales)
GAL.table = read.table("input/GAL_genes.xls", header=TRUE, sep="\t")
GAL.genes = GAL.table$Systematic[GAL.table$Repressor == "no"]
GAL.repressors = GAL.table$Systematic[GAL.table$Repressor == "yes"]

edgeR.results = data.frame(logFC=pooled.counts@logFC$edgeR, pvalue=pooled.counts@pvalue$edgeR,
					   DESeq.pvalue=pooled.counts@pvalue$DESeq)
edgeR.results$Gal = ifelse(rownames(pooled.counts@counts) %in% GAL.genes, "Pathway", ifelse(rownames(pooled.counts@counts) %in% GAL.repressors, "Repressor", "Significant"))
edgeR.results$Reads = rowSums(pooled.counts@counts)
edgeR.results$qvalue = pooled.counts@qvalue$edgeR

# code from http://stackoverflow.com/questions/11053899
reverseloglog_trans <- function(base = exp(1)) {
    trans <- function(x) log(-log(x, base), base)
    inv <- function(x) base^(-(base^x))
    trans_new(paste0("reverseloglog-", format(base)), trans, inv, 
              function(x) c(.1, 1e-05, 1e-25, 1e-100), 
              domain = c(1e-300, Inf))
}

library(ggplot2)
library(gridExtra)

edgeR.results$Gal = factor(edgeR.results$Gal,
                            levels=c("Pathway", "Repressor", "Significant"))
edgeR.results$Deletion = rownames(pooled.counts@counts)

edgeR.results$color = edgeR.results$Gal
levels(edgeR.results$color) = c(levels(edgeR.results$color), "Not Significant")
edgeR.results$color[edgeR.results$qvalue > .05] = "Not Significant"
edgeR.results$color = factor(edgeR.results$color,
                             levels=c("Not Significant", "Significant", "Pathway", "Repressor"))

color.scale = scale_colour_manual(name="color", values=c("grey", "black", "red", "blue"))

edgeR.results = edgeR.results[order(edgeR.results$color), ]
volcano.plot = ggplot(edgeR.results, aes(x=qvalue, y=logFC, col=color)) + scale_x_continuous(trans=reverseloglog_trans(10)) + geom_point(size=2) + ylab("Log2(YPGal/YPD) FC") + xlab("EdgeR p-value") + theme(legend.direction = "horizontal", legend.position = "top") + color.scale + theme(axis.text.x = element_text(angle = 30, hjust = 1))
pval.hist = ggplot2::qplot(edgeR.results$pvalue, binwidth=.03, xlim=0:1, xlab="EdgeR p-value", ylab="Frequency")
logFC.hist = ggplot2::qplot(edgeR.results$logFC, binwidth=.25) + xlab("Log2(YPGal/YPD) Fold Change") + ylab("Frequency")


## ----figure_2, dependson="volcano_plot", include=FALSE, dev=c('tiff', 'jpeg'), dev.args=list(bg="white"), fig.width=6,fig.height=6, dpi=500----
library(gridExtra)

g.a = textGrob(.5, unit(1,"npc") - unit(1,"line"), label="a)")
g.b = textGrob(.5, unit(1,"npc") - unit(1,"line"), label="b)")

smear.plot = ggplot(edgeR.results, aes(Reads, logFC, col=color)) + scale_x_log10() + geom_point() + theme(legend.position="none") + color.scale + ylab("Log2(YPGal/YPD) FC") + xlab("Reads per Strain")
grid.arrange(g.a, volcano.plot, g.b, smear.plot, nrow=2, widths=c(.03, .97))


## ----figure_3, dependson=c("go_term_graph", "barseq_processed_graph"), results=FALSE, include=FALSE, dev=c('tiff', 'jpeg'), dev.args=list(bg="white"), fig.width=7,fig.height=7, dpi=500----
p2 = base.plot + geom_line(aes(y=MSE.fit)) + xlim(0, 25e6) + ylab("MSE Compared to Full Experiment")

p4 = base.plot + geom_line(aes(y=fdr.fit)) + xlim(0, 25e6) + ylab("FDR Relative to Full Experiment")

g_legend <- function(a.gplot) {
	tmp <- ggplot_gtable(ggplot_build(a.gplot))
	leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
	legend <- tmp$grobs[[leg]]
	return(legend)
}

legend<-g_legend(p1 + theme(legend.position="top"))

## using grid.arrange for convenience
## could also manually push viewports
library(grid)
library(gridExtra)

g.a = textGrob(.5, unit(1,"npc") - unit(1,"line"), label="a)")
g.b = textGrob(.5, unit(1,"npc") - unit(1,"line"), label="b)")
g.c = textGrob(.5, unit(1,"npc") - unit(1,"line"), label="c)")
g.d = textGrob(.5, unit(1,"npc") - unit(1,"line"), label="d)")
g.n = textGrob(.5, unit(1,"npc") - unit(1,"line"), label="")

f3 = grid.arrange(legend,
	arrangeGrob(g.a, p1 + theme(legend.position="none"), 
                          #g.n, legend,
                          g.b, p2 + theme(legend.position="none"),
                          g.c, p3 + theme(legend.position="none"),
                          g.d, p4 + theme(legend.position="none"),
                          widths=c(.03, .47, .03, .47), nrow=2),
                 nrow=2, heights=c(.1, .9))

print(f3)


## ----figure_4, dependson="by_logFC", include=FALSE, dev=c('tiff', 'jpeg'), dev.args=list(bg="white"), fig.width=5.5,fig.height=5.5, dpi=500----
library(ggplot2)
print(ggplot(percents.logFC.report, aes(depth, percent.fit, lty=variable)) + geom_line(col="red") +
                scale_x_continuous(breaks=seq(0, 16e6, 2e6), limits=c(0, 16e6)) +
                geom_hline(aes(yintercept=.9), lty=2, col="black") +
                scale_y_continuous(breaks=seq(0, 1, .1)) +
                theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
                scale_linetype_manual(name = "Fold Change Threshold", values=c(2, 3, 1)) +
                xlab("Read Depth Per Condition") + ylab("Proportion of Genes Identified As Significant"))


## ----figure_S1, dependson="import_data", include=FALSE-------------------
library(ggplot2)

ggplot(data.frame(x=rowSums(original@counts)), aes(x=x)) + geom_histogram() + xlab("Reads per strain") + scale_x_log10() + ylab("Frequency")


## ----figure_S2, dependson="import_data", fig.height=5, fig.width=7.5, results="hide", include=FALSE----
library(ggplot2)
nonmax = filter.counts(original, max.reads=10000000)
up_totals = rowSums(nonmax@counts[, nonmax@design$tag == 1])
dn_totals = rowSums(nonmax@counts[, nonmax@design$tag == 2])
totals = data.frame(Up=up_totals, Down=dn_totals)
ggplot(totals, aes(x=Up, y=Down)) + geom_point()


## ----prepare_DGE---------------------------------------------------------
library(edgeR)

#unfiltered = combine.counts(combine.counts(original, "tag"), "technical")
#unfiltered = filter.counts(unfiltered, min.reads=100, max.reads=1e6)
#unfiltered = unfiltered[unfiltered@design$condition != "t0"]

maincounts = pooled.counts[pooled.counts@design$condition != "t0"]

dge = DGEList(maincounts@counts[rowSums(maincounts@counts) >= 100, ], group=maincounts@design$condition)
dge = estimateTagwiseDisp(dge)
dge = estimateCommonDisp(dge)


## ----figure_S3, dependson="percentvar", include=FALSE--------------------
library(ggplot2)

print(ggplot(var.melt, aes(quantile, value, col=Method, lty=variable)) + geom_line() +
      xlab("Percentile Threshold") + ylab("Percent of Variance Explained") +
      ylim(0, 100) + scale_x_continuous(breaks=seq(0, .75, .25), lim=c(0, .75)) + scale_colour_discrete(name = "Source"))


## ----figure_S4, dependson="prepare_DGE", include=FALSE-------------------
#par(mfrow=c(2, 1))
plotMeanVar(dge, show.raw.vars=TRUE, show.tagwise.vars=TRUE, xlab="Mean abundance level (log10 scale)", ylab="Variance of abundance level (log10 scale)")
#plotBCV(dge)


## ----figure_S5, dependson="DESeq", include=FALSE-------------------------
library(ggplot2)
ggplot(data=pooled.counts@pvalue, aes(x=edgeR, y=DESeq)) + geom_point()


## ----figure_S6, fig.width=10, dependson="barseq_processed_graph", include=FALSE----
data.2.reps = subset(summaries, summaries$tech.reps == 2)
base.plot + geom_point(data=data.2.reps, aes(y=significant, col=bio.reps), size=1) + geom_line(data=data.2.reps) + xlim(0, 25e6)


## ----table_dir-----------------------------------------------------------
dir.create("tables", showWarnings = FALSE)


## ----table_S3, dependson="import_data"-----------------------------------
bio.names = c("A", "B", "C", "D", "A", "B", "C", "D", "A", "B")
design.names = paste(experiment.design$condition, bio.names[experiment.design$biological], rep(1:2, 10, each=2), c("UP", "DN"), sep="_")
matrix.table = as.data.frame(experiment.counts)
colnames(matrix.table) = design.names

library(org.Sc.sgd.db)
nameindex = as.list(org.Sc.sgdGENENAME)

incl = rowSums(experiment.counts) > 100 & rowSums(experiment.counts) < 1000000

matrix.table = cbind(Strain=rownames(matrix.table), Gene=as.character(nameindex[rownames(matrix.table)]), Included=incl, matrix.table)
matrix.table = matrix.table[order(-matrix.table$Included), ]

write.table(matrix.table, file="tables/table_S3.xls", sep="\t", quote=FALSE, row.names=FALSE)


## ----table_S4, dependson="volcano_plot"----------------------------------
# rearrange rows and columns
edgeR.results = edgeR.results[order(edgeR.results$pvalue), ]

library(org.Sc.sgd.db)
nameindex = as.list(org.Sc.sgdGENENAME)
edgeR.results$Gene = as.character(nameindex[edgeR.results$Deletion])

edgeR.results$Gene[is.na(edgeR.results$Gene)] = edgeR.results$Deletion[is.na(edgeR.results$Gene)]

# reorder, remove and rename some columns
library(qvalue)
edgeR.results$DESeq.qvalue = qvalue(edgeR.results$DESeq.pvalue)$qvalue
edgeR.results = edgeR.results[c("Deletion", "Gene", "pvalue", "qvalue", "DESeq.pvalue", "DESeq.qvalue", "Gal", "Reads", "logFC")]
colnames(edgeR.results) = c("Deletion", "Gene", "edgeR pvalue", "edgeR qvalue", "DESeq pvalue", "DESeq qvalue", "Gal", "Reads", "logFC")
write.table(edgeR.results, file="tables/table_S4.xls", quote=FALSE, row.names=FALSE, sep="\t")


## ----table_S5, dependson="wilcoxon_test"---------------------------------
library(qvalue)

ids = names(gsc.terms)
terms = sapply(gsc.terms, Term)
definitions = sapply(gsc.terms, Definition)
go.term.table = data.frame(GOID=ids, Term=terms, Definition=definitions, Number=colSums(gene.contained), pvalue=wilcoxon.pvalues, qvalue=qvalue(wilcoxon.pvalues)$qvalue)
go.term.table = go.term.table[order(go.term.table$qvalue), ]

write.table(go.term.table, file="tables/table_S5.xls", quote=FALSE, row.names=FALSE, sep="\t")


## ----save_data-----------------------------------------------------------
counts = counts.obj

save.counts = counts[counts@design$condition != "t0"]

save.counts@design$biological = rep(1:4, 2, each=2)
save.counts@design$technical = rep(1:2, 8)

Full = combine.counts(save.counts, c("technical"))

bio.2 = apply(combn(4, 2), 2, function(pair) {
	combine.counts(save.counts[save.counts@design$biological %in% pair], c("technical"))
})
names(bio.2) = apply(combn(4, 2), 2, function(pair) paste("Bio.2.", pair[1], ".", pair[2], sep=""))

bio.3 = apply(combn(4, 3), 2, function(pair) {
	combine.counts(save.counts[save.counts@design$biological %in% pair], c("technical"))
})
names(bio.3) = apply(combn(4, 3), 2, function(pair) paste("Bio.3.", pair[1], ".", pair[2], ".", pair[3], sep=""))


tech.1 = save.counts[save.counts@design$technical == 1]
tech.2 = save.counts[save.counts@design$technical == 2]

subsets = c(list(Full=Full, tech.1=tech.1, tech.2=tech.2), bio.2, bio.3)

# smaller subsets

combinations.3 = combn(4, 3)

for (b in 1:NCOL(combinations.3)) {
	for (tech in 1:2) {
		bio.set = combinations.3[, b]
		sub = list(save.counts[save.counts@design$biological %in% bio.set & save.counts@design$technical == tech])
		names(sub) = paste("Bio.31", paste(bio.set, collapse="."), tech, sep=".")
		subsets = c(subsets, sub)
	}
}

combinations.2 = combn(4, 2)

for (b in 1:NCOL(combinations.2)) {
	for (tech in 1:2) {
		bio.set = combinations.2[, b]
		sub = list(save.counts[save.counts@design$biological %in% bio.set & save.counts@design$technical == tech])
		names(sub) = paste("Bio.21", paste(bio.set, collapse="."), tech, sep=".")
		subsets = c(subsets, sub)
	}
}

# turn each subset into a matrix with column names that are conditions
subsets = lapply(subsets, function(s) {
	ret = s@counts
	colnames(ret) = s@design$condition
	ret
})

save(subsets, file="BarSeq_subsets.Rdata")



