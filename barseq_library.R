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
        #TODO: initialize from file/files
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
