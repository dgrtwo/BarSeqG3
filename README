This document contains instructions for reproducing the manuscript:

David G. Robinson, Wei Chen, John D. Storey and David Gresham. "Design and
Analysis of Bar-seq Experiments." G3. 10.1534/g3.113.008565

We use knitr (http://yihui.name/knitr/) to generate the manuscript from the
input data, and compile the manuscript with LaTeX.

SETUP
================

This manuscript requires R at least version 2.15 (recommended 3.0) to reproduce,
as well as pdflatex to compile. You'll also need to install the
following packages from CRAN and BioConductor:

    install.packages(c("plyr", "data.table", "reshape", "ggplot2", "gplots",
                       "colorRamps", "gridExtra", "xtable"))

    source("http://bioconductor.org/biocLite.R")
    biocLite(c("GSEABase", "org.Sc.sgd.db", "GO.db", "edgeR", "DESeq",
               "qvalue"))

You also have to install the eigenR2 package from here:

http://www.genomine.org/eigenr2/

REPRODUCTION
================

To reproduce the manuscript, run the following lines on your terminal:

Rscript -e "library(knitr); knit('Robinson_2013.Rnw')"   # reproduce manuscript
pdflatex Robinson_2013.tex                               # compile LaTeX

In a few minutes it should generate the full PDF of the manuscript, along with
all figures except 1A and all tables except Tables 1 and 2.

SESSION
================

The manuscript was compiled with the following sessionInfo():

R version 3.0.1 (2013-05-16)
Platform: x86_64-apple-darwin10.8.0 (64-bit)

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] grid      parallel  stats     graphics  grDevices utils     datasets 
[8] methods   base     

other attached packages:
 [1] qvalue_1.34.0        gridExtra_0.9.1      scales_0.2.3        
 [4] colorRamps_2.3       gplots_2.11.3        MASS_7.3-28         
 [7] KernSmooth_2.23-10   caTools_1.14         gdata_2.13.2        
[10] gtools_3.0.0         eigenR2_1.0          xtable_1.7-1        
[13] reshape_0.8.4        plyr_1.8             ggplot2_0.9.3.1.99  
[16] GO.db_2.9.0          org.Sc.sgd.db_2.9.1  RSQLite_0.11.4      
[19] DBI_0.2-7            GSEABase_1.22.0      graph_1.38.3        
[22] annotate_1.38.0      AnnotationDbi_1.22.6 DESeq_1.12.1        
[25] lattice_0.20-23      locfit_1.5-9.1       Biobase_2.20.1      
[28] BiocGenerics_0.6.0   edgeR_3.2.4          limma_3.16.7        
[31] data.table_1.8.8     knitr_1.4.1         

loaded via a namespace (and not attached):
 [1] bitops_1.0-6       colorspace_1.2-2   dichromat_2.0-0    digest_0.6.3      
 [5] evaluate_0.4.7     formatR_0.9        genefilter_1.42.0  geneplotter_1.38.0
 [9] gtable_0.1.2       IRanges_1.18.3     labeling_0.2       munsell_0.4.2     
[13] proto_0.3-10       RColorBrewer_1.0-5 reshape2_1.2.2     splines_3.0.1     
[17] stats4_3.0.1       stringr_0.6.2      survival_2.37-4    tcltk_3.0.1       
[21] tools_3.0.1        XML_3.95-0.2      
