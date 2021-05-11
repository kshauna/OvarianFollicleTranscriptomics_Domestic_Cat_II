#  Early preantral follicles of the domestic cat express a gonadotropic and steroidogenic signalling potential
#### Kehoe, S.<sup>1</sup>, Jewgenow, K.<sup>1</sup>, Johnston, P.R.<sup>2</sup>, Braun, B.C.<sup>1</sup> 
#### <sup>1</sup> Department of Reproduction Biology, Leibniz-Institute for Zoo and Wildlife Research, Alfred-Kowalke-Straße 17, 10315 Berlin, Germany, kehoe@izw-berlin.de 

#### Previously see https://github.com/kshauna/OvarianFollicleTranscriptomics-DomesticCat
 
| Author of repository: | Author of script: | Date:   |
| --------------- | --------------- | --------------- |
| Shauna Kehoe<sup>1</sup> | Paul R. Johnston<sup>2</sup> |  May 03, 2021  |

<sup>1</sup> Department of Reproduction Biology, Leibniz-Institute for Zoo and Wildlife Research (IZW), Alfred-Kowalke-Straße 17, 10315 Berlin, Germany, corresponding author: kehoe@izw-berlin.de 

<sup>2</sup> Berlin Center for Genomics in Biodiversity Research BeGenDiv, Königin-Luise-Straße 6-8, D-14195, Berlin, Germany; Leibniz-Institute of Freshwater Ecology and Inland Fisheries, Müggelseedamm 310, 12587, Berlin, Germany; Freie Universität Berlin, Institut für Biologie, Königin-Luise-Straße 1-3, 14195, Berlin, Germany


## Figure 2: Gene expression of gonadotropin receptors, steroidogenic enzymes and a transporter, and steroid receptor genes in the early preantral follicles of the domestic cat.
### and
## Supplementary Figure 1: Gene expression of gonadotropin receptors, steroidogenic enzymes and a transporter, and steroid receptor genes in the early preantral follicles of the domestic cat. 

#### plotting genes

```
library(DESeq2)
library(dplyr)
library(tidyr)
library(MASS)
library(visreg)
library(cowplot)
library(ggplot2)
```
Accesssing the normalised counts

```
ddc <- counts(dds, normalized = TRUE) %>% as.data.frame
```
Carrying over the gene names
```
ddc$gene_id <- rownames(ddc)
```
Reformatting from wide to long
```
ddc_long <- gather(ddc, key = library, value = count, -gene_id)
```
Adding the sample data
```
ddc_long_meta <- inner_join(ddc_long, data.frame(dds@colData, library = dds@colData[,1]))
```

Taking the normalised counts back out, fitting a GLM, then extracting the coefficients and intervals for plotting. 

```
glm.nb(count ~ type, data = subset(ddc_long_meta, gene_id %in% c("ENSFCAG00000000001.5"))) %>%
  visreg
gene_fit <- glm.nb(count ~ type, data = subset(ddc_long_meta, gene_id %in% c("ENSFCAG00000000001.5"))) %>%
  visreg(plot = FALSE) %$% fit %>% dplyr::select(-count)
# visreg returns logged coefficients/bounds so compute the exponentials
gene_fit$visregFit <- exp(gene_fit$visregFit)
gene_fit$visregUpr <- exp(gene_fit$visregUpr)
gene_fit$visregLwr <- exp(gene_fit$visregLwr)

```

Plotting with ggplot2:

```
ggplot(gene_fit, aes(x = type, y = visregFit)) +
  scale_y_log10(breaks = c(1, 10, 100, 1000, 10000, 100000), limits = c(0.1, 100000)) +
  geom_pointrange(aes(ymin = visregLwr, ymax = visregUpr), size = 1) +
  xlab("Follicle type") + ylab("Gene expression (normalised count)")
```

## Supplementary File 2: The normalised gene count per early preantral follicle sample from the domestic cat. 

Extracting the y-axis values (normalised transcript counts) from Figure 2 and Supplemental Figure 1 for each sample type (PrF, PF, and SF): 
The counts are present as a matrix in a slot of the `DESeqDataSet` S4 object:

```R
head(dds@assays$data$counts)
```
```
                     Sample_1_S1 Sample_2_S2 Sample_3_S3 Sample_4_S4
ENSFCAG00000000001.5         143         176           0        1293
ENSFCAG00000000007.5         630        2288           3        5600
ENSFCAG00000000015.5          30          89           0          40
ENSFCAG00000000022.5         692        1175          46         948
ENSFCAG00000000023.4        1256         706          21         322
ENSFCAG00000000024.5        2513        1973          88         766
                     Sample_5_S5 Sample_6_S6 Sample_7_S7 Sample_8_S8
ENSFCAG00000000001.5         102         374         589         112
ENSFCAG00000000007.5        2057        2669        2882         322
ENSFCAG00000000015.5          57          58          16          17
ENSFCAG00000000022.5         759         706         721         159
ENSFCAG00000000023.4         545         686         670         228
ENSFCAG00000000024.5        1427        1155        2164         475
                     Sample_9_S9
ENSFCAG00000000001.5         129
ENSFCAG00000000007.5         781
ENSFCAG00000000015.5           3
ENSFCAG00000000022.5         197
ENSFCAG00000000023.4         209
ENSFCAG00000000024.5         637
```

Accessing with the `counts` function:

```R
head(counts(dds, normalized = FALSE))
```
```
                     Sample_1_S1 Sample_2_S2 Sample_3_S3 Sample_4_S4
ENSFCAG00000000001.5         143         176           0        1293
ENSFCAG00000000007.5         630        2288           3        5600
ENSFCAG00000000015.5          30          89           0          40
ENSFCAG00000000022.5         692        1175          46         948
ENSFCAG00000000023.4        1256         706          21         322
ENSFCAG00000000024.5        2513        1973          88         766
                     Sample_5_S5 Sample_6_S6 Sample_7_S7 Sample_8_S8
ENSFCAG00000000001.5         102         374         589         112
ENSFCAG00000000007.5        2057        2669        2882         322
ENSFCAG00000000015.5          57          58          16          17
ENSFCAG00000000022.5         759         706         721         159
ENSFCAG00000000023.4         545         686         670         228
ENSFCAG00000000024.5        1427        1155        2164         475
                     Sample_9_S9
ENSFCAG00000000001.5         129
ENSFCAG00000000007.5         781
ENSFCAG00000000015.5           3
ENSFCAG00000000022.5         197
ENSFCAG00000000023.4         209
ENSFCAG00000000024.5         637
```

Normalised:

```R
head(counts(dds, normalized = TRUE))
```
```
                     Sample_1_S1 Sample_2_S2 Sample_3_S3 Sample_4_S4
ENSFCAG00000000001.5   116.04889    95.67859     0.00000   448.24405
ENSFCAG00000000007.5   506.14547  1234.06625    54.54146  1935.50110
ENSFCAG00000000015.5    23.79297    47.92357     0.00000    14.14661
ENSFCAG00000000022.5   559.53148   638.51692   786.57956   330.82780
ENSFCAG00000000023.4  1003.31522   382.37526   352.22692   114.19296
ENSFCAG00000000024.5  1862.85185  1042.89987  1299.16458   302.24317
                     Sample_5_S5 Sample_6_S6 Sample_7_S7 Sample_8_S8
ENSFCAG00000000001.5     58.7268   201.43228  281.502162   223.30534
ENSFCAG00000000007.5   1185.3566  1429.91481 1369.461759   637.62999
ENSFCAG00000000015.5     33.4190    30.62059    7.680903    34.22374
ENSFCAG00000000022.5    440.8851   380.60610  344.792302   317.76698
ENSFCAG00000000023.4    320.6252   364.85721  323.029021   460.64544
ENSFCAG00000000024.5    895.0796   566.98705 1111.323391  1015.48956
                     Sample_9_S9
ENSFCAG00000000001.5  163.566419
ENSFCAG00000000007.5  987.365098
ENSFCAG00000000015.5    3.804836
ENSFCAG00000000022.5  251.041827
ENSFCAG00000000023.4  266.457613
ENSFCAG00000000024.5  825.534020
```

Retrieving counts for a specific gene with the `plotCounts` function:

```R
plotCounts(dds, gene = "ENSFCAG00000000001.5", intgroup = "type", returnData = TRUE, normalized = TRUE)
```
```
                count type
Sample_1_S1 116.54889    A
Sample_2_S2  96.17859    A
Sample_3_S3   0.50000    A
Sample_4_S4 448.74405    B
Sample_5_S5  59.22680    B
Sample_6_S6 201.93228    B
Sample_7_S7 282.00216    C
Sample_8_S8 223.80534    C
Sample_9_S9 164.06642    C
```

#### sessionInfo()
```
> sessionInfo()
R version 3.4.4 (2018-03-15)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 18.04.3 LTS

Matrix products: default
BLAS: /usr/lib/x86_64-linux-gnu/openblas/libblas.so.3
LAPACK: /usr/lib/x86_64-linux-gnu/libopenblasp-r0.2.20.so

locale:
 [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C               LC_TIME=en_GB.UTF-8        LC_COLLATE=en_GB.UTF-8    
 [5] LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_GB.UTF-8    LC_PAPER=en_GB.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] biomaRt_2.34.2

loaded via a namespace (and not attached):
  [1] bitops_1.0-6               matrixStats_0.55.0         bit64_0.9-7                httr_1.4.1                
  [5] RColorBrewer_1.1-2         progress_1.2.2             GenomeInfoDb_1.14.0        tools_3.4.4               
  [9] backports_1.1.5            R6_2.4.1                   rpart_4.1-15               Hmisc_4.3-0               
 [13] DBI_1.0.0                  lazyeval_0.2.2             BiocGenerics_0.24.0        colorspace_1.4-1          
 [17] nnet_7.3-12                withr_2.1.2                prettyunits_1.0.2          tidyselect_0.2.5          
 [21] gridExtra_2.3              DESeq2_1.18.1              bit_1.1-14                 compiler_3.4.4            
 [25] cli_1.1.0                  Biobase_2.38.0             htmlTable_1.13.3           DelayedArray_0.4.1        
 [29] scales_1.1.0               checkmate_1.9.4            genefilter_1.60.0          stringr_1.4.0             
 [33] digest_0.6.23              foreign_0.8-72             DOSE_3.4.0                 XVector_0.18.0            
 [37] base64enc_0.1-3            pkgconfig_2.0.3            htmltools_0.4.0            sessioninfo_1.1.1         
 [41] htmlwidgets_1.5.1          rlang_0.4.2                rstudioapi_0.10            RSQLite_2.1.4             
 [45] BiocParallel_1.12.0        acepack_1.4.1              GOSemSim_2.4.1             dplyr_0.8.3               
 [49] RCurl_1.95-4.12            magrittr_1.5               GO.db_3.5.0                GenomeInfoDbData_1.0.0    
 [53] Formula_1.2-3              Matrix_1.2-18              Rcpp_1.0.3                 munsell_0.5.0             
 [57] S4Vectors_0.16.0           lifecycle_0.1.0            stringi_1.5.3              SummarizedExperiment_1.8.1
 [61] zlibbioc_1.24.0            plyr_1.8.4                 qvalue_2.10.0              grid_3.4.4                
 [65] blob_1.2.0                 parallel_3.4.4             DO.db_2.9                  crayon_1.3.4              
 [69] lattice_0.20-38            splines_3.4.4              annotate_1.56.2            hms_0.5.2                 
 [73] locfit_1.5-9.1             zeallot_0.1.0              knitr_1.26                 pillar_1.4.2              
 [77] fgsea_1.4.1                igraph_1.2.4.2             GenomicRanges_1.30.3       geneplotter_1.56.0        
 [81] reshape2_1.4.3             stats4_3.4.4               fastmatch_1.1-0            XML_3.98-1.20             
 [85] glue_1.3.1                 latticeExtra_0.6-28        data.table_1.12.6          BiocManager_1.30.10       
 [89] vctrs_0.2.0                gtable_0.3.0               purrr_0.3.3                assertthat_0.2.1          
 [93] ggplot2_3.2.1              xfun_0.11                  xtable_1.8-4               survival_3.1-8            
 [97] tibble_2.1.3               rvcheck_0.1.7              AnnotationDbi_1.40.0       memoise_1.1.0             
[101] IRanges_2.12.0             cluster_2.1.0       
```
#### Packages mentioned in the manuscript:

**cowplot** https://cran.r-project.org/web/packages/cowplot/index.html Author: Claus O. Wilke ORCID iD [aut, cre]

**DESeq2** https://bioconductor.org/packages/release/bioc/html/DESeq2.html Love MI, Huber W, Anders S (2014). “Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2.” Genome Biology, 15, 550. doi: 10.1186/s13059-014-0550-8.

**dplyr** https://cran.r-project.org/web/packages/dplyr/index.html Author:	Hadley Wickham ORCID iD [aut, cre], Romain François ORCID iD [aut], Lionel Henry [aut], Kirill Müller ORCID iD [aut], RStudio [cph, fnd]

**ggplot2** https://cran.r-project.org/web/packages/ggplot2/index.html Author: Hadley Wickham [aut, cre], Winston Chang [aut], Lionel Henry [aut], Thomas Lin Pedersen [aut], Kohske Takahashi [aut], Claus Wilke [aut], Kara Woo [aut], Hiroaki Yutani [aut], Dewey Dunnington [aut], RStudio [cph]

**IHW** http://bioconductor.org/packages/release/bioc/html/IHW.html Ignatiadis N, Klaus B, Zaugg J, Huber W (2016). “Data-driven hypothesis weighting increases detection power in genome-scale multiple testing.” Nature Methods. doi: 10.1038/nmeth.3885 & Ignatiadis N, Huber W (2017). “Covariate-powered weighted multiple testing with false discovery rate control.” arXiv. doi: arXiv:1701.05179.

**MASS** https://cran.r-project.org/web/packages/MASS/index.html Author: Brian Ripley [aut, cre, cph], Bill Venables [ctb], Douglas M. Bates [ctb], Kurt Hornik [trl] (partial port ca 1998), Albrecht Gebhardt [trl] (partial port ca 1998), David Firth [ctb]

**tidyr** https://cran.r-project.org/web/packages/tidyr/index.html Author: Hadley Wickham [aut, cre], Lionel Henry [aut], RStudio [cph]

**visreg** https://cran.r-project.org/web/packages/visreg/visreg.pdf Author: Patrick Breheny, Woodrow Burchett     
