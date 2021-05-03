# Analysis workflow
## Previously see https://github.com/kshauna/OvarianFollicleTranscriptomics-DomesticCat

>We're interested to know if transcript count data can be extracted from two sections: import and create sample info and plotting genes. For the first section we would like to extract the data which summarizes transcript abundance, counts, and length prior to the create DESeqDataSet section. My supervisor thinks that there may be transcripts in this data that are not included in the DESeq2 output - as she phrased it: does DESeq2 deal with all transcripts independent if they are DEGs or not? In my opinion, I would assume that the significant and non-significant DESeq2 gene lists are all the transcripts found full stop but I could be wrong! 


I'm not really sure what they mean. But all the transcripts in `tx.fa` seem to belong to a gene:

```sh
grep "^>" tx.fa -c
```
```
42556
```
```sh
grep "^>" tx.fa | grep "gene:" -c
```
```
42556
```

The transcript-level quantification is used to make the gene-level differential expression more sensitive/accurate. By the time the `salmon` quantification results are read into `R`, they are already summarized at the gene level. The `txi` list contains 3 matrices (`abundance`, `counts`, and `length`) and a vector (`countsFromAbundance`).

```R
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
str(txi)
```
```
List of 4
 $ abundance          : num [1:26003, 1:9] 1.482 6.185 0.679 7.649 26.476 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : chr [1:26003] "ENSFCAG00000000001.5" "ENSFCAG00000000007.5" "ENSFCAG00000000015.5" "ENSFCA
G00000000022.5" ...
  .. ..$ : chr [1:9] "Sample_1_S1" "Sample_2_S2" "Sample_3_S3" "Sample_4_S4" ...
 $ counts             : num [1:26003, 1:9] 143 630 30 692 1256 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : chr [1:26003] "ENSFCAG00000000001.5" "ENSFCAG00000000007.5" "ENSFCAG00000000015.5" "ENSFCA
G00000000022.5" ...
  .. ..$ : chr [1:9] "Sample_1_S1" "Sample_2_S2" "Sample_3_S3" "Sample_4_S4" ...
 $ length             : num [1:26003, 1:9] 2313 2443 1060 2170 1138 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : chr [1:26003] "ENSFCAG00000000001.5" "ENSFCAG00000000007.5" "ENSFCAG00000000015.5" "ENSFCA
G00000000022.5" ...
  .. ..$ : chr [1:9] "Sample_1_S1" "Sample_2_S2" "Sample_3_S3" "Sample_4_S4" ...
 $ countsFromAbundance: chr "no"
```

```R
head(txi$abundance)
```
```
                     Sample_1_S1 Sample_2_S2 Sample_3_S3 Sample_4_S4
ENSFCAG00000000001.5    1.482485    1.271409    0.000000    5.663440
ENSFCAG00000000007.5    6.185080   15.686640    0.532664   23.392618
ENSFCAG00000000015.5    0.678950    1.422530    0.000000    0.399262
ENSFCAG00000000022.5    7.648640    9.079310    8.593270    4.472780
ENSFCAG00000000023.4   26.476400   10.496200    7.428500    2.980420
ENSFCAG00000000024.5  280.334000  163.253000  156.250000   44.985200
                     Sample_5_S5 Sample_6_S6 Sample_7_S7 Sample_8_S8
ENSFCAG00000000001.5   0.7836445    2.606769    3.627981    2.981844
ENSFCAG00000000007.5  15.1304680   17.701200   16.883139    8.144680
ENSFCAG00000000015.5   0.9961330    0.885168    0.221123    1.020830
ENSFCAG00000000022.5   6.2953190    5.270570    4.754994    4.540500
ENSFCAG00000000023.4   8.8379800    9.753650    8.599950   12.706400
ENSFCAG00000000024.5 140.7000000   86.435800  168.723000  159.738000
                     Sample_9_S9
ENSFCAG00000000001.5    1.997554
ENSFCAG00000000007.5   11.534600
ENSFCAG00000000015.5    0.103796
ENSFCAG00000000022.5    3.280650
ENSFCAG00000000023.4    6.722110
ENSFCAG00000000024.5  118.764000
```

The results contain the same genes as the matrices:

```R
all(rownames(txi$abundance) == rownames(res_AB))
```
```
[1] TRUE
```
```R
all(rownames(txi$counts) == rownames(res_AB))
```
```
[1] TRUE
```
```R
all(rownames(txi$length) == rownames(res_AB))
```
```
[1] TRUE
```

So all genes are in there (and all the transcripts belonged to a gene). If they are interested in the transcript abundance/counts then those are in the `salmon` output directories:

```sh
head Sample_1_S1/quant.sf
```
```
Name    Length  EffectiveLength TPM     NumReads
ENSFCAT00000050682.1    300     125.391 0       0
ENSFCAT00000064547.1    324     139.607 2.73975 15.95
ENSFCAT00000061608.1    324     139.607 1.12081 6.525
ENSFCAT00000054676.1    285     117.225 0       0
ENSFCAT00000054867.1    324     139.607 1.12081 6.525
ENSFCAT00000061928.1    285     117.225 0       0
ENSFCAT00000056337.1    285     117.225 0       0
ENSFCAT00000027343.3    423     212.461 1.34608 11.9259
ENSFCAT00000032528.3    423     212.461 3.84593 34.0741
```

>For the second section, we would like to extract the y-axis values from the dot plots so for each sample type (PrF, PF, and SF) we would like to extract the normalised transcript counts. Do you know how this can be done?

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

It's easier to access them with the `counts` function:

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

Or normalized:

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

If you just want the counts for a specific gene then they can also be accessed with the `plotCounts` function:

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
## plotting genes
```
library(DESeq2)
library(dplyr)
library(tidyr)
library(MASS)
library(visreg)
library(cowplot)
# get the normalized counts out
ddc <- counts(dds, normalized = TRUE) %>% as.data.frame
# carry over the gene names
ddc$gene_id <- rownames(ddc)
# reformat from wide to long
ddc_long <- gather(ddc, key = library, value = count, -gene_id)
# add the sample data
ddc_long_meta <- inner_join(ddc_long, data.frame(dds@colData, library = dds@colData[,1]))
# This is a bit backwards because it involves taking the normalized counts back out, fitting a GLM, 
# and then extracting the coefficients and intervals for plotting. 
# DESeq has already fitted GLMs for each gene but the data aren't in a very useful state.
# plot a given gene
glm.nb(count ~ type, data = subset(ddc_long_meta, gene_id %in% c("ENSFCAG00000028290.3"))) %>%
  visreg
bmp15_fit <- glm.nb(count ~ type, data = subset(ddc_long_meta, gene_id %in% c("ENSFCAG00000028290.3"))) %>%
  visreg(plot = FALSE) %$% fit %>% dplyr::select(-count)
# visreg returns logged coefficients/bounds so compute the exponentials
bmp15_fit$visregFit <- exp(bmp15_fit$visregFit)
bmp15_fit$visregUpr <- exp(bmp15_fit$visregUpr)
bmp15_fit$visregLwr <- exp(bmp15_fit$visregLwr)
# 
ggplot(bmp15_fit, aes(x = type, y = visregFit)) +
  scale_y_log10(breaks = c(1, 10, 100, 1000, 10000, 100000), limits = c(0.1, 100000)) +
  geom_pointrange(aes(ymin = visregLwr, ymax = visregUpr), size = 1) +
  xlab("Follicle type") + ylab("Gene expression (normalized count)")
# With n=3 it's probably worth putting the actual data points on too.
# DESeq2 has already fitted glms. So it would be possible to use those values instead. 
# The coefficients it provides aren't in a very useful format:
coef(dds, SE = FALSE) %>% subset(., rownames(.) == "ENSFCAG00000028290.3")
# So if you want to use them directly from DESeq, I guess the easiest thing to do would be to remove the intercept. i.e. run it again without an intercept:
dds2 <- dds
dds2@design <- ~ type + 0
dds2 <- DESeq(dds2)
# This would give more useful coefficients:
coef(dds2, SE = FALSE) %>% subset(., rownames(.) == "ENSFCAG00000028290.3")
# but they are log2, so:
2^coef(dds2, SE = FALSE) %>% subset(., rownames(.) == "ENSFCAG00000028290.3") %>% t
# These are almost identical to the coefficients from glm.nb:
bmp15_fit
# The intervals could will be a bit different for some genes (e.g. because of the shrinkage of the dispersion parameter in DESeq). 
# For "ENSFCAG00000028290.3", they are basically the same.
ce <- coef(dds2, SE = FALSE) %>% subset(., rownames(.) == "ENSFCAG00000028290.3") %>% t %>% as.data.frame
cee <- coef(dds2, SE = TRUE) %>% subset(., rownames(.) == "ENSFCAG00000028290.3") %>% t %>% as.data.frame
bmp15_fd <- data.frame(type = rownames(ce), coef = ce[,1], se = cee[,1])
bmp15_fd$lower <- bmp15_fd$coef - (2*bmp15_fd$se)
bmp15_fd$upper <- bmp15_fd$coef + (2*bmp15_fd$se)
# From DESeq:
bmp15_fd %>% dplyr::select(-se) %>% mutate_at(c("coef", "lower", "upper"), function(x) (2^x))
# From basic glm:
bmp15_fit
```
## Functional Annotation Clustering
## DAVID
```
* Input ENTREZ IDs and select for Felis catus
* DAVID accepts: ENSEMBL_GENE_ID, ENSEMBL_TRANSCRIPT_ID, ENTREZ_GENE_ID
* Input DEG ENTREZs that were created from the 'no version ENSFCAGs' with biomaRt, down- and up-regulated were not input seperately
* Deselect defaults and select  GOTERM_BP_FAT,  GOTERM_CC_FAT, GOTERM_MF_FAT, and KEGG pathways
* Select Functional Annotation Clustering
* Options: 
        * Kappa Similarity: Similarity Term Overlap 3 and Similarity Threshold 0.60
        * Classification: Initial Group Membership 3 Final Group Membership 3 Multiple Linkage Threshold 0.50
        * Enrichment Thresholds (EASE) 0.2
        * Display: Fold Change, Benjamini, and FDR
        * Outputs: Annotation Cluster No., Enrichment Score, Count, P_Value, Benjamini, Fold Change, FDR.
```
## MetaScape
Metascape's current architecture does not support domestic cat thus, we converted Felis Catus ENTREZ IDs into human orthologs, and then proceeded with Metascape analysis. 

* go to https://biit.cs.ut.ee/gprofiler/gorth.cgi
* use g:Orth tool
* convert those identifiers from your species to human (notice, your IDs should be space separated, comma does not work) 
* then you can take the converted ENSG (Ensembl Gene IDs) IDs for Metascape analysis
* Metascape's 'Express analysis' was selected

## sessionInfo()
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
# Cited packages
**bioMart** https://bioconductor.org/packages/release/bioc/html/biomaRt.html Durinck S, Spellman P, Birney E, Huber W (2009). “Mapping identifiers for the integration of genomic datasets with the R/Bioconductor package biomaRt.” Nature Protocols, 4, 1184–1191 & Durinck S, Moreau Y, Kasprzyk A, Davis S, De Moor B, Brazma A, Huber W (2005). “BioMart and Bioconductor: a powerful link between biological databases and microarray data analysis.” Bioinformatics, 21, 3439–3440.

**cowplot** https://cran.r-project.org/web/packages/cowplot/index.html Author: Claus O. Wilke ORCID iD [aut, cre]

**DESeq2** https://bioconductor.org/packages/release/bioc/html/DESeq2.html Love MI, Huber W, Anders S (2014). “Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2.” Genome Biology, 15, 550. doi: 10.1186/s13059-014-0550-8.

**dplyr** https://cran.r-project.org/web/packages/dplyr/index.html Author:	Hadley Wickham ORCID iD [aut, cre], Romain François ORCID iD [aut], Lionel Henry [aut], Kirill Müller ORCID iD [aut], RStudio [cph, fnd]

**FastQC** https://www.bioinformatics.babraham.ac.uk/projects/fastqc/ Andrews, S. (2010). FastQC:  A Quality Control Tool for High Throughput Sequence Data [Online]. 

**ggplot2** https://cran.r-project.org/web/packages/ggplot2/index.html Author: Hadley Wickham [aut, cre], Winston Chang [aut], Lionel Henry [aut], Thomas Lin Pedersen [aut], Kohske Takahashi [aut], Claus Wilke [aut], Kara Woo [aut], Hiroaki Yutani [aut], Dewey Dunnington [aut], RStudio [cph]

**hexbin** https://cran.r-project.org/web/packages/hexbin/hexbin.pdf Author: Dan Carr <dcarr@voxel.galaxy.gmu.edu>, ported by Nicholas
Lewin-Koh and Martin Maechler <maechler@stat.math.ethz.ch>, contains copies of lattice functions written by Deepayan Sarkar <deepayan.sarkar@r-project.org>

**IHW** http://bioconductor.org/packages/release/bioc/html/IHW.html Ignatiadis N, Klaus B, Zaugg J, Huber W (2016). “Data-driven hypothesis weighting increases detection power in genome-scale multiple testing.” Nature Methods. doi: 10.1038/nmeth.3885 & Ignatiadis N, Huber W (2017). “Covariate-powered weighted multiple testing with false discovery rate control.” arXiv. doi: arXiv:1701.05179.

**MASS** https://cran.r-project.org/web/packages/MASS/index.html Author: Brian Ripley [aut, cre, cph], Bill Venables [ctb], Douglas M. Bates [ctb], Kurt Hornik [trl] (partial port ca 1998), Albrecht Gebhardt [trl] (partial port ca 1998), David Firth [ctb]

**pheatmap** https://cran.r-project.org/web/packages/pheatmap/index.html Author: Raivo Kolde

**RColorBrewer** https://cran.r-project.org/web/packages/RColorBrewer/index.html Author: Erich Neuwirth [aut, cre]

**readr** https://cran.r-project.org/web/packages/readr/index.html Author: Hadley Wickham [aut], Jim Hester [aut, cre], Romain Francois [aut], R Core Team [ctb] (Date time code adapted from R), RStudio [cph, fnd], Jukka Jylänki [ctb, cph] (grisu3 implementation), Mikkel Jørgensen [ctb, cph] (grisu3 implementation)

**Salmon** https://combine-lab.github.io/salmon/ Patro, R., Duggal, G., Love, M. I., Irizarry, R. A., & Kingsford, C. (2017). Salmon provides fast and bias-aware quantification of transcript expression. Nature Methods.

**tidyr** https://cran.r-project.org/web/packages/tidyr/index.html Author: Hadley Wickham [aut, cre], Lionel Henry [aut], RStudio [cph]

**tximport** https://bioconductor.org/packages/release/bioc/html/tximport.html Soneson C, Love MI, Robinson MD (2015). “Differential analyses for RNA-seq: transcript-level estimates improve gene-level inferences.” F1000Research, 4. doi: 10.12688/f1000research.7563.1.

**visreg** https://cran.r-project.org/web/packages/visreg/visreg.pdf Author: Patrick Breheny, Woodrow Burchett
      
