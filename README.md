## Recurrent implantation failure endometrium transcriptomics

The repository containes *renv* project with code and materials for the study:

*A note on the importance of considering menstrual cycle progression when interpreting endometrium gene expression data: an example for recurrent implantation failure issue*, by Pavel Deryabin and Aleksandra Borodkina
Mechanisms of Cellular Senescence Research Laboratory
Institute of Cytology of the Russian Academy of Science
194064 Tikhoretsky av. 4, St. Petersburg, Russia

Analysis of bulk microarray datasets is fully reproducible, analysis of bulk and single cell RNA sequencing datasets can be reproduced from the stage of processed expression matrices

```r
writeLines(capture.output(sessionInfo()), 'sessionInfo.txt')
```

```r
R version 4.5.0 (2025-04-11)
Platform: x86_64-pc-linux-gnu
Running under: Ubuntu 24.04.2 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.12.0 
LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.12.0  LAPACK version 3.12.0

locale:
 [1] LC_CTYPE=ru_RU.UTF-8       LC_NUMERIC=C               LC_TIME=ru_RU.UTF-8        LC_COLLATE=ru_RU.UTF-8    
 [5] LC_MONETARY=ru_RU.UTF-8    LC_MESSAGES=ru_RU.UTF-8    LC_PAPER=ru_RU.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=ru_RU.UTF-8 LC_IDENTIFICATION=C       

time zone: Europe/Moscow
tzcode source: system (glibc)

attached base packages:
[1] stats4    splines   stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] pals_1.10                   lemon_0.5.0                 ggVennDiagram_1.5.2         RColorBrewer_1.1-3         
 [5] patchwork_1.3.2             ggExtra_0.10.1              lubridate_1.9.4             forcats_1.0.0              
 [9] purrr_1.0.4                 readr_2.1.5                 tidyr_1.3.1                 tibble_3.3.0               
[13] ggplot2_4.0.1               tidyverse_2.0.0             dplyr_1.1.4                 stringr_1.5.1              
[17] msigdbr_24.1.0              clusterProfiler_4.16.0      UCell_2.13.1                harmony_1.2.3              
[21] Rcpp_1.0.12.4               DoubletFinder_2.0.4         SeuratWrappers_0.4.0        Seurat_5.3.0               
[25] SeuratObject_5.1.0          sp_2.2-0                    tximport_1.36.0             gprofiler2_0.2.3           
[29] biomaRt_2.64.0              DESeq2_1.48.1               SummarizedExperiment_1.38.1 MatrixGenerics_1.20.0      
[33] matrixStats_1.5.0           GenomicRanges_1.60.0        GenomeInfoDb_1.44.0         tximeta_1.26.1             
[37] org.Hs.eg.db_3.21.0         AnnotationDbi_1.70.0        IRanges_2.42.0              S4Vectors_0.46.0           
[41] limma_3.64.1                endest_0.1.1                preprocessCore_1.70.0       sva_3.56.0                 
[45] BiocParallel_1.42.1         genefilter_1.90.0           mgcv_1.9-3                  nlme_3.1-168               
[49] psych_2.5.3                 GEOquery_2.76.0             Biobase_2.68.0              BiocGenerics_0.54.0        
[53] generics_0.1.4             

loaded via a namespace (and not attached):
  [1] R.methodsS3_1.8.2           dichromat_2.0-0.1           progress_1.2.3              goftest_1.2-3              
  [5] Biostrings_2.76.0           vctrs_0.6.5                 ggtangle_0.0.6              spatstat.random_3.4-1      
  [9] digest_0.6.37               png_0.1-8                   ggrepel_0.9.6               deldir_2.0-4               
 [13] parallelly_1.45.0           MASS_7.3-65                 reshape2_1.4.4              httpuv_1.6.16              
 [17] qvalue_2.40.0               withr_3.0.2                 xfun_0.52                   ggfun_0.1.8                
 [21] survival_3.8-3              memoise_2.0.1               gson_0.1.0                  tidytree_0.4.6             
 [25] zoo_1.8-14                  pbapply_1.7-2               R.oo_1.27.1                 prettyunits_1.2.0          
 [29] KEGGREST_1.48.0             promises_1.3.3              httr_1.4.7                  restfulr_0.0.15            
 [33] globals_0.18.0              fitdistrplus_1.2-2          rstudioapi_0.17.1           UCSC.utils_1.4.0           
 [37] miniUI_0.1.2                DOSE_4.2.0                  babelgene_22.9              curl_6.3.0                 
 [41] fields_16.3.1               polyclip_1.10-7             GenomeInfoDbData_1.2.14     SparseArray_1.8.0          
 [45] xtable_1.8-4                evaluate_1.0.3              S4Arrays_1.8.1              BiocFileCache_2.16.0       
 [49] hms_1.1.3                   irlba_2.3.5.1               colorspace_2.1-1            filelock_1.0.3             
 [53] ROCR_1.0-11                 reticulate_1.42.0           spatstat.data_3.1-6         magrittr_2.0.3             
 [57] lmtest_0.9-40               later_1.4.2                 ggtree_3.16.0               lattice_0.22-7             
 [61] mapproj_1.2.12              spatstat.geom_3.4-1         future.apply_1.20.0         scattermore_1.2            
 [65] XML_3.99-0.18               cowplot_1.1.3               RcppAnnoy_0.0.22            pillar_1.10.2              
 [69] compiler_4.5.0              RSpectra_0.16-2             stringi_1.8.7               tensor_1.5                 
 [73] GenomicAlignments_1.44.0    plyr_1.8.9                  crayon_1.5.3                abind_1.4-8                
 [77] BiocIO_1.18.0               gridGraphics_0.5-1          locfit_1.5-9.12             bit_4.6.0                  
 [81] fastmatch_1.1-6             codetools_0.2-20            plotly_4.10.4               mime_0.13                  
 [85] fastDummies_1.7.5           dbplyr_2.5.0                knitr_1.50                  blob_1.2.4                 
 [89] BiocVersion_3.21.1          AnnotationFilter_1.32.0     fs_1.6.6                    listenv_0.9.1              
 [93] ggplotify_0.1.2             Matrix_1.7-3                statmod_1.5.0               tzdb_0.5.0                 
 [97] pkgconfig_2.0.3             tools_4.5.0                 cachem_1.1.0                RSQLite_2.4.1              
[101] viridisLite_0.4.2           DBI_1.2.3                   fastmap_1.2.0               scales_1.4.0               
[105] grid_4.5.0                  ica_1.0-3                   Rsamtools_2.24.0            AnnotationHub_3.16.0       
[109] BiocManager_1.30.26         dotCall64_1.2               RANN_2.6.2                  farver_2.1.2               
[113] yaml_2.3.10                 rtracklayer_1.68.0          cli_3.6.5                   txdbmaker_1.4.1            
[117] lifecycle_1.0.4             uwot_0.2.3                  annotate_1.86.0             timechange_0.3.0           
[121] gtable_0.3.6                rjson_0.2.23                ggridges_0.5.6              progressr_0.15.1           
[125] parallel_4.5.0              ape_5.8-1                   jsonlite_2.0.0              edgeR_4.6.2                
[129] RcppHNSW_0.6.0              bitops_1.0-9                bit64_4.6.0-1               assertthat_0.2.1           
[133] Rtsne_0.17                  yulab.utils_0.2.0           spatstat.utils_3.1-4        BiocNeighbors_2.2.0        
[137] GOSemSim_2.34.0             spatstat.univar_3.1-3       R.utils_2.13.0              lazyeval_0.2.2             
[141] shiny_1.10.0                htmltools_0.5.8.1           enrichplot_1.28.2           GO.db_3.21.0               
[145] sctransform_0.4.2           rappdirs_0.3.3              ensembldb_2.32.0            glue_1.8.0                 
[149] spam_2.11-1                 httr2_1.1.2                 XVector_0.48.0              RCurl_1.98-1.17            
[153] treeio_1.32.0               mnormt_2.1.1                gridExtra_2.3               igraph_2.1.4               
[157] R6_2.6.1                    SingleCellExperiment_1.30.1 GenomicFeatures_1.60.0      cluster_2.1.8.1            
[161] aplot_0.2.6                 DelayedArray_0.34.1         tidyselect_1.2.1            ProtGenerics_1.40.0        
[165] maps_3.4.3                  xml2_1.3.8                  future_1.58.0               rsvd_1.0.5                 
[169] KernSmooth_2.23-26          S7_0.2.0                    data.table_1.17.6           htmlwidgets_1.6.4          
[173] fgsea_1.34.0                rlang_1.1.6                 spatstat.sparse_3.1-0       spatstat.explore_3.4-3     
[177] remotes_2.5.0               rentrez_1.2.4              
```