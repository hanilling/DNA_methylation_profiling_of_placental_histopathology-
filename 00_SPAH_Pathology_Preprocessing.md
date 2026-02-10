---
title: "SPAH Pathology Project Processing"
author: "Hannah Illing"
date: "2026 February 09"
output: 
  html_document:
    code_folding: show
    keep_md: yes
    toc: true  
    toc_depth: 5
    toc_float: 
      collapsed: true 
      smooth_scroll: true
---

## 0.0 Set-up

### 0.1 Packages


``` r
library(here)
library(tidyverse)
library(wateRmelon)
library(readxl)
library(minfi)

#load the probe information
probeInfo <- as.data.frame(cbind(IlluminaHumanMethylationEPICanno.ilm10b5.hg38::Manifest,             
                                 IlluminaHumanMethylationEPICanno.ilm10b5.hg38::Locations,              
                                 IlluminaHumanMethylationEPICanno.ilm10b5.hg38::Other)) 
probeInfo$probeID <- rownames(probeInfo)

chrXprobes <- probeInfo %>% filter(chr == "chrX")
chrYprobes <- probeInfo %>% filter(chr == "chrY")
autoProbes <- probeInfo %>% filter(!(chr %in% c("chrX", "chrY")))

set.seed(4815) #set seed for script so we get the same result every time it is run
```

## 1.0 Read in the Data 

Compilation of sample sheet, metadata and clinical information files is not shown in this script. Epivariable estimation and quality control checks were performed in another script. 

In the SPAH dataset there were 520 samples in total run on the Illumina EPIC v1.0 array. Of those 520 samples there were 509 unique placentas sampled and 11 duplicates. Prior to this script 10 samples were removed during the initial QC checks for this dataset: 8 samples were run as a preliminary test and used a different sampling strategy, thus were removed; 2 samples failed numerous QC checks (low fluorescence intensity, Illumina control probes, >1% of probes failing). At the end of this script I remove an additional 17 samples that failed other checks (n=5 sex-mismatches, n=1 likely CPM17, n=8 different sampling strategy, n=3 duplicates).


``` r
metadata_SPAH <- readRDS(here::here("C. Processing/02_Outputs/00_Metadata/SPAH_updated_metadata.rds"))
dim(metadata_SPAH) #507 x 95
```

```
## [1] 507  94
```

``` r
array_qc_ss <- readRDS(here::here("C. Processing/02_Outputs/00_Metadata/SPAH_updated_array_qc_ss.rds")) 
dim(array_qc_ss) #510 x 51, includes 3 replicate samples
```

```
## [1] 510  51
```

``` r
#check the ids match
table(metadata_SPAH$ID %in% array_qc_ss$Sample_ID)
```

```
## 
## TRUE 
##  507
```

``` r
metadata_SPAH_ss <- array_qc_ss %>% left_join(metadata_SPAH, by=c("Sample_ID"="ID"))
```

## 2.0 Normalization 


``` r
rgset_raw <- readRDS("//fs/teams/RobinsonLab/ROBLAB6 InfiniumSequenom/EPIC Batch Analysis Code and QC/EPICv1_Batch10_MILLER_code_qc/02_Output/00_Miller_ArrayQC/rgset_all_samples.rds")
dim(rgset_raw) #1051943     520
```

```
## [1] 1051943     520
```

``` r
rgset_raw <- rgset_raw[,metadata_SPAH_ss$rgset_colnames]
dim(rgset_raw) # x 510
```

```
## [1] 1051943     510
```

``` r
mset_noob <- preprocessNoob(rgset_raw)
```

```
## Loading required package: IlluminaHumanMethylationEPICmanifest
```

```
## Loading required package: IlluminaHumanMethylationEPICanno.ilm10b4.hg19
```

```
## 
## Attaching package: 'IlluminaHumanMethylationEPICanno.ilm10b4.hg19'
```

```
## The following objects are masked from 'package:IlluminaHumanMethylation450kanno.ilmn12.hg19':
## 
##     Islands.UCSC, Locations, Manifest, Other, SNPs.132CommonSingle,
##     SNPs.135CommonSingle, SNPs.137CommonSingle, SNPs.138CommonSingle,
##     SNPs.141CommonSingle, SNPs.142CommonSingle, SNPs.144CommonSingle,
##     SNPs.146CommonSingle, SNPs.147CommonSingle, SNPs.Illumina
```

``` r
mset_dasnoob <- wateRmelon::dasen(mset_noob)
betas_dasnoob <- getBeta(mset_dasnoob)
```

## 3.0 Probe Filtering

### 3.1 SNP (rs) probes


``` r
rs <- getSnpBeta(rgset_raw) #get the beta values for the 59 rs SNP probes
dim (rs) #59
```

```
## [1]  59 510
```

``` r
table(rownames(rs) %in% rownames(betas_dasnoob)) #FALSE, indicates they have been removed
```

```
## 
## FALSE 
##    59
```

### 3.2 Filter to b5 Manifest


``` r
#check the chip numbers to see if we should use the b4 or b5 manifest
table(metadata_SPAH_ss$Sentrix_ID > 204220330001) #all are true so use b5 manifest
```

```
## 
## TRUE 
##  510
```

``` r
#remove the bad b4 probes labeled in the MFG_Change_Flagged column
probeInfo_b5 <- probeInfo %>% filter(!MFG_Change_Flagged ==TRUE)
dim(probeInfo_b5) #864869  
```

```
## [1] 864869     44
```

``` r
table(rownames(betas_dasnoob) %in% probeInfo_b5$probeID)
```

```
## 
##  FALSE   TRUE 
##   1369 864869
```

``` r
b5_probes <- rownames(betas_dasnoob) %in% probeInfo_b5$probeID
betas_dasnoob_filt <- betas_dasnoob[b5_probes,] 
dim(betas_dasnoob_filt) #864869
```

```
## [1] 864869    510
```

### 3.3 Polymorphic & cross-hybridizing probes

The EPIC zhou annotation object for this script was downloaded April 30, 2024. We selected the hg38 manifest. Downloaded from [here](https://zwdzwd.github.io/InfiniumAnnotation) 


``` r
zhou_mask <- read.delim(here::here("B. Data/B. Reference Data/EPIC.hg38.mask.tsv.gz")) #newer masking info
dim(zhou_mask)#865918
```

```
## [1] 865918     12
```

``` r
zhou_anno <- read.delim(here::here("B. Data/B. Reference Data/EPIC.hg38.manifest.tsv.gz")) #newer annotation info
dim(zhou_anno) #866553 
```

```
## [1] 866553     28
```

``` r
# Determine how many probes will be removed if we remove everything from the MASK_general column
table(zhou_mask$MASK_general) #105454 probes will be removed with MASK_general
```

```
## 
##  FALSE   TRUE 
## 760464 105454
```

``` r
zhou_maskgeneral <- zhou_mask %>% filter(MASK_general) #save the MASK_general probes in a new object

dim(betas_dasnoob_filt) #
```

```
## [1] 864869    510
```

``` r
table(rownames(betas_dasnoob_filt) %in% zhou_maskgeneral$probeID) 
```

```
## 
##  FALSE   TRUE 
## 759566 105303
```

``` r
betas_dasnoob_filt <- betas_dasnoob_filt[!(rownames(betas_dasnoob_filt) %in% zhou_maskgeneral$probeID),]
dim(betas_dasnoob_filt) #
```

```
## [1] 759566    510
```

### 3.4 DetP and Beadcount

We remove poor quality probes with a detection P-value > 0.01 in more than 5% of samples. 


``` r
# get the detection P (detP) values from the raw rgset.  
detP <- minfi::detectionP(rgset_raw) 
detP <- detP[rownames(detP) %in% probeInfo_b5$probeID,] #filter to b5 probes
dim(detP) #864869
```

```
## [1] 864869    510
```

``` r
# rather than actually sex-stratifying for the detP probes I am going to change 
# all the female Y chromosome values to 0 so that they aren't considered for detP removal 
metadata_SPAH_ss_f <- metadata_SPAH_ss %>% filter(Sex_Predic == "XX")

f_detP <- detP #make an object where we will modifty the detP value in females 
f_detP[rownames(f_detP) %in% chrYprobes$probeID, colnames(f_detP) %in% metadata_SPAH_ss_f$rgset_colnames] <- 0
femaleYdetP <- f_detP[rownames(f_detP) %in% chrYprobes$probeID, colnames(f_detP) %in% metadata_SPAH_ss_f$rgset_colnames]
all(femaleYdetP == 0) #check this step worked
```

```
## [1] TRUE
```

``` r
failed_detP <- f_detP > 0.01 #get a matrix of how many probes failed detP

#get the beadcount values from the raw rgset
beadcount <- wateRmelon::beadcount(rgset_raw)
head(colnames(beadcount))#check the colnames of the beadcount object. An X is added
```

```
## [1] "X207179230009_R01C01" "X207179230009_R02C01" "X207179230009_R03C01"
## [4] "X207179230009_R04C01" "X207179230009_R05C01" "X207179230009_R06C01"
```

``` r
colnames(beadcount) <- gsub("X","",colnames(beadcount))
head(colnames(beadcount)) #fixed :)
```

```
## [1] "207179230009_R01C01" "207179230009_R02C01" "207179230009_R03C01"
## [4] "207179230009_R04C01" "207179230009_R05C01" "207179230009_R06C01"
```

``` r
beadcount <- beadcount [rownames(beadcount) %in% probeInfo_b5$probeID,] #filter to b5 probes
dim(beadcount) #864869
```

```
## [1] 864869    510
```

``` r
identical(rownames(betas_dasnoob), rownames(beadcount)) #check
```

```
## [1] FALSE
```

``` r
betas_raw <- preprocessRaw(rgset_raw)

# create an logical matrix with which probes failed a check (true = fail)
# sum to see if a probe failed any checks (false + false = 0)
failed_probes <- failed_detP + is.na(beadcount) + is.na(betas_raw)
```

```
## Warning in failed_detP + is.na(beadcount) + is.na(betas_raw): longer object
## length is not a multiple of shorter object length
```

``` r
failed_probes <- failed_probes > 0 # >0 means failed any of the checks
table(rowSums(failed_probes) >= (nrow(metadata_SPAH)*0.05)) #12737 probes to remove
```

```
## 
##  FALSE   TRUE 
## 852132  12737
```

``` r
# index them for removal
badprobes <- data.frame(probebad = rowSums(failed_probes) >= nrow(metadata_SPAH)*0.05) %>% filter(probebad == TRUE)
dim(badprobes) #12737
```

```
## [1] 12737     1
```

``` r
## Remove the probes from the betas object
dim(betas_dasnoob_filt)
```

```
## [1] 759566    510
```

``` r
betas_dasnoob_filt <- betas_dasnoob_filt[!(row.names(betas_dasnoob_filt) %in% rownames(badprobes)),]
dim(betas_dasnoob_filt) #
```

```
## [1] 748756    510
```

### 3.5 X and Y Specific Poor Probes


``` r
#load additional X and Y annotation info
chrXY_ext <- read.csv(here::here("B. Data/B. Reference Data/Inkster_extended_XY_probe_annotation.csv"),header = T, fileEncoding="UTF-8-BOM")

bad_xy_probes <- chrXY_ext %>% filter(MASK_X_transposed_region == TRUE |
                                        MASK_AnyRepetitive == TRUE)
dim(bad_xy_probes) #726
```

```
## [1] 726  20
```

``` r
table(rownames(betas_dasnoob_filt) %in% bad_xy_probes$probeID) #475
```

```
## 
##  FALSE   TRUE 
## 748281    475
```

``` r
#subset the methylation object
betas_dasnoob_filt <- betas_dasnoob_filt[!rownames(betas_dasnoob_filt) %in% bad_xy_probes$probeID,]
dim(betas_dasnoob_filt)  #748281
```

```
## [1] 748281    510
```

## 4.0 Remove Samples that failed QC checks


``` r
remove_samps <- read_xlsx(here::here("C. Processing/02_Outputs/01_QC_Preprocessing/SPAH_remove_samples_list.xlsx")) 
table(metadata_SPAH_ss$Sample_ID %in% remove_samps$Sample_ID)
```

```
## 
## FALSE  TRUE 
##   496    14
```

``` r
metadata_SPAH_ss_filt <- metadata_SPAH_ss %>% filter(!Sample_ID %in% remove_samps$Sample_ID)
dim(metadata_SPAH_ss_filt) #496
```

```
## [1] 496 144
```

``` r
#remove remaining replicate samples
metadata_SPAH_ss_filt <- metadata_SPAH_ss_filt %>% filter(!grepl("-redo",NUSeq_Core_Facility_ID))
dim(metadata_SPAH_ss_filt) #493
```

```
## [1] 493 144
```

``` r
#subset the methylation object
betas_dasnoob_filt <- betas_dasnoob_filt[,colnames(betas_dasnoob_filt) %in% metadata_SPAH_ss_filt$rgset_colnames]
dim(betas_dasnoob_filt)
```

```
## [1] 748281    493
```

``` r
saveRDS(betas_dasnoob_filt, here::here("D. Pathology Project/02_Outputs/AA. Cleaned/SPAH_dasen_noob.rds"))
saveRDS(metadata_SPAH_ss_filt, here::here("D. Pathology Project/02_Outputs/AA. Cleaned/SPAH_metadata_ss_filtered.rds"))
```

## sessionInfo()


``` r
sessionInfo()
```

```
## R version 4.3.1 (2023-06-16 ucrt)
## Platform: x86_64-w64-mingw32/x64 (64-bit)
## Running under: Windows Server 2016 x64 (build 14393)
## 
## Matrix products: default
## 
## 
## locale:
## [1] LC_COLLATE=English_Canada.1252  LC_CTYPE=English_Canada.1252   
## [3] LC_MONETARY=English_Canada.1252 LC_NUMERIC=C                   
## [5] LC_TIME=English_Canada.1252    
## 
## time zone: America/Vancouver
## tzcode source: internal
## 
## attached base packages:
## [1] parallel  stats4    stats     graphics  grDevices utils     datasets 
## [8] methods   base     
## 
## other attached packages:
##  [1] IlluminaHumanMethylationEPICanno.ilm10b4.hg19_0.6.0
##  [2] IlluminaHumanMethylationEPICmanifest_0.3.0         
##  [3] readxl_1.4.5                                       
##  [4] wateRmelon_2.6.0                                   
##  [5] illuminaio_0.42.0                                  
##  [6] IlluminaHumanMethylation450kanno.ilmn12.hg19_0.6.1 
##  [7] ROC_1.76.0                                         
##  [8] lumi_2.52.0                                        
##  [9] methylumi_2.46.0                                   
## [10] minfi_1.46.0                                       
## [11] bumphunter_1.42.0                                  
## [12] locfit_1.5-9.8                                     
## [13] iterators_1.0.14                                   
## [14] foreach_1.5.2                                      
## [15] Biostrings_2.68.1                                  
## [16] XVector_0.40.0                                     
## [17] SummarizedExperiment_1.30.2                        
## [18] MatrixGenerics_1.12.3                              
## [19] FDb.InfiniumMethylation.hg19_2.2.0                 
## [20] org.Hs.eg.db_3.17.0                                
## [21] TxDb.Hsapiens.UCSC.hg19.knownGene_3.2.2            
## [22] GenomicFeatures_1.52.1                             
## [23] AnnotationDbi_1.62.2                               
## [24] GenomicRanges_1.52.1                               
## [25] GenomeInfoDb_1.36.1                                
## [26] IRanges_2.34.1                                     
## [27] S4Vectors_0.38.1                                   
## [28] reshape2_1.4.4                                     
## [29] scales_1.4.0                                       
## [30] matrixStats_1.0.0                                  
## [31] limma_3.56.2                                       
## [32] Biobase_2.60.0                                     
## [33] BiocGenerics_0.46.0                                
## [34] lubridate_1.9.4                                    
## [35] forcats_1.0.0                                      
## [36] stringr_1.5.0                                      
## [37] dplyr_1.1.4                                        
## [38] purrr_1.0.2                                        
## [39] readr_2.1.4                                        
## [40] tidyr_1.3.0                                        
## [41] tibble_3.2.1                                       
## [42] ggplot2_3.5.2                                      
## [43] tidyverse_2.0.0                                    
## [44] here_1.0.1                                         
## 
## loaded via a namespace (and not attached):
##   [1] splines_4.3.1                                      
##   [2] BiocIO_1.10.0                                      
##   [3] bitops_1.0-7                                       
##   [4] filelock_1.0.2                                     
##   [5] cellranger_1.1.0                                   
##   [6] preprocessCore_1.62.1                              
##   [7] XML_3.99-0.14                                      
##   [8] lifecycle_1.0.4                                    
##   [9] rprojroot_2.0.3                                    
##  [10] lattice_0.21-8                                     
##  [11] MASS_7.3-60                                        
##  [12] base64_2.0.1                                       
##  [13] scrime_1.3.5                                       
##  [14] magrittr_2.0.3                                     
##  [15] sass_0.4.10                                        
##  [16] rmarkdown_2.29                                     
##  [17] jquerylib_0.1.4                                    
##  [18] yaml_2.3.10                                        
##  [19] doRNG_1.8.6                                        
##  [20] askpass_1.1                                        
##  [21] DBI_1.1.3                                          
##  [22] RColorBrewer_1.1-3                                 
##  [23] abind_1.4-5                                        
##  [24] zlibbioc_1.46.0                                    
##  [25] quadprog_1.5-8                                     
##  [26] RCurl_1.98-1.12                                    
##  [27] rappdirs_0.3.3                                     
##  [28] GenomeInfoDbData_1.2.10                            
##  [29] genefilter_1.82.1                                  
##  [30] IlluminaHumanMethylationEPICanno.ilm10b5.hg38_0.0.1
##  [31] annotate_1.78.0                                    
##  [32] DelayedMatrixStats_1.22.5                          
##  [33] codetools_0.2-19                                   
##  [34] DelayedArray_0.26.7                                
##  [35] xml2_1.3.8                                         
##  [36] tidyselect_1.2.1                                   
##  [37] farver_2.1.2                                       
##  [38] beanplot_1.3.1                                     
##  [39] BiocFileCache_2.8.0                                
##  [40] GenomicAlignments_1.36.0                           
##  [41] jsonlite_2.0.0                                     
##  [42] multtest_2.56.0                                    
##  [43] survival_3.5-7                                     
##  [44] tools_4.3.1                                        
##  [45] progress_1.2.2                                     
##  [46] Rcpp_1.0.14                                        
##  [47] glue_1.8.0                                         
##  [48] xfun_0.52                                          
##  [49] mgcv_1.9-0                                         
##  [50] HDF5Array_1.28.1                                   
##  [51] withr_3.0.2                                        
##  [52] BiocManager_1.30.22                                
##  [53] fastmap_1.2.0                                      
##  [54] rhdf5filters_1.12.1                                
##  [55] fansi_1.0.6                                        
##  [56] openssl_2.1.0                                      
##  [57] digest_0.6.37                                      
##  [58] timechange_0.3.0                                   
##  [59] R6_2.5.1                                           
##  [60] colorspace_2.1-1                                   
##  [61] dichromat_2.0-0.1                                  
##  [62] biomaRt_2.56.1                                     
##  [63] RSQLite_2.3.1                                      
##  [64] utf8_1.2.5                                         
##  [65] generics_0.1.3                                     
##  [66] data.table_1.17.6                                  
##  [67] rtracklayer_1.60.1                                 
##  [68] prettyunits_1.1.1                                  
##  [69] httr_1.4.7                                         
##  [70] S4Arrays_1.0.5                                     
##  [71] pkgconfig_2.0.3                                    
##  [72] gtable_0.3.6                                       
##  [73] blob_1.2.4                                         
##  [74] siggenes_1.74.0                                    
##  [75] htmltools_0.5.8.1                                  
##  [76] png_0.1-8                                          
##  [77] knitr_1.50                                         
##  [78] rstudioapi_0.17.1                                  
##  [79] tzdb_0.4.0                                         
##  [80] rjson_0.2.21                                       
##  [81] nlme_3.1-163                                       
##  [82] curl_5.0.2                                         
##  [83] cachem_1.1.0                                       
##  [84] rhdf5_2.44.0                                       
##  [85] KernSmooth_2.23-22                                 
##  [86] restfulr_0.0.15                                    
##  [87] GEOquery_2.68.0                                    
##  [88] pillar_1.9.0                                       
##  [89] grid_4.3.1                                         
##  [90] reshape_0.8.9                                      
##  [91] vctrs_0.6.5                                        
##  [92] dbplyr_2.3.3                                       
##  [93] xtable_1.8-4                                       
##  [94] evaluate_1.0.4                                     
##  [95] cli_3.6.3                                          
##  [96] compiler_4.3.1                                     
##  [97] Rsamtools_2.16.0                                   
##  [98] rlang_1.1.4                                        
##  [99] crayon_1.5.2                                       
## [100] rngtools_1.5.2                                     
## [101] nor1mix_1.3-0                                      
## [102] mclust_6.0.0                                       
## [103] affy_1.78.2                                        
## [104] plyr_1.8.9                                         
## [105] stringi_1.8.7                                      
## [106] BiocParallel_1.34.2                                
## [107] nleqslv_3.3.4                                      
## [108] Matrix_1.6-1                                       
## [109] hms_1.1.3                                          
## [110] sparseMatrixStats_1.12.2                           
## [111] bit64_4.0.5                                        
## [112] Rhdf5lib_1.22.0                                    
## [113] KEGGREST_1.40.0                                    
## [114] memoise_2.0.1                                      
## [115] affyio_1.70.0                                      
## [116] bslib_0.9.0                                        
## [117] bit_4.0.5
```
