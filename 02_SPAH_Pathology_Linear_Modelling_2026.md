---
title: "SPAH Pathology Project Linear Modelling"
author: "Hannah Illing"
date: "2026 February 11"
output: 
  html_document:
    keep_md: yes
    code_folding: show
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
library(readxl)
library(xlsx)
library(limma)
library(lumi)
library(biobroom)
library(cowplot)
library(irlba)
library(planet)
library(ggrepel)
library(plomics)
library(missMethyl)
library(bacon)
library(ggpubr)

#load Illumina b5 manifest
probeInfo <- readRDS(here::here("E. Sex and Environment Project [abandon]/02_Outputs/02_Variability/hg38_b5_EPIC_probeinfo.rds"))
chrXprobes <- probeInfo %>% filter(chr == "chrX")
chrYprobes <- probeInfo %>% filter(chr == "chrY")
autoProbes <- probeInfo %>% filter(!(chr %in% c("chrX", "chrY")))
remove(probeInfo)

set.seed(4815)

theme_set(theme_minimal())
theme_update(
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(), 
  axis.line = element_line(colour = "black"),
  axis.title.x = element_text(size = 16),
  axis.title.y = element_text(size = 16),
  axis.text = element_text(size = 14),
  strip.text= element_text(size = 16, color = "black"),
  legend.position="none"
)
```

### 0.2 Load Data

1.    `metadata` - demographic and epiphenotyping variables for all **493** samples we want to include in our analysis. It includes all of the variables calculated in my pathology epivariables script.
2.    `betas_dasnoob_filt` - dasen+noob normalized beta values from the QC and pre-processing script. Includes methylation data for **506** placentas. Needs to be filtered down to match the `metadata` object. Column names in this object match the `rgset_colnames` variable in the `metadata` object


``` r
## 1) Load the metadata
metadata <- readRDS(here::here("D. Pathology Project/02_Outputs/AA. Cleaned/SPAH_metadata_ss_filtered.rds"))
dim(metadata) #493 x 162
```

```
## [1] 493 144
```

``` r
## 2) load the processed betas
betas_dasnoob_filt <- readRDS(here::here("D. Pathology Project/02_Outputs/AA. Cleaned/SPAH_dasen_noob.rds"))
dim(betas_dasnoob_filt) #748007 x 506
```

```
## [1] 748281    493
```

``` r
#filter to the correct samples
betas_dasnoob_filt <- betas_dasnoob_filt[,metadata$rgset_colnames]
dim(betas_dasnoob_filt) #748007 x 493
```

```
## [1] 748281    493
```

``` r
all.equal(colnames(betas_dasnoob_filt), metadata$rgset_colnames) #check that the meta data and methylation data match
```

```
## [1] TRUE
```

## 1.0 Linear Modelling Set-Up

### 1.1 Function Set-up

The following functions can be used to run linear models looking for associations with MVM presence/absence, perform multiple test correction by false discovery rate (FDR) and output a table of how many CpGs are significantly associated with your co-variate of interest at different effect sizes (delta B).


``` r
linear_modelling <- function(betas,mvals,model){
  #run model on m values & calculate empirical Bayesian statistics
  fit<- lmFit(mvals, model) %>% eBayes()
  td <- biobroom::tidy.MArrayLM(fit) # throws a tibble warning, still works

  #calculate delta beta
  td <- td %>% mutate(delB =
       (biobroom::tidy.MArrayLM(lmFit(betas, model) %>% 
        eBayes()))$estimate)
  
  #multiple test correct
  td <- td %>% 
    dplyr::rename(probeID = gene) %>% #biobroom assumes tested genes, not CpGs, rename
    group_by(term) %>%
    mutate(fdr = p.adjust(p.value, method="fdr")) %>%
    ungroup() %>%
    as.data.frame() 
  return(td)
}

lm_significant_hits <- function(td){
  #significant hits 
  sig_hits <- matrix(c(
    nrow(td %>% filter(fdr<0.05)),
    nrow(td %>% filter(fdr<0.05 & abs(delB) > 0.05)),
    nrow(td %>% filter(fdr<0.05 & abs(delB) > 0.1)),
    nrow(td %>% filter(fdr<0.01)),
    nrow(td %>% filter(fdr<0.01 & abs(delB) > 0.05)),
    nrow(td %>% filter(fdr<0.01 & abs(delB) > 0.1)))
    ,ncol=2)
  
  colnames(sig_hits) <- c("FDR < 0.05", "FDR < 0.01")
  rownames(sig_hits) <- c("delB > 0", "delB > 0.05", "delB > 0.1")
  sig_hits <-as.table(sig_hits)
  return(sig_hits)
}
```

### 1.2 Object Set-up


``` r
## 1) Make the metadata object for modelling
# 1a) Check how many samples have NAs or EOPE that would be removed for modelling
table(metadata$PE_GA_state) #7 EOPE samples
```

```
## 
## EOPE LOPE nPTB  nTB 
##    7   42   40  402
```

``` r
colSums(is.na(metadata %>%
                dplyr::select(pcr_pre_eclampsia,
                          HarmonizedPath, #9 samples
                          ga_continuous)))
```

```
## pcr_pre_eclampsia    HarmonizedPath     ga_continuous 
##                 2                 9                 1
```

``` r
# 1b) Select the metadata variables and remove missing data
metadata_naomit <- metadata %>%
  dplyr::select(Sex_Predic,
                PE_GA_state,
                HarmonizedPath,
                mvm_di,
                mvm_score_3cat,
                FETAL_VASC_PATH,
                Fetal_vasc_path_3cat,
                CHRONIC_INFLAMMATION,
                chronic_inflammation_3cat,
                ACUTE_INFLAMMATION,
                acute_inflammation_3cat,
                ga_continuous,
                Trophoblasts,
                Syncytiotrophoblast,
                Stromal,
                Endothelial,
                Hofbauer,
                nRBC,
                Prob_African,
                Prob_Asian,
                Prob_European,
                rgset_colnames) %>%
 #filter(!PE_GA_state == "EOPE") %>%
  na.omit() 

dim(metadata_naomit) #483 samples 
```

```
## [1] 483  22
```

``` r
## 3) select probes 
## 4) convert to M-values

# autosomes
betas_auto <- na.omit(betas_dasnoob_filt [rownames(betas_dasnoob_filt) %in% autoProbes$probeID, metadata_naomit$rgset_colnames])
mvals_auto <- lumi::beta2m(betas_auto)
dim(mvals_auto) #732102 x 483
```

```
## [1] 732102    483
```

## 2.0 MVM

### 2.1 MVM vs. No MVM


``` r
table(metadata_naomit$mvm_score_3cat)
```

```
## 
##   0   1   2 
## 340  84  59
```

``` r
## MVM modelling
mod_mvm_auto_basic <- model.matrix(~ mvm_di, metadata_naomit,row.names=T)
mod_mvm_auto <- model.matrix(~ mvm_di +
                               Sex_Predic + 
                               Prob_European + 
                               Prob_African + 
                               ga_continuous, 
                             metadata_naomit,row.names=T)

mod_mvm_cells <- model.matrix(~ mvm_di + 
                                Sex_Predic + 
                                Prob_European + 
                                Prob_African + 
                                Syncytiotrophoblast + 
                                Stromal + 
                                Endothelial + 
                                nRBC + 
                                ga_continuous, 
                              metadata_naomit,row.names=T)

head(mod_mvm_cells)
```

```
##   (Intercept) mvm_di1 Sex_PredicXY Prob_European Prob_African
## 1           1       1            1     0.9944981 0.0019729407
## 2           1       1            0     0.8665399 0.0087901395
## 3           1       1            1     0.9949846 0.0010767692
## 4           1       0            0     0.9949141 0.0008266452
## 5           1       1            1     0.7087023 0.0108689808
## 6           1       0            1     0.9949257 0.0012142192
##   Syncytiotrophoblast     Stromal Endothelial       nRBC ga_continuous
## 1           0.7721641 0.017105452  0.07351688 0.01571232         30.71
## 2           0.8344046 0.019953485  0.07704861 0.01673799         39.29
## 3           0.9196496 0.012145654  0.04171552 0.01755493         39.43
## 4           0.8725699 0.020448455  0.06373648 0.02670345         41.29
## 5           0.9094798 0.006426069  0.06016573 0.01993460         39.00
## 6           0.8239425 0.042196156  0.05986966 0.01807088         38.00
```

``` r
#run models
lm_mvm_auto_basic <-linear_modelling (betas_auto,mvals_auto,mod_mvm_auto_basic) %>% filter (term == "mvm_di1")
```

```
## Warning: `tbl_df()` was deprecated in dplyr 1.0.0.
## i Please use `tibble::as_tibble()` instead.
## i The deprecated feature was likely used in the biobroom package.
##   Please report the issue at <https://github.com/StoreyLab/biobroom/issues>.
## This warning is displayed once every 8 hours.
## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
## generated.
```

``` r
lm_mvm_auto <-linear_modelling (betas_auto,mvals_auto,mod_mvm_auto) %>% filter (term == "mvm_di1")
lm_mvm_cells <-linear_modelling (betas_auto,mvals_auto,mod_mvm_cells) %>% filter (term == "mvm_di1")

#are there any significant hits? 
lm_significant_hits(lm_mvm_auto_basic) 
```

```
##             FDR < 0.05 FDR < 0.01
## delB > 0         21533       5700
## delB > 0.05         27         18
## delB > 0.1           0          0
```

``` r
lm_significant_hits(lm_mvm_auto) 
```

```
##             FDR < 0.05 FDR < 0.01
## delB > 0            20          0
## delB > 0.05          0          0
## delB > 0.1           0          0
```

``` r
lm_significant_hits(lm_mvm_cells)
```

```
##             FDR < 0.05 FDR < 0.01
## delB > 0            30          0
## delB > 0.05          0          0
## delB > 0.1           0          0
```

``` r
#calculate inflation using the bacon package
inflation(bacon(teststatistics = lm_mvm_auto_basic$statistic)) #1.189345 
```

```
##  sigma.0 
## 1.189345
```

``` r
inflation(bacon(teststatistics = lm_mvm_auto$statistic)) #0.935061
```

```
##   sigma.0 
## 0.9472515
```

``` r
inflation(bacon(teststatistics = lm_mvm_cells$statistic)) #1.089234 
```

```
##  sigma.0 
## 1.080608
```

### 2.2  High-grade MVM vs. No MVM

Here I will run linear regressions to see which CpGs are differentially methylated in association with MVM in a comparison between samples with high grade MVM pathology (n=59) and samples with no MVM (n=340).


``` r
table(metadata_naomit$mvm_score_3cat) #59 with high grade MVM, 340 with no MVM
```

```
## 
##   0   1   2 
## 340  84  59
```

``` r
mvm_grade_naomit <- metadata_naomit %>% filter(!mvm_score_3cat == 1)
dim(mvm_grade_naomit) #399 (340+59)
```

```
## [1] 399  22
```

``` r
betas_auto_mvm <- betas_auto[,mvm_grade_naomit$rgset_colnames]
dim(betas_auto_mvm) #x 399
```

```
## [1] 732102    399
```

``` r
mvals_auto_mvm <- mvals_auto[,mvm_grade_naomit$rgset_colnames]
dim(mvals_auto_mvm) #x 399
```

```
## [1] 732102    399
```

``` r
#set-up models
mod_HGMVM_auto_basic <- model.matrix(~ mvm_di, mvm_grade_naomit,row.names=T)
mod_HGMVM_auto <- model.matrix(~ mvm_di + 
                                 Sex_Predic + 
                                 Prob_European + 
                                 Prob_African + 
                                 ga_continuous, 
                               mvm_grade_naomit,row.names=T)
mod_HGMVM_auto_cells <- model.matrix(~ mvm_di + 
                                       Sex_Predic + 
                                       Prob_European + 
                                       Prob_African + 
                                       ga_continuous + 
                                       Syncytiotrophoblast + 
                                       Endothelial + 
                                       Stromal + 
                                       nRBC, 
                                     mvm_grade_naomit,row.names=T)

head(mod_HGMVM_auto_cells)
```

```
##   (Intercept) mvm_di1 Sex_PredicXY Prob_European Prob_African ga_continuous
## 1           1       1            1     0.9944981 0.0019729407         30.71
## 2           1       1            0     0.8665399 0.0087901395         39.29
## 3           1       1            1     0.9949846 0.0010767692         39.43
## 4           1       0            0     0.9949141 0.0008266452         41.29
## 5           1       0            1     0.9949257 0.0012142192         38.00
## 6           1       0            1     0.9937753 0.0010489338         37.71
##   Syncytiotrophoblast Endothelial    Stromal       nRBC
## 1           0.7721641  0.07351688 0.01710545 0.01571232
## 2           0.8344046  0.07704861 0.01995348 0.01673799
## 3           0.9196496  0.04171552 0.01214565 0.01755493
## 4           0.8725699  0.06373648 0.02044846 0.02670345
## 5           0.8239425  0.05986966 0.04219616 0.01807088
## 6           0.9091273  0.03837857 0.03036592 0.01329400
```

``` r
#run models
lm_HGMVM_auto_basic <-linear_modelling (betas_auto_mvm,mvals_auto_mvm,mod_HGMVM_auto_basic) %>% filter (term == "mvm_di1")
lm_HGMVM_auto <-linear_modelling (betas_auto_mvm,mvals_auto_mvm,mod_HGMVM_auto) %>% filter (term == "mvm_di1")
lm_HGMVM_auto_cells <-linear_modelling (betas_auto_mvm,mvals_auto_mvm,mod_HGMVM_auto_cells) %>% filter (term == "mvm_di1")

#are there any significant hits? 
lm_significant_hits(lm_HGMVM_auto_basic) # 1286
```

```
##             FDR < 0.05 FDR < 0.01
## delB > 0         82360      26877
## delB > 0.05       1286        910
## delB > 0.1           7          5
```

``` r
lm_significant_hits(lm_HGMVM_auto) # 281 + sex + ancestry + gestational age
```

```
##             FDR < 0.05 FDR < 0.01
## delB > 0          5833        453
## delB > 0.05        281         70
## delB > 0.1           1          1
```

``` r
lm_significant_hits(lm_HGMVM_auto_cells) # 17 + sex + ancestry + gestational age + cells
```

```
##             FDR < 0.05 FDR < 0.01
## delB > 0           284         14
## delB > 0.05         17          2
## delB > 0.1           1          0
```

``` r
#calculate inflation using the bacon package
inflation(bacon(teststatistics = lm_HGMVM_auto_basic$statistic)) #1.166218
```

```
##  sigma.0 
## 1.166218
```

``` r
inflation(bacon(teststatistics = lm_HGMVM_auto$statistic)) #1.071452
```

```
##  sigma.0 
## 1.071452
```

``` r
inflation(bacon(teststatistics = lm_HGMVM_auto_cells$statistic)) #1.052073 
```

```
##  sigma.0 
## 1.052073
```

### 2.3 Sensitivity Analysis


``` r
#repeat but remove the cases with eoPE
mvm_grade_noeope_naomit <- mvm_grade_naomit %>% filter(!PE_GA_state == "EOPE")

table(mvm_grade_noeope_naomit$mvm_score_3cat) #52 with high grade MVM, 340 with no MVM
```

```
## 
##   0   1   2 
## 340   0  52
```

``` r
betas_auto_mvm_noeope <- betas_auto[,mvm_grade_noeope_naomit$rgset_colnames]
dim(betas_auto_mvm_noeope) #x 392
```

```
## [1] 732102    392
```

``` r
mvals_auto_mvm_noeope <- mvals_auto[,mvm_grade_noeope_naomit$rgset_colnames]
dim(mvals_auto_mvm_noeope) #x 392
```

```
## [1] 732102    392
```

``` r
#set-up models
mod_HGMVM_noeope_auto_basic <- model.matrix(~ mvm_di, mvm_grade_noeope_naomit,row.names=T)
mod_HGMVM_noeope_auto <- model.matrix(~ mvm_di + 
                                 Sex_Predic + 
                                 Prob_European + 
                                 Prob_African + 
                                 ga_continuous, 
                               mvm_grade_noeope_naomit,row.names=T)
mod_HGMVM_noeope_auto_cells <- model.matrix(~ mvm_di + 
                                       Sex_Predic + 
                                       Prob_European + 
                                       Prob_African + 
                                       ga_continuous + 
                                       Syncytiotrophoblast + 
                                       Endothelial + 
                                       Stromal + 
                                       nRBC, 
                                     mvm_grade_noeope_naomit,row.names=T)

head(mod_HGMVM_noeope_auto_cells)
```

```
##   (Intercept) mvm_di1 Sex_PredicXY Prob_European Prob_African ga_continuous
## 1           1       1            1     0.9944981 0.0019729407         30.71
## 2           1       1            0     0.8665399 0.0087901395         39.29
## 3           1       1            1     0.9949846 0.0010767692         39.43
## 4           1       0            0     0.9949141 0.0008266452         41.29
## 5           1       0            1     0.9949257 0.0012142192         38.00
## 6           1       0            1     0.9937753 0.0010489338         37.71
##   Syncytiotrophoblast Endothelial    Stromal       nRBC
## 1           0.7721641  0.07351688 0.01710545 0.01571232
## 2           0.8344046  0.07704861 0.01995348 0.01673799
## 3           0.9196496  0.04171552 0.01214565 0.01755493
## 4           0.8725699  0.06373648 0.02044846 0.02670345
## 5           0.8239425  0.05986966 0.04219616 0.01807088
## 6           0.9091273  0.03837857 0.03036592 0.01329400
```

``` r
#run models
lm_HGMVM_noeope_auto_basic <-linear_modelling (betas_auto_mvm_noeope,
                                               mvals_auto_mvm_noeope,
                                               mod_HGMVM_noeope_auto_basic) %>% filter (term == "mvm_di1")
lm_HGMVM_noeope_auto <-linear_modelling (betas_auto_mvm_noeope,
                                         mvals_auto_mvm_noeope,
                                         mod_HGMVM_noeope_auto) %>% filter (term == "mvm_di1")
lm_HGMVM_noeope_auto_cells <-linear_modelling (betas_auto_mvm_noeope,
                                               mvals_auto_mvm_noeope,
                                               mod_HGMVM_noeope_auto_cells) %>% filter (term == "mvm_di1")

#are there any significant hits? 
lm_significant_hits(lm_HGMVM_noeope_auto_basic) # 213
```

```
##             FDR < 0.05 FDR < 0.01
## delB > 0          8706       1515
## delB > 0.05        213         86
## delB > 0.1           1          1
```

``` r
lm_significant_hits(lm_HGMVM_noeope_auto) # 100 + sex + ancestry + gestational age
```

```
##             FDR < 0.05 FDR < 0.01
## delB > 0          1436         98
## delB > 0.05        100         14
## delB > 0.1           1          1
```

``` r
lm_significant_hits(lm_HGMVM_noeope_auto_cells) # 10 + sex + ancestry + gestational age + cells
```

```
##             FDR < 0.05 FDR < 0.01
## delB > 0           198          3
## delB > 0.05         10          0
## delB > 0.1           1          0
```

``` r
#calculate inflation using the bacon package
inflation(bacon(teststatistics = lm_HGMVM_noeope_auto_basic$statistic)) #1.010836
```

```
##  sigma.0 
## 1.014356
```

``` r
inflation(bacon(teststatistics = lm_HGMVM_noeope_auto$statistic)) #1.018541
```

```
##  sigma.0 
## 1.017186
```

``` r
inflation(bacon(teststatistics = lm_HGMVM_noeope_auto_cells$statistic)) #1.042864
```

```
##  sigma.0 
## 1.042941
```


``` r
comp_mvm_models <- lm_HGMVM_noeope_auto %>%
  dplyr::select(probeID, delB, fdr) %>%
  dplyr::rename(delB_noeope = delB, fdr_noeope = fdr) %>%
  left_join((lm_HGMVM_auto %>%
  dplyr::select(probeID, delB, fdr))) 
```

```
## Joining with `by = join_by(probeID)`
```

``` r
sup_fig_sens_a <- comp_mvm_models %>%
  ggplot(aes(y=delB, x=delB_noeope))+
  labs(y = expression(paste(Delta*beta, " HG MVM - No MVM")),
       x = expression(paste(Delta*beta, " HG MVM - No MVM (no EOPE)")))+
  geom_point()+
  stat_cor()
  
sup_fig_sens_b <- comp_mvm_models %>%
  ggplot(aes(x=abs(delB-delB_noeope)))+
  geom_density()+
  labs(x = expression(paste("|",Delta*Delta*beta, "| With/Without EOPE samples")))

comp_mvm_models %>%
  filter(fdr<0.05 | fdr_noeope <0.05)%>%
  ggplot(aes(x=fdr, y=fdr_noeope))+
  geom_point()
```

![](02_SPAH_Pathology_Linear_Modelling_2026_files/figure-html/compare-mvm-models-1.png)<!-- -->

``` r
sup_fig_sens_c <- lm_HGMVM_noeope_auto %>%
  ggplot(aes(x=delB, y=-log10(fdr), color=case_when(
    fdr < 0.05 & abs(delB) > 0.05 ~ "DMC",
    .default = "Not Significant"))) +
  geom_point(alpha=0.6)+
  geom_hline(yintercept = -log10(0.05), color="grey22", linetype="dashed")+
  geom_vline(xintercept=c(-0.05,0.05), color="grey22", linetype="dashed")+
  scale_color_manual(values=c("DMC" ="#033F63","Not Significant" ="grey"))+
  labs(x=expression(paste(Delta*beta, " HG MVM - No MVM")),y="-log(FDR)",color="",
       title = "DNAme ~ HG MVM + GA + Sex + Ancestry",
       subtitle = "100 hits")+
  annotate(geom="text", x=-0.1, y=1.2, col="black", size = 4, label="FDR < 0.05")+
  theme(plot.title = element_text(hjust=0.5, size = 14, face="italic"),
        plot.subtitle = element_text(hjust=0.5, size = 13, color ="#033F63", face="bold"))+
  ylim(0,4.25)
```

### 2.4 Hits & GO Analysis


``` r
lm_hg_mvm <- lm_HGMVM_auto_basic%>%
  dplyr::select(probeID, delB, fdr) %>% dplyr::rename(Simple_delB = delB, Simple_fdr = fdr) %>%
  left_join((
    lm_HGMVM_auto %>% 
      dplyr::select(probeID, delB, fdr) %>% dplyr::rename(Main_delB = delB, Main_fdr = fdr)), by="probeID") %>%
  left_join((
    lm_HGMVM_auto_cells %>% 
      dplyr::select(probeID, delB, fdr) %>% dplyr::rename(Cells_delB = delB, Cells_fdr = fdr)), by="probeID") %>%
  left_join((
    autoProbes %>%
      dplyr::select(probeID, chr, Start_hg38, UCSC_RefGene_Name, Regulatory_Region)), by="probeID")

saveRDS(lm_hg_mvm, here::here("D. Pathology Project/02_Outputs/AA. Cleaned/LM objects/lm_hg_mvm.rds"))

lm_hg_mvm_hits <- lm_hg_mvm %>%
  filter(Main_fdr < 0.05 & abs(Main_delB)>0.05)
dim(lm_hg_mvm_hits) #281 probes 
```

```
## [1] 281  11
```

``` r
hg_mvm_hits_cells <- (lm_HGMVM_auto_cells %>% filter(fdr < 0.05 & abs(delB) > 0.05))$probeID
table(hg_mvm_hits_cells %in% lm_hg_mvm_hits$probeID) #all 17 cells hits overlap with the 281 main hits
```

```
## 
## TRUE 
##   17
```

``` r
# write.xlsx(lm_hg_mvm_hits, here::here("D. Pathology Project/02_Outputs/AA. Cleaned/Supplementary File 1.xlsx"),sheetName="HG MVM", append=T,row.names=F)

eopred_cpgs <- read.csv(here::here("B. Data/B. Reference Data/eopred_cpgs.csv"))

table(lm_hg_mvm_hits$probeID %in% eopred_cpgs$ï..probeID) #5/281
```

```
## 
## FALSE  TRUE 
##   276     5
```

``` r
table(hg_mvm_hits_cells %in% eopred_cpgs$ï..probeID) #1/17
```

```
## 
## FALSE  TRUE 
##    16     1
```

``` r
table(lm_hg_mvm_hits$Regulatory_Region) #mostly enhancers (53.7%)
```

```
## 
##       Dual   Enhancer  Gene body Intergenic   Promoter 
##         37        151         62         21         10
```

``` r
table(lm_hg_mvm_hits$Main_delB>0) #240/281 (85.4%)
```

```
## 
## FALSE  TRUE 
##   240    41
```

``` r
#Do the hits overlap the cell-type DMCs?
#Load Victor's Third Trimester Cell DMCs
cells_dmcs <- read.csv("//fs/teams/RobinsonLab/Victor/Projects/NIH - cells/outs/2_4_cell_dmcs_third.csv")
table(lm_hg_mvm_hits$probeID  %in% cells_dmcs$cpg) #155/281
```

```
## 
## FALSE  TRUE 
##   126   155
```

``` r
go_hg_mvm<- gometh(sig.cpg = lm_hg_mvm_hits$probeID ,
                    all.cpg = rownames(betas_dasnoob_filt),
                    array.type = "EPIC")
```

```
## All input CpGs are used for testing.
```

``` r
go_hg_mvm<- go_hg_mvm %>% mutate(GOterm = rownames(go_hg_mvm)) %>% as_tibble %>% arrange(FDR)
go_hg_mvm %>% filter(FDR < 0.05) %>% nrow() #0 significantly associated terms
```

```
## [1] 0
```

``` r
go_hg_mvm[1:10,]
```

```
## # A tibble: 10 x 7
##    ONTOLOGY TERM                                      N    DE  P.DE   FDR GOterm
##    <chr>    <chr>                                 <dbl> <dbl> <dbl> <dbl> <chr> 
##  1 BP       mitochondrial genome maintenance         30     0 1.00      1 GO:00~
##  2 BP       reproduction                           1456     9 0.853     1 GO:00~
##  3 MF       alpha-1,6-mannosyltransferase activi~     2     0 1         1 GO:00~
##  4 MF       trans-hexaprenyltranstransferase act~     2     0 1.00      1 GO:00~
##  5 BP       single strand break repair               11     0 1.00      1 GO:00~
##  6 MF       single-stranded DNA endodeoxyribonuc~    10     0 1         1 GO:00~
##  7 CC       phosphopyruvate hydratase complex         4     0 1.00      1 GO:00~
##  8 MF       lactase activity                          1     0 1.00      1 GO:00~
##  9 BP       alpha-glucoside transport                 2     0 1.00      1 GO:00~
## 10 BP       regulation of DNA recombination         133     1 0.614     1 GO:00~
```

``` r
#Do the hits overlap with Sam's 599 EOPE-associated CpGs?
eope_hits <- read_xlsx("//fs/teams/RobinsonLab/Hannah/05_Projects/EOPE Review Paper/Data/Wilson_2017_PE_hits.xlsx")
dim(eope_hits) #599
```

```
## [1] 599  11
```

``` r
table(lm_hg_mvm_hits$probeID  %in% eope_hits$Probe) #25/281 (8.9%)
```

```
## 
## FALSE  TRUE 
##   256    25
```

### 2.5 Plots


``` r
#fix genes with duplicates 
#fix genes with duplicates in gene name
lm_hg_mvm$UCSC_RefGene_Name <- str_extract(lm_hg_mvm$UCSC_RefGene_Name, "[^;]+")

fig_2_a <- lm_hg_mvm %>%
  ggplot(aes(x=Main_delB, y=-log10(Main_fdr), color=case_when(
    Main_fdr < 0.05 & abs(Main_delB) > 0.05 ~ "DMC",
    .default = "Not Significant"))) +
  geom_point(alpha=0.6)+
  geom_hline(yintercept = -log10(0.05), color="grey22", linetype="dashed")+
  geom_vline(xintercept=c(-0.05,0.05), color="grey22", linetype="dashed")+
  scale_color_manual(values=c("DMC" ="#033F63","Not Significant" ="grey"))+
  labs(x=expression(paste(Delta*beta, " HG MVM - No MVM")),y="-log(FDR)",color="",
       title = "DNAme ~ HG MVM + GA + Sex + Ancestry",
       subtitle = "281 hits")+
  annotate(geom="text", x=-0.11, y=1.2, col="black", size = 4, label="FDR < 0.05")+
  geom_text_repel(data=(lm_hg_mvm %>%
                          filter(!is.na(UCSC_RefGene_Name))%>%
                          filter(Main_delB < -0.08 & Main_fdr < 0.03 | Main_delB >0.07 & Main_fdr < 0.03)),
                  aes(x=Main_delB, y=-log10(Main_fdr), label=UCSC_RefGene_Name),
                      nudge_x = c(0,0.005,0,0,0,0,0.005,-0.01),
                      nudge_y = c(0.1,-0.5,0,0,0.3,0,-0.1,0),
  show.legend=F, max.overlaps = Inf, size=3, min.segment.length=0)+
  theme(plot.title = element_text(hjust=0.5, size = 14, face="italic"),
        plot.subtitle = element_text(hjust=0.5, size = 13, color ="#033F63", face="bold"))+
  ylim(0,5)

fig_2_b <- lm_hg_mvm %>%
  ggplot(aes(x=Cells_delB, y=-log10(Cells_fdr), color=case_when(
    Cells_fdr < 0.05 & abs(Cells_delB) > 0.05 ~ "DMC",
    .default = "Not Significant"))) +
  geom_point(alpha=0.6)+
  geom_hline(yintercept = -log10(0.05), color="grey22", linetype="dashed")+
  geom_vline(xintercept=c(-0.05,0.05), color="grey22", linetype="dashed")+
  scale_color_manual(values=c("DMC" ="#033F63","Not Significant" ="grey"))+
  labs(x=expression(paste(Delta*beta, " HG MVM - No MVM")),y="-log(FDR)",color="",
       title = "DNAme ~ HG MVM + GA + Sex + Ancestry + Cells",
       subtitle = "17 hits")+
  annotate(geom="text", x=0.1, y=1.2, col="black", size = 4, label="FDR < 0.05")+
  geom_text_repel(data=(lm_hg_mvm %>%
                          filter(!is.na(UCSC_RefGene_Name))%>%
                          filter(Cells_delB < -0.05 & Cells_fdr < 0.05)),
                  aes(x=Cells_delB, y=-log10(Cells_fdr), label=UCSC_RefGene_Name),
  show.legend=F, max.overlaps = Inf, size=3, min.segment.length=0, box.padding=0.5)+
  theme(plot.title = element_text(hjust=0.5, size = 14, face="italic"),
        plot.subtitle = element_text(hjust=0.5, size = 13, color ="#033F63", face="bold"))+
  ylim(0,5)
```

#### 2.5.1 Hits Plots


``` r
table(grepl("FLNB",autoProbes$UCSC_RefGene_Name))
```

```
## 
##  FALSE   TRUE 
## 846130    102
```

``` r
flnb_probes <- autoProbes %>% filter(chr == "chr3" & as.numeric(Start_hg38) > 58008300 & as.numeric(Start_hg38) < 58012000)
nrow(flnb_probes) #10
```

```
## [1] 10
```

``` r
betas_flnb <- betas_dasnoob_filt[rownames(betas_dasnoob_filt) %in% flnb_probes$probeID,, drop=F] %>%
    as.data.frame() %>%
    rownames_to_column(var="probeID")%>%
    pivot_longer(cols=!probeID,
    names_to = "rgset_colnames",
    values_to = "beta") %>%
    left_join(metadata, by="rgset_colnames") %>% #combine with metadata 
    left_join(flnb_probes, by = "probeID") %>%
  mutate(DMC = ifelse(probeID %in% eopred_cpgs$ï..probeID, "Yes","No"))

mvm_hits_flnb_plot <- betas_flnb %>%
  filter(!is.na(mvm_di))%>%
  ggplot(aes(x=reorder(probeID,Start_hg38), y=beta, color=mvm_score_3cat, fill=mvm_score_3cat))+
  geom_boxplot(alpha=0.4,outlier.shape=NA, position=position_dodge(0.6), width=0.4)+
  geom_point(alpha=0.5, position=position_jitterdodge(dodge.width=0.6,jitter.width=0.1), size=2.25,aes(shape=DMC))+
  ylim(0,1)+
  scale_color_manual(values = c("darkgrey","#87b8d6","#033F63"),
                     labels= c("0" = "No MVM","1" = "Low Grade","2" = "High Grade"))+
  scale_fill_manual(values = c("darkgrey","#87b8d6","#033F63"),
                     labels= c("0" = "No MVM","1" = "Low Grade","2" = "High Grade"))+
  scale_shape_manual(values=c(3,19)) +
  theme(legend.position = "right",
        plot.title = element_text(face="italic", hjust=0.5))+
  labs(x="", y= expression(paste("DNAme (",beta,")")), fill = "MVM Pathology",
       color = "MVM Pathology",
       title="FLNB", shape="eoPRED CpG")
mvm_hits_flnb_plot
```

![](02_SPAH_Pathology_Linear_Modelling_2026_files/figure-html/hgmvm-hits-plots-a-1.png)<!-- -->

``` r
ggsave("Figure S4.png", plot = mvm_hits_flnb_plot, path = here::here("D. Pathology Project/02_Outputs/AA. Cleaned/"), width = 13, height = 5, device = "png")
```


``` r
krt_probes <- autoProbes %>% filter(chr == "chr17" & as.numeric(Start_hg38) > 41522055 & as.numeric(Start_hg38) < 41522987)
nrow(krt_probes) #8
```

```
## [1] 8
```

``` r
betas_krt <- betas_dasnoob_filt[rownames(betas_dasnoob_filt) %in% krt_probes$probeID,, drop=F] %>%
    as.data.frame() %>%
    rownames_to_column(var="probeID")%>%
    pivot_longer(cols=!probeID,
    names_to = "rgset_colnames",
    values_to = "beta") %>%
    left_join(metadata, by="rgset_colnames") %>% #combine with metadata 
    left_join(krt_probes, by = "probeID") %>%
  mutate(DMC = ifelse(probeID %in% lm_hg_mvm_hits$probeID , "Yes","No"))

mvm_hits_krt15_plot <- betas_krt %>%
  filter(!is.na(mvm_di))%>%
  ggplot(aes(x=reorder(probeID,Start_hg38), y=beta, color=mvm_score_3cat, fill=mvm_score_3cat))+
  geom_boxplot(alpha=0.4,outlier.shape=NA, position=position_dodge(0.6), width=0.4)+
  geom_point(alpha=0.5, position=position_jitterdodge(dodge.width=0.6,jitter.width=0.1), size=2.25,aes(shape=DMC))+
  ylim(0,1)+
  scale_color_manual(values = c("darkgrey","#87b8d6","#033F63"),
                     labels= c("0" = "No MVM","1" = "Low Grade","2" = "High Grade"))+
  scale_fill_manual(values = c("darkgrey","#87b8d6","#033F63"),
                     labels= c("0" = "No MVM","1" = "Low Grade","2" = "High Grade"))+
  scale_shape_manual(values=c(3,19)) +
  theme(legend.position = "right",
        plot.title = element_text(face="italic", hjust=0.5))+
  labs(x="", y= expression(paste("DNAme (",beta,")")), fill = "MVM Pathology",
       color = "MVM Pathology",
       title="KRT15-KRT19 Intergenic Region", shape="DMP")

mvm_hits_krt15_plot
```

![](02_SPAH_Pathology_Linear_Modelling_2026_files/figure-html/hgmvm-hits-plots-b-1.png)<!-- -->

``` r
ggsave("krt15_krt19_region.png", plot = mvm_hits_krt15_plot, path = here::here("D. Pathology Project/02_Outputs/AA. Cleaned/"), width = 12, height = 5, device = "png")
```

#### 2.5.2 Cells Plots


``` r
## Load Victor's cell-type data
#metadata
pDat <- readRDS('//fs/teams/RobinsonLab/Victor/Projects/NIH - cells/data/main/interim/2_3_pDat_contam.rds')

#remove bad samples and filter to tissues of interest
pDat_filt <- pDat %>% 
     filter(maternal_contamination_norm_flip < 0.35,
            !Sample_Name %in% c('PM364_hofb_cs', 'PL293_v_R2', 'PM366_vc_R2', 'P131_hofb_cs')) %>%
  filter(Tissue %in% c("Endothelial","Hofbauer","Stromal","Syncytiotrophoblast","Trophoblasts","Villi")) 
pDat_filt[pDat_filt == "Trophoblasts"] <- "Cytotrophoblasts"
pDat_filt <- pDat_filt %>% mutate(Tissue = factor(Tissue, levels=  c("Endothelial","Hofbauer","Stromal","Syncytiotrophoblast","Cytotrophoblasts","Villi")))
table(pDat_filt$Tissue)
```

```
## 
##         Endothelial            Hofbauer             Stromal Syncytiotrophoblast 
##                  27                  21                  28                   5 
##    Cytotrophoblasts               Villi 
##                  24                  46
```

``` r
#only the raw rgset is available, need to normalize
rgset <- readRDS('//fs/teams/RobinsonLab/Victor/Projects/NIH - cells/data/main/interim/0_1_rgset_raw.rds')
mset_noob <- preprocessNoob(rgset)
```

```
## Loading required package: IlluminaHumanMethylationEPICmanifest
```

``` r
betas_noob <- getBeta(mset_noob)

saveRDS(betas_noob, here::here("D. Pathology Project/02_Outputs/AA. Cleaned/betas_noob_cell_sorted.rds"))

betas_gene <- t(betas_noob[rownames(betas_noob) %in% lm_hg_mvm_hits$probeID ,
                            colnames(betas_noob) %in% pDat_filt$Sentrix]) %>%
  as.data.frame() %>%
  rownames_to_column(var="Sentrix") 

betas_gene <- betas_gene %>%
  # reshape into longer format
  pivot_longer(cols = -Sentrix,
               names_to = 'probeID', 
               values_to = 'Beta')  %>%
  # add tissue and trimester info
  left_join(pDat_filt %>% select(Sentrix, Trimester, Tissue), by = 'Sentrix') %>%
  # calculate mean and sd for each cpg for each group
  group_by(Tissue, Trimester, probeID) %>%
  summarize(mean = mean(Beta),
            sd = sd(Beta)) %>%
  mutate(lower = mean-sd, upper = mean+sd) %>%
  filter(Trimester == "Third")
```

```
## `summarise()` has grouped output by 'Tissue', 'Trimester'. You can override
## using the `.groups` argument.
```

``` r
labs <- c("\u0394\u03B2 < 0", "\u0394\u03B2 > 0")
names(labs) <- c(FALSE,TRUE)

## see how the DMPs look in the individual cell types

sup_cells_mvm_a <- betas_gene %>%
  mutate(lower = ifelse(lower<0,0,lower),
         upper = ifelse(upper>1,1,upper))%>% #fix error
  filter(!Tissue == "Villi") %>%
  left_join(lm_hg_mvm_hits, by="probeID") %>%
  mutate(Vars = Main_delB >0) %>%
  filter(Vars == FALSE)%>%
  ggplot() +
  geom_linerange(alpha = 0.5, size = 1,
                 aes(x = reorder(probeID,Main_delB), ymin =lower, ymax = upper, color = Tissue)) +
  geom_point(alpha = 1, aes(x = probeID, y = mean, color = Tissue), size = 2) +
  theme(plot.title = element_text(hjust=0.5, size = 16, face="bold"),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=12),
        legend.position="none",
        strip.background = element_rect(color="black", fill="grey88")) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_x_discrete(expand = c(0.01,0.01))+
  scale_color_manual(values= c('#6A1B9A', '#1565C0','#388E3C',    '#E64A19',
                               '#FBC02D', '#C62828')) +
  labs(y = "DNA Methylation", x = "", color = '')+
  facet_wrap(vars(Vars), labeller = labeller(Vars = labs))
```

```
## Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
## i Please use `linewidth` instead.
## This warning is displayed once every 8 hours.
## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
## generated.
```

``` r
sup_cells_mvm_b <- betas_gene %>%
  mutate(lower = ifelse(lower<0,0,lower),
         upper = ifelse(upper>1,1,upper))%>% #fix error
  filter(!Tissue == "Villi") %>%
  left_join(lm_hg_mvm_hits, by="probeID") %>%
  mutate(Vars = Main_delB >0) %>%
  filter(Vars == TRUE)%>%
  ggplot() +
  geom_linerange(alpha = 0.5, size = 1,
                 aes(x = reorder(probeID,Main_delB), ymin =lower, ymax = upper, color = Tissue)) +
  geom_point(alpha = 1, aes(x = probeID, y = mean, color = Tissue), size = 2) +
  theme(plot.title = element_text(hjust=0.5, size = 16, face="bold"),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        legend.position="right",
        strip.background = element_rect(color="black", fill="grey88")) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_x_discrete(expand = c(0.01,0.01))+
  scale_color_manual(values= c('#6A1B9A', '#1565C0','#388E3C',    '#E64A19',
                               '#FBC02D', '#C62828')) +
  labs(y="",x = "", color = '')+
  facet_wrap(vars(Vars), labeller = labeller(Vars = labs))

sup_cells_mvm <- plot_grid(sup_cells_mvm_a, NULL, sup_cells_mvm_b,
                           nrow=1,
                           rel_widths=c(1,-0.03,0.7))
```

```
## Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
## conversion failure on 'Î”Î² < 0' in 'mbcsToSbcs': dot substituted for <ce>
```

```
## Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
## conversion failure on 'Î”Î² < 0' in 'mbcsToSbcs': dot substituted for <94>
```

```
## Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
## conversion failure on 'Î”Î² < 0' in 'mbcsToSbcs': dot substituted for <ce>
```

```
## Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
## conversion failure on 'Î”Î² < 0' in 'mbcsToSbcs': dot substituted for <94>
```

```
## Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
## conversion failure on 'Î”Î² < 0' in 'mbcsToSbcs': dot substituted for <ce>
```

```
## Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
## conversion failure on 'Î”Î² < 0' in 'mbcsToSbcs': dot substituted for <94>
```

```
## Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
## conversion failure on 'Î”Î² < 0' in 'mbcsToSbcs': dot substituted for <ce>
```

```
## Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
## conversion failure on 'Î”Î² < 0' in 'mbcsToSbcs': dot substituted for <94>
```

```
## Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
## conversion failure on 'Î”Î² > 0' in 'mbcsToSbcs': dot substituted for <ce>
```

```
## Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
## conversion failure on 'Î”Î² > 0' in 'mbcsToSbcs': dot substituted for <94>
```

```
## Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
## conversion failure on 'Î”Î² > 0' in 'mbcsToSbcs': dot substituted for <ce>
```

```
## Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
## conversion failure on 'Î”Î² > 0' in 'mbcsToSbcs': dot substituted for <94>
```

```
## Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
## conversion failure on 'Î”Î² > 0' in 'mbcsToSbcs': dot substituted for <ce>
```

```
## Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
## conversion failure on 'Î”Î² > 0' in 'mbcsToSbcs': dot substituted for <94>
```

```
## Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
## conversion failure on 'Î”Î² > 0' in 'mbcsToSbcs': dot substituted for <ce>
```

```
## Warning in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y, :
## conversion failure on 'Î”Î² > 0' in 'mbcsToSbcs': dot substituted for <94>
```

``` r
sup_cells_mvm
```

![](02_SPAH_Pathology_Linear_Modelling_2026_files/figure-html/cells-plots-mvm-1.png)<!-- -->

## 3.0 FVM

### 3.1 FVM vs. No FVM


``` r
#note that one sample is mistakenly labelled as FVM in the FETAL_VASC_PATH column
#Use the FETAL_VASC_PATH column to create an updated FETAL_VASC_PATH variable
table(metadata_naomit$FETAL_VASC_PATH, metadata_naomit$Fetal_vasc_path_3cat)
```

```
##    
##       0   1   2
##   0 320   0   0
##   1   1 125  37
```

``` r
metadata_naomit <- metadata_naomit %>% 
  mutate(FETAL_VASC_PATH = as.factor(ifelse(Fetal_vasc_path_3cat==0,0,1)))

#set-up models
mod_FVM_auto_basic <- model.matrix(~FETAL_VASC_PATH, metadata_naomit,row.names=T)

mod_FVM_auto <- model.matrix(~FETAL_VASC_PATH+
                              ga_continuous+
                              Sex_Predic+
                              Prob_African+
                              Prob_European,
                            metadata_naomit,
                            row.names=T)

mod_FVM_auto_cells <- model.matrix(~FETAL_VASC_PATH+
                              ga_continuous+
                              Sex_Predic+
                              Prob_African+
                              Prob_European +
                              Syncytiotrophoblast+
                              Stromal+
                              Endothelial+
                              nRBC,
                            metadata_naomit,
                            row.names=T)

head(mod_FVM_auto_cells)
```

```
##   (Intercept) FETAL_VASC_PATH1 ga_continuous Sex_PredicXY Prob_African
## 1           1                0         30.71            1 0.0019729407
## 2           1                0         39.29            0 0.0087901395
## 3           1                0         39.43            1 0.0010767692
## 4           1                1         41.29            0 0.0008266452
## 5           1                0         39.00            1 0.0108689808
## 6           1                0         38.00            1 0.0012142192
##   Prob_European Syncytiotrophoblast     Stromal Endothelial       nRBC
## 1     0.9944981           0.7721641 0.017105452  0.07351688 0.01571232
## 2     0.8665399           0.8344046 0.019953485  0.07704861 0.01673799
## 3     0.9949846           0.9196496 0.012145654  0.04171552 0.01755493
## 4     0.9949141           0.8725699 0.020448455  0.06373648 0.02670345
## 5     0.7087023           0.9094798 0.006426069  0.06016573 0.01993460
## 6     0.9949257           0.8239425 0.042196156  0.05986966 0.01807088
```

``` r
#run regression
lm_FVM_auto_basic <-linear_modelling (betas_auto,mvals_auto,mod_FVM_auto_basic) %>% filter (term == "FETAL_VASC_PATH1")
lm_FVM_auto <-linear_modelling (betas_auto,mvals_auto,mod_FVM_auto) %>% filter (term == "FETAL_VASC_PATH1")
lm_FVM_auto_cells <-linear_modelling (betas_auto,mvals_auto,mod_FVM_auto_cells) %>% filter (term == "FETAL_VASC_PATH1")

#are there any significant hits?
lm_significant_hits(lm_FVM_auto_basic)
```

```
##             FDR < 0.05 FDR < 0.01
## delB > 0          7350       1725
## delB > 0.05          1          0
## delB > 0.1           0          0
```

``` r
lm_significant_hits(lm_FVM_auto)
```

```
##             FDR < 0.05 FDR < 0.01
## delB > 0          3626        996
## delB > 0.05          1          0
## delB > 0.1           0          0
```

``` r
lm_significant_hits(lm_FVM_auto_cells)
```

```
##             FDR < 0.05 FDR < 0.01
## delB > 0            21          9
## delB > 0.05          0          0
## delB > 0.1           0          0
```

``` r
lm_FVM_auto %>% filter(fdr<0.05 &abs(delB)>0.05) %>%
  arrange(probeID) %>%
  left_join(autoProbes, by="probeID")
```

```
##      probeID             term   estimate statistic     p.value      lod
## 1 cg25916282 FETAL_VASC_PATH1 -0.3605682 -3.936841 9.47733e-05 1.176151
##          delB       fdr AddressA AddressB
## 1 -0.05598493 0.0300753 51683281  6799916
##                                            ProbeSeqA
## 1 TACCCAAATTAATATAATAAATCATCATAAATAATTCCTAAAAAAATACA
##                                            ProbeSeqB Type NextBase Color  chr
## 1 TACCCGAATTAATATAATAAATCATCGTAAATAATTCCTAAAAAAATACG    I        A   Red chr2
##         pos strand
## 1 177023760      -
##                                                                                                               Forward_Sequence
## 1 CAGTGTGACAATTGCCCGGGTTGGTGTGATAAATCATCGTAAGTAATTCCTGAAAGGGTG[CG]AGACTGTTGGGGGCCGGGCGAGGACTGTAAATCTTTCCGGTTTATTGCTCTATGAACATA
##                                            SourceSeq UCSC_RefGene_Name
## 1 CGCACCCTTTCAGGAATTACTTACGATGATTTATCACACCAACCCGGGCA                  
##   UCSC_RefGene_Accession UCSC_RefGene_Group Phantom4_Enhancers
## 1                                                             
##   Phantom5_Enhancers DMR X450k_Enhancer            HMM_Island
## 1                                       2:176731887-176732046
##   Regulatory_Feature_Name Regulatory_Feature_Group GencodeBasicV12_NAME
## 1                                                                 HOXD3
##   GencodeBasicV12_Accession GencodeBasicV12_Group GencodeCompV12_NAME
## 1         ENST00000468418.2                 5'UTR         HOXD3;HOXD3
##              GencodeCompV12_Accession GencodeCompV12_Group
## 1 ENST00000432796.2;ENST00000468418.2          5'UTR;5'UTR
##   DNase_Hypersensitivity_NAME DNase_Hypersensitivity_Evidence_Count
## 1    chr2:177023365-177023995                                     3
##   OpenChromatin_NAME OpenChromatin_Evidence_Count TFBS_NAME TFBS_Evidence_Count
## 1                                                                              
##   Methyl27_Loci Methyl450_Loci Random_Loci MFG_Change_Flagged CHR_hg38
## 1                         TRUE                          FALSE     chr2
##   Start_hg38  End_hg38 Strand_hg38 Regulatory_Region Relation_to_Island
## 1  176159031 176159033           +          Promoter            N_Shore
##               Islands_Name
## 1 chr2:177024501-177025692
```

``` r
#calculate inflation using the bacon package
inflation(bacon(teststatistics = lm_FVM_auto_basic$statistic)) #1.074454 
```

```
##  sigma.0 
## 1.067284
```

``` r
inflation(bacon(teststatistics = lm_FVM_auto$statistic)) #1.040359   
```

```
##  sigma.0 
## 1.088041
```

``` r
inflation(bacon(teststatistics = lm_FVM_auto_cells$statistic)) #0.9973838 
```

```
##   sigma.0 
## 0.9991745
```

### 3.2 High-grade FVM vs. No FVM


``` r
#remove samples with low grade FVM
fvm_grade_naomit <- metadata_naomit %>% filter(!Fetal_vasc_path_3cat == 1) %>%
  mutate(Fetal_vasc_path_3cat = as.numeric(Fetal_vasc_path_3cat))
table(fvm_grade_naomit$Fetal_vasc_path_3cat) #37 HG FVM, 321 No FVM
```

```
## 
##   1   3 
## 321  37
```

``` r
betas_auto_fvm <- betas_auto[,fvm_grade_naomit$rgset_colnames]
mvals_auto_fvm <- mvals_auto[,fvm_grade_naomit$rgset_colnames]
dim(mvals_auto_fvm) #x358
```

```
## [1] 732102    358
```

``` r
#set-up models
mod_HGFVM_auto_basic <- model.matrix(~FETAL_VASC_PATH,fvm_grade_naomit,row.names=T)

mod_HGFVM_auto <- model.matrix(~FETAL_VASC_PATH+
                              ga_continuous+
                              Sex_Predic+
                              Prob_African+
                              Prob_European,
                            fvm_grade_naomit,
                            row.names=T)

mod_HGFVM_auto_cells <- model.matrix(~FETAL_VASC_PATH+
                              ga_continuous+
                              Sex_Predic+
                              Prob_African+
                              Prob_European +
                              Syncytiotrophoblast+
                              Stromal+
                              Endothelial+
                              nRBC,
                            fvm_grade_naomit,
                            row.names=T)

head(mod_HGFVM_auto_cells)
```

```
##   (Intercept) FETAL_VASC_PATH1 ga_continuous Sex_PredicXY Prob_African
## 1           1                0         30.71            1  0.001972941
## 2           1                0         39.29            0  0.008790140
## 3           1                0         39.43            1  0.001076769
## 4           1                0         39.00            1  0.010868981
## 5           1                0         38.00            1  0.001214219
## 6           1                0         37.71            1  0.001048934
##   Prob_European Syncytiotrophoblast     Stromal Endothelial       nRBC
## 1     0.9944981           0.7721641 0.017105452  0.07351688 0.01571232
## 2     0.8665399           0.8344046 0.019953485  0.07704861 0.01673799
## 3     0.9949846           0.9196496 0.012145654  0.04171552 0.01755493
## 4     0.7087023           0.9094798 0.006426069  0.06016573 0.01993460
## 5     0.9949257           0.8239425 0.042196156  0.05986966 0.01807088
## 6     0.9937753           0.9091273 0.030365924  0.03837857 0.01329400
```

``` r
#run regression
lm_HGFVM_auto_basic <-linear_modelling (betas_auto_fvm,mvals_auto_fvm,mod_HGFVM_auto_basic) %>% filter (term == "FETAL_VASC_PATH1")
lm_HGFVM_auto <-linear_modelling (betas_auto_fvm,mvals_auto_fvm,mod_HGFVM_auto) %>% filter (term == "FETAL_VASC_PATH1")
lm_HGFVM_auto_cells <-linear_modelling (betas_auto_fvm,mvals_auto_fvm,mod_HGFVM_auto_cells) %>% filter (term == "FETAL_VASC_PATH1")

#are there any significant hits?
lm_significant_hits(lm_HGFVM_auto_basic)
```

```
##             FDR < 0.05 FDR < 0.01
## delB > 0          2245        239
## delB > 0.05        135         14
## delB > 0.1           1          1
```

``` r
lm_significant_hits(lm_HGFVM_auto)
```

```
##             FDR < 0.05 FDR < 0.01
## delB > 0           931         63
## delB > 0.05         56          3
## delB > 0.1           0          0
```

``` r
lm_significant_hits(lm_HGFVM_auto_cells)
```

```
##             FDR < 0.05 FDR < 0.01
## delB > 0             0          0
## delB > 0.05          0          0
## delB > 0.1           0          0
```

``` r
#calculate inflation using the bacon package
inflation(bacon(teststatistics = lm_HGFVM_auto_basic$statistic)) #1.054842  
```

```
## sigma.0 
## 1.06594
```

``` r
inflation(bacon(teststatistics = lm_HGFVM_auto$statistic)) #1.046303 
```

```
##  sigma.0 
## 1.046072
```

``` r
inflation(bacon(teststatistics = lm_HGFVM_auto_cells$statistic)) #0.9479838 
```

```
##   sigma.0 
## 0.9491475
```

#### 3.2.1 Hits


``` r
fvm_hits <- lm_FVM_auto %>% filter(fdr <0.05 & abs(delB) > 0.05)
hg_fvm_hits <- lm_HGFVM_auto %>% filter(fdr <0.05 & abs(delB) > 0.05) 
dim(hg_fvm_hits) #56
```

```
## [1] 56  8
```

``` r
table(hg_fvm_hits$delB>0)
```

```
## 
## FALSE  TRUE 
##    22    34
```

``` r
table(hg_fvm_hits$probeID %in% fvm_hits$probeID) #No
```

```
## 
## FALSE 
##    56
```

``` r
#Do the hits overlap eoPRED CpGs?
table(hg_fvm_hits$probeID %in% eopred_cpgs$ï..probeID) #0
```

```
## 
## FALSE 
##    56
```

``` r
#Do the hits overlap with MVM CpGs?
table(hg_fvm_hits$probeID %in% lm_hg_mvm_hits$probeID ) #0
```

```
## 
## FALSE 
##    56
```

``` r
table(fvm_hits$probeID %in% lm_hg_mvm_hits$probeID )  #0
```

```
## 
## FALSE 
##     1
```

``` r
## Are the hits cell specific DMCs?
table(hg_fvm_hits$probeID %in% cells_dmcs$cpg) #48/56
```

```
## 
## FALSE  TRUE 
##     8    48
```

``` r
## GO analysis
go_hg_fvm<- gometh(sig.cpg = hg_fvm_hits$probeID,
                    all.cpg = rownames(betas_auto),
                    array.type = "EPIC")
```

```
## All input CpGs are used for testing.
```

``` r
go_hg_fvm<- go_hg_fvm %>% mutate(GOterm = rownames(go_hg_fvm)) %>% as_tibble %>% arrange(FDR)
go_hg_fvm %>% filter(FDR < 0.05) %>% nrow() #0 significantly associated terms
```

```
## [1] 0
```

``` r
go_hg_fvm[1:10,]
```

```
## # A tibble: 10 x 7
##    ONTOLOGY TERM                                      N    DE  P.DE   FDR GOterm
##    <chr>    <chr>                                 <dbl> <dbl> <dbl> <dbl> <chr> 
##  1 BP       mitochondrial genome maintenance         30     0 1.00      1 GO:00~
##  2 BP       reproduction                           1400     3 0.610     1 GO:00~
##  3 MF       alpha-1,6-mannosyltransferase activi~     2     0 1.00      1 GO:00~
##  4 MF       trans-hexaprenyltranstransferase act~     2     0 1.00      1 GO:00~
##  5 BP       single strand break repair               11     0 1         1 GO:00~
##  6 MF       single-stranded DNA endodeoxyribonuc~    10     0 1         1 GO:00~
##  7 CC       phosphopyruvate hydratase complex         4     0 1         1 GO:00~
##  8 MF       lactase activity                          1     0 1         1 GO:00~
##  9 BP       alpha-glucoside transport                 2     0 1         1 GO:00~
## 10 BP       regulation of DNA recombination         126     0 1         1 GO:00~
```

``` r
lm_hg_fvm <- lm_HGFVM_auto_basic %>%
  dplyr::select(probeID, delB, fdr) %>% dplyr::rename(Simple_delB = delB, Simple_fdr = fdr) %>%
  left_join((
    lm_HGFVM_auto %>% 
      dplyr::select(probeID, delB, fdr) %>% dplyr::rename(Main_delB = delB, Main_fdr = fdr)), by="probeID") %>%
  left_join((
    lm_HGFVM_auto_cells%>% 
      dplyr::select(probeID, delB, fdr) %>% dplyr::rename(Cells_delB = delB, Cells_fdr = fdr)), by="probeID") %>%
  left_join((
    autoProbes %>%
      dplyr::select(probeID, chr, Start_hg38, UCSC_RefGene_Name, Regulatory_Region)
  ), by="probeID")

saveRDS(lm_hg_fvm, here::here("D. Pathology Project/02_Outputs/AA. Cleaned/LM objects/lm_hg_fvm.rds"))

lm_hg_fvm_hits <- lm_hg_fvm %>%
  filter(Main_fdr<0.05 & abs(Main_delB) >0.05)
dim(lm_hg_fvm_hits) #56 probes
```

```
## [1] 56 11
```

``` r
table(lm_hg_fvm_hits$Regulatory_Region)
```

```
## 
##       Dual   Enhancer  Gene body Intergenic   Promoter 
##         11         26         14          1          4
```

``` r
# write.xlsx(lm_hg_fvm_hits, here::here("D. Pathology Project/02_Outputs/AA. Cleaned/Supplementary File 1.xlsx"),sheetName="HG FVM", append=T,row.names=F)
```

#### 3.2.2 Plots


``` r
#fix genes with duplicates in gene name
lm_hg_fvm$UCSC_RefGene_Name <- str_extract(lm_hg_fvm$UCSC_RefGene_Name, "[^;]+")

fig_2_c <- lm_hg_fvm %>%
  ggplot(aes(x=Main_delB, y=-log10(Main_fdr), color=case_when(
    Main_fdr < 0.05 & abs(Main_delB) > 0.05 ~ "DMP",
    .default = "Not Significant"))) +
  geom_point(alpha=0.6)+
  geom_hline(yintercept = -log10(0.05), color="grey22", linetype="dashed")+
  geom_vline(xintercept=c(-0.05,0.05), color="grey22", linetype="dashed")+
  scale_color_manual(values=c("DMP" ="#3b6b4a","Not Significant" ="grey"))+
  labs(x=expression(paste(Delta*beta, " HG FVM - No FVM")),y="-log(FDR)",color="",
       title = "DNAme ~ HG FVM + GA + Sex + Ancestry",
       subtitle = "56 hits")+
  theme(plot.title = element_text(hjust=0.5, size = 14, face="italic"),
        plot.subtitle = element_text(hjust=0.5, size = 13, color ="#3b6b4a", face="bold"),
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16))+
  annotate(geom="text", x=-0.1, y=1.2, col="black", size = 4, label="FDR < 0.05")+
  geom_text_repel(data=(lm_hg_fvm %>%
                          filter(!is.na(UCSC_RefGene_Name)) %>%
                          filter(abs(Main_delB) > 0.065 & Main_fdr < 0.03|
                                                   abs(Main_delB) >0.05 & Main_fdr<0.01)),
                  aes(x=Main_delB, y=-log10(Main_fdr), label=UCSC_RefGene_Name),
  show.legend=F, max.overlaps = Inf, size=3, min.segment.length=0.2)+
  ylim(0,3.4)

fig_2_d <- lm_hg_fvm %>%
  ggplot(aes(x=Cells_delB, y=-log10(Cells_fdr), color=case_when(
    Cells_fdr < 0.05 & abs(Cells_delB) > 0.05 ~ "DMP",
    .default = "Not Significant"))) +
  geom_point(alpha=0.6)+
  geom_hline(yintercept = -log10(0.05), color="grey22", linetype="dashed")+
  geom_vline(xintercept=c(-0.05,0.05), color="grey22", linetype="dashed")+
  scale_color_manual(values=c("DMP" ="#3b6b4a","Not Significant" ="grey"))+
  labs(x=expression(paste(Delta*beta, " HG FVM - No FVM")),y="-log(FDR)",color="",
       title = "DNAme ~ HG FVM + GA + Sex + Ancestry + Cells",
       subtitle = "0 hits")+
  theme(plot.title = element_text(hjust=0.5, size = 14, face="italic"),
        plot.subtitle = element_text(hjust=0.5, size = 13, color ="#3b6b4a", face="bold"),
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16))+
  annotate(geom="text", x=0.075, y=1.2, col="black", size = 4, label="FDR < 0.05")+
  ylim(0,3.4)
```

#### 3.2.3 Cells Plots


``` r
betas_gene2 <- t(betas_noob[rownames(betas_noob) %in% lm_hg_fvm_hits$probeID,
                            colnames(betas_noob) %in% pDat_filt$Sentrix]) %>%
  as.data.frame() %>%
  rownames_to_column(var="Sentrix") 

betas_gene2 <- betas_gene2 %>%
  # reshape into longer format
  pivot_longer(cols = -Sentrix,
               names_to = 'probeID', 
               values_to = 'Beta')  %>%
  # add tissue and trimester info
  left_join(pDat_filt %>% select(Sentrix, Trimester, Tissue), by = 'Sentrix') %>%
  # calculate mean and sd for each cpg for each group
  group_by(Tissue, Trimester, probeID) %>%
  summarize(mean = mean(Beta),
            sd = sd(Beta)) %>%
  mutate(lower = mean-sd, upper = mean+sd) %>%
  filter(Trimester == "Third")
```

```
## `summarise()` has grouped output by 'Tissue', 'Trimester'. You can override
## using the `.groups` argument.
```

``` r
## see how the DMPs look in the individual cell types
betas_gene2 %>%
  filter(probeID %in% lm_hg_fvm_hits$probeID)%>%
  left_join(lm_hg_fvm_hits, by="probeID") %>%
  mutate(Vars = Main_delB >0) %>%
  ggplot(aes(y=mean, x=Tissue, color=Tissue))+
  geom_point()+
  geom_boxplot(alpha=0.5,outlier.shape=NA)+
  theme(legend.position="none")+
  facet_wrap(vars(Vars), labeller = labeller(Vars = labs))+
  labs(y="Mean Beta")
```

![](02_SPAH_Pathology_Linear_Modelling_2026_files/figure-html/cells-plots-fvm-1.png)<!-- -->

``` r
betas_gene2 %>%
  filter(probeID %in% lm_hg_fvm_hits$probeID)%>%
  left_join(lm_hg_fvm_hits, by="probeID") %>%
  mutate(Vars = Main_delB >0) %>%
  ggplot(aes(y=mean, x=Tissue, color=probeID, group=probeID))+
  geom_point()+
  geom_line()+
  theme(legend.position="none")+
  facet_wrap(vars(Vars), labeller = labeller(Vars = labs))+
  labs(y="Mean Beta")
```

![](02_SPAH_Pathology_Linear_Modelling_2026_files/figure-html/cells-plots-fvm-2.png)<!-- -->

``` r
## plot all of the 54 DMPs by individual cell-type
sup_cells_fvm_a <- betas_gene2 %>%
  filter(!Tissue == "Villi") %>%
  left_join(lm_hg_fvm_hits, by="probeID") %>%
  mutate(Vars = Main_delB >0) %>%
  filter(Vars == FALSE) %>%
  ggplot() +
  geom_linerange(alpha = 0.5, size = 1,
                 aes(x = probeID, ymin =lower, ymax = upper, color = Tissue)) +
  geom_point(alpha = 1, aes(x = reorder(probeID,Main_delB), y = mean, color = Tissue), size = 2) +
  theme(plot.title = element_text(hjust=0.5, size = 16, face="bold"),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=12),
        legend.position="none",
        strip.background = element_rect(color="black", fill="grey88")) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_x_discrete(expand = c(0.01,0.01))+
  scale_color_manual(values= c('#6A1B9A', '#1565C0','#388E3C',    '#E64A19',
                               '#FBC02D', '#C62828')) +
  labs(y = "DNA Methylation", x = "", color = '')+
  facet_wrap(vars(Vars), labeller = labeller(Vars = labs), nrow=1)

sup_cells_fvm_b <- betas_gene2 %>%
  filter(!Tissue == "Villi") %>%
  left_join(lm_hg_fvm_hits, by="probeID") %>%
  mutate(Vars = Main_delB >0) %>%
  filter(Vars == TRUE) %>%
  ggplot() +
  geom_linerange(alpha = 0.5, size = 1,
                 aes(x = probeID, ymin =lower, ymax = upper, color = Tissue)) +
  geom_point(alpha = 1, aes(x = reorder(probeID,Main_delB), y = mean, color = Tissue), size = 2) +
  theme(plot.title = element_text(hjust=0.5, size = 16, face="bold"),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        legend.position="right",
        strip.background = element_rect(color="black", fill="grey88")) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_x_discrete(expand = c(0.01,0.01))+
  scale_color_manual(values= c('#6A1B9A', '#1565C0','#388E3C',    '#E64A19',
                               '#FBC02D', '#C62828')) +
  labs(y = "", x = "", color = "")+
  facet_wrap(vars(Vars), labeller = labeller(Vars = labs), nrow=1)

sup_cells_fvm <- plot_grid(sup_cells_fvm_a, NULL, sup_cells_fvm_b, nrow=1, rel_widths=c(0.6,-0.03,1))
sup_cells_fvm
```

![](02_SPAH_Pathology_Linear_Modelling_2026_files/figure-html/cells-plots-fvm-3.png)<!-- -->

## 4.0 AI

### 4.1 AI vs. No AI


``` r
table(metadata_naomit$acute_inflammation_3cat)
```

```
## 
##   0   1   2 
## 216 202  65
```

``` r
#set-up models
mod_AI_auto_basic <- model.matrix(~ACUTE_INFLAMMATION,metadata_naomit,row.names=T)

mod_AI_auto <- model.matrix(~ACUTE_INFLAMMATION+
                              ga_continuous+
                              Sex_Predic+
                              Prob_African+
                              Prob_European,
                            metadata_naomit,
                            row.names=T)

mod_AI_auto_cells <- model.matrix(~ACUTE_INFLAMMATION+
                              ga_continuous+
                              Sex_Predic+
                              Prob_African+
                              Prob_European +
                              Syncytiotrophoblast+
                              Stromal+
                              Endothelial+
                              nRBC,
                            metadata_naomit,
                            row.names=T)

head(mod_AI_auto_cells)
```

```
##   (Intercept) ACUTE_INFLAMMATION1 ga_continuous Sex_PredicXY Prob_African
## 1           1                   1         30.71            1 0.0019729407
## 2           1                   1         39.29            0 0.0087901395
## 3           1                   0         39.43            1 0.0010767692
## 4           1                   1         41.29            0 0.0008266452
## 5           1                   1         39.00            1 0.0108689808
## 6           1                   0         38.00            1 0.0012142192
##   Prob_European Syncytiotrophoblast     Stromal Endothelial       nRBC
## 1     0.9944981           0.7721641 0.017105452  0.07351688 0.01571232
## 2     0.8665399           0.8344046 0.019953485  0.07704861 0.01673799
## 3     0.9949846           0.9196496 0.012145654  0.04171552 0.01755493
## 4     0.9949141           0.8725699 0.020448455  0.06373648 0.02670345
## 5     0.7087023           0.9094798 0.006426069  0.06016573 0.01993460
## 6     0.9949257           0.8239425 0.042196156  0.05986966 0.01807088
```

``` r
#run regression
lm_AI_auto_basic <-linear_modelling (betas_auto,mvals_auto,mod_AI_auto_basic) %>% filter (term == "ACUTE_INFLAMMATION1")
lm_AI_auto <-linear_modelling (betas_auto,mvals_auto,mod_AI_auto) %>% filter (term == "ACUTE_INFLAMMATION1")
lm_AI_auto_cells <-linear_modelling (betas_auto,mvals_auto,mod_AI_auto_cells) %>% filter (term == "ACUTE_INFLAMMATION1")

#are there any significant hits?
lm_significant_hits(lm_AI_auto_basic)
```

```
##             FDR < 0.05 FDR < 0.01
## delB > 0             0          0
## delB > 0.05          0          0
## delB > 0.1           0          0
```

``` r
lm_significant_hits(lm_AI_auto)
```

```
##             FDR < 0.05 FDR < 0.01
## delB > 0             0          0
## delB > 0.05          0          0
## delB > 0.1           0          0
```

``` r
lm_significant_hits(lm_AI_auto_cells)
```

```
##             FDR < 0.05 FDR < 0.01
## delB > 0             0          0
## delB > 0.05          0          0
## delB > 0.1           0          0
```

``` r
#calculate inflation using the bacon package
inflation(bacon(teststatistics = lm_AI_auto_basic$statistic)) #0.976918
```

```
##   sigma.0 
## 0.9769249
```

``` r
inflation(bacon(teststatistics = lm_AI_auto$statistic)) #0.9073235 
```

```
##   sigma.0 
## 0.9073364
```

``` r
inflation(bacon(teststatistics = lm_AI_auto_cells$statistic)) #0.9587326  
```

```
##  sigma.0 
## 0.958732
```

### 4.2 High-grade AI vs. No AI


``` r
#remove samples with low grade AI
ai_grade_naomit <- metadata_naomit %>% filter(!acute_inflammation_3cat == 1)
table(ai_grade_naomit$acute_inflammation_3cat) #216 No AI, 65 HG AI
```

```
## 
##   0   2 
## 216  65
```

``` r
betas_auto_ai <- betas_auto[,ai_grade_naomit$rgset_colnames]
mvals_auto_ai <- mvals_auto[,ai_grade_naomit$rgset_colnames]
dim(mvals_auto_ai) #x281
```

```
## [1] 732102    281
```

``` r
#set-up models
mod_HGAI_auto_basic <- model.matrix(~ACUTE_INFLAMMATION,ai_grade_naomit,row.names=T)

mod_HGAI_auto <- model.matrix(~ACUTE_INFLAMMATION+
                              ga_continuous+
                              Sex_Predic+
                              Prob_African+
                              Prob_European,
                            ai_grade_naomit,
                            row.names=T)

mod_HGAI_auto_cells <- model.matrix(~ACUTE_INFLAMMATION+
                              ga_continuous+
                              Sex_Predic+
                              Prob_African+
                              Prob_European +
                              Syncytiotrophoblast+
                              Stromal+
                              Endothelial+
                              nRBC,
                            ai_grade_naomit,
                            row.names=T)

head(mod_HGAI_auto_cells)
```

```
##   (Intercept) ACUTE_INFLAMMATION1 ga_continuous Sex_PredicXY Prob_African
## 1           1                   1         39.29            0  0.008790140
## 2           1                   0         39.43            1  0.001076769
## 3           1                   0         38.00            1  0.001214219
## 4           1                   0         37.29            0  0.002362143
## 5           1                   1         40.71            1  0.983179884
## 6           1                   0         33.14            0  0.007721731
##   Prob_European Syncytiotrophoblast    Stromal Endothelial       nRBC
## 1    0.86653994           0.8344046 0.01995348 0.077048609 0.01673799
## 2    0.99498462           0.9196496 0.01214565 0.041715517 0.01755493
## 3    0.99492574           0.8239425 0.04219616 0.059869660 0.01807088
## 4    0.99300050           0.6885865 0.05808867 0.077162778 0.01645350
## 5    0.01153420           0.9767056 0.00000000 0.009968167 0.01332622
## 6    0.05328532           0.9503356 0.00000000 0.018286540 0.02808383
```

``` r
#run regression
lm_HGAI_auto_basic <-linear_modelling (betas_auto_ai,mvals_auto_ai,mod_HGAI_auto_basic) %>% filter (term == "ACUTE_INFLAMMATION1")
lm_HGAI_auto <-linear_modelling (betas_auto_ai,mvals_auto_ai,mod_HGAI_auto) %>% filter (term == "ACUTE_INFLAMMATION1")
lm_HGAI_auto_cells <-linear_modelling (betas_auto_ai,mvals_auto_ai,mod_HGAI_auto_cells) %>% filter (term == "ACUTE_INFLAMMATION1")

#are there any significant hits?
lm_significant_hits(lm_HGAI_auto_basic)
```

```
##             FDR < 0.05 FDR < 0.01
## delB > 0             0          0
## delB > 0.05          0          0
## delB > 0.1           0          0
```

``` r
lm_significant_hits(lm_HGAI_auto)
```

```
##             FDR < 0.05 FDR < 0.01
## delB > 0             0          0
## delB > 0.05          0          0
## delB > 0.1           0          0
```

``` r
lm_significant_hits(lm_HGAI_auto_cells)
```

```
##             FDR < 0.05 FDR < 0.01
## delB > 0             0          0
## delB > 0.05          0          0
## delB > 0.1           0          0
```

``` r
#calculate inflation using the bacon package
inflation(bacon(teststatistics = lm_HGAI_auto_basic$statistic)) #0.9357896  
```

```
## sigma.0 
## 0.93528
```

``` r
inflation(bacon(teststatistics = lm_HGAI_auto$statistic)) #0.9879485
```

```
##  sigma.0 
## 0.987971
```

``` r
inflation(bacon(teststatistics = lm_HGAI_auto_cells$statistic)) #1.032015  
```

```
##  sigma.0 
## 1.032012
```

## 5.0 CI

### 5.1 CI vs. No CI


``` r
table(metadata_naomit$chronic_inflammation_3cat)
```

```
## 
##   0   1   2 
## 217 161 105
```

``` r
#set-up models
mod_CI_auto_basic <- model.matrix(~CHRONIC_INFLAMMATION,metadata_naomit,row.names=T)

mod_CI_auto <- model.matrix(~CHRONIC_INFLAMMATION+
                              ga_continuous+
                              Sex_Predic+
                              Prob_African+
                              Prob_European,
                            metadata_naomit,
                            row.names=T)

mod_CI_auto_cells <- model.matrix(~CHRONIC_INFLAMMATION+
                              ga_continuous+
                              Sex_Predic+
                              Prob_African+
                              Prob_European +
                              Syncytiotrophoblast+
                              Stromal+
                              Endothelial+
                              nRBC,
                            metadata_naomit,
                            row.names=T)

head(mod_CI_auto_cells)
```

```
##   (Intercept) CHRONIC_INFLAMMATION1 ga_continuous Sex_PredicXY Prob_African
## 1           1                     0         30.71            1 0.0019729407
## 2           1                     1         39.29            0 0.0087901395
## 3           1                     0         39.43            1 0.0010767692
## 4           1                     0         41.29            0 0.0008266452
## 5           1                     0         39.00            1 0.0108689808
## 6           1                     0         38.00            1 0.0012142192
##   Prob_European Syncytiotrophoblast     Stromal Endothelial       nRBC
## 1     0.9944981           0.7721641 0.017105452  0.07351688 0.01571232
## 2     0.8665399           0.8344046 0.019953485  0.07704861 0.01673799
## 3     0.9949846           0.9196496 0.012145654  0.04171552 0.01755493
## 4     0.9949141           0.8725699 0.020448455  0.06373648 0.02670345
## 5     0.7087023           0.9094798 0.006426069  0.06016573 0.01993460
## 6     0.9949257           0.8239425 0.042196156  0.05986966 0.01807088
```

``` r
#run regression
lm_CI_auto_basic <-linear_modelling (betas_auto,mvals_auto,mod_CI_auto_basic) %>% filter (term == "CHRONIC_INFLAMMATION1")
lm_CI_auto <-linear_modelling (betas_auto,mvals_auto,mod_CI_auto) %>% filter (term == "CHRONIC_INFLAMMATION1")
lm_CI_auto_cells <-linear_modelling (betas_auto,mvals_auto,mod_CI_auto_cells) %>% filter (term == "CHRONIC_INFLAMMATION1")

#are there any significant hits?
lm_significant_hits(lm_CI_auto_basic)
```

```
##             FDR < 0.05 FDR < 0.01
## delB > 0             0          0
## delB > 0.05          0          0
## delB > 0.1           0          0
```

``` r
lm_significant_hits(lm_CI_auto)
```

```
##             FDR < 0.05 FDR < 0.01
## delB > 0             0          0
## delB > 0.05          0          0
## delB > 0.1           0          0
```

``` r
lm_significant_hits(lm_CI_auto_cells)
```

```
##             FDR < 0.05 FDR < 0.01
## delB > 0             0          0
## delB > 0.05          0          0
## delB > 0.1           0          0
```

``` r
#calculate inflation using the bacon package
inflation(bacon(teststatistics = lm_CI_auto_basic$statistic)) #0.9158924
```

```
##   sigma.0 
## 0.9159065
```

``` r
inflation(bacon(teststatistics = lm_CI_auto$statistic)) #0.9259814  
```

```
##   sigma.0 
## 0.9259692
```

``` r
inflation(bacon(teststatistics = lm_CI_auto_cells$statistic)) #0.9280374  
```

```
##   sigma.0 
## 0.9280388
```

### 5.2 High-grade CI vs. No CI


``` r
#remove samples with low grade CI
ci_grade_naomit <- metadata_naomit %>% filter(!chronic_inflammation_3cat == 1)
table(ci_grade_naomit$chronic_inflammation_3cat) #217 No CI, 105 HG CI
```

```
## 
##   0   2 
## 217 105
```

``` r
betas_auto_ci <- betas_auto[,ci_grade_naomit$rgset_colnames]
mvals_auto_ci <- mvals_auto[,ci_grade_naomit$rgset_colnames]
dim(mvals_auto_ci) #x322
```

```
## [1] 732102    322
```

``` r
#set-up models
mod_HGCI_auto_basic <- model.matrix(~CHRONIC_INFLAMMATION,ci_grade_naomit,row.names=T)

mod_HGCI_auto <- model.matrix(~CHRONIC_INFLAMMATION+
                              ga_continuous+
                              Sex_Predic+
                              Prob_African+
                              Prob_European,
                            ci_grade_naomit,
                            row.names=T)

mod_HGCI_auto_cells <- model.matrix(~CHRONIC_INFLAMMATION+
                              ga_continuous+
                              Sex_Predic+
                              Prob_African+
                              Prob_European +
                              Syncytiotrophoblast+
                              Stromal+
                              Endothelial+
                              nRBC,
                            ci_grade_naomit,
                            row.names=T)

head(mod_HGCI_auto_cells)
```

```
##   (Intercept) CHRONIC_INFLAMMATION1 ga_continuous Sex_PredicXY Prob_African
## 1           1                     0         30.71            1 0.0019729407
## 2           1                     0         39.43            1 0.0010767692
## 3           1                     0         41.29            0 0.0008266452
## 4           1                     0         39.00            1 0.0108689808
## 5           1                     0         38.00            1 0.0012142192
## 6           1                     0         37.71            1 0.0010489338
##   Prob_European Syncytiotrophoblast     Stromal Endothelial       nRBC
## 1     0.9944981           0.7721641 0.017105452  0.07351688 0.01571232
## 2     0.9949846           0.9196496 0.012145654  0.04171552 0.01755493
## 3     0.9949141           0.8725699 0.020448455  0.06373648 0.02670345
## 4     0.7087023           0.9094798 0.006426069  0.06016573 0.01993460
## 5     0.9949257           0.8239425 0.042196156  0.05986966 0.01807088
## 6     0.9937753           0.9091273 0.030365924  0.03837857 0.01329400
```

``` r
#run regression
lm_HGCI_auto_basic <-linear_modelling (betas_auto_ci,mvals_auto_ci,mod_HGCI_auto_basic)%>% filter (term == "CHRONIC_INFLAMMATION1")
lm_HGCI_auto <-linear_modelling (betas_auto_ci,mvals_auto_ci,mod_HGCI_auto)%>% filter (term == "CHRONIC_INFLAMMATION1")
lm_HGCI_auto_cells <-linear_modelling (betas_auto_ci,mvals_auto_ci,mod_HGCI_auto_cells)%>% filter (term == "CHRONIC_INFLAMMATION1")

#are there any significant hits?
lm_significant_hits(lm_HGCI_auto_basic)
```

```
##             FDR < 0.05 FDR < 0.01
## delB > 0             0          0
## delB > 0.05          0          0
## delB > 0.1           0          0
```

``` r
lm_significant_hits(lm_HGCI_auto)
```

```
##             FDR < 0.05 FDR < 0.01
## delB > 0             0          0
## delB > 0.05          0          0
## delB > 0.1           0          0
```

``` r
lm_significant_hits(lm_HGCI_auto_cells)
```

```
##             FDR < 0.05 FDR < 0.01
## delB > 0             0          0
## delB > 0.05          0          0
## delB > 0.1           0          0
```

``` r
#calculate inflation using the bacon package
inflation(bacon(teststatistics = lm_HGCI_auto_basic$statistic)) #0.96176 
```

```
##   sigma.0 
## 0.9617684
```

``` r
inflation(bacon(teststatistics = lm_HGCI_auto$statistic)) #0.971272 
```

```
##   sigma.0 
## 0.9712896
```

``` r
inflation(bacon(teststatistics = lm_HGCI_auto_cells$statistic)) #1.024292
```

```
##  sigma.0 
## 1.024317
```

## 6.0 Manuscript Figures

### 6.1 Figure 2.0


``` r
fig_2 <- plot_grid(fig_2_a, NULL,fig_2_b, 
                   NULL, NULL, NULL,
                   fig_2_c, NULL, fig_2_d,
                   labels = c("A.","","B.","","","","C.","","D."), nrow = 3, ncol=3,
                   rel_widths = c(1,0.1,1), rel_heights = c(1,0.1,1))

fig_2
```

![](02_SPAH_Pathology_Linear_Modelling_2026_files/figure-html/fig-4-1.png)<!-- -->

``` r
ggsave("Figure 2.png", plot = fig_2, path = here::here("D. Pathology Project/02_Outputs/AA. Cleaned/"), width = 11, height =9, device = "png")
```

### 6.2 Supplementary Figure 2


``` r
qq_plot <- function(df){
  
  plot <- df %>%
    dplyr::select(p.value) %>%
    arrange(p.value) %>%
    mutate(Exp.p.value = row_number()/nrow(df)) %>% #calculate based on uniform distribution 
    ggplot(aes(x=-log10(Exp.p.value), y=-log10(p.value)))+
    geom_point()+
    labs(x="Expected -log10(p-value)", y="Observed -log10(p-value)")+
    geom_abline(slope=1, intercept =0)+
    annotate(geom="text", x=1, y=5, col="black", size = 4,
             label=paste0("\u03BB = ", round(inflation(bacon(teststatistics = df$statistic)),3)))
  return(plot)
  
}

title_ai <- ggdraw() + draw_label("AI",fontface = 'bold', size=12,hjust = 0.5) 
title_ci <- ggdraw() + draw_label("CI",fontface = 'bold', size=12,hjust = 0.5)
title_fvm <- ggdraw() + draw_label("FVM",fontface = 'bold', size=12,hjust = 0.5)
title_mvm <- ggdraw() + draw_label("MVM",fontface = 'bold', size=12,hjust = 0.5)

title_hgai <- ggdraw() + draw_label("HG AI",fontface = 'bold', size=12,hjust = 0.5) 
title_hgci <- ggdraw() + draw_label("HG CI",fontface = 'bold', size=12,hjust = 0.5)
title_hgfvm <- ggdraw() + draw_label("HG FVM",fontface = 'bold', size=12,hjust = 0.5)
title_hgmvm <- ggdraw() + draw_label("HG MVM",fontface = 'bold', size=12,hjust = 0.5)

qq_all <- plot_grid(
  
  title_ai,
  title_ci,
  title_fvm,
  title_mvm,
  
  qq_plot(lm_AI_auto_basic),
  qq_plot(lm_CI_auto_basic),
  qq_plot(lm_FVM_auto_basic),
  qq_plot(lm_mvm_auto_basic),
  
  qq_plot(lm_AI_auto),
  qq_plot(lm_CI_auto),
  qq_plot(lm_FVM_auto),
  qq_plot(lm_mvm_auto),
  
  qq_plot(lm_AI_auto_cells),
  qq_plot(lm_CI_auto_cells),
  qq_plot(lm_FVM_auto_cells),
  qq_plot(lm_mvm_cells),
  
  title_hgai,
  title_hgci,
  title_hgfvm,
  title_hgmvm,
  
  qq_plot(lm_HGAI_auto_basic),
  qq_plot(lm_HGCI_auto_basic),
  qq_plot(lm_HGFVM_auto_basic),
  qq_plot(lm_HGMVM_auto_basic),
  
  qq_plot(lm_HGAI_auto),
  qq_plot(lm_HGCI_auto),
  qq_plot(lm_HGFVM_auto),
  qq_plot(lm_HGMVM_auto),
  
  qq_plot(lm_HGAI_auto_cells),
  qq_plot(lm_HGCI_auto_cells),
  qq_plot(lm_HGFVM_auto_cells),
  qq_plot(lm_HGMVM_auto_cells),
  
  nrow=8, ncol=4,
  rel_heights = c(0.2,1,1,1,0.2,1,1,1)
)

ggsave("Figure S2.png", plot = qq_all, path = here::here("D. Pathology Project/02_Outputs/AA. Cleaned/"), width = 15, height =19, device = "png")
```

### 6.3 Supplementary Figure 3


``` r
sup_fig_sens <- plot_grid(sup_fig_sens_a, sup_fig_sens_b, sup_fig_sens_c,
                          nrow=2, ncol=2, 
                          labels = c("A.","B.", "C."))

ggsave("Figure S3.png", plot = sup_fig_sens, path = here::here("D. Pathology Project/02_Outputs/AA. Cleaned/"), width = 13, height = 8, device = "png")
```

### 6.4 Supplementary Figure 5


``` r
title_mvm <- ggdraw() + 
  draw_label("High-grade MVM DMPs",
    fontface = 'bold',
    size=16,
    hjust = 0.5) 

title_fvm <- ggdraw() + 
  draw_label("High-grade FVM DMPs",
    fontface = 'bold',
    size=16,
    hjust = 0.5) 

sup_cells_figure <- plot_grid(title_mvm, sup_cells_mvm, 
                              title_fvm, sup_cells_fvm,
                              nrow=4, rel_heights=c(0.1,1,0.1,1),
                              labels=c("A.","","B.",""))

ggsave("Figure S6.png", plot = sup_cells_figure, path = here::here("D. Pathology Project/02_Outputs/AA. Cleaned/"), width = 16, height =12, device = "png")
```

## Save Outputs


``` r
saveRDS(lm_mvm_auto_basic, here::here("D. Pathology Project/02_Outputs/AA. Cleaned/LM objects/lm_mvm_auto_basic.rds"))
saveRDS(lm_mvm_auto, here::here("D. Pathology Project/02_Outputs/AA. Cleaned/LM objects/lm_mvm_auto.rds"))
saveRDS(lm_mvm_cells, here::here("D. Pathology Project/02_Outputs/AA. Cleaned/LM objects/lm_mvm_cells.rds"))

saveRDS(lm_HGMVM_auto_basic, here::here("D. Pathology Project/02_Outputs/AA. Cleaned/LM objects/lm_HGMVM_auto_basic.rds"))
saveRDS(lm_HGMVM_auto, here::here("D. Pathology Project/02_Outputs/AA. Cleaned/LM objects/lm_HGMVM_auto.rds"))
saveRDS(lm_HGMVM_auto_cells, here::here("D. Pathology Project/02_Outputs/AA. Cleaned/LM objects/lm_HGMVM_auto_cells.rds"))

saveRDS(lm_FVM_auto_basic, here::here("D. Pathology Project/02_Outputs/AA. Cleaned/LM objects/lm_FVM_auto_basic.rds"))
saveRDS(lm_FVM_auto, here::here("D. Pathology Project/02_Outputs/AA. Cleaned/LM objects/lm_FVM_auto.rds"))
saveRDS(lm_FVM_auto_cells, here::here("D. Pathology Project/02_Outputs/AA. Cleaned/LM objects/lm_FVM_auto_cells.rds"))

saveRDS(lm_HGFVM_auto_basic, here::here("D. Pathology Project/02_Outputs/AA. Cleaned/LM objects/lm_HGFVM_auto_basic.rds"))
saveRDS(lm_HGFVM_auto, here::here("D. Pathology Project/02_Outputs/AA. Cleaned/LM objects/lm_HGFVM_auto.rds"))
saveRDS(lm_HGFVM_auto_cells, here::here("D. Pathology Project/02_Outputs/AA. Cleaned/LM objects/lm_HGFVM_auto_cells.rds"))

saveRDS(lm_AI_auto_basic, here::here("D. Pathology Project/02_Outputs/AA. Cleaned/LM objects/lm_AI_auto_basic.rds"))
saveRDS(lm_AI_auto, here::here("D. Pathology Project/02_Outputs/AA. Cleaned/LM objects/lm_AI_auto.rds"))
saveRDS(lm_AI_auto_cells, here::here("D. Pathology Project/02_Outputs/AA. Cleaned/LM objects/lm_AI_auto_cells.rds"))

saveRDS(lm_HGAI_auto_basic, here::here("D. Pathology Project/02_Outputs/AA. Cleaned/LM objects/lm_HGAI_auto_basic.rds"))
saveRDS(lm_HGAI_auto, here::here("D. Pathology Project/02_Outputs/AA. Cleaned/LM objects/lm_HGAI_auto.rds"))
saveRDS(lm_HGAI_auto_cells, here::here("D. Pathology Project/02_Outputs/AA. Cleaned/LM objects/lm_HGAI_auto_cells.rds"))

saveRDS(lm_CI_auto_basic, here::here("D. Pathology Project/02_Outputs/AA. Cleaned/LM objects/lm_CI_auto_basic.rds"))
saveRDS(lm_CI_auto, here::here("D. Pathology Project/02_Outputs/AA. Cleaned/LM objects/lm_CI_auto.rds"))
saveRDS(lm_CI_auto_cells, here::here("D. Pathology Project/02_Outputs/AA. Cleaned/LM objects/lm_CI_auto_cells.rds"))

saveRDS(lm_HGCI_auto_basic, here::here("D. Pathology Project/02_Outputs/AA. Cleaned/LM objects/lm_HGCI_auto_basic.rds"))
saveRDS(lm_HGCI_auto, here::here("D. Pathology Project/02_Outputs/AA. Cleaned/LM objects/lm_HGCI_auto.rds"))
saveRDS(lm_HGCI_auto_cells, here::here("D. Pathology Project/02_Outputs/AA. Cleaned/LM objects/lm_HGCI_auto_cells.rds"))
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
##  [1] IlluminaHumanMethylationEPICmanifest_0.3.0         
##  [2] ggpubr_0.6.0                                       
##  [3] bacon_1.28.0                                       
##  [4] ellipse_0.5.0                                      
##  [5] BiocParallel_1.34.2                                
##  [6] missMethyl_1.34.0                                  
##  [7] IlluminaHumanMethylationEPICanno.ilm10b4.hg19_0.6.0
##  [8] IlluminaHumanMethylation450kanno.ilmn12.hg19_0.6.1 
##  [9] minfi_1.46.0                                       
## [10] bumphunter_1.42.0                                  
## [11] locfit_1.5-9.8                                     
## [12] iterators_1.0.14                                   
## [13] foreach_1.5.2                                      
## [14] Biostrings_2.68.1                                  
## [15] XVector_0.40.0                                     
## [16] SummarizedExperiment_1.30.2                        
## [17] MatrixGenerics_1.12.3                              
## [18] matrixStats_1.0.0                                  
## [19] GenomicRanges_1.52.1                               
## [20] GenomeInfoDb_1.36.1                                
## [21] IRanges_2.34.1                                     
## [22] S4Vectors_0.38.1                                   
## [23] plomics_0.2.0                                      
## [24] ggrepel_0.9.3                                      
## [25] planet_1.8.0                                       
## [26] irlba_2.3.5.1                                      
## [27] Matrix_1.6-1                                       
## [28] cowplot_1.1.3                                      
## [29] biobroom_1.32.0                                    
## [30] broom_1.0.5                                        
## [31] lumi_2.52.0                                        
## [32] Biobase_2.60.0                                     
## [33] BiocGenerics_0.46.0                                
## [34] limma_3.56.2                                       
## [35] xlsx_0.6.5                                         
## [36] readxl_1.4.5                                       
## [37] lubridate_1.9.4                                    
## [38] forcats_1.0.0                                      
## [39] stringr_1.5.0                                      
## [40] dplyr_1.1.4                                        
## [41] purrr_1.0.2                                        
## [42] readr_2.1.4                                        
## [43] tidyr_1.3.0                                        
## [44] tibble_3.2.1                                       
## [45] ggplot2_3.5.2                                      
## [46] tidyverse_2.0.0                                    
## [47] here_1.0.1                                         
## 
## loaded via a namespace (and not attached):
##   [1] splines_4.3.1             BiocIO_1.10.0            
##   [3] bitops_1.0-7              filelock_1.0.2           
##   [5] BiasedUrn_2.0.11          cellranger_1.1.0         
##   [7] preprocessCore_1.62.1     methylumi_2.46.0         
##   [9] XML_3.99-0.14             lifecycle_1.0.4          
##  [11] rstatix_0.7.2             rprojroot_2.0.3          
##  [13] lattice_0.21-8            MASS_7.3-60              
##  [15] base64_2.0.1              scrime_1.3.5             
##  [17] backports_1.5.0           magrittr_2.0.3           
##  [19] sass_0.4.10               rmarkdown_2.29           
##  [21] jquerylib_0.1.4           yaml_2.3.10              
##  [23] doRNG_1.8.6               askpass_1.1              
##  [25] DBI_1.1.3                 RColorBrewer_1.1-3       
##  [27] abind_1.4-5               zlibbioc_1.46.0          
##  [29] quadprog_1.5-8            RCurl_1.98-1.12          
##  [31] xlsxjars_0.6.1            rappdirs_0.3.3           
##  [33] GenomeInfoDbData_1.2.10   genefilter_1.82.1        
##  [35] annotate_1.78.0           DelayedMatrixStats_1.22.5
##  [37] codetools_0.2-19          DelayedArray_0.26.7      
##  [39] xml2_1.3.8                tidyselect_1.2.1         
##  [41] farver_2.1.2              beanplot_1.3.1           
##  [43] BiocFileCache_2.8.0       illuminaio_0.42.0        
##  [45] GenomicAlignments_1.36.0  jsonlite_2.0.0           
##  [47] multtest_2.56.0           survival_3.5-7           
##  [49] systemfonts_1.2.3         tools_4.3.1              
##  [51] progress_1.2.2            ragg_1.2.5               
##  [53] Rcpp_1.0.14               glue_1.8.0               
##  [55] xfun_0.52                 mgcv_1.9-0               
##  [57] HDF5Array_1.28.1          withr_3.0.2              
##  [59] BiocManager_1.30.22       fastmap_1.2.0            
##  [61] rhdf5filters_1.12.1       fansi_1.0.6              
##  [63] openssl_2.1.0             digest_0.6.37            
##  [65] timechange_0.3.0          R6_2.5.1                 
##  [67] textshaping_1.0.1         colorspace_2.1-1         
##  [69] GO.db_3.17.0              dichromat_2.0-0.1        
##  [71] biomaRt_2.56.1            RSQLite_2.3.1            
##  [73] utf8_1.2.5                generics_0.1.3           
##  [75] data.table_1.17.6         rtracklayer_1.60.1       
##  [77] prettyunits_1.1.1         httr_1.4.7               
##  [79] S4Arrays_1.0.5            pkgconfig_2.0.3          
##  [81] rJava_1.0-6               gtable_0.3.6             
##  [83] blob_1.2.4                siggenes_1.74.0          
##  [85] htmltools_0.5.8.1         carData_3.0-5            
##  [87] scales_1.4.0              png_0.1-8                
##  [89] knitr_1.50                rstudioapi_0.17.1        
##  [91] tzdb_0.4.0                rjson_0.2.21             
##  [93] nlme_3.1-163              curl_5.0.2               
##  [95] org.Hs.eg.db_3.17.0       cachem_1.1.0             
##  [97] rhdf5_2.44.0              KernSmooth_2.23-22       
##  [99] AnnotationDbi_1.62.2      restfulr_0.0.15          
## [101] GEOquery_2.68.0           pillar_1.9.0             
## [103] grid_4.3.1                reshape_0.8.9            
## [105] vctrs_0.6.5               car_3.1-2                
## [107] dbplyr_2.3.3              xtable_1.8-4             
## [109] evaluate_1.0.4            GenomicFeatures_1.52.1   
## [111] cli_3.6.3                 compiler_4.3.1           
## [113] Rsamtools_2.16.0          rlang_1.1.4              
## [115] crayon_1.5.2              rngtools_1.5.2           
## [117] ggsignif_0.6.4            labeling_0.4.3           
## [119] nor1mix_1.3-0             mclust_6.0.0             
## [121] affy_1.78.2               plyr_1.8.9               
## [123] stringi_1.8.7             nleqslv_3.3.4            
## [125] hms_1.1.3                 sparseMatrixStats_1.12.2 
## [127] bit64_4.0.5               Rhdf5lib_1.22.0          
## [129] statmod_1.5.0             KEGGREST_1.40.0          
## [131] memoise_2.0.1             affyio_1.70.0            
## [133] bslib_0.9.0               bit_4.0.5
```
