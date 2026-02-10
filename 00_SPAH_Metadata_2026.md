---
title: "SPAH Cohort Metadata Updates"
author: "Hannah Illing"
date: "05 February, 2026"
output: 
 html_document:
    keep_md: yes
    toc: true
    toc_depth: 5
    toc_float:
      collapsed: true 
      smooth_scroll: true
editor_options: 
  chunk_output_type: inline
  
---

## 0.0 Introduction

In this script I load the raw clinical, demographic, and array run ('metadata') data for all placentas in the **S**tress **P**regnancy **a**nd **H**ealth (SPAH) study. I combine all of the metadata into one object and calculate some additional variables. 

NOTE: if you are using the data from GSE307289 then much of this data will be already available and formatted

## 1.0 Set-up

### 1.1 Load Packages 


``` r
library(here)
library(tidyverse)
library(readxl)
library(minfi)
library(wateRmelon)
library(tableone)
library(plomics)
library(planet)
library(EpiDISH)
library(ggpubr)
library(ExperimentHub) 
library(mixOmics) 

set.seed(4815) #set seed 
```

### 1.2 Load Data

In this section we load all of the files with sample information:

1. `pData_SPAH` - this is the metadata  This file contains information for all **605** pregnant mothers part of the SPAH cohort. Each individual has a unique identifier in the `ID` column, which matches the `Sample_ID` column in the `array_qc_ss` file.
2. `array_qc_ss` - this is the file with the information on how the samples were loaded as well as the calculated QC/methylation values created in the first QC script. This contains information for all **520** placentas we have methylation data for. The methylation arrays were conducted by the Northwestern University Core facility and each placenta was given a unique `NUSeq_Core_Facility_ID`. These IDs don't all match the SPAH study unique identifiers in the `pData_SPAH` file, so in the previous script a new column called `Sample_ID` was manually added to the `array_qc_ss` which matches the `pData_SPAH` `ID` column. This file has the column `rgset_colnames` which are the column names of the methylation data object and allows us to match the samples with the methylation data.
3. `delivery_locations` - The delivery locations for the SPAH cohort were missing from the initial meta data file. The file `delivery_location.9.7.23.csv` contains the locations from where the babies were delivered/where the placentas were collected.
4. `drug_info` - In the original meta data file there was information about smoking (cigarettes) and alcohol consumption pre- and during gestation. However, there is more information about the consumption of other drugs during gestation in the file `drug_info.csv`. There is a separate codebook for the `drug_info.csv` file in the same location we read the file in from.
5. `histodata` - the full metadata file for the histology variables for all **575** placentas in the SPAH study that were looked at by a pathologist.


``` r
# 1. Read in the clinical and demographic data
pData_SPAH <- read.csv("//fs/teams/RobinsonLab/ROBLAB6 InfiniumSequenom/EPIC Batch Analysis Code and QC/EPICv1_Batch10_MILLER_code_qc/03_MetaData/SPAH II UBC Data 2.8.24.csv", header = T, fileEncoding="UTF-8-BOM")
dim(pData_SPAH) #605 x 60, includes samples we don't have methylation data for
```

```
## [1] 605  60
```

``` r
# 2. Read in the sample sheet (updated in the first QC script) 
array_qc_ss <- read.csv("//fs/teams/RobinsonLab/ROBLAB6 InfiniumSequenom/EPIC Batch Analysis Code and QC/EPICv1_Batch10_MILLER_code_qc/02_Output/00_Miller_ArrayQC/allplates_allcalcs_samplesheet.csv",header = T) 
dim(array_qc_ss) #520 x 41, all of the samples we have methylation data for
```

```
## [1] 520  41
```

``` r
# 3. Read in the delivery locations data
delivery_locations <- read_excel("//fs/teams/RobinsonLab/ROBLAB6 InfiniumSequenom/EPIC Batch Analysis Code and QC/EPICv1_Batch10_MILLER_code_qc/03_MetaData/UBC_SPAHdata delivery location 9.7.23.xlsx")
```

```
## New names:
## * `` -> `...3`
```

``` r
dim(delivery_locations) #605 x 21 x 4
```

```
## [1] 605   4
```

``` r
head(delivery_locations) #note two columns have no info
```

```
## # A tibble: 6 x 4
##      ID DeliveryLocation ...3  Legend                                           
##   <dbl>            <dbl> <lgl> <chr>                                            
## 1  2001                0 NA    1 = Evanston Hospital and 0 = Other( locations w~
## 2  2071                0 NA    <NA>                                             
## 3  2285                0 NA    <NA>                                             
## 4  2302                0 NA    <NA>                                             
## 5  2365                0 NA    <NA>                                             
## 6  2372                0 NA    <NA>
```

``` r
# 4. Read in the drug and alcohol use data
drug_info <- read.csv("//fs/teams/RobinsonLab/ROBLAB6 InfiniumSequenom/EPIC Batch Analysis Code and QC/EPICv1_Batch10_MILLER_code_qc/03_MetaData/drug_info.csv")
dim(drug_info) #605 x 21
```

```
## [1] 605  21
```

``` r
colnames(drug_info) #note that some columns overlap with the pData_SPAH object
```

```
##  [1] "ID"                  "hp1.v1"              "hp2.v1"             
##  [4] "hp3.v1"              "hp5.v1"              "hp6.v1"             
##  [7] "hp7.v1"              "hp8.v1"              "hp10.v1"            
## [10] "hp4.v1"              "hp9.v1"              "hp11.v1"            
## [13] "hp12.v1"             "hp13.v1"             "hp14.v1"            
## [16] "hp15.v1"             "hp16.v1"             "hp.drinks.prior.v1" 
## [19] "hp.drinks.during.v1" "any_SU_during"       "alc_smoke_before"
```

``` r
# 5. Read in the extra histology data
histodata <- read.csv(here::here("B. Data/A. Cohort Data/histology_cleaned_11.15.23.csv"))
dim(histodata) #575 x 98
```

```
## [1] 575  98
```

### 1.3 Fix Data Sheets

Fix some of the formatting and naming in the data objects. Some of the variable names and classes when loaded from `.csv` files are incorrect. To make it easier in the rest of the script when we need to use these variables for calculations or for analysis we will correct the variable classes and names all at once at the beginning of the script.


``` r
pData_SPAH$ID <- sub("^","SPAH_", pData_SPAH$ID)
head(pData_SPAH$ID)
```

```
## [1] "SPAH_2001" "SPAH_2002" "SPAH_2003" "SPAH_2004" "SPAH_2005" "SPAH_2006"
```

``` r
array_qc_ss$Sample_ID <- sub("^","SPAH_", array_qc_ss$Sample_ID)
head(array_qc_ss$Sample_ID)
```

```
## [1] "SPAH_2578" "SPAH_2379" "SPAH_2159" "SPAH_2062" "SPAH_2170" "SPAH_2295"
```

``` r
delivery_locations$ID <- sub ("^", "SPAH_", delivery_locations$ID)
head(delivery_locations)
```

```
## # A tibble: 6 x 4
##   ID        DeliveryLocation ...3  Legend                                       
##   <chr>                <dbl> <lgl> <chr>                                        
## 1 SPAH_2001                0 NA    1 = Evanston Hospital and 0 = Other( locatio~
## 2 SPAH_2071                0 NA    <NA>                                         
## 3 SPAH_2285                0 NA    <NA>                                         
## 4 SPAH_2302                0 NA    <NA>                                         
## 5 SPAH_2365                0 NA    <NA>                                         
## 6 SPAH_2372                0 NA    <NA>
```

``` r
drug_info$ID <- sub ("^", "SPAH_", drug_info$ID)
head(drug_info)
```

```
##          ID hp1.v1 hp2.v1 hp3.v1 hp5.v1 hp6.v1 hp7.v1 hp8.v1 hp10.v1 hp4.v1
## 1 SPAH_2001      0     NA     NA     NA     NA      1      1       0      0
## 2 SPAH_2002      1      8      2     NA     NA      0      0       0      0
## 3 SPAH_2003      1      2      6     NA     NA      1      8       0      0
## 4 SPAH_2004      0     NA     NA     NA     NA      0      0       0      0
## 5 SPAH_2005      1      1      2     NA     NA      0      0       0      0
## 6 SPAH_2006      1      8      2     NA     NA      0      0       0      0
##   hp9.v1 hp11.v1 hp12.v1 hp13.v1 hp14.v1 hp15.v1 hp16.v1 hp.drinks.prior.v1
## 1      0       0       0       0       0       0       0                  0
## 2      0       0       0       0       0       0       0                  1
## 3      0       0       0       0       0       0       0                 12
## 4      0       0       0       0       0       0       0                  0
## 5      0       0       0       0       0       0       0                  2
## 6      0       0       0       0       0       0       0                  1
##   hp.drinks.during.v1 any_SU_during alc_smoke_before
## 1                   0             0                1
## 2                   0             0                1
## 3                   0             0                2
## 4                   0             0                0
## 5                   0             0                1
## 6                   0             0                1
```

``` r
histodata$ï..ID <- sub("^","SPAH_", histodata$ï..ID)
head(histodata)
```

```
##       ï..ID Placental_wt Placenta_SGA_AGA_LGA Wt_Pctile Bilobed Accessory_lobe
## 1 SPAH_2001          570                  AGA        10       0              0
## 2 SPAH_2002          420                  AGA         3       0              0
## 3 SPAH_2003          496                  AGA         5       0              0
## 4 SPAH_2004          432                  AGA         3       0              0
## 5 SPAH_2005          842                  LGA        11       0              0
## 6 SPAH_2007          570                  AGA         9       0              0
##   CORD_ABNORMALITY Decreased_Coiling Hyper_Coiling Velamentous_Insertion
## 1                0                 0             0                     0
## 2                0                 0             0                     0
## 3                0                 0             0                     0
## 4                0                 0             0                     0
## 5                0                 0             0                     0
## 6                1                 0             1                     0
##   Marginal_Cord_Insertion Single_Umbilical_artery Cord_knot Abnormal_length
## 1                       0                       0         0               0
## 2                       0                       0         0               0
## 3                       0                       0         0               0
## 4                       0                       0         0               0
## 5                       0                       0         0               0
## 6                       0                       0         0               0
##   Furcate_Insertion ACUTE_INFLAMMATION Maternal_inflamm Maternal_stage
## 1                 0                  0                0              0
## 2                 0                  1                1              2
## 3                 0                  0                0              0
## 4                 0                  0                0              0
## 5                 0                  1                1              1
## 6                 0                  1                1              1
##   Maternal_stage_3cat Maternal_high_stage Fetal_inflamm Fetal_stage
## 1                   0                   0             0           0
## 2                   2                   1             0           0
## 3                   0                   0             0           0
## 4                   0                   0             0           0
## 5                   1                   0             1           1
## 6                   1                   0             1           2
##   Fetal_stage_3cat Fetal_high_stage Funisitis Granular_fetal Periph_funisitis
## 1                0                0         0              0                0
## 2                0                0         0              0                0
## 3                0                0         0              0                0
## 4                0                0         0              0                0
## 5                1                0         0              1                0
## 6                2                1         1              4                0
##   acute_inflammation_3cat Villous_edema CHRONIC_INFLAMMATION
## 1                       0             0                    1
## 2                       2             0                    1
## 3                       0             0                    1
## 4                       0             0                    1
## 5                       1             0                    0
## 6                       2             0                    1
##   Chronic_inflammation_sum Chronic_villitis Chronic_villitis_di
## 1                        1                0                   0
## 2                        1                0                   0
## 3                        2                0                   0
## 4                        2                0                   0
## 5                        0                0                   0
## 6                        1                0                   0
##   Chronic_basal_villitis CDPC Chronic_chorioamnio Chronic_chorionitis
## 1                      1    0                   0                   0
## 2                      0    1                   0                   0
## 3                      0    1                   0                   0
## 4                      1    1                   0                   0
## 5                      0    0                   0                   0
## 6                      0    1                   0                   0
##   Chronic_marg_decid Chronic_decid_perivasc Chronic_intervillositis
## 1                  0                      0                       0
## 2                  0                      0                       0
## 3                  1                      0                       0
## 4                  0                      0                       0
## 5                  0                      0                       0
## 6                  0                      0                       0
##   Eosin_Tcell_vasculitis chronic_inflammation_comp chronic_inflammation_3cat
## 1                      0                         1                         1
## 2                      0                         1                         1
## 3                      0                         2                         2
## 4                      0                         2                         2
## 5                      0                         0                         0
## 6                      0                         1                         1
##   FETAL_VASC_PATH Fetal_vasc_thrombi FVT_Multifocal Thrombi_chorionic_vessel
## 1               1                  0              0                        0
## 2               0                  0              0                        0
## 3               0                  0              0                        0
## 4               1                  1              1                        0
## 5               1                  0              0                        0
## 6               1                  0              0                        0
##   TCV_Multifocal Thrombi_stem_villous_vessel TSV_Multifocal
## 1              0                           0              0
## 2              0                           0              0
## 3              0                           0              0
## 4              0                           1              1
## 5              0                           0              0
## 6              0                           0              0
##   Thrombi_umbilical_vessel TUV_Multifocal Avascular_villi AV_Multifocal
## 1                        0              0               1             1
## 2                        0              0               0             0
## 3                        0              0               0             0
## 4                        0              0               1             1
## 5                        0              0               1             0
## 6                        0              0               1             0
##   Fetal_thrombotic_vasculopathy Myonecrosis Fetal_vasc_path_3cat Vessel_path
## 1                             1           0                    2           1
## 2                             0           0                    0           0
## 3                             0           0                    0           0
## 4                             0           0                    1           0
## 5                             0           0                    1           0
## 6                             0           0                    1           0
##   FN_AA Musc_BP_arterioles MHMA Basal_Dec_vasc_thromb Villous_changes Infarct
## 1     0                  1    0                     0               0       0
## 2     0                  0    0                     0               0       0
## 3     0                  0    0                     0               0       0
## 4     0                  0    0                     0               0       0
## 5     0                  0    0                     0               0       0
## 6     0                  0    0                     0               1       1
##   Infarct_multiple Increased_SK Villous_AG Increased_PVF DVH_STV MVM_lesion_sum
## 1                0            0          0             0       0              1
## 2                0            0          0             0       0              0
## 3                0            0          0             0       0              0
## 4                0            0          0             0       0              0
## 5                0            0          0             0       0              0
## 6                0            0          0             0       0              1
##   MVM_score mvm_score_3cat mvm_di EV_OF_ABRUPTION
## 1         1              0      0               0
## 2         0              0      0               0
## 3         0              0      0               0
## 4         0              0      0               0
## 5         0              0      0               0
## 6         1              0      0               0
##   BP_acute_intradecid_hemmorhage RP_blood_hematoma RP_blood_hematoma_hemosid
## 1                              0                 0                         0
## 2                              0                 0                         0
## 3                              0                 0                         0
## 4                              0                 0                         0
## 5                              0                 0                         0
## 6                              0                 0                         0
##   RP_blood_hematoma_infarct Focal_acute_IV_hem OLD_MBP_BLEED Memhem
## 1                         0                  0             0      0
## 2                         0                  0             0      0
## 3                         0                  0             0      0
## 4                         0                  0             0      0
## 5                         0                  0             0      0
## 6                         0                  0             0      0
##   BP_hemosid_depo Chorioam_hemosiderosis SC_IV_hem_thrombi Amnion_Nodosum
## 1               0                      0                 1              0
## 2               0                      0                 1              0
## 3               0                      0                 0              0
## 4               0                      0                 0              0
## 5               0                      0                 0              0
## 6               0                      0                 0              0
##   Meconium_histiocytosis Villous_dysmaturity Villous_chorangiosis
## 1                      0                   0                    0
## 2                      0                   0                    0
## 3                      0                   0                    0
## 4                      0                   0                    0
## 5                      0                   0                    1
## 6                      1                   0                    0
##   Massive_PVF_deposition BPmyo other_chronic_path Phenotype Phenotype_Group
## 1                      0     1                  0        cF              3b
## 2                      0     1                  0        Ac              6b
## 3                      0     0                  0         C              4a
## 4                      0     0                  0        Cf              4a
## 5                      0     1                  1        af              6a
## 6                      0     0                  0       Acf              7c
##                                                                                                            phenotype_label
## 1     High-grade fetal vascular pathology with low-grade chronic inflammation and/or low-grade maternal vascular pathology
## 2                             High-grade acute inflammation and low-grade fetal vascular pathology or chronic inflammation
## 3                                     High-grade chronic inflammation (with or without low-grade fetal vascular pathology)
## 4                                     High-grade chronic inflammation (with or without low-grade fetal vascular pathology)
## 5                Low-grade fetal vascular pathology or chronic inflammation (with or without low-grade acute inflammation)
## 6 High-grade acute inflammation and low-grade maternal vascular pathology or any two or more low-grade chronic pathologies
##          Comments Covid_supplement SPAH_supplement Placenta_SGAvAGAvLGA
## 1              cF                1               1                    2
## 2                                1               0                    2
## 3                                1               1                    2
## 4                                1               0                    2
## 5              af                1               1                    3
## 6 infarct, 0.6 cm                1               0                    2
```

``` r
delivery_locations <- delivery_locations %>% dplyr::select(!c("...3", "Legend"))
dim(delivery_locations) 
```

```
## [1] 605   2
```

``` r
array_qc_ss <- array_qc_ss %>%
  dplyr::rename(Prob_European = Prob_Caucasian,
                Cytotrophoblast = Trophoblasts,
                Predicted_genetic_ancestry = Predicted_ethnicity)%>%
  mutate(Sentrix_ID = as.character(Sentrix_ID))

array_qc_ss[array_qc_ss == "Caucasian"] <-"European"
array_qc_ss[array_qc_ss == "Female"] <-"XX"
array_qc_ss[array_qc_ss == "Male"] <-"XY"

table(grepl("#NULL!",pData_SPAH)) 
```

```
## 
## FALSE  TRUE 
##     9    51
```

``` r
pData_SPAH <- replace(pData_SPAH, pData_SPAH == "#NULL!", NA)
table(grepl("#NULL!",pData_SPAH)) 
```

```
## 
## FALSE 
##    60
```

``` r
pData_SPAH <- pData_SPAH %>% 
  mutate(ID = as.factor(ID)) %>% 
  mutate_if(is.character, as.numeric) %>%
  mutate(ID = as.character(ID)) %>% 
  mutate(pcr_baby_sex = as.character(pcr_baby_sex),
         pcr_mode_conception = as.factor(pcr_mode_conception),
         diabetes_pregestational = as.factor(diabetes_pregestational),
         parity = as.factor(parity),
         pcr_gravida = as.factor(pcr_gravida),
         ACUTE_INFLAMMATION = as.factor(ACUTE_INFLAMMATION),
         CHRONIC_INFLAMMATION = as.factor(CHRONIC_INFLAMMATION),
         Chronic_villitis = as.factor(Chronic_villitis),
         FETAL_VASC_PATH = as.factor(FETAL_VASC_PATH),
         mvm_di =as.factor(mvm_di),
         Fetal_vasc_path_3cat = as.factor(Fetal_vasc_path_3cat),
         MVM_lesion_sum =as.factor(MVM_lesion_sum),
         mvm_score_3cat = as.factor(mvm_score_3cat),
         new_pcr_hypertension_gestation = as.factor(new_pcr_hypertension_gestation),
         pcr_pre_eclampsia = as.factor(pcr_pre_eclampsia),
         pcr_peripart_infect_chorio = as.factor(pcr_peripart_infect_chorio))

sapply(pData_SPAH,class) 
```

```
##                             ID                            Age 
##                    "character"                      "numeric" 
##             Race_AmIndALNative                     Race_Asian 
##                      "integer"                      "integer" 
##                     Race_Black            Race_HINativePacIsl 
##                      "integer"                      "integer" 
##                     Race_white                     Race_Other 
##                      "integer"                      "integer" 
##               Ethnicity_HispYN             ResourcesComposite 
##                      "integer"                      "numeric" 
##              PrestigeComposite          DisadvantageComposite 
##                      "numeric"                      "numeric" 
##                         EvntCT                         DiffCT 
##                      "numeric"                      "numeric" 
##                         EvntTH                         DiffTH 
##                      "numeric"                      "numeric" 
##                         BMI.v1                 prepreg_weight 
##                      "numeric"                      "numeric" 
##                    prepreg_BMI                    epds.sum.v1 
##                      "numeric"                      "numeric" 
##                         hp7.v1                         hp9.v1 
##                      "numeric"                      "numeric" 
##             hp.drinks.prior.v1            hp.drinks.during.v1 
##                      "numeric"                      "numeric" 
##                  ga_continuous                   pcr_baby_sex 
##                      "numeric"                    "character" 
##                    pcr_gravida                         parity 
##                       "factor"                       "factor" 
##            pcr_mode_conception              pcr_pre_eclampsia 
##                       "factor"                       "factor" 
## new_pcr_hypertension_gestation                            HDP 
##                       "factor"                      "numeric" 
##                       pcr_sptd                  pcr_sptd_sptl 
##                      "numeric"                      "numeric" 
##                 pcr_sptd_pprom               pcr_diabetes_any 
##                      "numeric"                      "numeric" 
##        diabetes_pregestational           diabetes_gestational 
##                       "factor"                      "numeric" 
##                   vag_delivery                      C_section 
##                      "numeric"                      "numeric" 
##             operative_delivery           pcr_postp_hemorrhage 
##                      "numeric"                      "numeric" 
##       pcr_peripartum_infection     pcr_peripart_infect_chorio 
##                      "numeric"                       "factor" 
##                   pcr_birth_wt               pcr_birth_length 
##                      "numeric"                      "numeric" 
##                  pcr_head_circ                   Placental_wt 
##                      "numeric"                      "numeric" 
##           Placenta_SGAvAGAvLGA                      Wt_Pctile 
##                      "numeric"                      "numeric" 
##             ACUTE_INFLAMMATION           CHRONIC_INFLAMMATION 
##                       "factor"                       "factor" 
##       Chronic_inflammation_sum               Chronic_villitis 
##                      "numeric"                       "factor" 
##      chronic_inflammation_3cat                FETAL_VASC_PATH 
##                      "numeric"                       "factor" 
##           Fetal_vasc_path_3cat                 MVM_lesion_sum 
##                       "factor"                       "factor" 
##                 mvm_score_3cat                         mvm_di 
##                       "factor"                       "factor"
```

### 1.4 Remove Poor Quality Samples

Based on the initial Array QC script, several samples should be removed:

  * All of the test plate samples - they had irregular beta value densities. These samples had duplicates that were run on other plates that looked more normal.
  * Sample `SPAH_2175` - failed most of the array QC tests
  * Sample `SPAH_2374` - failed the bead count threshold (more than 1% of probes had bead count <3)


``` r
#get the list of samples that should be removed based on initial QC of the cohort
remove_samples_list <- read.csv("//fs/teams/RobinsonLab/ROBLAB6 InfiniumSequenom/EPIC Batch Analysis Code and QC/EPICv1_Batch10_MILLER_code_qc/02_Output/00_Miller_ArrayQC/sample_removal_list.csv", header = T)
dim(remove_samples_list) #10
```

```
## [1] 10  1
```

``` r
remove_samples_list$Sample_ID <- sub("^","SPAH_", remove_samples_list$Sample_ID)
remove_samples_list$Sample_ID
```

```
##  [1] "SPAH_2175" "SPAH_2374" "SPAH_2002" "SPAH_2003" "SPAH_2004" "SPAH_2005"
##  [7] "SPAH_2007" "SPAH_2008" "SPAH_2009" "SPAH_2010"
```

``` r
table(array_qc_ss$Sample_ID %in% remove_samples_list$Sample_ID)
```

```
## 
## FALSE  TRUE 
##   502    18
```

``` r
array_qc_ss_remove <- array_qc_ss %>% filter(array_qc_ss$Sample_ID %in% remove_samples_list$Sample_ID)
dim(array_qc_ss_remove) 
```

```
## [1] 18 41
```

``` r
array_qc_ss_remove <- array_qc_ss_remove %>% arrange(Sentrix_ID)
array_qc_ss_remove <- array_qc_ss_remove %>% distinct(Sample_ID, .keep_all = TRUE) 
dim(array_qc_ss_remove) 
```

```
## [1] 10 41
```

``` r
array_qc_ss <- array_qc_ss %>% filter(!rgset_colnames %in% array_qc_ss_remove$rgset_colnames) 
dim(array_qc_ss) #510
```

```
## [1] 510  41
```

### 1.5 Combine Data Sheets


``` r
pData_SPAH <- pData_SPAH %>% filter(pData_SPAH$ID %in% array_qc_ss$Sample_ID)
dim(pData_SPAH) 
```

```
## [1] 507  60
```

``` r
metadata_SPAH <- pData_SPAH %>% left_join(delivery_locations, by="ID")
metadata_SPAH <- metadata_SPAH %>% dplyr::select(-c(hp.drinks.prior.v1, hp.drinks.during.v1, hp7.v1, hp9.v1))
metadata_SPAH <- metadata_SPAH %>% left_join(drug_info, by = "ID")
metadata_SPAH <- metadata_SPAH %>% left_join(histodata %>%
  dplyr::select(ï..ID, Maternal_inflamm, Maternal_stage_3cat, Fetal_inflamm, Fetal_stage_3cat, acute_inflammation_3cat, Villous_dysmaturity, Placenta_SGA_AGA_LGA), by = c("ID" = "ï..ID")
)
dim(metadata_SPAH)
```

```
## [1] 507  84
```

## 2.0 Calculate New Variables

### 2.1 Harmonize Variables

#### 2.1.1 Mode of Delivery


``` r
metadata_SPAH <- metadata_SPAH %>% mutate(mode_delivery = as.factor(case_when(
  vag_delivery == 1~ "vaginal delivery",
  C_section == 1~ "C section",
  operative_delivery == 1~ "operative delivery")))

table(metadata_SPAH$mode_delivery)
```

```
## 
##          C section operative delivery   vaginal delivery 
##                152                 13                338
```

#### 2.1.2 Preeclampsia Subgroup


``` r
#Classify Preeclampsia as EOPE or LOPE and other placentas as 
#term (nTB), preterm (nPTB), or very preterm (nvPTB)
metadata_SPAH <- metadata_SPAH %>% mutate(PE_GA_state = factor(case_when(
  is.na(ga_continuous) ~ "NA",
  is.na(pcr_pre_eclampsia) ~ NA,
  pcr_pre_eclampsia == 1 & ga_continuous < 34 ~ "EOPE",
  pcr_pre_eclampsia == 1 & ga_continuous >= 34 ~ "LOPE",
  #pcr_pre_eclampsia == 0 & ga_continuous < 34 ~ "nvPTB",
  pcr_pre_eclampsia == 0 & ga_continuous  < 37 ~ "nPTB",
  pcr_pre_eclampsia == 0 & ga_continuous >= 37 ~ "nTB"),
  levels = c("EOPE", "LOPE","nPTB","nTB")))

table(metadata_SPAH$PE_GA_state)
```

```
## 
## EOPE LOPE nPTB  nTB 
##    7   46   40  412
```

#### 2.1.3 Pathology Class


``` r
metadata_SPAH <- metadata_SPAH %>%
  mutate(AI = case_when(
      ACUTE_INFLAMMATION == 1 ~"AI",
      ACUTE_INFLAMMATION == 0 ~ NA),
    CI = case_when(
      CHRONIC_INFLAMMATION == 1 ~"CI",
      CHRONIC_INFLAMMATION == 0 ~ NA),
    FVM = case_when(
      FETAL_VASC_PATH == 1 ~"FVM",
      FETAL_VASC_PATH == 0 ~ NA),
    MVM = case_when(
      mvm_di == 1 ~"MVM",
      mvm_di == 0 ~ NA)) %>%
  unite("HarmonizedPath", "AI":"MVM", sep = ", ", na.rm = TRUE, remove = TRUE)

metadata_SPAH$HarmonizedPath[metadata_SPAH$HarmonizedPath == ""] <- "No Pathology" 

metadata_SPAH[is.na(metadata_SPAH$ACUTE_INFLAMMATION), "HarmonizedPath"] <- "NA"

# reorder the levels to have the exclusive classes first
metadata_SPAH <- metadata_SPAH %>% mutate(HarmonizedPath = factor(.$HarmonizedPath, levels=c("No Pathology",
  "AI",
  "CI",
  "FVM",
  "MVM",
  "AI, CI",
  "AI, CI, FVM",
  "AI, CI, MVM",
  "AI, FVM",
  "AI, FVM, MVM",
  "AI, MVM",
  "CI, FVM",
  "CI, FVM, MVM",
  "CI, MVM",
  "FVM, MVM",
  "AI, CI, FVM, MVM")))

table(metadata_SPAH$HarmonizedPath)
```

```
## 
##     No Pathology               AI               CI              FVM 
##               46               73               55               15 
##              MVM           AI, CI      AI, CI, FVM      AI, CI, MVM 
##               23               54               47               30 
##          AI, FVM     AI, FVM, MVM          AI, MVM          CI, FVM 
##               30               10               19               33 
##     CI, FVM, MVM          CI, MVM         FVM, MVM AI, CI, FVM, MVM 
##               14               29                9               11
```

#### 2.1.4 Maternal Self-Reported Race


``` r
metadata_SPAH <- metadata_SPAH %>%
  mutate(AmIndALNative = case_when(
    Race_AmIndALNative == 1 ~"American Indian or Alaskan Native",
    Race_AmIndALNative == 0 ~ NA),
    Asian = case_when(
      Race_Asian == 1 ~"Asian",
      Race_Asian == 0 ~ NA),
    Black = case_when(
      Race_Black == 1 ~"Black/African American",
      Race_Black == 0 ~ NA),
    HINativePacIsl = case_when(
      Race_HINativePacIsl == 1 ~"Hawaiian Native or Pacific Islander",
      Race_HINativePacIsl == 0 ~ NA),
    White = case_when(
      Race_white == 1 ~"White",
      Race_white == 0 ~ NA),
    Other = case_when(
      Race_Other == 1 ~"Other Race",
      Race_Other == 0 ~ NA))%>% 
  unite("HarmonizedRace", "AmIndALNative":"Other", sep = ", ", na.rm = TRUE, remove = TRUE)

table(metadata_SPAH$HarmonizedRace)
```

```
## 
##                                                   American Indian or Alaskan Native 
##                                                                                   1 
## American Indian or Alaskan Native, Asian, Black/African American, White, Other Race 
##                                                                                   1 
##                    American Indian or Alaskan Native, Black/African American, White 
##                                                                                   1 
##                                            American Indian or Alaskan Native, White 
##                                                                                   3 
##                                                                               Asian 
##                                                                                  49 
##      Asian, Black/African American, Hawaiian Native or Pacific Islander, Other Race 
##                                                                                   1 
##                                                                        Asian, White 
##                                                                                   5 
##                                                              Black/African American 
##                                                                                  76 
##                                                  Black/African American, Other Race 
##                                                                                   2 
##                                                       Black/African American, White 
##                                                                                   7 
##                                                 Hawaiian Native or Pacific Islander 
##                                                                                   4 
##                                                                          Other Race 
##                                                                                  43 
##                                                                               White 
##                                                                                 314
```

### 2.2 Placental Efficiency

In this section we calculate the placental efficiency residuals from the placental weight and birth weigh [Christians et al. (2018)](https://www.sciencedirect.com/science/article/pii/S014340041830170X?via%3Dihub). 


``` r
metadata_SPAH_pe <- metadata_SPAH %>%
  filter(!is.na(pcr_birth_wt))%>%
  filter(!is.na(Placental_wt)) %>%
  mutate(Placental_efficiency = residuals(lm(pcr_birth_wt ~
       Placental_wt+ 
       ga_continuous + 
       pcr_baby_sex, 
       pData_SPAH)))

summary(metadata_SPAH_pe$Placental_efficiency)
```

```
##     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
## -1051.53  -218.75   -11.05     0.00   200.60  1653.30
```

``` r
metadata_SPAH_pe <- metadata_SPAH_pe %>%
  dplyr::select(c(Placental_efficiency, ID))

metadata_SPAH <- metadata_SPAH %>% 
  left_join(metadata_SPAH_pe, by= "ID")
```

### 2.3 Birth Weight SD

Here I calculate the birthweight standard deviation and SGA/AGA/LGA calls using the [Fenton 2013 birth chart reference](https://ucalgary.ca/resource/preterm-growth-chart/preterm-growth-chart).

The file `SizeatbirthcalculatorFenton2013.xlsx` was downloaded from the Univeristy of Calgary's site on July 29th, 2024. I used this spreadsheet's calculations to create the code below and the dataset I read in as a reference. 


``` r
bw_scores <- read.csv(here::here("B. Data/B. Reference Data/BW_SD_data_Fenton_2013.csv"), header = T, fileEncoding="UTF-8-BOM")
colnames(bw_scores)
```

```
##  [1] "GA_weeks_completed" "L_wt"               "M_wt"              
##  [4] "S_wt"               "Per_10"             "Per_90"            
##  [7] "Sex"                "S3Pos"              "S23Pos"            
## [10] "S3Neg"              "S23Neg"
```

``` r
#round the GA to # of weeks completed, to match the Fenton calculator
metadata_SPAH <- metadata_SPAH %>% mutate (GA_wk_completed = floor(ga_continuous))

bw_sd_Fenton <- function (x, y, z){ #x = GA, y = BW, z = Sex
  if (is.na(x)| is.na (y)){
    BW_sd <- data.frame(
      Bw_call = NA,
      Bw_sd_Fenton = NA)
    
  } else{
    if(z == 2){ 
      temp <- bw_scores %>% filter (GA_weeks_completed == x & Sex ==1) #1 is male in Fenton data
    }else{
      temp <- bw_scores %>% filter (GA_weeks_completed == x & Sex ==0)} #0 is female in Fenton data
    
    if(((y/temp$M_wt)^(temp$L_wt)-1)/(temp$L_wt*temp$S_wt)> 3){
      bw_sd <- 3+ (y-temp$S3Pos)/temp$S23Pos
    } else if(((y/temp$M_wt)^(temp$L_wt)-1)/(temp$L_wt*temp$S_wt)< -3){
      bw_sd <- -3+ (y-temp$S3Neg)/temp$S23Neg
    } else {
      bw_sd <- (((y/temp$M_wt)^(temp$L_wt)-1)/(temp$L_wt*temp$S_wt))}
    
    if(y < temp$Per_10){ #Per_10 = 10th percentile weight
      bw_call <- "SGA"
    } else if(y > temp$Per_90){ #Per_90 = 90th percentile weight
      bw_call <- "LGA"
    } else{
      bw_call <- "AGA"}
    
    BW_sd <- data.frame(
      Bw_call = bw_call,
      Bw_sd_Fenton = bw_sd)
  }
  return(BW_sd)
}

BWSD <- data.frame(
  Bw_call = c(),
  Bw_sd_Fenton = c())

for (i in 1:nrow(metadata_SPAH)){
  BW_temp <- bw_sd_Fenton (metadata_SPAH[i, "GA_wk_completed"], metadata_SPAH[i, "pcr_birth_wt"], metadata_SPAH[i, "pcr_baby_sex"])
  BWSD<- BWSD %>% rbind(BW_temp)
}

metadata_SPAH <- metadata_SPAH %>%
  cbind(BWSD) %>%
  mutate(Bw_call = factor(Bw_call,
                          levels = c("SGA","AGA","LGA"))) 

table(metadata_SPAH$Bw_call)
```

```
## 
## SGA AGA LGA 
##  23 414  67
```


### 2.4 Epiphenotyping Variables


``` r
rgset_raw <- readRDS("//fs/teams/RobinsonLab/ROBLAB6 InfiniumSequenom/EPIC Batch Analysis Code and QC/EPICv1_Batch10_MILLER_code_qc/02_Output/00_Miller_ArrayQC/rgset_all_samples.rds")
dim(rgset_raw) #1051943 x 520
```

```
## [1] 1051943     520
```

``` r
rgset_raw <- rgset_raw[,array_qc_ss$rgset_colnames] 
dim(rgset_raw) 
```

```
## [1] 1051943     510
```

``` r
mset_noob <- preprocessNoob(rgset_raw) #noob normalize and background correct
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
betas_bmiq <- wateRmelon::BMIQ(mset_noob) #bmiq normalize
```

### 2.4.1 Update Cell-type Estimates

Calculate the cell-type estimates using the robust partial correlations (RPC, constrained) method.


``` r
data(plCellCpGsThird) 

#use the robust partial correlations method
epidish_RPC <- epidish(
    beta.m = betas_bmiq[rownames(plCellCpGsThird), ],
    ref.m = plCellCpGsThird,
    method = "RPC"
)

#rename the old variables
array_qc_ss <- array_qc_ss %>%
  dplyr::rename(
    Stromal.old = Stromal,
    Hofbauer.old = Hofbauer,
    Endothelial.old = Endothelial,
    nRBC.old = nRBC,
    Syncytiotrophoblast.old = Syncytiotrophoblast,
    Cytotrophoblast.old = Cytotrophoblast
  ) %>%
  left_join(epidish_RPC$estF %>% 
              as.data.frame() %>%
              rownames_to_column(., var="rgset_colnames"))%>%
  mutate(CTB_STB = Trophoblasts/Syncytiotrophoblast) #create the CTB/STB variable
```

```
## Joining with `by = join_by(rgset_colnames)`
```



### 2.4.2 eoPRED scores

Here I will calculate the eoPRED scores using Icíar Fernández Boyano's function.
 


``` r
#' @title predictPreeclampsia
#'
#' @description Uses 45 CpGs to predict early preeclampsia (PE delivered before or at 34 weeks of gestation)
#' on placental DNA methylation microarray data.
#'
#' @details Assigns the class labels "early-PE" or "normotensive" to each sample
#' and returns a class probability.
#'
#' # It is recommended that users apply beta-mixture quantile normalization (BMIQ) to their data
#' prior to prediction. This was the normalization method used on the training data.
#'
#' @param betas matrix or array of methylation values on the beta scale (0, 1),
#' where the variables are arranged in rows, and samples in columns.
#' 
#' @param ... feeds into outersect function
#'
#' @return produces a list with components detailed in the `mixOmics::predict` R documentation
#'
#' @examples
#'
#' # To predict early preeclampsia on 450k/850k samples
#'
#' # Load data
#' # data(peBetas)
#' # predictPreeclampsia(peBetas, dist = "max.dist")
#' 
#' @import mixOmics
#' @import ExperimentHub
#' @export predictPreeclampsia
#'

predictPreeclampsia <- function(betas, ...){
  
  # read in data to generate model
  eh = ExperimentHub()
  mod = eh[['EH8090']]
  trainCpGs = colnames(mod$X)
  
  # check that there are no NAs in the predictors (or if there are, how many)
  peCpGs = mixOmics::selectVar(mod)$name
  pp <- intersect(colnames(betas), peCpGs)
  
  if(length(pp) < length(peCpGs)){
    stop(paste(
      "Only", length(pp), "out of 45 predictive CpGs present. All 45 predictive CpGs are needed to run the function."
    ))
  } else {
    message(paste(length(pp), "of 45 predictive CpGs present."))
    message("BMIQ normalization is recommended for best results. If choosing other method, it is recommended to compare results to predictions on BMIQ normalized data.")
  }
  
  # set up data for prediction
  
  # if input data is missing any of the cpgs present in the training data, this function
  # adds the ones that are missing as NAs
  # necessary for `mixOmics::predict` to work
  
  outersect = function(x, y) {
    sort(c(x[!x%in%y],
           y[!y%in%x]))
  }
  
  if(inherits(betas, 'matrix')){
  } else if (inherits(betas, 'array')) {
  } else {
    
    # throw an error
    print(paste0("Input data must be a matrix or an array"))
  }
  
  betasSubset <- betas[,colnames(betas) %in% trainCpGs]
  
  # order
  betasSubset <- betasSubset[drop=FALSE,, trainCpGs]
  
  if(all(colnames(betasSubset) == trainCpGs) == FALSE){
    stop()
  } else
    
    # predict
    out <- mixOmics:::predict.mixo_spls(mod, betasSubset)
  
  # get class probabilities
  CP <- out$predict[,,1]
  CP <- t(apply(as.matrix(CP), 1, function(data) exp(data)/sum(exp(data))))
  CP <- as.data.frame(CP) %>% tibble::rownames_to_column("Sample_ID")
  CP$PE_Status <- CP$comp1
  CP <- CP %>%
    dplyr::mutate(PE_Status = dplyr::case_when(EOPE > 0.55 ~ "EOPE",
                                                EOPE < 0.55 ~ "Normotensive"))
  
  return(CP)
}
```



``` r
#eoPred was designed to work on transposed, unfiltered, and BMIQ normalized beta values
eopred_cpgs <- read.csv(here::here("B. Data/B. Reference Data/eopred_cpgs.csv"))

table(eopred_cpgs$ï..probeID %in% rownames(betas_bmiq)) 
```

```
## 
## TRUE 
##   45
```

``` r
t_betas_bmiq <- t(betas_bmiq)

eoPRED_results <- predictPreeclampsia(t_betas_bmiq)
```

```
## Bioconductor version 3.17 (BiocManager 1.30.22), R 4.3.1 (2023-06-16 ucrt)
```

```
## Installing package(s) 'eoPredData'
```

```
## Warning: package 'eoPredData' is not available for Bioconductor version '3.17'
## 
## A version of this package for your version of R might be available elsewhere,
## see the ideas at
## https://cran.r-project.org/doc/manuals/r-patched/R-admin.html#Installing-packages
```

```
## Warning: unable to access index for repository https://bioconductor.org/packages/3.17/data/annotation/bin/windows/contrib/4.3:
##   cannot open URL 'https://bioconductor.org/packages/3.17/data/annotation/bin/windows/contrib/4.3/PACKAGES'
```

```
## Warning: unable to access index for repository https://bioconductor.org/packages/3.17/data/experiment/bin/windows/contrib/4.3:
##   cannot open URL 'https://bioconductor.org/packages/3.17/data/experiment/bin/windows/contrib/4.3/PACKAGES'
```

```
## Warning: unable to access index for repository https://bioconductor.org/packages/3.17/workflows/bin/windows/contrib/4.3:
##   cannot open URL 'https://bioconductor.org/packages/3.17/workflows/bin/windows/contrib/4.3/PACKAGES'
```

```
## Installation paths not writeable, unable to update packages
##   path: C:/Program Files/R/R-4.3.1/library
##   packages:
##     abind, ape, askpass, backports, base64, base64enc, bayesm, BH, BiasedUrn,
##     BiocManager, bit, bit64, bitops, blob, boot, brew, brio, broom, cachem,
##     Cairo, callr, car, carData, caret, caTools, class, clock, clue, cluster,
##     codetools, compositions, corrplot, crayon, credentials, crosstalk, curl,
##     DBI, dbplyr, DelayedMatrixStats, deldir, DEoptimR, Deriv, desc, DescTools,
##     diffobj, doBy, doRNG, dotCall64, downlit, DT, dtplyr, dunn.test, e1071,
##     emmeans, estimability, expm, extrafont, extrafontdb, FactoMineR, farver,
##     fastcluster, fastICA, fastmap, ff, fields, filelock, fontawesome, forcats,
##     foreign, fs, FSA, futile.logger, future, future.apply, gargle, gdata,
##     generics, GenomeInfoDb, GenomicFeatures, GenomicRanges, gert, GetoptLong,
##     ggdendro, ggplotify, ggrepel, ggsci, ggstats, ggupset, ggvenn, gh, gld,
##     glmnet, GlobalOptions, globals, glue, gmodels, googledrive, googlesheets4,
##     gower, gplots, gtable, gtools, GWASExactHW, hardhat, haven, here, hexbin,
##     highr, hms, htmlTable, htmlwidgets, httpuv, httr2, igraph, interp, ipred,
##     irlba, isoband, janitor, jpeg, jsonlite, KEGGREST, KernSmooth, khroma,
##     labeling, labelled, later, latex2exp, lattice, latticeExtra, lava, leaps,
##     listenv, lme4, lmodel2, lmom, locfit, logistf, magrittr, maps, markdown,
##     MatrixModels, matrixStats, mclust, mediation, mice, mime, miniUI, minqa,
##     minty, modeltools, munsell, mvtnorm, nleqslv, nlme, nloptr, nnet, nor1mix,
##     openssl, openxlsx, operator.tools, ordinal, parallelly, patchwork,
##     pbkrtest, permute, pheatmap, pillar, pkgbuild, pkgdown, pkgload, plotly,
##     plotrix, plyr, prettyunits, pROC, processx, prodlim, profvis, progress,
##     progressr, promises, proxy, ps, purrr, quantreg, R.oo, R.utils, R6, ragg,
##     rappdirs, RcppArmadillo, RcppEigen, RCurl, readODS, readr, readxl, recipes,
##     REDCapR, rematch, remotes, reprex, reshape, restfulr, Rhdf5lib, rJava,
##     rjson, RMariaDB, robustbase, roxygen2, rpart, rprojroot, rsample, RSQLite,
##     Rttf2pt1, rversions, rvest, S4Arrays, S4Vectors, sandwich, sass, scales,
##     scrime, selectr, sessioninfo, shape, shiny, showtext, slider, snakecase,
##     spam, SparseM, spatial, statmod, strawr, stringi, stringr, survey,
##     survival, svglite, sys, sysfonts, tensorA, testthat, tibble, tidyr,
##     tidyREDCap, tidyselect, timeDate, tzdb, ucminf, usethis, uuid, viridis,
##     viridisLite, waldo, warp, wesanderson, withr, writexl, xlsxjars, XML,
##     xopen, xts, yulab.utils, zip, zoo
```

```
## Old packages: 'ADMM', 'bayestestR', 'broom.helpers', 'bslib', 'cards', 'cardx',
##   'checkmate', 'circlize', 'cli', 'colorspace', 'commonmark', 'cowplot',
##   'cpp11', 'data.table', 'datawizard', 'dbscan', 'dendextend', 'devtools',
##   'digest', 'dplyr', 'evaluate', 'fansi', 'geomtextpath', 'GGally', 'ggExtra',
##   'ggfortify', 'ggplot2', 'ggpmisc', 'ggpp', 'ggpubr', 'ggridges', 'ggtern',
##   'ggthemes', 'gprofiler2', 'gt', 'gtsummary', 'Hmisc', 'htmltools', 'insight',
##   'jtools', 'knitr', 'labdsv', 'lifecycle', 'litedown', 'lubridate', 'maotai',
##   'mclustcomp', 'multcomp', 'parameters', 'pixmap', 'pracma', 'qtl',
##   'rbibutils', 'Rcpp', 'Rdimtools', 'Rdpack', 'reactable', 'reformulas',
##   'reshape2', 'rgl', 'rlang', 'rmarkdown', 'Rmpfr', 'rstatix', 'rstudioapi',
##   'rsvg', 'S7', 'segmented', 'shapes', 'shinyjs', 'srvyr', 'systemfonts',
##   'textshaping', 'TH.data', 'timechange', 'tinytex', 'TMB', 'unmarked', 'utf8',
##   'V8', 'vctrs', 'VennDiagram', 'VGAM', 'vroom', 'WGCNA', 'xfun', 'xml2',
##   'yaml'
```

```
## loading from cache
```

```
## 45 of 45 predictive CpGs present.
```

```
## BMIQ normalization is recommended for best results. If choosing other method, it is recommended to compare results to predictions on BMIQ normalized data.
```

``` r
array_qc_ss <- array_qc_ss %>% left_join(eoPRED_results, by= c("rgset_colnames"="Sample_ID"))
```

### 2.4.3 Epigenetic Age Acceleration

Calculated as the residuals of the linear regression between clinical gestational age and PlaNET gestational age adjusting for sex (`extAccel`) and cell estimates (`intAccel`).


``` r
table(metadata_SPAH$ID %in% array_qc_ss$Sample_ID)
```

```
## 
## TRUE 
##  507
```

``` r
metadata_eaa <- array_qc_ss %>%
  filter(!grepl("-redo",NUSeq_Core_Facility_ID)) %>% 
  left_join((metadata_SPAH %>% dplyr::select(ID, ga_continuous)),
             by = c("Sample_ID" = "ID")) %>%
  filter(!is.na(ga_continuous)) %>%
  filter(!is.na(Sex_Predic))

metadata_eaa <- metadata_eaa %>%
  #calculate extrinsic epigenetic age acceleration 
  mutate(extAccel = residuals(lm(GA_CPC ~
                                   ga_continuous+
                                   Sex_Predic,
                                 metadata_eaa)),
         #calculate intrinsic epigenetic age acceleration
         intAccel = residuals(lm(GA_CPC ~
                                   ga_continuous + 
                                   Syncytiotrophoblast+
                                   Stromal+
                                   Endothelial+
                                   nRBC+
                                   Sex_Predic,
                                 metadata_eaa)))


metadata_SPAH <- metadata_SPAH %>% left_join(
  metadata_eaa %>% dplyr::select(Sample_ID, extAccel, intAccel),
  by = c("ID" = "Sample_ID")
)
```

## 3.0 Save Outputs


``` r
write.csv(metadata_SPAH, here::here("C. Processing/02_Outputs/00_Metadata/SPAH_updated_metadata.csv"), row.names = FALSE)
write.csv (array_qc_ss, here::here("C. Processing/02_Outputs/00_Metadata/SPAH_updated_array_qc_ss.csv"), row.names = FALSE)

saveRDS(metadata_SPAH, here::here("C. Processing/02_Outputs/00_Metadata/SPAH_updated_metadata.rds"))
saveRDS(array_qc_ss, here::here("C. Processing/02_Outputs/00_Metadata/SPAH_updated_array_qc_ss.rds"))
saveRDS(betas_bmiq, here::here("C. Processing/02_Outputs/00_Metadata/SPAH_bmiq_noob_betas_510samps.rds"))
```



