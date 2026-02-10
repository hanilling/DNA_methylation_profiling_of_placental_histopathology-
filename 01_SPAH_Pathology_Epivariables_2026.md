---
title: "SPAH Placental Pathologies and Epivariables"
author: "Hannah Illing"
date: "2026 February 10"
output: 
  html_document:
    keep_md: yes
    code_folding: show
    toc: true  
    toc_depth: 4
    toc_float: 
      collapsed: true 
      smooth_scroll: true
---
  
## 0.0 Set-up

### 0.1 Load Packages 


``` r
library(here)
library(tidyverse)
library(ggpubr)
library(plomics)
library(cowplot)
library(rstatix)
library(jtools)
library(tableone)
library(Hmisc)
library(corrplot)
library(betareg)

set.seed(4815)
```

### 0.2 Style Guide


``` r
#set document style
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


#categories <- c("#5e356e","#de8e04","#35666e","grey33") 
categories <- c("#bf4e52","#631215","grey55","grey22") 
categories2 <- c("darkgrey","#D3AB04","#e05600","#3b6b4a","#033F63")
mvm_pal<- c("#87b8d6","#6099bd","#407799","#033F63","#02253B")
fvm_pal <- c("#7aa387","#5d966f","#3b6b4a","#114020")
ai_pal <- c("darkgrey","#FCDC5F","#D3AB04")
ci_pal <- c("darkgrey","#e38c56","#e05600")
```

### 0.3 Load Data


``` r
metadata <- readRDS(here::here("D. Pathology Project/02_Outputs/AA. Cleaned/SPAH_metadata_ss_filtered.rds"))
dim(metadata) #493 x 162
```

```
## [1] 493 144
```

``` r
#fix some variable classes
metadata <- metadata %>%
  mutate(Sex_Predic = as.factor(Sex_Predic),
         EOPE_class = as.factor(PE_GA_state == "EOPE"),
         LOPE_class = as.factor(PE_GA_state == "LOPE"),
         primi_gravida = as.factor(pcr_gravida == "1"),
         primi_para = as.factor(parity == "0"),
         Placenta_SGAvAGAvLGA = as.factor(Placenta_SGAvAGAvLGA),
         Feto_placental_wt = pcr_birth_wt/Placental_wt,
         pcr_diabetes_any = as.factor(pcr_diabetes_any),
         Race_AmIndALNative = as.factor(Race_AmIndALNative),
         Race_Asian = as.factor(Race_Asian),
         Race_Black = as.factor(Race_Black),
         Race_HINativePacIsl = as.factor(Race_HINativePacIsl), 
         Race_white = as.factor(Race_white),
         Race_Other = as.factor(Race_Other),
         Ethnicity_HispYN = as.factor(Ethnicity_HispYN),
         acute_inflammation_3cat = as.factor(acute_inflammation_3cat),
         chronic_inflammation_3cat = as.factor(chronic_inflammation_3cat)) %>% 
#note that one sample is mistakenly labelled as FVM in the FETAL_VASC_PATH column
#Use the FETAL_VASC_PATH column to create an updated FETAL_VASC_PATH variable
  mutate(FETAL_VASC_PATH = as.factor(ifelse(Fetal_vasc_path_3cat==0,0,1))) %>%
  mutate(diabetes_pregestational.new = ifelse(diabetes_gestational == 0,1,0),
         diabetes_gestational.new = ifelse(diabetes_pregestational == "0",1,0))
```
## 1.0 Tables 

### 1.1 Table 1


``` r
metadata_filt <- metadata %>% filter(!is.na(HarmonizedPath))
dim(metadata_filt) #9 samples removed
```

```
## [1] 484 151
```

``` r
metadata_filt <- metadata_filt %>%
  mutate(NoPath = ifelse(HarmonizedPath == "No Pathology","No Path","Path"))

myVars <- c("Sex_Predic", "ga_continuous", "Bw_call",
  "pcr_pre_eclampsia", "pcr_diabetes_any", "primi_para",
  "Age", "HarmonizedRace", "Ethnicity_HispYN","ACUTE_INFLAMMATION",
  "acute_inflammation_3cat", "CHRONIC_INFLAMMATION",
  "chronic_inflammation_3cat", "FETAL_VASC_PATH", "Fetal_vasc_path_3cat",
  "mvm_di", "mvm_score_3cat")

catVars <- c("Sex_Predic", "Bw_call",
  "pcr_pre_eclampsia", "pcr_diabetes_any", "primi_para",
  "HarmonizedRace", "Ethnicity_HispYN","ACUTE_INFLAMMATION",
  "acute_inflammation_3cat", "CHRONIC_INFLAMMATION",
  "chronic_inflammation_3cat", "FETAL_VASC_PATH", "Fetal_vasc_path_3cat",
  "mvm_di", "mvm_score_3cat")

CreateTableOne(vars = myVars, data = metadata_filt, factorVars = catVars)
```

```
##                                                                      
##                                                                       Overall      
##   n                                                                     484        
##   Sex_Predic = XY (%)                                                   271 (56.0) 
##   ga_continuous (mean (SD))                                           38.61 (2.03) 
##   Bw_call (%)                                                                      
##      SGA                                                                 23 ( 4.8) 
##      AGA                                                                397 (82.4) 
##      LGA                                                                 62 (12.9) 
##   pcr_pre_eclampsia = 1 (%)                                              49 (10.1) 
##   pcr_diabetes_any = 1 (%)                                               99 (20.5) 
##   primi_para = TRUE (%)                                                 181 (37.4) 
##   Age (mean (SD))                                                     33.37 (5.64) 
##   HarmonizedRace (%)                                                               
##      American Indian or Alaskan Native                                    1 ( 0.2) 
##      American Indian or Alaskan Native, Black/African American, White     1 ( 0.2) 
##      American Indian or Alaskan Native, White                             3 ( 0.6) 
##      Asian                                                               47 ( 9.7) 
##      Asian, White                                                         4 ( 0.8) 
##      Black/African American                                              75 (15.5) 
##      Black/African American, Other Race                                   2 ( 0.4) 
##      Black/African American, White                                        7 ( 1.4) 
##      Hawaiian Native or Pacific Islander                                  4 ( 0.8) 
##      Other Race                                                          43 ( 8.9) 
##      White                                                              297 (61.4) 
##   Ethnicity_HispYN = 1 (%)                                              119 (24.6) 
##   ACUTE_INFLAMMATION = 1 (%)                                            267 (55.2) 
##   acute_inflammation_3cat (%)                                                      
##      0                                                                  217 (44.8) 
##      1                                                                  202 (41.7) 
##      2                                                                   65 (13.4) 
##   CHRONIC_INFLAMMATION = 1 (%)                                          267 (55.2) 
##   chronic_inflammation_3cat (%)                                                    
##      0                                                                  217 (44.8) 
##      1                                                                  161 (33.3) 
##      2                                                                  106 (21.9) 
##   FETAL_VASC_PATH = 1 (%)                                               162 (33.5) 
##   Fetal_vasc_path_3cat (%)                                                         
##      0                                                                  322 (66.5) 
##      1                                                                  125 (25.8) 
##      2                                                                   37 ( 7.6) 
##   mvm_di = 1 (%)                                                        144 (29.8) 
##   mvm_score_3cat (%)                                                               
##      0                                                                  340 (70.2) 
##      1                                                                   84 (17.4) 
##      2                                                                   60 (12.4)
```

## 2.0 Figure 1


``` r
metadata_filter <- metadata %>% 
  filter(!is.na(HarmonizedPath))%>%
  dplyr::rename(CTB = Trophoblasts,
                  STB = Syncytiotrophoblast) %>%
    pivot_longer(cols=c(CTB,
                        STB,
                        Hofbauer,
                        Endothelial,
                        nRBC,
                        Stromal),
                 names_to = "Cell_Type", values_to="Estimate") %>%
    mutate(Cell_Type = factor(Cell_Type,
                              levels=c("nRBC", "Hofbauer", "Endothelial",
                                       "Stromal", "CTB", "STB")))

fig_1_a <- metadata_filter %>%
  ggplot(aes(x = Cell_Type, y = Estimate, fill = mvm_di))+
  geom_boxplot(alpha=0.5, outlier.shape = NA)+
  scale_fill_manual(values = c("darkgrey", "#147ab8"), labels=c("No","Yes"))+
  labs(x="Cell Estimate", y="Proportion of Sample", fill="MVM")+
  theme(axis.text = element_text(size=15, color="black"),
        strip.text = element_text(
        size = 15, color = "black", face = "bold"),
        legend.position="bottom",
        legend.text = element_text(size=13),
        legend.title = element_text(size=15))+
  ylim(0,1)

fig_1_b <- metadata_filter %>%
  filter(!mvm_score_3cat == 1) %>%
  ggplot(aes(x = Cell_Type, y = Estimate, fill = mvm_score_3cat))+
  geom_boxplot(alpha=0.5, outlier.shape = NA)+
  scale_color_manual(values = c("darkgrey", "#147ab8"), guide="none")+
  scale_fill_manual(values = c("darkgrey", "#147ab8"), labels=c("No MVM","High-grade MVM"))+
  labs(x="Cell Estimate", y="Proportion of Sample", fill="MVM")+
  theme(axis.text = element_text(size=15, color="black"),
        strip.text = element_text(
        size = 15, color = "black", face = "bold"),
        legend.position="bottom",
        legend.text = element_text(size=13),
        legend.title = element_text(size=15))+
  annotate(geom="text",label = "***", x=2,y=0.96, size=7)+
  annotate(geom="text",label = "**", x=4,y=0.96, size=7)+
  ylim(0,1)


fig_1_c <- metadata_filter %>%
  ggplot(aes(x = Cell_Type, y = Estimate, fill = FETAL_VASC_PATH))+
  geom_boxplot(alpha=0.5, outlier.shape = NA)+
  scale_fill_manual(values = c("darkgrey", "#2ea353"), labels=c("No","Yes"))+
  labs(x="Cell Estimate", y="Proportion of Sample", fill="FVM")+
  theme(axis.text = element_text(size=15, color="black"),
        strip.text = element_text(
        size = 15, color = "black", face = "bold"),
        legend.position="bottom",
        legend.text = element_text(size=13),
        legend.title = element_text(size=15))+
  annotate(geom="text",label = "***", x=2,y=0.96, size=7)+
  ylim(0,1)

fig_1_d <- metadata_filter %>%
  filter(!Fetal_vasc_path_3cat == 1) %>%
  ggplot(aes(x = Cell_Type, y = Estimate, fill = Fetal_vasc_path_3cat))+
  geom_boxplot(alpha=0.5, outlier.shape = NA)+
  scale_fill_manual(values = c("darkgrey", "#2ea353"), labels=c("No","High-grade"))+
  labs(x="Cell Estimate", y="Proportion of Sample", fill="FVM")+
  theme(axis.text = element_text(size=15, color="black"),
        strip.text = element_text(
        size = 15, color = "black", face = "bold"),
        legend.position="bottom",
        legend.text = element_text(size=13),
        legend.title = element_text(size=15))+
  annotate(geom="text",label = "*", x=2,y=0.96, size=7)+
  ylim(0,1)


fig_1 <- plot_grid(fig_1_a, NULL, fig_1_b, 
                   NULL, NULL, NULL,
                fig_1_c, NULL, fig_1_d,
                   labels = c("A.","","B.","","","","C.","","D."), nrow = 3, ncol=3,
                   rel_widths = c(1,0.1,1), rel_heights = c(1,0.1,1))

fig_1
```

![](01_SPAH_Pathology_Epivariables_2026_files/figure-html/figure-1-1.png)<!-- -->

``` r
ggsave("Figure 1_2026.png", plot = fig_1, path = here::here("D. Pathology Project/02_Outputs/AA. Cleaned/"), width = 15, height =9, device = "png")
```

## 3.0 Beta-Regressions


``` r
metadata_naomit <- metadata_filt %>%
  filter(!PE_GA_state == "EOPE") %>%
  filter(!is.na(ga_continuous))

## Individual Cell Types by AI
#AI ~ GA + STB
fit_ai_stb <- betareg(Syncytiotrophoblast ~ ACUTE_INFLAMMATION + ga_continuous, data = metadata_naomit)
summary(fit_ai_stb)
```

```
## 
## Call:
## betareg(formula = Syncytiotrophoblast ~ ACUTE_INFLAMMATION + ga_continuous, 
##     data = metadata_naomit)
## 
## Quantile residuals:
##     Min      1Q  Median      3Q     Max 
## -3.0232 -0.6800 -0.0980  0.6116  2.6168 
## 
## Coefficients (mean model with logit link):
##                     Estimate Std. Error z value Pr(>|z|)    
## (Intercept)         -0.58366    0.62449  -0.935 0.349983    
## ACUTE_INFLAMMATION1  0.00925    0.05910   0.157 0.875630    
## ga_continuous        0.05469    0.01630   3.355 0.000793 ***
## 
## Phi coefficients (precision model with identity link):
##       Estimate Std. Error z value Pr(>|z|)    
## (phi)  14.7606     0.9474   15.58   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
## 
## Type of estimator: ML (maximum likelihood)
## Log-likelihood: 477.8 on 4 Df
## Pseudo R-squared: 0.02042
## Number of iterations: 10 (BFGS) + 3 (Fisher scoring)
```

``` r
#AI ~ GA + CTB
fit_ai_ctb <- betareg(Trophoblasts ~ ACUTE_INFLAMMATION + ga_continuous, data = metadata_naomit)
```

```
## Loading required namespace: numDeriv
```

``` r
summary(fit_ai_ctb)
```

```
## 
## Call:
## betareg(formula = Trophoblasts ~ ACUTE_INFLAMMATION + ga_continuous, 
##     data = metadata_naomit)
## 
## Randomized quantile residuals:
##     Min      1Q  Median      3Q     Max 
## -2.4151 -0.7035 -0.0169  0.6549  2.5994 
## 
## Coefficients (mu model with logit link):
##                      Estimate Std. Error z value Pr(>|z|)    
## (Intercept)          2.052268   0.578541   3.547 0.000389 ***
## ACUTE_INFLAMMATION1 -0.008373   0.065315  -0.128 0.897990    
## ga_continuous       -0.111159   0.015464  -7.188 6.55e-13 ***
## 
## Phi coefficients (phi model with identity link):
##       Estimate Std. Error z value Pr(>|z|)    
## (phi)   38.343      5.595   6.853 7.24e-12 ***
## 
## Exceedence parameter (extended-support xbetax model):
##         Estimate Std. Error z value Pr(>|z|)    
## Log(nu)  -2.6978     0.1356   -19.9   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
## 
## Exceedence parameter nu: 0.0674
## Type of estimator: ML (maximum likelihood)
## Log-likelihood: 358.1 on 5 Df
## Number of iterations in BFGS optimization: 33
```

``` r
#AI ~ GA + Endothelial
fit_ai_end <- betareg(Endothelial ~ ACUTE_INFLAMMATION + ga_continuous, data = metadata_naomit)
summary(fit_ai_end)
```

```
## 
## Call:
## betareg(formula = Endothelial ~ ACUTE_INFLAMMATION + ga_continuous, data = metadata_naomit)
## 
## Quantile residuals:
##     Min      1Q  Median      3Q     Max 
## -4.8432 -0.4230  0.2064  0.6316  2.1279 
## 
## Coefficients (mean model with logit link):
##                      Estimate Std. Error z value Pr(>|z|)    
## (Intercept)         -2.674579   0.490534  -5.452 4.97e-08 ***
## ACUTE_INFLAMMATION1 -0.013629   0.044462  -0.307    0.759    
## ga_continuous       -0.001132   0.012774  -0.089    0.929    
## 
## Phi coefficients (precision model with identity link):
##       Estimate Std. Error z value Pr(>|z|)    
## (phi)   69.921      4.585   15.25   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
## 
## Type of estimator: ML (maximum likelihood)
## Log-likelihood:  1050 on 4 Df
## Pseudo R-squared: 0.0001994
## Number of iterations: 19 (BFGS) + 2 (Fisher scoring)
```

``` r
#AI ~ GA + Stromal
fit_ai_str <- betareg(Stromal ~ ACUTE_INFLAMMATION + ga_continuous, data = metadata_naomit)
summary(fit_ai_str)
```

```
## 
## Call:
## betareg(formula = Stromal ~ ACUTE_INFLAMMATION + ga_continuous, data = metadata_naomit)
## 
## Randomized quantile residuals:
##     Min      1Q  Median      3Q     Max 
## -3.0387 -0.6257  0.0243  0.6490  2.8190 
## 
## Coefficients (mu model with logit link):
##                     Estimate Std. Error z value Pr(>|z|)    
## (Intercept)         -3.73892    0.72569  -5.152 2.57e-07 ***
## ACUTE_INFLAMMATION1  0.04346    0.06097   0.713    0.476    
## ga_continuous        0.02067    0.01871   1.105    0.269    
## 
## Phi coefficients (phi model with identity link):
##       Estimate Std. Error z value Pr(>|z|)    
## (phi)   65.460      7.521   8.703   <2e-16 ***
## 
## Exceedence parameter (extended-support xbetax model):
##         Estimate Std. Error z value Pr(>|z|)    
## Log(nu)  -3.8414     0.1279  -30.03   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
## 
## Exceedence parameter nu: 0.0215
## Type of estimator: ML (maximum likelihood)
## Log-likelihood: 744.8 on 5 Df
## Number of iterations in BFGS optimization: 16
```

``` r
#AI ~ GA + nRBC
fit_ai_rbc <- betareg(nRBC ~ ACUTE_INFLAMMATION + ga_continuous, data = metadata_naomit)
summary(fit_ai_rbc)
```

```
## 
## Call:
## betareg(formula = nRBC ~ ACUTE_INFLAMMATION + ga_continuous, data = metadata_naomit)
## 
## Randomized quantile residuals:
##     Min      1Q  Median      3Q     Max 
## -3.0544 -0.5572  0.0087  0.5199  5.9441 
## 
## Coefficients (mu model with logit link):
##                      Estimate Std. Error z value Pr(>|z|)    
## (Intercept)         -4.108978   0.563298  -7.295    3e-13 ***
## ACUTE_INFLAMMATION1 -0.048046   0.049358  -0.973    0.330    
## ga_continuous        0.003985   0.014635   0.272    0.785    
## 
## Phi coefficients (phi model with identity link):
##       Estimate Std. Error z value Pr(>|z|)    
## (phi)   207.50      18.02   11.51   <2e-16 ***
## 
## Exceedence parameter (extended-support xbetax model):
##         Estimate Std. Error z value Pr(>|z|)    
## Log(nu)  -5.7767     0.1517  -38.08   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
## 
## Exceedence parameter nu: 0.0031
## Type of estimator: ML (maximum likelihood)
## Log-likelihood:  1469 on 5 Df
## Number of iterations in BFGS optimization: 28
```

``` r
#AI ~ GA + Hofbauer
fit_ai_hbc <- betareg(Hofbauer ~ ACUTE_INFLAMMATION + ga_continuous, data = metadata_naomit)
summary(fit_ai_hbc)
```

```
## 
## Call:
## betareg(formula = Hofbauer ~ ACUTE_INFLAMMATION + ga_continuous, data = metadata_naomit)
## 
## Randomized quantile residuals:
##     Min      1Q  Median      3Q     Max 
## -2.5390 -0.6332  0.0141  0.6479  3.4372 
## 
## Coefficients (mu model with logit link):
##                     Estimate Std. Error z value Pr(>|z|)    
## (Intercept)         -7.12247    0.90107  -7.904 2.69e-15 ***
## ACUTE_INFLAMMATION1 -0.08252    0.06719  -1.228 0.219402    
## ga_continuous        0.07695    0.02310   3.331 0.000865 ***
## 
## Phi coefficients (phi model with identity link):
##       Estimate Std. Error z value Pr(>|z|)    
## (phi)   195.43      24.36   8.024 1.03e-15 ***
## 
## Exceedence parameter (extended-support xbetax model):
##         Estimate Std. Error z value Pr(>|z|)    
## Log(nu)  -4.8645     0.1255  -38.77   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
## 
## Exceedence parameter nu: 0.0077
## Type of estimator: ML (maximum likelihood)
## Log-likelihood:  1050 on 5 Df
## Number of iterations in BFGS optimization: 26
```

``` r
## Individual Cell Types by High-Grade AI
ai_grade <- metadata_naomit %>%
  filter(!acute_inflammation_3cat == 1) 
table(ai_grade$acute_inflammation_3cat)
```

```
## 
##   0   1   2 
## 211   0  64
```

``` r
#AI ~ GA + STB
fit_hgai_stb <- betareg(Syncytiotrophoblast ~ ACUTE_INFLAMMATION + ga_continuous, data = ai_grade)
summary(fit_hgai_stb)
```

```
## 
## Call:
## betareg(formula = Syncytiotrophoblast ~ ACUTE_INFLAMMATION + ga_continuous, 
##     data = ai_grade)
## 
## Quantile residuals:
##     Min      1Q  Median      3Q     Max 
## -2.8764 -0.7243 -0.1184  0.6553  2.5888 
## 
## Coefficients (mean model with logit link):
##                      Estimate Std. Error z value Pr(>|z|)   
## (Intercept)         -0.863021   0.732607  -1.178  0.23879   
## ACUTE_INFLAMMATION1 -0.007549   0.091884  -0.082  0.93452   
## ga_continuous        0.061897   0.019144   3.233  0.00122 **
## 
## Phi coefficients (precision model with identity link):
##       Estimate Std. Error z value Pr(>|z|)    
## (phi)   14.356      1.211   11.85   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
## 
## Type of estimator: ML (maximum likelihood)
## Log-likelihood: 270.9 on 4 Df
## Pseudo R-squared: 0.03179
## Number of iterations: 10 (BFGS) + 3 (Fisher scoring)
```

``` r
#AI ~ GA + CTB
fit_hgai_ctb <- betareg(Trophoblasts ~ ACUTE_INFLAMMATION + ga_continuous, data = ai_grade)
summary(fit_hgai_ctb)
```

```
## 
## Call:
## betareg(formula = Trophoblasts ~ ACUTE_INFLAMMATION + ga_continuous, 
##     data = ai_grade)
## 
## Randomized quantile residuals:
##     Min      1Q  Median      3Q     Max 
## -2.9323 -0.6952  0.0106  0.6770  2.6263 
## 
## Coefficients (mu model with logit link):
##                     Estimate Std. Error z value Pr(>|z|)    
## (Intercept)          1.82067    0.64725   2.813  0.00491 ** 
## ACUTE_INFLAMMATION1  0.05325    0.09469   0.562  0.57389    
## ga_continuous       -0.10435    0.01733  -6.020 1.75e-09 ***
## 
## Phi coefficients (phi model with identity link):
##       Estimate Std. Error z value Pr(>|z|)    
## (phi)   42.482      7.875   5.394 6.88e-08 ***
## 
## Exceedence parameter (extended-support xbetax model):
##         Estimate Std. Error z value Pr(>|z|)    
## Log(nu)  -2.6299     0.1647  -15.97   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
## 
## Exceedence parameter nu: 0.0721
## Type of estimator: ML (maximum likelihood)
## Log-likelihood: 210.6 on 5 Df
## Number of iterations in BFGS optimization: 23
```

``` r
#AI ~ GA + Endothelial
fit_hgai_end <- betareg(Endothelial ~ ACUTE_INFLAMMATION + ga_continuous, data = ai_grade)
summary(fit_hgai_end)
```

```
## 
## Call:
## betareg(formula = Endothelial ~ ACUTE_INFLAMMATION + ga_continuous, data = ai_grade)
## 
## Quantile residuals:
##     Min      1Q  Median      3Q     Max 
## -3.6480 -0.4005  0.2068  0.6518  1.8998 
## 
## Coefficients (mean model with logit link):
##                      Estimate Std. Error z value Pr(>|z|)    
## (Intercept)         -2.627624   0.583829  -4.501 6.77e-06 ***
## ACUTE_INFLAMMATION1 -0.077586   0.070327  -1.103     0.27    
## ga_continuous       -0.002303   0.015214  -0.151     0.88    
## 
## Phi coefficients (precision model with identity link):
##       Estimate Std. Error z value Pr(>|z|)    
## (phi)   68.574      5.919   11.59   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
## 
## Type of estimator: ML (maximum likelihood)
## Log-likelihood: 605.8 on 4 Df
## Pseudo R-squared: 0.004113
## Number of iterations: 20 (BFGS) + 4 (Fisher scoring)
```

``` r
#AI ~ GA + Stromal
fit_hgai_str <- betareg(Stromal ~ ACUTE_INFLAMMATION + ga_continuous, data = ai_grade)
summary(fit_hgai_str)
```

```
## 
## Call:
## betareg(formula = Stromal ~ ACUTE_INFLAMMATION + ga_continuous, data = ai_grade)
## 
## Randomized quantile residuals:
##     Min      1Q  Median      3Q     Max 
## -2.7909 -0.7161  0.0303  0.7412  2.4751 
## 
## Coefficients (mu model with logit link):
##                      Estimate Std. Error z value Pr(>|z|)    
## (Intercept)         -2.975020   0.896710  -3.318 0.000908 ***
## ACUTE_INFLAMMATION1 -0.019558   0.102990  -0.190 0.849386    
## ga_continuous       -0.000206   0.023138  -0.009 0.992896    
## 
## Phi coefficients (phi model with identity link):
##       Estimate Std. Error z value Pr(>|z|)    
## (phi)    54.21       8.81   6.153  7.6e-10 ***
## 
## Exceedence parameter (extended-support xbetax model):
##         Estimate Std. Error z value Pr(>|z|)    
## Log(nu)   -3.988      0.204  -19.55   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
## 
## Exceedence parameter nu: 0.0185
## Type of estimator: ML (maximum likelihood)
## Log-likelihood: 425.4 on 5 Df
## Number of iterations in BFGS optimization: 25
```

``` r
#AI ~ GA + nRBC
fit_hgai_rbc <- betareg(nRBC ~ ACUTE_INFLAMMATION + ga_continuous, data = ai_grade)
summary(fit_hgai_rbc)
```

```
## 
## Call:
## betareg(formula = nRBC ~ ACUTE_INFLAMMATION + ga_continuous, data = ai_grade)
## 
## Randomized quantile residuals:
##     Min      1Q  Median      3Q     Max 
## -3.4270 -0.6102  0.0256  0.4450  5.3354 
## 
## Coefficients (mu model with logit link):
##                      Estimate Std. Error z value Pr(>|z|)    
## (Intercept)         -4.077557   0.707258  -5.765 8.15e-09 ***
## ACUTE_INFLAMMATION1  0.045438   0.079794   0.569    0.569    
## ga_continuous        0.002205   0.018390   0.120    0.905    
## 
## Phi coefficients (phi model with identity link):
##       Estimate Std. Error z value Pr(>|z|)    
## (phi)   163.71      16.96   9.654   <2e-16 ***
## 
## Exceedence parameter (extended-support xbetax model):
##         Estimate Std. Error z value Pr(>|z|)    
## Log(nu)  -6.2247     0.2589  -24.05   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
## 
## Exceedence parameter nu: 0.002
## Type of estimator: ML (maximum likelihood)
## Log-likelihood: 851.2 on 5 Df
## Number of iterations in BFGS optimization: 12
```

``` r
#AI ~ GA + Hofbauer
fit_hgai_hbc <- betareg(Hofbauer ~ ACUTE_INFLAMMATION + ga_continuous, data = ai_grade)
summary(fit_hgai_hbc)
```

```
## 
## Call:
## betareg(formula = Hofbauer ~ ACUTE_INFLAMMATION + ga_continuous, data = ai_grade)
## 
## Randomized quantile residuals:
##     Min      1Q  Median      3Q     Max 
## -2.7787 -0.6361  0.0216  0.6990  3.1938 
## 
## Coefficients (mu model with logit link):
##                     Estimate Std. Error z value Pr(>|z|)    
## (Intercept)         -6.34439    1.07627  -5.895 3.75e-09 ***
## ACUTE_INFLAMMATION1 -0.06579    0.10898  -0.604    0.546    
## ga_continuous        0.05647    0.02763   2.044    0.041 *  
## 
## Phi coefficients (phi model with identity link):
##       Estimate Std. Error z value Pr(>|z|)    
## (phi)   167.20      27.45    6.09 1.13e-09 ***
## 
## Exceedence parameter (extended-support xbetax model):
##         Estimate Std. Error z value Pr(>|z|)    
## Log(nu)  -4.9025     0.1784  -27.47   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
## 
## Exceedence parameter nu: 0.0074
## Type of estimator: ML (maximum likelihood)
## Log-likelihood: 596.9 on 5 Df
## Number of iterations in BFGS optimization: 25
```

``` r
## Individual Cell Types by MVM
#MVM ~ GA + STB
fit_mvm_stb <- betareg(Syncytiotrophoblast ~ mvm_di + ga_continuous, data = metadata_naomit)
summary(fit_mvm_stb)
```

```
## 
## Call:
## betareg(formula = Syncytiotrophoblast ~ mvm_di + ga_continuous, data = metadata_naomit)
## 
## Quantile residuals:
##     Min      1Q  Median      3Q     Max 
## -2.9282 -0.6839 -0.1313  0.5968  2.6602 
## 
## Coefficients (mean model with logit link):
##               Estimate Std. Error z value Pr(>|z|)    
## (Intercept)   -0.80980    0.63980  -1.266 0.205620    
## mvm_di1        0.07064    0.06566   1.076 0.281958    
## ga_continuous  0.06016    0.01642   3.664 0.000248 ***
## 
## Phi coefficients (precision model with identity link):
##       Estimate Std. Error z value Pr(>|z|)    
## (phi)  14.7987     0.9499   15.58   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
## 
## Type of estimator: ML (maximum likelihood)
## Log-likelihood: 478.4 on 4 Df
## Pseudo R-squared: 0.02216
## Number of iterations: 10 (BFGS) + 3 (Fisher scoring)
```

``` r
#MVM ~ GA + CTB
fit_mvm_ctb <- betareg(Trophoblasts ~ mvm_di + ga_continuous, data = metadata_naomit)
summary(fit_mvm_ctb)
```

```
## 
## Call:
## betareg(formula = Trophoblasts ~ mvm_di + ga_continuous, data = metadata_naomit)
## 
## Randomized quantile residuals:
##     Min      1Q  Median      3Q     Max 
## -3.6882 -0.7236  0.0200  0.6699  2.5605 
## 
## Coefficients (mu model with logit link):
##               Estimate Std. Error z value Pr(>|z|)    
## (Intercept)    2.34103    0.61010   3.837 0.000124 ***
## mvm_di1       -0.07475    0.07342  -1.018 0.308595    
## ga_continuous -0.11816    0.01591  -7.428  1.1e-13 ***
## 
## Phi coefficients (phi model with identity link):
##       Estimate Std. Error z value Pr(>|z|)    
## (phi)   38.478      5.619   6.848 7.49e-12 ***
## 
## Exceedence parameter (extended-support xbetax model):
##         Estimate Std. Error z value Pr(>|z|)    
## Log(nu)  -2.6931     0.1352  -19.93   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
## 
## Exceedence parameter nu: 0.0677
## Type of estimator: ML (maximum likelihood)
## Log-likelihood: 358.5 on 5 Df
## Number of iterations in BFGS optimization: 23
```

``` r
#MVM ~ GA + Endothelial
fit_mvm_end <- betareg(Endothelial ~ mvm_di + ga_continuous, data = metadata_naomit)
summary(fit_mvm_end)
```

```
## 
## Call:
## betareg(formula = Endothelial ~ mvm_di + ga_continuous, data = metadata_naomit)
## 
## Quantile residuals:
##     Min      1Q  Median      3Q     Max 
## -4.8712 -0.4206  0.2058  0.6243  2.1084 
## 
## Coefficients (mean model with logit link):
##                Estimate Std. Error z value Pr(>|z|)    
## (Intercept)   -2.643911   0.500062  -5.287 1.24e-07 ***
## mvm_di1       -0.003129   0.049105  -0.064    0.949    
## ga_continuous -0.002096   0.012811  -0.164    0.870    
## 
## Phi coefficients (precision model with identity link):
##       Estimate Std. Error z value Pr(>|z|)    
## (phi)   69.907      4.584   15.25   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
## 
## Type of estimator: ML (maximum likelihood)
## Log-likelihood:  1050 on 4 Df
## Pseudo R-squared: 4.613e-05
## Number of iterations: 19 (BFGS) + 2 (Fisher scoring)
```

``` r
#MVM ~ GA + Stromal
fit_mvm_str <- betareg(Stromal ~ mvm_di + ga_continuous, data = metadata_naomit)
summary(fit_mvm_str)
```

```
## 
## Call:
## betareg(formula = Stromal ~ mvm_di + ga_continuous, data = metadata_naomit)
## 
## Randomized quantile residuals:
##     Min      1Q  Median      3Q     Max 
## -3.0894 -0.6287  0.0104  0.6935  3.0141 
## 
## Coefficients (mu model with logit link):
##               Estimate Std. Error z value Pr(>|z|)    
## (Intercept)   -3.44710    0.73993  -4.659 3.18e-06 ***
## mvm_di1       -0.12717    0.06915  -1.839   0.0659 .  
## ga_continuous  0.01470    0.01882   0.781   0.4347    
## 
## Phi coefficients (phi model with identity link):
##       Estimate Std. Error z value Pr(>|z|)    
## (phi)   66.146      7.595   8.709   <2e-16 ***
## 
## Exceedence parameter (extended-support xbetax model):
##         Estimate Std. Error z value Pr(>|z|)    
## Log(nu)  -3.8384     0.1272  -30.18   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
## 
## Exceedence parameter nu: 0.0215
## Type of estimator: ML (maximum likelihood)
## Log-likelihood: 746.3 on 5 Df
## Number of iterations in BFGS optimization: 16
```

``` r
#MVM ~ GA + nRBC
fit_mvm_rbc <- betareg(nRBC ~ mvm_di + ga_continuous, data = metadata_naomit)
summary(fit_mvm_rbc)
```

```
## 
## Call:
## betareg(formula = nRBC ~ mvm_di + ga_continuous, data = metadata_naomit)
## 
## Randomized quantile residuals:
##     Min      1Q  Median      3Q     Max 
## -2.8823 -0.5473  0.0315  0.5199  5.9364 
## 
## Coefficients (mu model with logit link):
##                Estimate Std. Error z value Pr(>|z|)    
## (Intercept)   -4.048784   0.569134  -7.114 1.13e-12 ***
## mvm_di1        0.026459   0.054613   0.484    0.628    
## ga_continuous  0.001564   0.014577   0.107    0.915    
## 
## Phi coefficients (phi model with identity link):
##       Estimate Std. Error z value Pr(>|z|)    
## (phi)   203.94      17.59   11.59   <2e-16 ***
## 
## Exceedence parameter (extended-support xbetax model):
##         Estimate Std. Error z value Pr(>|z|)    
## Log(nu)  -5.7831     0.1515  -38.17   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
## 
## Exceedence parameter nu: 0.0031
## Type of estimator: ML (maximum likelihood)
## Log-likelihood:  1469 on 5 Df
## Number of iterations in BFGS optimization: 31
```

``` r
#MVM ~ GA + Hofbauer
fit_mvm_hbc <- betareg(Hofbauer ~ mvm_di + ga_continuous, data = metadata_naomit)
summary(fit_mvm_hbc)
```

```
## 
## Call:
## betareg(formula = Hofbauer ~ mvm_di + ga_continuous, data = metadata_naomit)
## 
## Randomized quantile residuals:
##     Min      1Q  Median      3Q     Max 
## -3.6403 -0.6307 -0.0105  0.6686  3.4287 
## 
## Coefficients (mu model with logit link):
##               Estimate Std. Error z value Pr(>|z|)    
## (Intercept)   -6.65008    0.89845  -7.402 1.34e-13 ***
## mvm_di1       -0.11716    0.07701  -1.521  0.12817    
## ga_continuous  0.06434    0.02275   2.828  0.00469 ** 
## 
## Phi coefficients (phi model with identity link):
##       Estimate Std. Error z value Pr(>|z|)    
## (phi)   194.55      24.27   8.017 1.09e-15 ***
## 
## Exceedence parameter (extended-support xbetax model):
##         Estimate Std. Error z value Pr(>|z|)    
## Log(nu)  -4.8718     0.1263  -38.58   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
## 
## Exceedence parameter nu: 0.0077
## Type of estimator: ML (maximum likelihood)
## Log-likelihood:  1051 on 5 Df
## Number of iterations in BFGS optimization: 25
```

``` r
## Individual Cell Types by High-grade MVM
mvm_grade <- metadata_naomit %>%
  filter(!mvm_score_3cat == 1) 

#MVM ~ GA + STB
fit_hgmvm_stb <- betareg(Syncytiotrophoblast ~ mvm_di + ga_continuous, data = mvm_grade)
summary(fit_hgmvm_stb)
```

```
## 
## Call:
## betareg(formula = Syncytiotrophoblast ~ mvm_di + ga_continuous, data = mvm_grade)
## 
## Quantile residuals:
##     Min      1Q  Median      3Q     Max 
## -2.6567 -0.6964 -0.1443  0.6040  2.7454 
## 
## Coefficients (mean model with logit link):
##               Estimate Std. Error z value Pr(>|z|)    
## (Intercept)   -1.79710    0.75829  -2.370   0.0178 *  
## mvm_di1        0.20717    0.09812   2.111   0.0347 *  
## ga_continuous  0.08560    0.01948   4.394 1.11e-05 ***
## 
## Phi coefficients (precision model with identity link):
##       Estimate Std. Error z value Pr(>|z|)    
## (phi)   15.149      1.072   14.13   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
## 
## Type of estimator: ML (maximum likelihood)
## Log-likelihood: 400.2 on 4 Df
## Pseudo R-squared: 0.04081
## Number of iterations: 10 (BFGS) + 2 (Fisher scoring)
```

``` r
#MVM ~ GA + CTB
fit_hgmvm_ctb <- betareg(Trophoblasts ~ mvm_di + ga_continuous, data = mvm_grade)
summary(fit_hgmvm_ctb)
```

```
## 
## Call:
## betareg(formula = Trophoblasts ~ mvm_di + ga_continuous, data = mvm_grade)
## 
## Randomized quantile residuals:
##     Min      1Q  Median      3Q     Max 
## -2.8725 -0.7512  0.0062  0.6432  2.6042 
## 
## Coefficients (mu model with logit link):
##               Estimate Std. Error z value Pr(>|z|)    
## (Intercept)    3.76775    0.72625   5.188 2.13e-07 ***
## mvm_di1       -0.28658    0.11073  -2.588  0.00965 ** 
## ga_continuous -0.15375    0.01899  -8.095 5.70e-16 ***
## 
## Phi coefficients (phi model with identity link):
##       Estimate Std. Error z value Pr(>|z|)    
## (phi)   42.465      7.135   5.951 2.66e-09 ***
## 
## Exceedence parameter (extended-support xbetax model):
##         Estimate Std. Error z value Pr(>|z|)    
## Log(nu)  -2.5608     0.1427  -17.95   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
## 
## Exceedence parameter nu: 0.0772
## Type of estimator: ML (maximum likelihood)
## Log-likelihood: 273.2 on 5 Df
## Number of iterations in BFGS optimization: 22
```

``` r
#MVM ~ GA + Endothelial
fit_hgmvm_end <- betareg(Endothelial ~ mvm_di + ga_continuous, data = mvm_grade)
summary(fit_hgmvm_end)
```

```
## 
## Call:
## betareg(formula = Endothelial ~ mvm_di + ga_continuous, data = mvm_grade)
## 
## Quantile residuals:
##     Min      1Q  Median      3Q     Max 
## -4.8823 -0.4038  0.2079  0.6109  2.0120 
## 
## Coefficients (mean model with logit link):
##               Estimate Std. Error z value Pr(>|z|)    
## (Intercept)   -2.27133    0.60143  -3.777 0.000159 ***
## mvm_di1       -0.05981    0.07429  -0.805 0.420718    
## ga_continuous -0.01158    0.01542  -0.751 0.452512    
## 
## Phi coefficients (precision model with identity link):
##       Estimate Std. Error z value Pr(>|z|)    
## (phi)   68.227      4.932   13.83   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
## 
## Type of estimator: ML (maximum likelihood)
## Log-likelihood: 861.6 on 4 Df
## Pseudo R-squared: 0.001912
## Number of iterations: 10 (BFGS) + 3 (Fisher scoring)
```

``` r
#MVM ~ GA + Stromal
fit_hgmvm_str <- betareg(Stromal ~ mvm_di + ga_continuous, data = mvm_grade)
summary(fit_hgmvm_str)
```

```
## 
## Call:
## betareg(formula = Stromal ~ mvm_di + ga_continuous, data = mvm_grade)
## 
## Randomized quantile residuals:
##     Min      1Q  Median      3Q     Max 
## -3.5126 -0.6574 -0.0170  0.6611  2.7768 
## 
## Coefficients (mu model with logit link):
##               Estimate Std. Error z value Pr(>|z|)   
## (Intercept)   -2.45047    0.79436  -3.085  0.00204 **
## mvm_di1       -0.28607    0.10016  -2.856  0.00429 **
## ga_continuous -0.01030    0.02022  -0.509  0.61056   
## 
## Phi coefficients (phi model with identity link):
##       Estimate Std. Error z value Pr(>|z|)    
## (phi)    83.48      11.53   7.242 4.42e-13 ***
## 
## Exceedence parameter (extended-support xbetax model):
##         Estimate Std. Error z value Pr(>|z|)    
## Log(nu)  -3.7623     0.1326  -28.37   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
## 
## Exceedence parameter nu: 0.0232
## Type of estimator: ML (maximum likelihood)
## Log-likelihood: 640.6 on 5 Df
## Number of iterations in BFGS optimization: 23
```

``` r
#MVM ~ GA + nRBC
fit_hgmvm_rbc <- betareg(nRBC ~ mvm_di + ga_continuous, data = mvm_grade)
summary(fit_hgmvm_rbc)
```

```
## 
## Call:
## betareg(formula = nRBC ~ mvm_di + ga_continuous, data = mvm_grade)
## 
## Randomized quantile residuals:
##     Min      1Q  Median      3Q     Max 
## -3.2600 -0.5488  0.0192  0.4956  5.9582 
## 
## Coefficients (mu model with logit link):
##               Estimate Std. Error z value Pr(>|z|)    
## (Intercept)   -4.45925    0.68121  -6.546 5.91e-11 ***
## mvm_di1        0.10222    0.07754   1.318    0.187    
## ga_continuous  0.01132    0.01743   0.650    0.516    
## 
## Phi coefficients (phi model with identity link):
##       Estimate Std. Error z value Pr(>|z|)    
## (phi)    201.9       18.8   10.74   <2e-16 ***
## 
## Exceedence parameter (extended-support xbetax model):
##         Estimate Std. Error z value Pr(>|z|)    
## Log(nu)  -5.9819     0.1871  -31.96   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
## 
## Exceedence parameter nu: 0.0025
## Type of estimator: ML (maximum likelihood)
## Log-likelihood:  1233 on 5 Df
## Number of iterations in BFGS optimization: 30
```

``` r
#MVM ~ GA + Hofbauer
fit_hgmvm_hbc <- betareg(Hofbauer ~ mvm_di + ga_continuous, data = mvm_grade)
summary(fit_hgai_hbc)
```

```
## 
## Call:
## betareg(formula = Hofbauer ~ ACUTE_INFLAMMATION + ga_continuous, data = ai_grade)
## 
## Randomized quantile residuals:
##     Min      1Q  Median      3Q     Max 
## -2.4310 -0.6234  0.0216  0.6990  3.1938 
## 
## Coefficients (mu model with logit link):
##                     Estimate Std. Error z value Pr(>|z|)    
## (Intercept)         -6.34439    1.07627  -5.895 3.75e-09 ***
## ACUTE_INFLAMMATION1 -0.06579    0.10898  -0.604    0.546    
## ga_continuous        0.05647    0.02763   2.044    0.041 *  
## 
## Phi coefficients (phi model with identity link):
##       Estimate Std. Error z value Pr(>|z|)    
## (phi)   167.20      27.45    6.09 1.13e-09 ***
## 
## Exceedence parameter (extended-support xbetax model):
##         Estimate Std. Error z value Pr(>|z|)    
## Log(nu)  -4.9025     0.1784  -27.47   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
## 
## Exceedence parameter nu: 0.0074
## Type of estimator: ML (maximum likelihood)
## Log-likelihood: 596.9 on 5 Df
## Number of iterations in BFGS optimization: 25
```

``` r
newdat <- data.frame(
  mvm_di = factor(c(0, 1)),
  ga_continuous = mean(mvm_grade$ga_continuous, na.rm = TRUE)
)

pred <- predict(fit_hgmvm_stb, newdat, type = "response")
diff(pred)
```

```
##          2 
## 0.02838372
```

## 4.0 Supplementary Figures

### 4.1 Supplementary Figure 1

Does epigenetic age acceleration age differ by pathology? 


``` r
eaa_fig_a <- metadata %>%
  filter(!is.na(ACUTE_INFLAMMATION))%>%
  ggplot(aes(x=acute_inflammation_3cat, y=extAccel, color=acute_inflammation_3cat, fill=acute_inflammation_3cat))+
  stat_compare_means(aes(group=acute_inflammation_3cat), label.x = 0.95)+
  geom_point(alpha=0.5, position=position_jitterdodge(jitter.width=0.2),size=3)+
  geom_boxplot(alpha=0.5, outlier.shape = NA)+
  scale_color_manual(values=ai_pal)+
  scale_fill_manual(values=ai_pal)+
  labs(y="Extrinsic EAA (weeks)",x="AI grade")+
  scale_x_discrete(labels = c("0" = "No AI", "1" = "Low-Grade AI", "2" = "High-Grade AI"))+
  ylim(-2.7,3.5)

eaa_fig_b <- metadata %>%
  filter(!is.na(ACUTE_INFLAMMATION))%>%
  ggplot(aes(x=acute_inflammation_3cat, y=intAccel, color=acute_inflammation_3cat, fill=acute_inflammation_3cat))+
  stat_compare_means(aes(group=acute_inflammation_3cat), label.x = 0.95)+
  geom_point(alpha=0.5, position=position_jitterdodge(jitter.width=0.2),size=3)+
  geom_boxplot(alpha=0.5, outlier.shape = NA)+
  scale_color_manual(values=ai_pal)+
  scale_fill_manual(values=ai_pal)+
  labs(y="Intrinsic EAA (weeks)",x="AI grade")+
  scale_x_discrete(labels = c("0" = "No AI", "1" = "Low-Grade AI", "2" = "High-Grade AI"))+
  ylim(-2.7,3.5)

eaa_fig_c <- metadata %>%
  filter(!is.na(CHRONIC_INFLAMMATION))%>%
  ggplot(aes(x=chronic_inflammation_3cat, y=extAccel, color=chronic_inflammation_3cat, fill=chronic_inflammation_3cat))+
  stat_compare_means(aes(group=chronic_inflammation_3cat), label.x = 0.95)+
  geom_point(alpha=0.5, position=position_jitterdodge(jitter.width=0.2),size=3)+
  geom_boxplot(alpha=0.5, outlier.shape = NA)+
  scale_color_manual(values=ci_pal)+
  scale_fill_manual(values=ci_pal)+
  labs(y="Extrinsic EAA (weeks)",x="CI grade")+
  scale_x_discrete(labels = c("0" = "No CI", "1" = "Low-Grade CI", "2" = "High-Grade CI"))+
  ylim(-2.7,3.5)

eaa_fig_d <- metadata %>%
  filter(!is.na(CHRONIC_INFLAMMATION))%>%
  ggplot(aes(x=chronic_inflammation_3cat, y=intAccel, color=chronic_inflammation_3cat, fill=chronic_inflammation_3cat))+
  stat_compare_means(aes(group=chronic_inflammation_3cat), label.x = 0.95)+
  geom_point(alpha=0.5, position=position_jitterdodge(jitter.width=0.2),size=3)+
  geom_boxplot(alpha=0.5, outlier.shape = NA)+
  scale_color_manual(values=ci_pal)+
  scale_fill_manual(values=ci_pal)+
  labs(y="Intrinsic EAA (weeks)",x="CI grade")+
  scale_x_discrete(labels = c("0" = "No CI", "1" = "Low-Grade CI", "2" = "High-Grade CI"))+
  ylim(-2.7,3.5)

eaa_fig_e <- metadata %>%
  filter(!is.na(FETAL_VASC_PATH))%>%
  ggplot(aes(x=Fetal_vasc_path_3cat, y=extAccel, color=Fetal_vasc_path_3cat, fill=Fetal_vasc_path_3cat))+
  stat_compare_means(aes(group=Fetal_vasc_path_3cat), label.x = 0.95)+
  geom_point(alpha=0.5, position=position_jitterdodge(jitter.width=0.2),size=3)+
  geom_boxplot(alpha=0.5, outlier.shape = NA)+
  scale_color_manual(values=c("darkgrey",fvm_pal[1],fvm_pal[3]))+
  scale_fill_manual(values=c("darkgrey",fvm_pal[1],fvm_pal[3]))+
  labs(y="Extrinsic EAA (weeks)",x="FVM grade")+
  scale_x_discrete(labels = c("0" = "No FVM", "1" = "Low-Grade FVM", "2" = "High-Grade FVM"))+
  ylim(-2.7,3.5)

eaa_fig_f <- metadata %>%
  filter(!is.na(FETAL_VASC_PATH))%>%
  ggplot(aes(x=Fetal_vasc_path_3cat, y=intAccel, color=Fetal_vasc_path_3cat, fill=Fetal_vasc_path_3cat))+
  stat_compare_means(aes(group=Fetal_vasc_path_3cat), label.x = 0.95)+
  geom_point(alpha=0.5, position=position_jitterdodge(jitter.width=0.2),size=3)+
  geom_boxplot(alpha=0.5, outlier.shape = NA)+
  scale_color_manual(values=c("darkgrey",fvm_pal[1],fvm_pal[3]))+
  scale_fill_manual(values=c("darkgrey",fvm_pal[1],fvm_pal[3]))+
  labs(y="Intrinsic EAA (weeks)",x="FVM grade")+
  scale_x_discrete(labels = c("0" = "No FVM", "1" = "Low-Grade FVM", "2" = "High-Grade FVM"))+
  ylim(-2.7,3.5)

eaa_fig_g <- metadata %>%
  filter(!is.na(mvm_di))%>%
  ggplot(aes(x=mvm_score_3cat, y=extAccel, color = mvm_score_3cat, fill= mvm_score_3cat))+
  geom_point(alpha=0.5, position=position_jitterdodge(jitter.width=0.2),size=3)+
  geom_boxplot(alpha=0.5, outlier.shape = NA)+
  scale_fill_manual(values = c(categories2[1],mvm_pal[3], mvm_pal[4]))+
  scale_color_manual(values = c(categories2[1],mvm_pal[3],mvm_pal[4]))+
  labs(x="MVM Grade", y="Extrinsic EAA (weeks)")+
  stat_compare_means(aes(group=mvm_score_3cat), label.x = 0.95)+
  scale_x_discrete(labels = c("0" = "No MVM", "1" = "Low-Grade MVM", "2" = "High-Grade MVM"))+
  ylim(-2.7,3.5)

eaa_fig_h <- metadata %>%
  filter(!is.na(mvm_di))%>%
  ggplot(aes(x=mvm_score_3cat, y=intAccel, color = mvm_score_3cat, fill= mvm_score_3cat))+
  geom_point(alpha=0.5, position=position_jitterdodge(jitter.width=0.2),size=3)+
  geom_boxplot(alpha=0.5, outlier.shape = NA)+
  scale_fill_manual(values = c(categories2[1],mvm_pal[3], mvm_pal[4]))+
  scale_color_manual(values = c(categories2[1],mvm_pal[3],mvm_pal[4]))+
  labs(x="MVM Grade", y="Intrinsic EAA (weeks)")+
  stat_compare_means(aes(group=mvm_score_3cat), label.x = 0.95)+
  scale_x_discrete(labels = c("0" = "No MVM", "1" = "Low-Grade MVM", "2" = "High-Grade MVM"))+
  ylim(-2.7,3.5)

eaa_fig <- plot_grid(eaa_fig_a, eaa_fig_b, eaa_fig_c, eaa_fig_d,
                       eaa_fig_e, eaa_fig_f, eaa_fig_g, eaa_fig_h,
                   labels = c("A.","B.","C.", "D.", "E.","F.", "G.","H."), 
                   nrow = 4, ncol=2)

ggsave("Figure S1.png", plot = eaa_fig, path = here::here("D. Pathology Project/02_Outputs/AA. Cleaned/"), width = 12, height = 14, device = "png")
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
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
##  [1] betareg_3.2-4   corrplot_0.92   Hmisc_5.2-3     tableone_0.13.2
##  [5] jtools_2.3.0    rstatix_0.7.2   cowplot_1.1.3   plomics_0.2.0  
##  [9] ggpubr_0.6.0    lubridate_1.9.4 forcats_1.0.0   stringr_1.5.0  
## [13] dplyr_1.1.4     purrr_1.0.2     readr_2.1.4     tidyr_1.3.0    
## [17] tibble_3.2.1    ggplot2_3.5.2   tidyverse_2.0.0 here_1.0.1     
## 
## loaded via a namespace (and not attached):
##  [1] DBI_1.1.3           gridExtra_2.3       sandwich_3.1-1     
##  [4] rlang_1.1.4         magrittr_2.0.3      furrr_0.3.1        
##  [7] e1071_1.7-13        compiler_4.3.1      flexmix_2.3-20     
## [10] systemfonts_1.2.3   vctrs_0.6.5         pkgconfig_2.0.3    
## [13] fastmap_1.2.0       backports_1.5.0     labeling_0.4.3     
## [16] pander_0.6.6        utf8_1.2.5          rmarkdown_2.29     
## [19] tzdb_0.4.0          haven_2.5.3         ragg_1.2.5         
## [22] xfun_0.52           modeltools_0.2-24   cachem_1.1.0       
## [25] labelled_2.12.0     jsonlite_2.0.0      broom_1.0.5        
## [28] parallel_4.3.1      cluster_2.1.4       R6_2.5.1           
## [31] bslib_0.9.0         stringi_1.8.7       RColorBrewer_1.1-3 
## [34] parallelly_1.37.0   car_3.1-2           rpart_4.1.19       
## [37] numDeriv_2016.8-1.1 lmtest_0.9-40       jquerylib_0.1.4    
## [40] Rcpp_1.0.14         knitr_1.50          zoo_1.8-12         
## [43] base64enc_0.1-3     Matrix_1.6-1        splines_4.3.1      
## [46] nnet_7.3-19         timechange_0.3.0    tidyselect_1.2.1   
## [49] rstudioapi_0.17.1   dichromat_2.0-0.1   abind_1.4-5        
## [52] yaml_2.3.10         codetools_0.2-19    listenv_0.9.1      
## [55] lattice_0.21-8      withr_3.0.2         evaluate_1.0.4     
## [58] foreign_0.8-84      future_1.33.1       survival_3.5-7     
## [61] proxy_0.4-27        survey_4.2-1        pillar_1.9.0       
## [64] carData_3.0-5       checkmate_2.3.2     stats4_4.3.1       
## [67] generics_0.1.3      rprojroot_2.0.3     hms_1.1.3          
## [70] scales_1.4.0        globals_0.16.2      class_7.3-22       
## [73] glue_1.8.0          tools_4.3.1         data.table_1.17.6  
## [76] ggsignif_0.6.4      grid_4.3.1          mitools_2.4        
## [79] colorspace_2.1-1    nlme_3.1-163        htmlTable_2.4.3    
## [82] Formula_1.2-5       cli_3.6.3           textshaping_1.0.1  
## [85] fansi_1.0.6         gtable_0.3.6        broom.mixed_0.2.9.6
## [88] sass_0.4.10         digest_0.6.37       htmlwidgets_1.6.4  
## [91] farver_2.1.2        htmltools_0.5.8.1   lifecycle_1.0.4    
## [94] statmod_1.5.0
```
