# DNA methylation profiling of placental histopathology
We evaluated the association of DNA methylation with the four major classes of placental pathology: maternal vascular malperfusion (MVM), fetal vascular malperfusion (FVM), acute and chronic inflammation (AI and CI).

The Stress Pregnancy and Health (SPAH) dataset was used. It is available at GSE307289. Additional pathology information is available upon request to the authors. 

00_SPAH_Metadata.rmd - combining metadata files, updating variable names, calculating new variables
01_SPAH_Pathology_Preprocessing.rmd - probe filtering, normalization, removing samples that don't meet quality control criteria
02_SPAH_Pathology_Epivariables.rmd - associations between placental pathology, clinical variables, and epivariables
03_SPAH_Pathology_Linear_Modelling.rmd - linear modelling for the four pathology classes
