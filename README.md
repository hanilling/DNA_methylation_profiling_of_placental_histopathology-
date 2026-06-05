# DNA methylation profiling of placental histopathology
In the paper **Characterizing placental dysfunction with DNA methylation profiling and placental histopathology** (Placenta, 2026; [PMID: 41689925](https://doi.org/10.1016/j.placenta.2026.02.008)) we evaluated the association of DNA methylation with the four major classes of placental pathology: maternal vascular malperfusion (MVM), fetal vascular malperfusion (FVM), acute and chronic inflammation (AI and CI).
<img width="1920" height="1080" alt="Graphical_Abstract" src="https://github.com/user-attachments/assets/be8f6cd2-2135-4685-b861-6a46f683e924" />

The **Stress Pregnancy and Health (SPAH)** dataset was used. It is available at [GSE307289](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE307289). Additional pathology information is available upon request to the authors. 

### Code Files:

- **00_SPAH_Metadata.rmd** - combining metadata files, updating variable names, calculating new variables
- **00_SPAH_Pathology_Preprocessing.rmd**- probe filtering, normalization, removing samples that don't meet quality control criteria
- **01_SPAH_Pathology_Epivariables.rmd** - associations between placental pathology, clinical variables, and epivariables
- **02_SPAH_Pathology_Linear_Modelling.rmd** - linear modelling for the four pathology classes
