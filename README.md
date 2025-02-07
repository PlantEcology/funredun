# funredun
Calculates community functional redundancy as subtraction of functional diversity (Rao's Q) from Simpson's D.

## Install
Installation in R requires <a href="https://cran.r-project.org/package=devtools">devtools</a> or <a href="https://cran.r-project.org/package=remotes">remotes</a>. Requires <a href="https://cran.r-project.org/package=vegan">vegan</a> and <a href="https://cran.r-project.org/package=usedist">usedist</a> packages to run.
```
devtools::install_github("PlantEcology/funredun")

OR

remotes::install_github("PlantEcology/funredun")
```

## Data Format
Function uses two data frames to calculate site functional redundancy spDat and funDat.

```
funredun(spDat, funDat, method = "gower", redund = TRUE, funDiv = FALSE)
```
Variable | Description 
------|-----
spDat | Data frame with rows as sites, columns as species, and elements as counts
funDat | Data frame with rows as species (same as spDat column names), columns as functional traits, elements as counts, measures, binary, etc.
method | Available options include "bray", "gower", and "altGower". See <a href="https://cran.r-project.org/web/packages/vegan/index.html">vegan::vegdist</a> for details. Default is Gower distance.
redund | Redundancy calculation as difference from Simpson's D (R = D - Q) or uniqueness (U = Q/D). Default is difference (TRUE).
funDiv | Functional Diversity as Rao's Q (Botta-Dukát 2005). Default is FALSE.

## Output
A data frame is produced with rows as sites (from spDat) and a column of functional redundancy (includes a column of functional diversity if option selected).

## References
Botta-Dukát, Z. 2005. Rao's quadratic entropy as a measure of functional diversity based on multiple traits. Journal of Vegetation Science 16:533-540.
