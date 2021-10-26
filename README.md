# funredun
Calculates community functional redundancy using species abundance and functional traits

## Install
Installation in R requires <a href="https://cran.r-project.org/package=devtools">devtools</a>. Requires <a href="https://cran.r-project.org/package=vegan">vegan</a> and <a href="https://cran.r-project.org/package=usedist">usedist</a> packages to run.
```
install_github("PlantEcology/funredun")
```

## Data Format
Function uses two data frames to calculate site functional redundancy spDat and funDat.

```
funredun(spDat, funDat, method = "bray", redund = TRUE)
```
Variable | Description 
------|-----
spDat | Data frame with rows as sites, columns as species, and elements as counts
funDat | Data frame with rows as species (same as spDat column names), columns as functional traits, elements as counts, measures, binary, etc.
method | Available options include "bray", "gower", and "altGower". See <a href="https://cran.r-project.org/web/packages/vegan/index.html">vegan::vegdist</a> for details. Default is Bray-Curtis dissimilarity.
redund | Redundancy calculation as difference from Simpson's D (R = D - Q) or uniqueness (U = Q/D). Default is difference (TRUE).

## Output
A data frames is produced with rows as sites (from spDat) and a column of functional redundancy. Elements are calculated as subtraction of functional diversity (Q) from Simpson's D. Q is sum of the products of Bray-Curtis dissimilarity between species i and j, proportion of species i at the site, and proportion of species j at the site.
