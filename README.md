# funredun
Calculates community functional redundancy using species abundance and functional traits

# Install
Installation in R requires <a href="https://cran.r-project.org/package=devtools">devtools</a>. Requires <a href="https://cran.r-project.org/package=vegan">vegan</a> and <a href="https://cran.r-project.org/package=usedist">usedist</a> packages to run.
```
install_github("PlantEcology/funredun")
```

## Data Format
Function uses two data frames to calculate site functional redundancy spDat and funDat.

spDat = Data frame with rows as sites, columns as species, and elements as counts

funDat = Data frame with rows as species (same as spDat column names), columns as traits, elements as counts, measures, binary, etc.

## Output
A data frames is produced with rows as sites (from spDat) and a column of functional redundancy. Elements are calculated as subtraction of functional diversity (Q) from Simpson's D. Q is sum of the products of Bray-Curtis dissimilarity between species i and j, proportion of species i at the site, and proportion of species j at the site.
