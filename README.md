
## RIFT package v2.0

We present a novel method for prioritization of rare variants within a
given set of variants after the set of variants is found to be
significant using aggregate testing methods. Building on the rich
outlier-detection statistical literature, we present a computationally
efficient approach to be applied following identification of a set of
variants that is agnostic to putative function. Our approach, which we
refer to as RIFT for Rare Variant Influential Filtering Tool, leverages
the influence of the variant on the aggregate test of association by
quantifying the change in the aggregate test when that variant is
removed. It is particularly well suited for rare and uncommon variants,
the most common applications of aggregate tests, but is applicable to
aggregate testing of variants of all frequencies. When applied to a
significant set of rare/uncommon variants, RIFT provides a scheme for
quantifying the contribution of an individual variant to the overall
association signal, while adjusting for covariates. This method also
provides a quantitative measure by which to rank variants for further
investigation and several visualizations to aid in evaluation of a
region of interest.

v2.0 (Updated September 22, 2021) - includes functions to apply randomForest (Breiman, 2001) and viRandomForest (Liu & Zhao, 2017) methods to a set of variants.  This requires having the packages randomForest and viRandomForest installed.  

* randomForest can be installed via CRAN (https://cran.r-project.org/web/packages/randomForest/index.html)
* viRandomForest can be downloaded and installed from https://zhaocenter.org/softwares.html

We have also included a file containing the code used to simulate datasets of rare variant genetic data.  The file simulate_rare_genotype_data.R requires haplotypes simulated by COSI, specification of parameters can be modified based on the desired simulation parameters (size of region, number of cases/controls, proportion variants under the alternative, etc.).  

### Installation

``` r
install.packages("devtools")
library(devtools)
install_github("rachelzoeb/RIFT")
```

### Example code to run original RIFT method

``` r
library(RIFT)
data(exampleData)
y <- exampleData$PHENOTYPE - 1
geno <- as.matrix(exampleData[,-c(1,2)])
cov <- as.matrix(exampleData$SEX)
results <- calculateDelta_SKAT0(y, geno.mat = geno, cov.mat = cov)
results <- calculateTukeyFences(results)
results <- calculateMAD(results)
results$ind.stats[results$ind.stats$tukey.mild,]
plotDeltaTukey(results)
```

### Example code to run RF and vi-RF methods

``` r
library(RIFT)
data(exampleData)
y <- exampleData$PHENOTYPE - 1
geno <- as.matrix(exampleData[,-c(1,2)])
results <- calculateRF(y, geno.mat = geno)
results$rf[results$rf$MeanDecreaseAccuracy.tukey,]
results$virf[results$virf$MeanDecreaseAccuracy.tukey,]
```

