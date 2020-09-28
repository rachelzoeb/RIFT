#' @title Example dataset for RIFT
#'
#' @description A data set with the genotypes for 2000 subjects (1000 cases and 1000 controls) including sex as a covariate
#'
#' @docType data
#'
#' @usage data(exampleData)
#'
#' @format A R data frame with 2000 rows and 45 variables
#' \describe{
#'   \item{PHENOTYPE}{Controls (1), Cases (2)}
#'   \item{SEX}{Males (1), Females (2)}
#'   \item{genotype data}{Numeric genotype matrix with each row as a different individual and each column as a separate gene/snp. Each genotype should be coded as 0, 1, 2}
#' }
#' @source <https://www.github.com/rachelzoeb/RIFT>
#'
#' @examples
#' data(exampleData)
#' y <- exampleData$PHENOTYPE - 1
#' geno <- as.matrix(exampleData[,-c(1,2)])
#' cov <- as.matrix(exampleData$SEX)
#' results <- calculateDelta_SKAT0(y, geno.mat = geno, cov.mat = cov)
#' results <- calculateTukeyFences(results)
#' results <- calculateMAD(results)
#' results$ind.stats[results$ind.stats$tukey.mild,]
#' plotDeltaTukey(results)
"exampleData"
