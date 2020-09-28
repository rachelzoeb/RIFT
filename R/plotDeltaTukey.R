#' Function to calculate delta chi-square scores based on SKAT-O (Step 1 of RIFT)
#'
#' This function takes in y (phenotype), x (genotype matrix) and cov (covariate matrix)
#'
#' @param results output from calculateDelta_SKAT0()
#' @return Plot of delta chi-square score by variant
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
#' @export

plotDeltaTukey <- function(results){

  plot(x = results$ind.stats$SNP_idx,
       y = results$ind.stats$delta,
       type ="n", xaxt='n', ann=FALSE)
  points(x = results$ind.stats$SNP_idx,
         y = results$ind.stats$delta,
         pch = 19, col = c("black", "red")[as.factor(results$ind.stats$tukey.mild)], xlab = "Variant", ylab = "Delta Chi-Square")
  abline(h = 0, col = "red")
  axis(1, at = results$ind.stats$SNP_idx, labels = results$ind.stats$SNP_excluded, las = 2)

}

