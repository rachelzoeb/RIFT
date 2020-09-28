#' Function to calculate MAD and call influential variants
#'
#' This function takes in the results list object from calculateDelta_SKAT0()
#'
#' @param results output from calculateDelta_SKAT0()
#' @param cutoff MAD cutoff for outliers, recommended cutoff of 3
#' @param constant MAD constant, recommended constant for normal distribution of 1.4826
#' @return A list object of results
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

# Function to calculate MAD
calculateMAD <- function(results, cutoff = 3, constant = 1.4826){

  # Extract full parameter
  full_delta <- results$ind.stats$delta[results$ind.stats$SNP_excluded == "NONE"]

  # Extract jackknife etsimates
  ind_delta <- results$ind.stats$delta[results$ind.stats$SNP_excluded != "NONE"]
  n <- length(ind_delta)

  # Calculate MAD
  results$mad$mad <- mad(ind_delta, constant = constant)
  results$mad$cutoff <- cutoff
  results$ind.stats$mad.dev <- c(0, (abs(ind_delta - median(ind_delta))/results$mad$mad))

  # Determine if deviation is larger than cutoff
  results$ind.stats$mad.outlier <- ifelse(results$ind.stats$mad.dev > cutoff & results$ind.stats$delta < 0, TRUE, FALSE)

  # Resort results
  results$ind.stats <- arrange(results$ind.stats, SNP_idx)

  # Return plot
  return(results)

}
