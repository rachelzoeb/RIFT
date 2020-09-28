#' Function to calculate Tukey Fences and call influential variants
#'
#' This function takes in the results list object from calculateDelta_SKAT0()
#'
#' @param results output from calculateDelta_SKAT0()
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

calculateTukeyFences <- function(results){

  # Extract full parameter
  full_delta <- results$ind.stats$delta[results$ind.stats$SNP_excluded == "NONE"]

  # Extract jackknife etsimates
  ind_delta <- results$ind.stats$delta[results$ind.stats$SNP_excluded != "NONE"]
  n <- length(ind_delta)

  # Calculate IQR and quantiles
  lowerq = quantile(ind_delta)[2]
  upperq = quantile(ind_delta)[4]
  iqr = upperq - lowerq

  # Compute bounds for inner fences
  mild.threshold.upper = (iqr * 1.5) + upperq
  mild.threshold.lower = lowerq - (iqr * 1.5)

  # Compute bounds for outer fences
  extreme.threshold.upper = (iqr * 3) + upperq
  extreme.threshold.lower = lowerq - (iqr * 3)

  # Save thresholds
  results$iqr <- list(iqr = iqr, lowerq = lowerq, upperq = upperq,
                      mild.threshold.lower = mild.threshold.lower, mild.threshold.upper = mild.threshold.upper,
                      extreme.threshold.lower = extreme.threshold.lower, extreme.threshold.upper = extreme.threshold.upper)

  # Identify outlier
  results$ind.stats$fence <- sapply(results$ind.stats$delta, function(x){
    if ((x < extreme.threshold.lower) | (x > extreme.threshold.upper)) {
      return("extreme")
    } else if ((x < mild.threshold.lower) | (x > mild.threshold.upper)) {
      return("mild")
    } else {
      return("none")
    }
  })

  # Make the fence variable a factor
  results$ind.stats$fence <- factor(results$ind.stats$fence, levels = c("extreme", "mild", "none"))

  # Make test vectors
  # results$ind.stats$tukey.mild <- ifelse(results$ind.stats$fence %in% c("mild", "extreme"), TRUE, FALSE)
  results$ind.stats$tukey.mild <- ifelse((results$ind.stats$fence %in% c("mild", "extreme")) & (results$ind.stats$delta < 0), TRUE, FALSE)
  results$ind.stats$tukey.extreme <- ifelse((results$ind.stats$fence %in% c("extreme")) & (results$ind.stats$delta < 0), TRUE, FALSE)

  # Return thresholds
  return(results)

}
