#' Function to apply the standard random forest to genotype data
#'
#' This function takes in y (phenotype), x (genotype matrix).  This function does not allow adjusting for covariates.  This function includes the Tukey fences within the function.
#'
#' @param y outcome vector
#' @param geno.mat genotype matrix, requires no missing data
#' @return A list object of RF and vi-RF results
#' @importFrom randomForest viRandomForests
#' @importFrom dplyr arrange
#' @examples
#' data(exampleData)
#' y <- exampleData$PHENOTYPE - 1
#' geno <- as.matrix(exampleData[,-c(1,2)])
#' results <- calculateRF(y, geno.mat = geno)
#' results$rf[results$rf$MeanDecreaseAccuracy.tukey,]
#' @export
calculateRF <- function(y, geno.mat){

  if (requireNamespace("randomForest", quietly = TRUE) & requireNamespace("viRandomForests", quietly = TRUE)) {

    ## Format input
    SNPS <- colnames(geno.mat)

    ## format outcome
    y <- factor(y, levels = c(0,1))
    vars <- c("y", SNPS)
    form <- formula("y ~ .")

    ## Run standard random forest function
    rf_default <- randomForest(form, data = data.frame(y, geno.mat), importance = TRUE, proximity = FALSE)
    rf_default <- data.frame(rf_default$importance)

    ## Run default viRandomForest function
    virf_default <- viRandomForests(form, data = data.frame(y, geno.mat), importance = TRUE, proximity = FALSE)
    virf_default <- data.frame(virf_default$importance)

    #####**********  Function to call IV using Tukey Mild criteria
    tukey_mild_outlier <- function(x){

      # Extract jackknife etsimates
      n <- length(x)

      # Calculate IQR and quantiles
      lowerq = quantile(x)[2]
      upperq = quantile(x)[4]
      iqr = upperq - lowerq

      # Compute bounds for inner fences
      mild.threshold.upper = (iqr * 1.5) + upperq
      mild.threshold.lower = lowerq - (iqr * 1.5)

      # Compute bounds for outer fences
      extreme.threshold.upper = (iqr * 3) + upperq
      extreme.threshold.lower = lowerq - (iqr * 3)

      # Identify outlier
      fence <- sapply(x, function(y){
        if ((y < extreme.threshold.lower) | (y > extreme.threshold.upper)) {
          return("extreme")
        } else if ((y < mild.threshold.lower) | (y > mild.threshold.upper)) {
          return("mild")
        } else {
          return("none")
        }
      })

      # Make the fence variable a factor
      fence <- factor(fence, levels = c("extreme", "mild", "none"))

      # Make test vectors
      return(ifelse((fence %in% c("mild", "extreme")), TRUE, FALSE))

    }

    ## Apply tukey to RF results
    rf_default$MeanDecreaseAccuracy.tukey <- tukey_mild_outlier(rf_default$MeanDecreaseAccuracy)
    rf_default$MeanDecreaseGini.tukey <- tukey_mild_outlier(rf_default$MeanDecreaseGini)

    ## Apply tukey to vi-RF results
    virf_default$MeanDecreaseAccuracy.tukey <- tukey_mild_outlier(virf_default$MeanDecreaseAccuracy)
    virf_default$MeanDecreaseGini.tukey <- tukey_mild_outlier(virf_default$MeanDecreaseGini)

    ## Return results
    return(list(rf = rf_default,
                virf = virf_default))

      } else {

    # provide warning to install packaged
    stop("Please install randomForest and viRandomForests to use this function")
  }


}
