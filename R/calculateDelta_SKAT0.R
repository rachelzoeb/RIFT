#' Function to calculate delta chi-square scores based on SKAT-O (Step 1 of RIFT)
#'
#' This function takes in y (phenotype), x (genotype matrix) and cov (covariate matrix)
#'
#' @param y outcome vector
#' @param geno.mat genotype matrix, requires no missing data
#' @param cov.mat covariate matrix, if applicable
#' @return A list object of results
#' @importFrom SKAT SKAT_Null_Model
#' @importFrom SKAT SKAT
#' @importFrom dplyr arrange
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
calculateDelta_SKAT0 <- function(y, geno.mat, cov.mat = NULL){

  ## Format input
  SNPS <- colnames(geno.mat)

  if (is.null(cov.mat) == TRUE){

    obj <- SKAT_Null_Model(y ~ 1, out_type = "D", Adjustment = FALSE)

  } else {

    obj <- SKAT_Null_Model(y ~ as.matrix(cov.mat), out_type = "D", Adjustment = FALSE)

  }

  ## Create dataframe to hold results
  res <- NULL

  ## Loop through to remove one SNPS and re-run SKAT, saving results in iterative_results dataframe
  for (idx in c(0:dim(geno.mat)[2])){

    ## Set up genotype data to remove 0 or 1 snps
    if (idx == 0){

      GENO_sub <- geno.mat
      SNP_excluded <- "NONE"

      ## Run SKATO to find optimal parameter
      out_skato <- SKAT(GENO_sub, obj, method = "optimal.adj")
      rho_est <- out_skato$param$rho_est
      pval_rho <- out_skato$param$minp

    } else {

      ## print which SNP has been calculating
      print(paste0("Computing p-value for SNP #", idx))
      GENO_sub <- geno.mat[,-idx]
      SNP_excluded <- SNPS[idx]

      ## Run SKATO with rho from entire data
      out_skato <- SKAT(GENO_sub, obj, r.corr = rho_est)
      pval_rho <- out_skato$p.value


    }


    ## Save results to dataframe
    res <- rbind(data.frame(SNP_idx = idx, SNP_excluded = SNP_excluded, p.value = pval_rho, rho_est = rho_est), res)

  }

  ## Calculate Chi-Square statistics
  print("Computing chi-square statistics")
  res$chisq <- sapply(res$p.value, FUN = function(x){return(qchisq(x, df = 1, lower.tail = FALSE))})

  ## Calculate delta Chi-Square statistics
  print("Computing delta chi-square statistics")
  res$delta <- sapply(res$chisq, FUN = function(x){return(x - res$chisq[res$SNP_excluded == "NONE"])})

  ## Resort results
  res <- arrange(res, SNP_idx)

  ## Return results
  return(list(ind.stats = res))

}
