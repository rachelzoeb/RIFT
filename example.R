#' Code to run RIFT from plink recode genotype file

# read in example dataset
raw_dat <- read.table("example_data.txt")

# Outcome should be binary, coded as 0/1
y <- raw_dat$PHENOTYPE - 1

# Create genotype matrix
geno <- as.matrix(raw_dat[,-1])

# Calculate delta chi-square
results <- calculateDelta_SKAT0(y, geno.mat = geno)

# Call influential variants using Tukey Fences
results <- calculateTukeyFences(results)

# Call influential variants using MAD
results <- calculateMAD(results)

# Print IV
results$ind.stats[results$ind.stats$tukey.mild,]

#' Code to run RF and viRF methods from plink recode genotype file
#' Must have randomForest and viRandomForests packages downloade
results <- calculateRF(y, geno.mat = geno)

# Print IV from RF methods
results$rf[results$rf$MeanDecreaseAccuracy.tukey,]
results$virf[results$virf$MeanDecreaseAccuracy.tukey,]
