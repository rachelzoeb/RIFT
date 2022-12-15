######################################################################################################
#     simulate_rare_genotype_data.R
#
#     This code is to simulate rare variant genotype data using simulated haplotypes from COSI
#     This code was set up on a cluster, will need to change the paths accordingly
#
######################################################################################################

########################################################################################
#### Load functions
########################################################################################

#--------------- sample rare variant region fromm haplotype data  ---------------------

## Function for randomly sample rare variant region from haplotype data
sampleRareGeneFromHaplotypes <- function(haps, pos, size.kb = 3){

  # Randomly select variant
  ### Subset to only observed rare variants
  all.maf <- pos$FREQ1
  rare.vars <- which((0 < all.maf) & (all.maf < 0.03))
  rare.pos <- pos$CHROM_POS[rare.vars]

  start_variant <- sample(rare.vars, 1)
  start_pos <- pos$CHROM_POS[start_variant]
  end_pos <- start_pos + size.kb*1000

  rare.vars.sub <- rare.vars[((start_pos <= rare.pos) & (rare.pos <= end_pos))]

  # Return matrix of 50 variants
  return(haps[,rare.vars.sub])

}


#--------------- generate genotype data ---------------------

## function for simulating pseudo sample from population
### parameter setting
#   - population: population's haplotype data
#   - sample.size: sample size of simulated samples
### output
#   - genotype: simulated genotype data

genotype <- function(population, sample.size){

  id <- c(1:nrow(population))
  h1 <- sample(id, size = sample.size, replace = T)
  h2 <- sample(id, size = sample.size, replace = T)
  genotype <- population[h1,]+population[h2,]
  row.names(genotype) <- NULL
  return(genotype)

}


#---------------------- generate phenotype -------------------

## function for generate phenotype
### phenotype: Y ~ Bernoulli(p), logit(p) = beta0 + beta*SNP
### parameter setting
#   - only.intercept: a logic value indicates whether using null model (T) or not (F)
#   - genotype.data: import genotype data (coding = 0, 1, 2) including column names (SNPs names)
#   - asso.variant.id: ID of causal SNPs.
#   - intercept: intercept of logistic regression
#   - odds.ratio: ORs of SNPs. The leangth must be equal to asso.CV+asso.RV or the length of asso.variant.id.
#           ORs of asso.CV are given first as specifying asso.CV & asso.RV.
#           ORs are given in order as specifying asso.variant.id.
#   - sample.size: sample size
### output
#   - asso.variants: ID of causal variants
#   - MAF: MAF for causal variants
#   - odds.ratio: OR for each causal variant
#   - prob: probability of affecting the disease
#   - disease.status: subjects' phenotype
#   - case.no: a table for presenting no. of cases and no. of controls

phenotype <- function(only.intercept = F, genotype.data, asso.variant.id, odds.ratio, intercept, sample.size){

  if(only.intercept == F){
    odds.ratio <- as.vector(odds.ratio)
    genotype.data <- as.matrix(genotype.data)
    asso.variants <- as.matrix(genotype.data[,asso.variant.id])

    # calculate MAF
    MAF <- colMeans(asso.variants)/2
    MAF_ALL <- colMeans(genotype.data)/2

    # generate phenotype
    y <- intercept + asso.variants%*%log(odds.ratio)
    p <- exp(y)/(1+exp(y))
    disease.status <- rbinom(length(p), size = 1, prob = p)

    # output
    res <- list()
    res$asso.variants <- colnames(asso.variants)
    res$MAF <- MAF
    res$MAF_ALL <- MAF_ALL
    res$odds.ratio <- odds.ratio
    res$prob <- p
    res$disease.status <- disease.status
    res$case.no <- table(disease.status)

  } else if(only.intercept == T){

    # generate phenotype
    y <- intercept
    p <- exp(y)/(1+exp(y))
    disease.status <- rbinom(sample.size, size = 1, prob = p)

    # output
    res <- list()
    res$prob <- p
    res$disease.status <- disease.status
    res$case.no <- table(disease.status)
  }

  return(res)
}

########################################################################################
#### Set up parameters
########################################################################################

## Set Seed from job index - this is cluster specific code
# inparms <- commandArgs(trailingOnly = TRUE)
# # the value gets pulled in as a vector of strings, so have to
# # change to numeric for use in the R script
# row_id <- as.numeric(inparms[1])

## Set up directories
cosi_data_directory <- ".../cosi/bestfit/"  ## this is where the COSI output is located
output_directory <- ".../simulated_data_dir/" # this is where simulated datasets will be stored

## Read in HAPS and POS file
haps <- read.table(file = paste0(cosi_data_directory, "out.hap-1"))
pos <- read.table(file = paste0(cosi_data_directory, "out.pos-1"), header = TRUE)

### Remove first two columns and convert to 1 for derived and 0 for ancestral
haps <- haps[,-c(1:2)]
haps <- apply(haps, 2, function(x){ifelse(x == 2, 0, x)})
## Format SNP name
#pos$SNP_name = paste0("SNP.", pos$SNP)
pos$SNP_name = paste0("C", pos$CHROM, ".", pos$CHROM_POS)
colnames(haps) <- pos$SNP_name

## Simulation number and seed
set.seed(1) # used for manuscript
num.regions <- 50
num.sims <- 500

## set up simulation parameters
## Set up parameters
prevalence = 0.05
aa = log(prevalence/(1-prevalence))
prop.assoc.variant = 0.10  ## proportion of variants under the alternative
cc = 0.4 ## this is the effect size parameter c
cs.no = 1000 ## number of cases
cn.no = 1000 ## number of controls

## start loop to simulate data
for (region in 1:num.regions){

  ## This simulation approach treats the population as the source of the maf and definition of rare, here we set MAF cutoff and region size
  haps.region <- sampleRareGeneFromHaplotypes(haps, pos, maf.cutoff = 0.03, size.kb = 3)

  ## Use above parameters to set up associated variants
  num.assoc.variant = round(dim(haps.region)[2]*prop.assoc.variant)
  MAF = colMeans(haps.region)
  select.SNP <- sample(colnames(haps.region), num.assoc.variant, replace = F)
  beta_k <- cc*abs(log10(MAF[select.SNP]))

  # save region information
  region_dat <- data.frame(SNP = colnames(haps.region),
                           MAF,
                           select.SNP = (colnames(haps.region) %in% select.SNP),
                           odds.ratio = exp(beta_k)[match(colnames(haps.region), select.SNP)])
  write.table(region_dat, file = paste0(output_directory,"region", region, ".info"))

  ## Set up simulation matrix
  for (i in 1:num.sims){

      ## Sample genotypes
      sim.genotype <- genotype(population = haps.region, sample.size = 40000)  ## This must be large enough to simulate a large enough samples of case to randomly draw from

      ## Generate phenotypes
      sim.phenotype <- phenotype(genotype.data = sim.genotype, asso.variant.id = select.SNP, odds.ratio = exp(beta_k), intercept = aa, sample.size = dim(sim.genotype)[1])

      ## Make equal # cases and controls
      cn <- which(sim.phenotype$disease.status == 0)
      cs <- which(sim.phenotype$disease.status == 1)
      cn.id <- sample(cn, cn.no, replace = FALSE)
      cs.id <- sample(cs, cs.no, replace = FALSE)
      disease.status <- sim.phenotype$disease.status[c(cn.id, cs.id)]
      raw_dat <- data.frame(cbind(PHENOTYPE = (sim.phenotype$disease.status[c(cn.id, cs.id)] + 1),
                                  sim.genotype[c(cn.id, cs.id),]))

      # Subset data to only those variants observed
      observed <- which(colSums(raw_dat) != 0)
      raw_dat <- raw_dat[,observed]

      ## save raw data
      write.table(raw_dat, file = paste0(output_directory,"region", region, "_sim", i, ".data"))

     }

}
