# RIFT
Rare Variant Influential Filtering Tool (RIFT)

We present a novel method for prioritization of rare variants within a given set of variants after the set of variants is found to be significant using aggregate testing methods. Building on the rich outlier-detection statistical literature, we present a computationally efficient approach to be applied following identification of a set of variants that is agnostic to putative function.  Our approach, which we refer to as RIFT for Rare Variant Influential Filtering Tool, leverages the influence of the variant on the aggregate test of association by quantifying the change in the aggregate test when that variant is removed. It is particularly well suited for rare and uncommon variants, the most common applications of aggregate tests, but is applicable to aggregate testing of variants of all frequencies. When applied to a significant set of rare/uncommon variants, RIFT provides a scheme for quantifying the contribution of an individual variant to the overall association signal, while adjusting for covariates. This method also provides a quantitative measure by which to rank variants for further investigation and several visualizations to aid in evaluation of a region of interest. 

To download and install
install.packages("devtools")
library(devtools)
install_github("rachelzoeb/RIFT")
