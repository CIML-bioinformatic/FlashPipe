## ####################################################
## This script aims to load and send back plot and data
## from the TCR trust4 output.
## ####################################################

## @knitr computeTrust4_tcr

# If no BCR previously, but TCR then we launch the rest of the code
if (!BCR && TCR){
  cat("\n  \n")
  cat("## Trust4 {.tabset .tab-fade} \n\n") 
}

