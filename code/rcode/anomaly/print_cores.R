# ------------------------------------------------------------------------------
### Libraries and functions
library(doParallel)
library(foreach)
## Set up parallel backend
detectCores()
makeCluster(4)
