# BiTSCR

## Requirement
R version at least 4.1.1 (2021-08-10)

## Input data files
1. Two gene node covariate datasets from the first species and second species respectively
2. An edge matrix that includes the orthology information between these two species

## Output data
The frequency matrix $\bar{M}$ indicate how frequent two genes are connected in the simulations.

## Tutorial and notes
- Please see the R markdown "testbitsc.Rmd" in the vignettes for examples to install, run, and visualize the output cluster. 
- Please allow for at least 10 minutes for the algorithm to generate result for the example. 
- Please make sure R can take at least 7 GB space to store intermediate matrix during calculation. 
- Please expect the warnings such as "Warning in do.call(.Call, args = dot_call_args) : only 10 eigenvalue(s) converged, less than k = 15", which are caused by the inconsistency of the $K_0$ and the actual number of clusters may be generated when $K_0$ is large. Note that when the warnings exist, the performance and accuracy of the implementation are not affected.
