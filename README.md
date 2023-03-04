# fupcompleteness

Collection of functions to evaluate follow-up completeness. The functions are dependend on the Icens package (1) and interval package (2) and these should be installed first by entering: 

```{r} 
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("Icens")
install.packages("interval")
library(interval)
```



1: R. Gentleman and Alain Vandal (2020). Icens: NPMLE for Censored and Truncated Data. R package, version 1.62.0.
2: Michael P. Fay, Pamela A. Shaw (2010). Exact and Asymptotic Weighted Logrank Tests for Interval, Censored Data: The interval R Package. Journal of Statistical Software, 36(2), 1-34. doi:10.18637/jss.v036.i02


