---
title: "The PercentMaxDiff user's guide"
author: "Scott Tyler"
date: "11/05/2020"
output: 
  prettydoc::html_pretty:
    theme: cayman
    highlight: vignette
vignette: >
  %\VignetteEncoding{UTF-8}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteIndexEntry{"The PercentMaxDiff user's guide"}
---



# Introduction

Percent Maximum Difference provides a metric that enables the measurement of how similar (PMD near 0) or different (PMD near 1) two batches are based just on clustering results. This can be viewed as an alternative to a Chi-squared type analysis based on a contingency table in which the batches are noted in columns, and the abundance of cells in each cluster is noted in rows. PMD advances on Chi-squared and Fishers exact tests because its results are actually not dependent on the number of clusters found. It does this by first calculating the hypothetical maximum possible asymmetry between your batches. (i.e.: if each batch was only comprised of clusters that only ever appeared in that batch and there were no clusters that appeared in more than one batch). Then it calculates the observed asymmetry in your real results, then returning the percentage of the observed asymmetry relative to the maximum possible asymmetry. This gives a lower bound of zero, where the relative abundance of each cluster was spot-on the same between all of the batches. This also gives an upper bound of one, in which the observed clustering results were fully asymmetric across batch; or in other words, each cluster was fully batch specific.


# Installation

`devtools::github_install('scottyler89/PercentMaxDiff')`

Although, this package has also been submitted to Bioconductor; we'll update this section once it's accpeted there!


# Examples
## Running PMD on batch/clustering results


```r
library("PercentMaxDiff")
## generate a vector that represents batch
batch <- rep(seq_len(2),each=100)
## generate a vector that represents clusters
clusters <- rep(rep(seq_len(2),each=50),2)
## run pmd
pmd_res <- pmd(batch, clusters)
```

You can also run it on more than two batches! It actually works in n-dimentions both for batches or clusters.


```r
## generate a vector that represents batch
batch <- rep(seq_len(3),each=100)
## generate a vector that represents clusters
clusters <- rep(rep(seq_len(2),each=50),3)
## run pmd
pmd_res <- pmd(batch, clusters)
```

## Looking at the PMD results

The resultant object is a list with several useful metrics that describe how similar batches are to each other based on their clustering results.

* _*cont_table*_: The contingency table of clusters (rows) and batches (columns)

* _*expected*_: The expected matrix, as with a Chi-square test.

*  _*pmd_null*_: Simulations of the null distribution of PMDs using the observed global percentage of cluster abundances across all batches with the observed batch sizes to match the input. Note that this will approach zero, but due to random sampling, typically will never actually get there. That's what makes having this null background useful.

* _*pmd_null_lambda*_: The lambda value of the null distributions Poisson fit. 

* _*p.value*_: Significance for whether or not batches are different in their cellular composition. This is determined through the generating *num_sim* null distributions, and emperically measuring the number of times the observed PMD was greater than the PMDs generated from the null distribution. Low p-values indicate that the batches are indeed different from each other.

* _*pmd*): The percent maximum difference (pmd) of the input dataset.



## Looking at significance

Part of the pmd function is running a null background. This mimics the numbers and relative abundnance of your input data, but under the null-hypothesis that your batches were not actually different. These simulations yield a vector of null pmd calculations. You can look at that like so:


```r
hist(pmd_res$pmd_null)
```

![plot of chunk hist_of_null](figure/hist_of_null-1.png)

Using this null background, an empirical p-value is calculated:



```r
pmd_res$p.value
```

```
## [1] 0.996
```

## comparing two different 'batch correction' or normalization approaches

Lets say you had done batch integration & clustering using two different algoirthms. There were three datasets - two were very similar, and one was very different. 



```r
batch <- rep(seq_len(3),each=1000)
clust<- c(rep(seq_len(4),each=250),rep(seq_len(4),each=250))
clust<- c(clust, c(rep(seq_len(4),each=250)+4))
```

Let's take a look at what that looks like

```
unique(cbind(batch, clust))
      batch clust
# [1,]     1     1
# [2,]     1     2
# [3,]     1     3
# [4,]     1     4
# [5,]     2     1
# [6,]     2     2
# [7,]     2     3
# [8,]     2     4
# [9,]     3     5
#[10,]     3     6
#[11,]     3     7
#[12,]     3     8

```

So this is a situation where we had three batches, 4 clusters that appeared in batches 1 and 2 evenly, but batch 3 had 4 clusters that were specific to it, and had none of the clusters that appear in the first two batches.

Let's run the PMD on it to quantify how similar these batches are overall.

```
## This will take a few minutes because of the 10000 simulations.
pmd_res <- pmd(batch, clust)
print(pmd_res$pmd)
print(pmd_res$p.value)
# > print(pmd_res$p.value)
# [1] 0
# > print(pmd_res$pmd)
# [1] 0.6666667
```

The PMD in this case is .6666. That's because 2/3rds of the cells come from clusters that are shared across batches, while the remaining 1/3rd is from batch 3 that was completely different from the others. The P-value is 0 (or in this case just more significant that the 10000 null simulations run in which there was no pattern by batch). 

Note that you could also run these assays on each batch pairwise if you want to do something like a pairwise post-hoc.



