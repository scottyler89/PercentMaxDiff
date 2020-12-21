## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----two_batches, message=FALSE-----------------------------------------------
library("PercentMaxDiff")
## generate a vector that represents batch
batch <- rep(seq_len(2),each=100)
## generate a vector that represents clusters
clusters <- rep(rep(seq_len(2),each=50),2)
## run pmd
pmd_res <- pmd(batch, clusters)

## ----three_batches, message=FALSE---------------------------------------------
## generate a vector that represents batch
batch <- rep(seq_len(3),each=100)
## generate a vector that represents clusters
clusters <- rep(rep(seq_len(2),each=50),3)
## run pmd
pmd_res <- pmd(batch, clusters)

## ----hist_of_null, message=FALSE----------------------------------------------
hist(pmd_res$pmd_null, breaks = 15, main = paste("lambda:", pmd_res$pmd_null_lambda))

## ----p_val, message=FALSE-----------------------------------------------------
pmd_res$p.value

## ----make_three_, message=FALSE-----------------------------------------------
batch <- rep(seq_len(3),each=1000)
clust<- c(rep(seq_len(4),each=250),rep(seq_len(4),each=250))
clust<- c(clust, c(rep(seq_len(4),each=250)+4))

## ----look_at_three, message=FALSE---------------------------------------------
unique(cbind(batch, clust))

## ----run_three_batches, message=FALSE-----------------------------------------
## This will take a few seconds because of the simulations.
pmd_res_3_perfect <- pmd(batch, clust)
pmd_res_3_perfect$pmd_raw
pmd_res_3_perfect$pmd
pmd_res_3_perfect$p.value

## ----run_three_batches_realistic, message=FALSE-------------------------------
## generate the random cluster labels based on the probability vectors fed into the first argument
batch_size <- 500
batch1_clusters <- get_random_sample_cluster(c(rep(.25,4),rep(0,4)),batch_size)
batch2_clusters <- get_random_sample_cluster(c(rep(.25,4),rep(0,4)),batch_size)
batch3_clusters <- get_random_sample_cluster(c(rep(0,4),rep(.25,4)),batch_size)
## make the batch labels
batch1_labels <- rep("batch1",batch_size)
batch2_labels <- rep("batch2",batch_size)
batch3_labels <- rep("batch3",batch_size)
## collate the final batch and cluster vectors
realistic_batch <- c(batch1_labels,
                     batch2_labels,
                     batch3_labels)
realistic_clust <- c(batch1_clusters,
                     batch2_clusters,
                     batch3_clusters)
# run the pmd function!
pmd_res_3_realistic <- pmd(realistic_batch, realistic_clust)
pmd_res_3_realistic$pmd_raw
pmd_res_3_realistic$pmd
pmd_res_3_realistic$p.value

## ----run_post_hoc, message=FALSE----------------------------------------------
## first using batch and cluster labels
batch_labels <- paste("batch",rep(seq_len(3),each=100))
cluster_labels <- paste("cluster",c(rep(rep(seq_len(2),each=50),2),rep(seq_len(2),each=50)+1))
pmd_res <- pmd(batch_labels, cluster_labels)
## This is what the cluster/batch contingency table looks like:
print(pmd_res$cont_table)
## Which makes the percent maximum difference ~ 1/3rd
print(pmd_res$pmd_raw)
## in a real-world sceanrio, the data wouldn't be this clean beacuse of noise via
## Poisson sampling so the corrected PMD is actually a bit lower:
print(pmd_res$pmd)
## Using the pmd_posthoc function, we'll be able to figure out exactly which batch(es) 
## are causing the 1/3rd asymmetry
pmd_pairwise <- pmd_posthoc(pmd_res)
print(pmd_pairwise)
## looking at the pmd_pairwise$pmd_raw_table, we can clearly see that batch 3 is half different from both batches 1 and 2.
print(pmd_pairwise$pmd_raw_table)
## similar to above, because in a real-world example the data wouldn't be so clean, the adjusted PMDs are a bit lower
print(pmd_pairwise$pmd_table)
## Now if we look at the pairwise p-values, we see that batch3 is the clear outlier and batches 1 and 2 are very similar
print(pmd_pairwise$p.value_table)
## The above are the nominal p-values. Below are the adjusted p-values (BH correction by default)
print(pmd_pairwise$p.value_table_adjusted)

