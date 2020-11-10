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
## batch 1 vs batch 2 
pmd_res_1v2 <- pmd(c(batch1_labels, batch2_labels),
               c(batch1_clusters, batch2_clusters))
pmd_res_1v2$pmd_raw
pmd_res_1v2$pmd
pmd_res_1v2$p.value

## batch 1 vs batch 3
pmd_res_1v3 <- pmd(c(batch1_labels, batch3_labels),
               c(batch1_clusters, batch3_clusters))
pmd_res_1v3$pmd_raw
pmd_res_1v3$pmd
pmd_res_1v3$p.value

## batch 2 vs batch 3
pmd_res_2v3 <- pmd(c(batch2_labels, batch3_labels),
               c(batch2_clusters, batch3_clusters))
pmd_res_2v3$pmd_raw
pmd_res_2v3$pmd
pmd_res_2v3$p.value

