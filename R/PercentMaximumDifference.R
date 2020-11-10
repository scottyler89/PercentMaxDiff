#' get_cont_table
#' @description \code{get_cont_table} Gets a congingency table that quantifies two factors against each other, assuming that they are in the same order.
#' @param x a vector that can be interpreted as a factor. Should be clusters
#' @param y a vector that can be interpreted as a factor. Should be batches
#' @return a contingency matrix with x quantified in rows, y quantified in columns
#' @examples
#'    x <- rep( c( rep(1, 50), rep(2, 50) ), 2)
#'    y <- c(rep(1, 100), rep(2, 100))
#'    cont_table <- get_cont_table(x, y)
#' @name get_cont_table
#' @export
get_cont_table<-function(x, y){
    ## x and y are integeric labels
    x<-as.integer(as.factor(x))
    y<-as.integer(as.factor(y))
    if (length(x) != length(y)){
        print("x and y were different lengths!")
        return()
    }
    cont_mat <- matrix(ncol=length(unique(y)), nrow = length(unique(x)))
    cont_mat[,] <- 0 
    for (i in seq(from = 1, to = length(x)) ) {
        cont_mat[x[i],y[i]] = cont_mat[x[i],y[i]] + 1
    }
    return(cont_mat)
}




#' get_percent_max_resid
#' @description \code{get_percent_max_resid} Takes in the observed and expected matrix (as with \code{chisq.test}) and uses it to calculate the PMD.
#' @param observed observed matrix (i.e. contingency table)
#' @param expected expected matrix under the null hypothesis that there are no differences
#' @return the percentage of asymmetry relative to the maximum possible asymmetry.
#' @examples
#'    dummy_mat <- matrix(c(500, 250, 100, 1000, 500, 200),nrow = 3, ncol=2)
#'    get_percent_max_resid(dummy_mat, chisq.test(dummy_mat)$expected)
#' @importFrom stats chisq.test
#' @name get_percent_max_resid
#' @export
get_percent_max_resid<-function(observed, expected){
    num_diff <- sum(abs(observed - expected))
    num_dataset <- colSums(observed)
    baseline_mat <- matrix(ncol=length(num_dataset), nrow=length(num_dataset))
    baseline_mat[,] <- 0
    for (i in seq(from = 1, to =length(num_dataset)) ) {
        baseline_mat[i, i]<-num_dataset[i]
    }
    base_chi <- suppressWarnings(chisq.test(baseline_mat))
    max_diff <- sum(abs(base_chi$expected - baseline_mat))
    return(num_diff/max_diff)
}




#' get_percent_max_diff_from_tables
#' @description \code{get_percent_max_diff_from_tables} Calculates PMD from two input dataframes. The first data frame
#'                 has cell IDs in first column, and batch in the second column. The second
#'                 dataframe has the Cell IDs and then the clusters.
#' @param groups1_table a data frame in which the first column is the cell IDs, and second column is the cluster that cell belongs to.
#' @param groups2_table a data frame of the same format as above, but for the second dataset.
#' @return pmd A list object containing all of the same results as the \code{pmd} function
#' @examples 
#'       groups1_table = data.frame(cell = paste("cell",1:10),
#'                                  batch = c(rep(1,5),rep(2,5)))
#'       groups2_table = data.frame(cell = paste("cell",1:10),
#'                                  clust = c(rep(1,5),rep(2,5)))
#'       get_percent_max_diff_from_tables(groups1_table, groups2_table)
#' @importFrom stats na.omit
#' @name get_percent_max_diff_from_tables
#' @export
get_percent_max_diff_from_tables<-function(groups1_table, groups2_table){
    ## groups1_table HAS to be the batch/dataset labels
    ## assumes the first column is the cell and second column is the group
    colnames(groups1_table)[1]<-"cell"
    colnames(groups2_table)[1]<-"cell"
    merged_group_table<-na.omit(merge(groups1_table,groups2_table,by="cell"))
    group1_labs<-as.factor(merged_group_table[,2])
    group2_labs<-as.factor(merged_group_table[,3])
    # print(unique(group1_labs))
    # print(unique(group2_labs))
    return(pmd(group1_labs,group2_labs))
}




#' get_random_sample_cluster
#' @description \code{get_random_sample_cluster} Takes in a probability vector & the number of samples to simulate; returns that number of samples, fitting the specified proportions.
#' @param prob_vect a vector of probabilities that sum up to 1, that represent the relative abundance of each clusters that we'll be sampling from.
#' @param num_samples the number of samples to simulate.
#' @return a vector of simulated cluster labels of length \code{num_samples} based on the probility vector that codes the relative abundance of each cluster.
#' @examples
#'    prob_vect <- c(.25, .25, .5)
#'    sim_clusters <- get_random_sample_cluster(prob_vect, 1000)
#' @importFrom stats runif
#' @name get_random_sample_cluster
#' @export
get_random_sample_cluster<-function(prob_vect, num_samples){
    cum_sum<-cumsum(prob_vect)
    temp_rand<-runif(num_samples)
    group_vect<-c()
    for (i in seq(from = 1, to = num_samples)) {
        cur_clust<-1
        while (temp_rand[i]>cum_sum[cur_clust]){
            cur_clust<-cur_clust+1
        }
        group_vect<-c(group_vect, cur_clust)
    }
    return(group_vect)
}




#' get_pmd_null_vect
#' @description \code{get_pmd_null_vect} Takes in a contingency matrix (rows: clusters x columns: batch), and gets a null distribution of expected PMDs if there was no asymmetry based on batch.
#' @param expected_mat A matrix that represents a contingency matrix of clusters (rows) and batch (cols).
#' @param num_sim for generating a null distribution in calculating a p-value, the number of simulations to run.
#' @return a vector of floating points that represent a null distribution that matches your input for percent maximum difference
#' @examples
#'    dummy_mat <- matrix(c(500, 250, 100, 1000, 500, 200),nrow = 3, ncol=2)
#'    pmd_null <- get_pmd_null_vect(dummy_mat)
#' @name get_pmd_null_vect
#' @export
get_pmd_null_vect<-function(expected_mat, num_sim = 10000){
    total_cells<-sum(expected_mat)
    num_clust<-dim(expected_mat)[1]
    num_batch<-dim(expected_mat)[2]
    ## number of cells from each batch
    cells_per_batch<-colSums(expected_mat)
    ## figure out the null distribution of the relative abundance of each cluster (in rows)
    cells_per_clust<-rowSums(expected_mat)
    ## the relative abundance of each cluster
    clust_probs<-cells_per_clust/sum(cells_per_clust)
    pmd_null_vect<-c()
    for (p in 1:num_sim){
        batch_vect<-c()
        clust_vect<-c()
        for (b in seq(from = 1, to = num_batch)){
            ## generate each batch's sampling of cells from the main probability vector
            batch_vect<-c(batch_vect,rep(b,cells_per_batch[b]))
            temp_clusts<-get_random_sample_cluster(clust_probs, cells_per_batch[b])
            clust_vect<-c(clust_vect, temp_clusts)
            # ## populate the matrix
            # for (cell in temp_clusts){
            #   temp_mat[cell,b]<-temp_mat[cell,b]+1
            # }
        }
        pmd_null_vect<-c(pmd_null_vect, get_percent_max_diff(batch_vect,clust_vect))
    }
    return(pmd_null_vect)
}




#' get_percent_max_diff
#' @description \code{get_percent_max_diff} Gets the Percent Maximum Difference (PMD) for clustering results relative to batch.
#' @param group1_labs a vector of labels for the batch that each sample (or cell for scRNAseq) came from. Can be a string or factor.
#' @param group2_labs  a vector of labels for each cluster that each sample belongs to. Can be a string or factor.
#' @return a floating point number of the percent maximum difference
#' @examples
#'    batch <- rep(seq_len(2),each=100)
#'    clusters <- rep(rep(seq_len(2),each=50),2)
#'    pmd_num <- get_percent_max_diff(batch, clusters)
#' @importFrom stats chisq.test
#' @name get_percent_max_diff
#' @export
get_percent_max_diff<-function(group1_labs,group2_labs){
    cont_table<-get_cont_table(as.factor(group2_labs),as.factor(group1_labs))
    chi_result <- suppressWarnings(chisq.test(cont_table))
    observed<-cont_table
    expected<-chi_result$expected
    temp_percent_max_difference<-get_percent_max_resid(observed,expected)
    return(temp_percent_max_difference)
}




#' Percent Maximum Difference from a contingency table
#' @description \code{pmd} Percent Maximum Difference provides a metric that enables the measurement
#'                         of how similar (PMD near 0) or different (PMD near 1) two batches are
#'                         based just on clustering results. This can be viewed as an alternative
#'                         to a Chi-squared type analysis based on a contingency table in which 
#'                         the batches are noted in columns, and the abundance of cells in each cluster
#'                         is noted in rows. PMD advances on Chi-squared and Fishers exact tests
#'                         because its results are actually not dependent on the number of clusters found.
#'                         It does this by first calculating the hypothetical maximum possible asymmetry 
#'                         between your batches. (i.e.: if each batch was only comprised of clusters 
#'                         that only ever appeared in that batch and there were no clusters that appeared 
#'                         in more than one batch). Then it calculates the observed asymmetry in your 
#'                         real results, then returning the percentage of the observed asymmetry
#'                         relative to the maximum possible asymmetry. This gives a lower bound of zero,
#'                         where the relative abundance of each cluster was spot-on the same between
#'                         all of the batches. This also gives an upper bound of one, in which 
#'                         the observed clustering results were fully asymmetric across batch; 
#'                         or in other words, each cluster was fully batch specific.
#' @param cont_table As an alternative to passing in batch and cluster labels, you can also just pass in the contingency table.
#'                   Note that the clusters MUST be in the rows, and the batches must be the columns of this matrix
#' @param num_sim for generating a null distribution in calculating a p-value, the number of simulations to run.
#' @return list object
#'    \enumerate{
#'    \item \code{cont_table} The contingency table of clusters (rows) and batches (columns)
#'    \item \code{chi} A list containing all of the results from a traditional Chi Sqr
#'           as with the chisq.test function.
#'    \item \code{pmd} The percent maximum difference (pmd) of the input dataset.
#'    \item \code{pmd_null} Simulations of the null distribution of PMDs using the 
#'           observed global percentage of cluster abundances across all batches 
#'           with the observed batch sizes to match the input.
#'           Note that this will approach zero, but due to random sampling, typically will
#'           never actually get there. That's what makes having this null background useful.
#'    \item \code{pmd_z} The Z-statistic for how different the observed PMD is from the 
#'           distribution of null PMDs. However, in many cases, the
#'           null background distribution of PMDs is a skewed Poisson or beta-like distribution
#'           so a Z-statistic is less meaningful in that case. The returned \code{p.value} 
#'           however is generated as an empirical p-value as measured against the simulations.
#'     \item \code{p.value} - Significance for whether or not batches are different in their cellular composition.
#'           This is determined through the generating \code{num_sim} null distributions, 
#'           and emperically measuring the number of times the observed PMD was greater 
#'           than the PMDs generated from the null distribution.
#'           Low p-values indicate that the batches are indeed different from each other.
#'    }
#' @examples
#'    ## generate 
#'    cont_table <- matrix(c(50,100,30,80,90,50),nrow = 3, ncol = 2)
#'    ## Note that the clusters are in rows, and the batches are in columns 
#'    ## of this matrix. This is important!
#'    pmd_res <- pmd_from_cont_table(cont_table)
#' @importFrom stats chisq.test sd
#' @importFrom MASS fitdistr
#' @name pmd_from_cont_table
#' @export
pmd_from_cont_table<-function(cont_table, num_sim = 10000){
    pmd_results<-list()
    chi_result <- suppressWarnings(chisq.test(cont_table))
    cur_pmd<-get_percent_max_resid(cont_table,chi_result$expected)
    pmd_results$pmd_null<-get_pmd_null_vect(chi_result$expected, num_sim = num_sim)
    pmd_results$cont_table<-cont_table
    pmd_results$pmd_raw<-cur_pmd
    pmd_results$chi<-chi_result
    pmd_results$p.value<-sum(cur_pmd<pmd_results$pmd_null)/num_sim
    temp_fit<-suppressWarnings(fitdistr(pmd_results$pmd_null, densfun="poisson"))
    lambda <- as.numeric(temp_fit$estimate)
    pmd_results$pmd_null_lambda <- lambda
    ## lambda is the center of mass of the null distribution.
    ## This is what controls both the slope and intercept of the PMD line.
    ## Here we'll correct for that
    null_mean<-mean(pmd_results$pmd_null)
    null_sd<-sd(pmd_results$pmd_null)
    pmd_z<-(cur_pmd-null_mean)/null_sd
    pmd_results$pmd_z<-pmd_z
    ## log the final pmd
    pmd_results$pmd <- (cur_pmd - lambda) / (1 - lambda)
    return(pmd_results)
}



#' Percent Maximum Difference
#' @description \code{pmd} Percent Maximum Difference provides a metric that enables the measurement
#'                         of how similar (PMD near 0) or different (PMD near 1) two batches are
#'                         based just on clustering results. This can be viewed as an alternative
#'                         to a Chi-squared type analysis based on a contingency table in which 
#'                         the batches are noted in columns, and the abundance of cells in each cluster
#'                         is noted in rows. PMD advances on Chi-squared and Fishers exact tests
#'                         because its results are actually not dependent on the number of clusters found.
#'                         It does this by first calculating the hypothetical maximum possible asymmetry 
#'                         between your batches. (i.e.: if each batch was only comprised of clusters 
#'                         that only ever appeared in that batch and there were no clusters that appeared 
#'                         in more than one batch). Then it calculates the observed asymmetry in your 
#'                         real results, then returning the percentage of the observed asymmetry
#'                         relative to the maximum possible asymmetry. This gives a lower bound of zero,
#'                         where the relative abundance of each cluster was spot-on the same between
#'                         all of the batches. This also gives an upper bound of one, in which 
#'                         the observed clustering results were fully asymmetric across batch; 
#'                         or in other words, each cluster was fully batch specific.
#' @param batch_labels a vector of labels for the batch that each sample (or cell for scRNAseq) came from. Can be a string or factor.
#' @param cluster_labels  a vector of labels for each cluster that each sample belongs to. Can be a string or factor.
#' @param num_sim for generating a null distribution in calculating a p-value, the number of simulations to run.
#' @return list object
#'    \enumerate{
#'    \item \code{cont_table} The contingency table of clusters (rows) and batches (columns)
#'    \item \code{chi} A list containing all of the results from a traditional Chi Sqr
#'           as with the chisq.test function.
#'    \item \code{pmd} The percent maximum difference (pmd) of the input dataset.
#'    \item \code{pmd_null} Simulations of the null distribution of PMDs using the 
#'           observed global percentage of cluster abundances across all batches 
#'           with the observed batch sizes to match the input.
#'           Note that this will approach zero, but due to random sampling, typically will
#'           never actually get there. That's what makes having this null background useful.
#'    \item \code{pmd_z} The Z-statistic for how different the observed PMD is from the 
#'           distribution of null PMDs. However, in many cases, the
#'           null background distribution of PMDs is a skewed Poisson or beta-like distribution
#'           so a Z-statistic is less meaningful in that case. The returned \code{p.value} 
#'           however is generated as an empirical p-value as measured against the simulations.
#'     \item \code{p.value} - Significance for whether or not batches are different in their cellular composition.
#'           This is determined through the generating \code{num_sim} null distributions, 
#'           and emperically measuring the number of times the observed PMD was greater 
#'           than the PMDs generated from the null distribution.
#'           Low p-values indicate that the batches are indeed different from each other.
#'    }
#' @examples
#'    ## first using batch and cluster labels
#'    batch <- as.factor(rep(seq_len(2),each=100))
#'    clusters <- as.factor(rep(rep(seq_len(2),each=50),2))
#'    pmd_res <- pmd(batch, clusters)
#' @importFrom stats chisq.test sd
#' @importFrom MASS fitdistr
#' @name pmd
#' @export
pmd<-function(batch_labels, cluster_labels, num_sim = 10000){
    cont_table<-get_cont_table(as.factor(cluster_labels),as.factor(batch_labels))
    pmd_results <- pmd_from_cont_table(cont_table, num_sim = num_sim)
    return(pmd_results)
}

