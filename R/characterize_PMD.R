#' characterize_pmd
#' @description \code{characterize_pmd} Run the PMD characterization over a range of power conditions.
#' @param num_cells_b1 The number of cells in batch 1
#' @param num_cells_b2 The number of cells in batch 2
#' @param num_clust_shared the number of clusters that are shared across the two batches
#' @param iters the number of iterations for each condition
#' @importFrom MASS fitdistr
#' @importFrom stats chisq.test
#' @return list object
#'    \enumerate{
#'    \item \code{out_df}  A dataframe that has the results of the characterization. Also returns a vector 
#'    \item \code{null_df} A dataframe with the vector of the PMDs generated from the null distribution of the input.
#'    }
#' @include PercentMaximumDifference.R
#' @examples
#'    characterize_pmd()
#' @name characterize_pmd
#' @export
characterize_pmd <- function(num_cells_b1 = 1000,
                            num_cells_b2 = 1000,
                            num_clust_shared = seq(from=0, to=10),
                            iters = 10) {
    print(paste("simulating",num_cells_b1, "and", num_cells_b2))
    max_clust <- max(num_clust_shared)
    total_cells <- c()
    num_clust <- c()
    dataset_sizes <- c()
    null_ds_names <- c()
    pmd_null_vect <- c()
    min_of_expected_mat <- c()
    lambda_estimate <- c()
    iter_vect <- c()
    pmd_raw_vect <- c()
    pmd_vect <- c()
    num_shared_vect <- c()
    for (iter in seq(from = 1, to = iters)) {
        for (num_shared in num_clust_shared){
            temp_ds_name <- paste(num_cells_b1, "vs", num_cells_b2, sep="_")
            clust_prob_1 <- rep(1,max_clust) / max_clust
            clust_prob_2 <- rep(1,max_clust) / max_clust
            clust_vect_1 <- get_random_sample_cluster(clust_prob_1, num_cells_b1)
            clust_vect_2 <- get_random_sample_cluster(clust_prob_2, num_cells_b2) + num_shared
            clusters <- c(clust_vect_1, clust_vect_2)
            batches <- c(rep(1,num_cells_b1),rep(2, num_cells_b2))
            temp_pmd <- get_percent_max_diff(batches, clusters)
            ## now get the null, fit the Poisson distribution & log the params
            cont_mat <- get_cont_table(as.factor(clusters), as.factor(batches))
            expected_mat <- chisq.test(cont_mat)$expected
            pmd_null <- get_pmd_null_vect(expected_mat, num_sim = 100)
            if (num_shared == 0) {
                null_ds_names <- c(null_ds_names, rep(temp_ds_name, length(pmd_null)) )
                pmd_null_vect <- c(pmd_null_vect, pmd_null)
            }
            temp_fit <- fitdistr(pmd_null, densfun="poisson")
            temp_lambda <- as.numeric(temp_fit$estimate)
            final_pmd <- (temp_pmd - temp_lambda) / (1 - temp_lambda)
            ## now log it all
            iter_vect <- c(iter_vect, iter)
            total_cells <- c(total_cells, length(clusters))
            min_of_expected_mat <- c(min_of_expected_mat, min(expected_mat))
            dataset_sizes <- c(dataset_sizes, temp_ds_name)
            num_clust <- c(num_clust, length(unique(clusters)))
            lambda_estimate <- c(lambda_estimate, temp_lambda)
            pmd_raw_vect <- c(pmd_raw_vect, temp_pmd)
            pmd_vect <- c(pmd_vect, final_pmd)
            num_shared_vect <- c(num_shared_vect, num_shared)
        }
    }
    out_df <- data.frame(dataset_sizes = dataset_sizes,
                        iter = iter_vect,
                        total_cells = total_cells,
                        total_number_of_clusters = num_clust,
                        b1_size = rep(num_cells_b1, length(total_cells)),
                        b2_size = rep(num_cells_b2, length(total_cells)),
                        batch_size_difference = rep(abs(num_cells_b1-num_cells_b2), length(total_cells)),
                        min_of_expected_mat = min_of_expected_mat,
                        null_distribution_lambda_estimate = lambda_estimate,
                        num_non_shared_clusters = num_shared_vect,
                        raw_pmd = pmd_raw_vect,
                        pmd = pmd_vect)
    null_df <- data.frame(dataset_sizes = null_ds_names,
                          pmd_null = pmd_null_vect)
    return_obj <- list()
    return_obj$out_df <- out_df
    return_obj$null_df <- null_df
    return(return_obj)
}



#' do_full_pmd_characterization
#' @description \code{do_full_pmd_characterization} Run the PMD characterization over a range of power conditions.
#' @param directory the directory for putting the output image. This isn't needed though;
#'                  we'll skip the plotting in that case.
#' @param b1_size_vect a vector denoting the size for the first batch for the given iteration 
#' @param b2_size_vect a vector denoting the size for the second batch for the given iteration
#' @include PercentMaximumDifference.R
#' @return Null
#' @importFrom ggplot2 ggplot geom_point geom_smooth aes
#' @importFrom grDevices dev.off png
#' @importFrom stats loess lm
#' @examples
#'    pmd_bench_table_combined <- do_full_pmd_characterization(b1_size_vect = c(250, 250),
#'                                                             b2_size_vect = c(500, 1000))
#' @name do_full_pmd_characterization
#' @export
do_full_pmd_characterization<-function(directory = '', 
                                        b1_size_vect = c(10000, 10000, 1000, 250, 250, 100), 
                                        b2_size_vect = c(10000, 1000, 1000, 2000, 250, 100)) {
    results_list <- list()
    null_list <- list()
    run_names <- c()
    for (i in seq(from = 1, to = length(b1_size_vect))) {
        temp_b1_size <- b1_size_vect[i]
        temp_b2_size <- b2_size_vect[i]
        temp_run_name <- paste(temp_b1_size, "vs", temp_b2_size,sep="_")
        run_names <- c(run_names, temp_run_name)
        temp_obj <- characterize_pmd(num_cells_b1 = temp_b1_size,
                                     num_cells_b2 = temp_b2_size)
        results_list[[temp_run_name]] <- temp_obj$out_df
        null_list[[temp_run_name]] <- temp_obj$null_df
    }
    out_df_names <- names(results_list[[run_names[1]]])
    pmd_bench_table_combined <- results_list[[run_names[1]]]
    null_names <- names(null_list[[run_names[1]]])
    pmd_null_combined <- null_list[[run_names[1]]]
    for (i in seq(from = 2, to = length(b1_size_vect))){
        pmd_bench_table_combined <- rbind(pmd_bench_table_combined,
                                          results_list[[run_names[i]]])
        pmd_null_combined <- rbind(pmd_null_combined,
                                   null_list[[run_names[i]]])
    }
    names(pmd_bench_table_combined) <- out_df_names
    names(pmd_null_combined) <- null_names
    print(pmd_null_combined)
    if (directory != ''){  
        outfile <- paste(directory, "pmd_characterization%01d.png", sep = '/')
        print(outfile)
        png(outfile,width=3000, height=2100, res=600)
        p <- ggplot(pmd_bench_table_combined, aes(x = num_non_shared_clusters, y = raw_pmd, color = dataset_sizes)) +
                    geom_point(shape = 1) +    # Use hollow circles
                    geom_smooth(method = lm)   # Add linear regression line 
                                     #  (by default includes 95% confidence region)
        print(p)
        p <- ggplot(pmd_null_combined, aes(pmd_null, fill = dataset_sizes)) +
                    geom_density(adjust = 3, alpha=.3)+
              ylab("Density") + xlab("Raw PMD") +
              xlim(0,1) +
              ggtitle("Distribution of raw PMD when no difference by batch")
        print(p)
        p <- ggplot(pmd_bench_table_combined, aes(x = total_number_of_clusters, y = null_distribution_lambda_estimate, color = dataset_sizes)) +
                    geom_point(shape = 1) +    # Use hollow circles
                    geom_smooth(method = lm)   # Add linear regression line 
                                     #  (by default includes 95% confidence region)
        print(p)
        p <- ggplot(pmd_bench_table_combined, aes(x = log2(min_of_expected_mat), y = null_distribution_lambda_estimate, color = dataset_sizes)) +
                    geom_point(shape = 1) +    # Use hollow circles
                    geom_smooth(method = loess, span = 50)   # Add linear regression line 
                                     #  (by default includes 95% confidence region)
        print(p)
        p <- ggplot(pmd_bench_table_combined, aes(x = num_non_shared_clusters, y = pmd, color = dataset_sizes)) +
                    geom_point(shape = 1) +    # Use hollow circles
                    geom_smooth(method = lm)   # Add linear regression line 
                                     #  (by default includes 95% confidence region)
        print(p)
        dev.off()
    }
    return(pmd_bench_table_combined)
}


