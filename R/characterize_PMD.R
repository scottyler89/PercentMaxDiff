
#' get_avg_entropy
#' @description \code{get_avg_entropy} Gets the average entropy of clusters by batch
#' @param cont_mat batches in columns, clusters in rows
#' @importFrom downsample downsample_mat
#' @importFrom entropy entropy
#' @return numeric object
#'    \enumerate{
#'    \item \code{avg_clust_entropy_by_batch}  The average Shannon entropy for clusters as they segregate by batch. Input batches are downsampled so that entropy can actually be computed fairly.
#'    }
#' @include PercentMaximumDifference.R
#' @examples
#'    get_avg_entropy(matrix())
#' @name characterize_pmd
#' @export
get_avg_entropy<-function(cont_mat){
    if (dim(cont_mat)[1]==1){
        ## to easily handle the special case in which we have a single cluster
        ## we'll just duplicate the one row. This one make a difference because
        ## the entropy across both rows 
        cont_mat<-rbind(cont_mat,cont_mat)
    }
    cont_mat<-downsample_mat(cont_mat, quiet=TRUE)
    avg_clust_entropy_by_batch<-mean(apply(cont_mat,1,entropy))
    return(avg_clust_entropy_by_batch)
}



#' characterize_pmd
#' @description \code{characterize_pmd} Run the PMD characterization over a range of power conditions.
#' @param num_cells_b1 The number of cells in batch 1
#' @param num_cells_b2 The number of cells in batch 2
#' @param num_clust_shared the number of clusters that are shared across the two batches
#' @param iters the number of iterations for each condition
#' @param pmd_2 A boolean argument from when we had thought about a standardized version of PMD that used squares/sqrts instead of absolute values. It worked much worse, so don't bother using this arg.
#' @importFrom MASS fitdistr
#' @importFrom rcompanion cramerV
#' @importFrom stats chisq.test
#' @importFrom vegan diversity
#' @importFrom entropy entropy
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
                            iters = 10,
                            pmd_2 = FALSE) {
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
    chi_vect <- c()
    cramers_v <- c()
    chi_neg_log_p_vect <- c()
    inv_simp_idx_vect <- c()
    shan_entrop_vect <- c()
    for (iter in seq(from = 1, to = iters)) {
        for (num_shared in num_clust_shared){
            temp_ds_name <- paste(num_cells_b1, "vs", num_cells_b2, sep="_")
            clust_prob_1 <- rep(1,max_clust) / max_clust
            clust_prob_2 <- rep(1,max_clust) / max_clust
            clust_vect_1 <- get_random_sample_cluster(clust_prob_1, num_cells_b1)
            clust_vect_2 <- get_random_sample_cluster(clust_prob_2, num_cells_b2) + num_shared
            clusters <- c(clust_vect_1, clust_vect_2)
            batches <- c(rep(1,num_cells_b1),rep(2, num_cells_b2))
            if (pmd_2){
                temp_pmd <- get_percent_max_diff_2(batches, clusters)
            } else {
                temp_pmd <- get_percent_max_diff(batches, clusters)
            }
            ## now get the null, fit the Poisson distribution & log the params
            cont_mat <- get_cont_table(as.factor(clusters), as.factor(batches))
            expected_mat <- get_expected(cont_mat)
            if (pmd_2){
                pmd_null <- get_pmd_null_vect_2(expected_mat, num_sim = 100)
            } else {
                pmd_null <- get_pmd_null_vect(expected_mat, num_sim = 100)
            }
            if (num_shared == 0) {
                null_ds_names <- c(null_ds_names, rep(temp_ds_name, length(pmd_null)) )
                pmd_null_vect <- c(pmd_null_vect, pmd_null)
            }
            ## get the chi result
            temp_chi_res <- chisq.test(cont_mat)
            chi_vect <-c(chi_vect, temp_chi_res$statistic)
            chi_neg_log_p_vect <- c(chi_neg_log_p_vect,-log10(temp_chi_res$p.value))
            ## get the lambda
            temp_fit <- fitdistr(pmd_null, densfun="poisson")
            temp_lambda <- as.numeric(temp_fit$estimate)
            final_pmd <- (temp_pmd - temp_lambda) / (1 - temp_lambda)
            ## get Cramers V
            temp_cramers_v <- as.numeric(cramerV(cont_mat, bias.correct = TRUE))
            cramers_v <- c(cramers_v, temp_cramers_v)
            ## get inverse simpson's index
            inv_simp_idx_vect<-c(inv_simp_idx_vect,mean(diversity(t(cont_mat),MARGIN=2,index="invsimpson")))
            ## get shannon entropy
            shan_entrop_vect<-c(shan_entrop_vect,get_avg_entropy(cont_mat))
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
                        pmd = pmd_vect,
                        chi_sq = chi_vect,
                        chi_neg_log_p = chi_neg_log_p_vect,
                        cramers_v = cramers_v,
                        inverse_simp = inv_simp_idx_vect,
                        shannon_entropy = shan_entrop_vect)
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
#' @param pmd_2 A boolean argument from when we had thought about a standardized version of PMD that used squares/sqrts instead of absolute values. It worked much worse, so don't bother using this arg.
#' @include PercentMaximumDifference.R
#' @return Null
#' @importFrom ggplot2 ggplot geom_point geom_smooth aes geom_density ylab xlab xlim ggtitle
#' @importFrom grDevices dev.off png
#' @importFrom stats loess lm
#' @examples
#'    pmd_bench_table_combined <- do_full_pmd_characterization(b1_size_vect = c(250, 250),
#'                                                             b2_size_vect = c(500, 1000))
#' @name do_full_pmd_characterization
#' @export
do_full_pmd_characterization<-function(directory = '', 
                                        b1_size_vect = c(10000, 10000, 1000, 250, 250, 100), 
                                        b2_size_vect = c(10000, 1000, 1000, 2000, 250, 100),
                                        pmd_2 = FALSE) {
    results_list <- list()
    null_list <- list()
    run_names <- c()
    for (i in seq(from = 1, to = length(b1_size_vect))) {
        temp_b1_size <- b1_size_vect[i]
        temp_b2_size <- b2_size_vect[i]
        temp_run_name <- paste(temp_b1_size, "vs", temp_b2_size,sep="_")
        run_names <- c(run_names, temp_run_name)
        temp_obj <- characterize_pmd(num_cells_b1 = temp_b1_size,
                                     num_cells_b2 = temp_b2_size,
                                     pmd_2 = pmd_2)
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
    if (directory != ''){
        dir.create(directory, showWarnings=F)
        outfile <- paste(directory, "pmd_characterization%01d.png", sep = '/')
        print(outfile)
        png(outfile,width=3000, height=2100, res=600)
        p <- ggplot(pmd_bench_table_combined, aes(x = num_non_shared_clusters, y = raw_pmd, color = dataset_sizes)) +
                    geom_point(shape = 1) +    # Use hollow circles
                    geom_smooth(lwd = 2, method = loess)   # Add linear regression line 
                                     #  (by default includes 95% confidence region)
        print(p)
        p <- ggplot(pmd_bench_table_combined, aes(x = num_non_shared_clusters, y = chi_sq, color = dataset_sizes)) +
            geom_point(shape = 1) +    # Use hollow circles
            geom_smooth(lwd = 2, method = loess)   # Add linear regression line 
                             #  (by default includes 95% confidence region)
        print(p)
        p <- ggplot(pmd_bench_table_combined, aes(x = num_non_shared_clusters, y = chi_neg_log_p, color = dataset_sizes)) +
            geom_point(shape = 1) +    # Use hollow circles
            geom_smooth(lwd = 2, method = loess)   # Add linear regression line 
                             #  (by default includes 95% confidence region)
        print(p)
        p <- ggplot(pmd_bench_table_combined, aes(x = num_non_shared_clusters, y = cramers_v, color = dataset_sizes)) +
            geom_point(shape = 1) +    # Use hollow circles
            geom_smooth(lwd = 2, method = loess)   # Add linear regression line 
                             #  (by default includes 95% confidence region)
        print(p)
        ##########################################
        p <- ggplot(pmd_bench_table_combined, aes(x = num_non_shared_clusters, y = inverse_simp, color = dataset_sizes)) +
            geom_point(shape = 1) +    # Use hollow circles
            geom_smooth(lwd = 2, method = loess)   # Add linear regression line 
                             #  (by default includes 95% confidence region)
        print(p)
        p <- ggplot(pmd_bench_table_combined, aes(x = num_non_shared_clusters, y = shannon_entropy, color = dataset_sizes)) +
            geom_point(shape = 1) +    # Use hollow circles
            geom_smooth(lwd = 2, method = loess)   # Add linear regression line 
                             #  (by default includes 95% confidence region)
        print(p)
        ##########################################
        p <- ggplot(pmd_null_combined, aes(pmd_null, fill = dataset_sizes)) +
                    geom_density(adjust = 3, alpha=.3)+
              ylab("Density") + xlab("Raw PMD") +
              xlim(0,1) +
              ggtitle("Distribution of raw PMD when no difference by batch")
        print(p)
        p <- ggplot(pmd_bench_table_combined, aes(x = total_number_of_clusters, y = null_distribution_lambda_estimate, color = dataset_sizes)) +
                    geom_point(shape = 1) +    # Use hollow circles
                    geom_smooth(lwd = 2, method = loess)   # Add linear regression line 
                                     #  (by default includes 95% confidence region)
        print(p)
        p <- ggplot(pmd_bench_table_combined, aes(x = log2(min_of_expected_mat), y = null_distribution_lambda_estimate, color = dataset_sizes)) +
                    geom_point(shape = 1) +    # Use hollow circles
                    geom_smooth(lwd = 2, method = loess, span = 50)   # Add loess regression
                                     #  (by default includes 95% confidence region)
        print(p)
        p <- ggplot(pmd_bench_table_combined, aes(x = num_non_shared_clusters, y = pmd, color = dataset_sizes)) +
                    geom_point(shape = 1) +    # Use hollow circles
                    geom_smooth(lwd = 2, method = loess)   # Add linear regression line 
                                     #  (by default includes 95% confidence region)
        print(p)
        dev.off()
    }
    return(pmd_bench_table_combined)
}



#' characterize_invariance_property
#' @description \code{characterize_invariance_property} Run the PMD characterization over a range of power conditions.
#' @param directory the directory for putting the output image. This isn't needed though;
#'                  we'll skip the plotting in that case.
#' @param num_clust_per_batch_range a vector of number of clusters (as integers)
#' @param b1_size_vect a vector denoting the size for the first batch for the given iteration 
#' @param b2_size_vect a vector denoting the size for the second batch for the given iteration
#' @param percent_difference a floating point betwee 0. and 1. that corresponds to the 
#'                           percentage difference between the two batches
#' @param iterations an integer that 
#' @param expand_b2_only a TRUE/FALSE boolean operator that dictates the mechanism of simulation in terms of what the 
#'                \code{num_clust_per_batch_range} argument operates on. When FALSE, both batches
#'                are expanded with increasing numbers of groups & what determins \code{percent_difference} is the number
#'                of shared clusters. When TRUE, batch1 has all its cells allocated to a single cluster & batch 2
#'                has the \code{percent_difference} percent of cells allocated to that same cluster as batch 1, and the rest
#'                of its cells are partitioned equally across the remaining clusters.
#' @include PercentMaximumDifference.R
#' @return Null
#' @importFrom ggplot2 ggplot geom_point geom_smooth aes geom_density ylab xlab xlim ggtitle scale_y_continuous theme element_text stat_smooth
#' @importFrom grDevices dev.off png
#' @importFrom stats loess lm
#' @importFrom scales hue_pal
#' @importFrom cowplot plot_grid
#' @importFrom rcompanion cramerV
#' @importFrom vegan diversity
#' @importFrom entropy entropy
#' @examples
#'    pmd_bench_table_combined <- characterize_invariance_property(b1_size_vect = c(250, 250),
#'                                                                b2_size_vect = c(500, 1000),
#'                                                                iterations = 3)
#' @name characterize_invariance_property
#' @export
characterize_invariance_property<-function(directory='',
                                           num_clust_per_batch_range = seq(4, 17, by = 4),
                                           b1_size_vect = c(10000, 10000, 1000, 250, 250, 100), 
                                           b2_size_vect = c(10000, 1000, 1000, 2000, 250, 100),
                                           percent_difference = c(0, .25, .5, .75, 1),
                                           iterations = 10,
                                           expand_b2_only = FALSE ){
    iter_vect <- c()
    full_b1_size_vect <- c()
    full_b2_size_vect <- c()
    percent_difference_vect <- c()
    original_order_vect <- c()
    total_cells <- c()
    num_clust <- c()
    dataset_sizes <- c()
    num_clust_per_batch_vect <- c()
    percent_difference_vect <- c()
    pmd_raw_vect <- c()
    pmd_vect <- c()
    chi_stat_vect <- c()
    chi_neg_log_p_vect <- c()
    cramers_v <- c()
    inv_simp_idx_vect <- c()
    shan_entrop_vect <- c()
    max_possible_neg_log_p <- -log10(as.numeric(noquote(unlist(format(.Machine)))['double.xmin']))
    for (num_clust_per_batch in num_clust_per_batch_range){
        print(num_clust_per_batch)
        for (pcnt_diff in percent_difference){
            for (i in seq(1,length(b1_size_vect))) {
                for (iter in seq(1,iterations)){
                    num_cells_b1 <- b1_size_vect[i]
                    num_cells_b2 <- b2_size_vect[i]
                    full_b1_size_vect <- c(full_b1_size_vect, num_cells_b1)
                    full_b2_size_vect <- c(full_b2_size_vect, num_cells_b2)
                    temp_ds_name <- paste(num_cells_b1, "vs", num_cells_b2, sep="_")
                    if (!(temp_ds_name %in% original_order_vect)){
                        original_order_vect <- c(original_order_vect, temp_ds_name)
                    }
                    dataset_sizes <- c(dataset_sizes, temp_ds_name)
                    prob_vect<-rep(1,num_clust_per_batch)/num_clust_per_batch
                    if (expand_b2_only) {
                        clust_vect_1 <- get_random_sample_cluster(1,num_cells_b1)
                        #clust_vect_1 <- get_random_sample_cluster(rep(1,4)/4,num_cells_b1)
                    } else {
                        clust_vect_1 <- get_random_sample_cluster(prob_vect,num_cells_b1)
                    }
                    if ( num_clust_per_batch == 1 ){
                        ## wanted to include 1, but this messes with any percent overlap that isn't 
                        ## either 0 or 1
                        shift <- as.integer(round(pcnt_diff))
                        temp_pcnt_diff <- shift
                    }
                    else {
                        if (expand_b2_only){
                            shift <- as.integer((4)*pcnt_diff)
                        } else {
                            shift <- as.integer((num_clust_per_batch)*pcnt_diff)
                        }
                        temp_pcnt_diff <- pcnt_diff
                    }
                    if (expand_b2_only){
                    #if (expand_b2_only || num_clust_per_batch == 1 ){
                        clust_vect_2 <- get_random_sample_cluster(c(1-pcnt_diff,## this is the percent overlap cluster
                                                                    rep(pcnt_diff/num_clust_per_batch,num_clust_per_batch)),
                                                                    num_cells_b2)
                        # clust_vect_2 <- get_random_sample_cluster(c(rep((1-pcnt_diff)/4,4),## this is the percent overlap cluster
                        #                      rep(pcnt_diff/num_clust_per_batch,num_clust_per_batch)),
                        #                      num_cells_b2)
                    } else {
                        clust_vect_2 <- get_random_sample_cluster(prob_vect,num_cells_b2)+shift
                    }
                    # actual_percent_diff <- 1-( length(intersect(unique(clust_vect_1),unique(clust_vect_2)))/length(union(unique(clust_vect_1),unique(clust_vect_2)) ) )
                    # print(paste("   shift clusters:", as.integer(num_clust_per_batch*pcnt_diff)))
                    # print(paste("    actual_percent_diff:",actual_percent_diff))
                    clusters <- c(clust_vect_1, clust_vect_2)
                    batches <- c(rep(1,num_cells_b1),rep(2, num_cells_b2))
                    temp_cont_table <- get_cont_table(clusters, batches)
                    temp_pmd_res <- pmd(batches, clusters)
                    temp_chi_res <- chisq.test(temp_cont_table)
                    temp_cramers_v <- as.numeric(cramerV(temp_cont_table, bias.correct = TRUE))
                    # print(paste("intended:",pcnt_diff,"raw_pmd:",temp_pmd_res$pmd_raw,"pmd:",temp_pmd_res$pmd,"chi-sq:",temp_chi_res$statistic))
                    # print(temp_pmd_res$cont_table)
                    ## log the results
                    iter_vect<-c(iter_vect, iter)
                    total_cells <- c(total_cells, length(clusters))
                    if (expand_b2_only){
                    #if (expand_b2_only || num_clust_per_batch == 1 ){
                        percent_difference_vect <- c(percent_difference_vect, pcnt_diff)
                    } else {
                        percent_difference_vect <- c(percent_difference_vect, temp_pcnt_diff)
                    }
                    num_clust <- c(num_clust, length(unique(clusters)))
                    num_clust_per_batch_vect <- c(num_clust_per_batch_vect, num_clust_per_batch)
                    pmd_raw_vect <- c(pmd_raw_vect, temp_pmd_res$pmd_raw)
                    pmd_vect <- c(pmd_vect, temp_pmd_res$pmd)
                    chi_stat_vect <- c(chi_stat_vect, temp_chi_res$statistic)
                    cramers_v <- c(cramers_v, temp_cramers_v)
                    chi_neg_log_p_vect <- c(chi_neg_log_p_vect, min(c(max_possible_neg_log_p,-log10(temp_chi_res$p.value))))
                    ## get inverse simpson's index
                    inv_simp_idx_vect<-c(inv_simp_idx_vect,mean(diversity(t(temp_cont_table),MARGIN=2,index="invsimpson")))
                    ## get shannon entropy
                    shan_entrop_vect<-c(shan_entrop_vect,get_avg_entropy(temp_cont_table))

                }
            }
        }
    }
    pmd_bench_table_combined <- data.frame(dataset_sizes = factor(dataset_sizes, levels = original_order_vect),
                                            PrcntDiffClusters = factor(percent_difference_vect, levels = sort(unique(percent_difference_vect))),
                                            percent_different_clusters_numeric = percent_difference_vect,
                                            NumberClusters = as.factor(num_clust_per_batch_vect),
                                            NumberOfClusters = num_clust_per_batch_vect,
                                            iter = iter_vect,
                                            total_cells = total_cells,
                                            total_number_of_clusters = num_clust,
                                            b1_size = full_b1_size_vect,
                                            b2_size = full_b2_size_vect,
                                            raw_pmd = pmd_raw_vect,
                                            pmd = pmd_vect,
                                            chi_sq = chi_stat_vect,
                                            chi_neg_log_p = chi_neg_log_p_vect,
                                            cramers_v = cramers_v,
                                            inverse_simp = inv_simp_idx_vect,
                                            shannon_entropy = shan_entrop_vect)
    if (directory != ''){ 
        if (!file.exists(directory)){
            dir.create(directory, showWarnings = FALSE)
        }
        outfile <- paste(directory, "pmd_characterization%01d.png", sep = '/')
        print(outfile)
        png(outfile,width=3000, height=2100, res=600)
        p1 <- ggplot(pmd_bench_table_combined, aes(x = NumberOfClusters, y = raw_pmd, color = PrcntDiffClusters, linetype = dataset_sizes)) +#color = dataset_sizes, linetype = percent_difference)
                    geom_point(shape = 1) +    # Use hollow circles
                    geom_smooth(lwd = 2, method = loess, formula = y ~ x, span = 8) +  # Add linear regression line 
                    #ylim(-0.1,1.05) +                #  (by default includes 95% confidence region)
                    scale_y_continuous(breaks = seq(0,1,by=.25)) +
                    ggtitle("Raw PMD") +
                    xlab("Number Clusters") +
                    ylab("Raw PMD") + 
                    theme(text = element_text(size = 18))
        print(p1)
        p2 <- ggplot(pmd_bench_table_combined, aes(x = NumberOfClusters, y = pmd, color = PrcntDiffClusters, linetype = dataset_sizes)) +
                    geom_point(shape = 1) +    # Use hollow circles
                    geom_smooth(lwd = 2, method = loess, formula = y ~ x, span = 8) +  # Add linear regression line 
                    #ylim(-0.1,1.1) +                #  (by default includes 95% confidence region)
                    scale_y_continuous(breaks = seq(0,1,by=.25)) +
                    ggtitle("PMD") +
                    xlab("Number Clusters") +
                    ylab("PMD") + 
                    theme(text = element_text(size = 18))
        print(p2)
        p3 <- ggplot(pmd_bench_table_combined, aes(x = NumberOfClusters, y = chi_sq, color = PrcntDiffClusters, linetype = dataset_sizes)) +
            geom_point(shape = 1) +    # Use hollow circles
            geom_smooth(lwd = 2, method = loess, formula = y ~ x, span = 8) +  # Add linear regression line 
            ggtitle("Chi-square Statistic") +
            xlab("Number Clusters") +
            ylab("Chi-square Stat") + 
            theme(text = element_text(size = 18))
        print(p3)
        p4 <- ggplot(pmd_bench_table_combined, aes(x = NumberOfClusters, y = chi_neg_log_p, color = PrcntDiffClusters, linetype = dataset_sizes)) +
            geom_point(shape = 1) +    # Use hollow circles
            geom_smooth(lwd = 2, method = loess, formula = y ~ x, span = 8) +  # Add linear regression line 
            ggtitle("Chi-square -log10(p-value)")  +
            xlab("Number Clusters") +
            ylab("Chi-square -log10(p-value)") + 
            theme(text = element_text(size = 18))
        print(p4)
        p42 <- ggplot(pmd_bench_table_combined, aes(x = NumberOfClusters, y = cramers_v, color = PrcntDiffClusters, linetype = dataset_sizes)) +
            geom_point(shape = 1) +    # Use hollow circles
            geom_smooth(lwd = 2, method = loess, formula = y ~ x, span = 8) +  # Add linear regression line 
            ggtitle("Bias Corrected Cramer's V")  +
            xlab("Number Clusters") +
            ylab("Cramer's V") + 
            theme(text = element_text(size = 18))
        print(p42)
        ##
        pInvSimp1 <- ggplot(pmd_bench_table_combined, aes(x = NumberOfClusters, y = inverse_simp, color = PrcntDiffClusters, linetype = dataset_sizes)) +
            geom_point(shape = 1) +    # Use hollow circles
            geom_smooth(lwd = 2, method = loess, formula = y ~ x, span = 8) +  # Add linear regression line 
            ggtitle("Inverse Simpsons Index")  +
            xlab("Number Clusters") +
            ylab("Inverse Simpsons Index") + 
            theme(text = element_text(size = 18))
        print(pInvSimp1)
        ##
        ##
        pShanEnt1 <- ggplot(pmd_bench_table_combined, aes(x = NumberOfClusters, y = shannon_entropy, color = PrcntDiffClusters, linetype = dataset_sizes)) +
            geom_point(shape = 1) +    # Use hollow circles
            geom_smooth(lwd = 2, method = loess, formula = y ~ x, span = 8) +  # Add linear regression line 
            ggtitle("Shannon Entropy")  +
            xlab("Number Clusters") +
            ylab("Shannon Entropy") + 
            theme(text = element_text(size = 18))
        print(pShanEnt1)
        ##
        #####################################################
        p5 <- ggplot(pmd_bench_table_combined, aes(x = percent_different_clusters_numeric, y = raw_pmd, color = NumberClusters, linetype = dataset_sizes)) +#color = dataset_sizes, linetype = percent_difference)#p5 <- ggplot(pmd_bench_table_combined, aes(x = percent_different_clusters, y = raw_pmd, color = as.factor(NumberOfClusters), linetype = dataset_sizes)) +
            geom_point(shape = 1) +    # Use hollow circles
            #geom_smooth(lwd = 2, method = loess, formula = y ~ x, span = 1e7) +  # Add linear regression line 
            #ylim(-0.1,1.05) +                #  (by default includes 95% confidence region)
            scale_y_continuous(breaks = seq(0,1,by=.25))  +
            xlab("Percent Different Clusters") +
            ylab("Raw PMD") + 
            theme(text = element_text(size = 18))
        ############################
        all_batches <- unique(pmd_bench_table_combined$NumberClusters)
        temp_num_colors <- length(all_batches)
        temp_colors <- hue_pal()(temp_num_colors)
        for (temp_col_group_index in seq(1,temp_num_colors)){
            temp_batch <- all_batches[temp_col_group_index]
            temp_color <- temp_colors[temp_col_group_index]
            temp_subset <- pmd_bench_table_combined[pmd_bench_table_combined$NumberClusters==temp_batch,]
            #print(temp_batch)
            #print(head(temp_subset))
            #p6 <- p6 + ggplot(temp_subset, aes(x = percent_different_clusters_numeric, y = pmd, color = temp_color, linetype = dataset_sizes))
            if (length(unique(temp_subset$percent_different_clusters_numeric)) == 2){
                p5 <- p5 + stat_smooth(lwd = 2, data = temp_subset, method = lm, formula = y ~ x, span=1e6, se=TRUE)   # Add linear regression line 
            } else {
                p5 <- p5 + stat_smooth(lwd = 2, data = temp_subset, method = lm, formula = y ~ poly(x,2), span=1e6, se=TRUE)
            }
        }
        ############################
        print(p5)
        p6 <- ggplot(pmd_bench_table_combined, aes(x = percent_different_clusters_numeric, y = pmd, color = NumberClusters, linetype = dataset_sizes)) +
            geom_point(shape = 1) +    # Use hollow circles
            #geom_smooth(lwd = 2, method = lm, formula = y ~ poly(x,2), span=1e6) +
            #geom_smooth(lwd = 2, method = "gam", formula = y ~ s(x, bs = "cs")) +
            #stat_smooth(lwd = 2, method = loess, formula = jitter(y, factor=1) ~ jitter(x, factor=1), span=1e6) +  # Add linear regression line 
            #stat_smooth(lwd = 2, method = lm, formula = y ~ poly(x,min((length(unique(x)) - 1),2)), span=1e6, se=FALSE) +  # Add linear regression line 
            #stat_smooth(lwd = 2, method = lm, formula = jitter(y, factor=.01*max(y)) ~ poly(jitter(x, factor=.01),1), span=1e6, se=FALSE) +  # Add linear regression line 
            #ylim(-0.1,1.1) +                #  (by default includes 95% confidence region)
            scale_y_continuous(breaks = seq(0,1,by=.25))   +
            xlab("Percent Different Clusters") +
            ylab("PMD") + 
            theme(text = element_text(size = 18))
        ############################
        all_batches <- unique(pmd_bench_table_combined$NumberClusters)
        temp_num_colors <- length(all_batches)
        temp_colors <- hue_pal()(temp_num_colors)
        for (temp_col_group_index in seq(1,temp_num_colors)){
            temp_batch <- all_batches[temp_col_group_index]
            temp_color <- temp_colors[temp_col_group_index]
            temp_subset <- pmd_bench_table_combined[pmd_bench_table_combined$NumberClusters==temp_batch,]
            #print(temp_batch)
            #print(head(temp_subset))
            #p6 <- p6 + ggplot(temp_subset, aes(x = percent_different_clusters_numeric, y = pmd, color = temp_color, linetype = dataset_sizes))
            if (length(unique(temp_subset$percent_different_clusters_numeric)) == 2){
                p6 <- p6 + stat_smooth(lwd = 2, data = temp_subset, method = lm, formula = y ~ x, span=1e6, se=TRUE)   # Add linear regression line 
            } else {
                p6 <- p6 + stat_smooth(lwd = 2, data = temp_subset, method = lm, formula = y ~ poly(x,2), span=1e6, se=TRUE)
            }
        }
        ############################
        print(p6)
        ############################
        p7 <- ggplot(pmd_bench_table_combined, aes(x = percent_different_clusters_numeric, y = chi_sq, color = NumberClusters, linetype = dataset_sizes)) +
            geom_point(shape = 1)   +
            xlab("Percent Different Clusters") +
            ylab("Chi-square Statistic")  + 
            theme(text = element_text(size = 18))   # Use hollow circles
            #stat_smooth(lwd = 2, method = lm, formula = jitter(y, factor=.01*sqrt(max(y))) ~ poly(jitter(x, factor=.01),2), span=1e6, se=FALSE)
            #geom_smooth(lwd = 2, method = loess, formula = y ~ x, span = 1e7)   # Add linear regression line 
                             #  (by default includes 95% confidence region)
        all_batches <- unique(pmd_bench_table_combined$NumberClusters)
        temp_num_colors <- length(all_batches)
        temp_colors <- hue_pal()(temp_num_colors)
        for (temp_col_group_index in seq(1,temp_num_colors)){
            temp_batch <- all_batches[temp_col_group_index]
            temp_color <- temp_colors[temp_col_group_index]
            temp_subset <- pmd_bench_table_combined[pmd_bench_table_combined$NumberClusters==temp_batch,]
            #print(temp_batch)
            #print(head(temp_subset))
            #p6 <- p6 + ggplot(temp_subset, aes(x = percent_different_clusters_numeric, y = pmd, color = temp_color, linetype = dataset_sizes))
            if (length(unique(temp_subset$percent_different_clusters_numeric)) == 2){
                p7 <- p7 + stat_smooth(lwd = 2, data = temp_subset, method = lm, formula = y ~ x, span=1e6, se=TRUE)   # Add linear regression line 
            } else {
                p7 <- p7 + stat_smooth(lwd = 2, data = temp_subset, method = lm, formula = y ~ poly(x,2), span=1e6, se=TRUE)
            }
        }
        ############################
        print(p7)
        ############################
        p8 <- ggplot(pmd_bench_table_combined, aes(x = percent_different_clusters_numeric, y = chi_neg_log_p, color = NumberClusters, linetype = dataset_sizes)) +
            geom_point(shape = 1)   +
            xlab("Percent Different Clusters") +
            ylab("Chi-square -log10(p-value)")  + 
            theme(text = element_text(size = 18))   # Use hollow circles
            #stat_smooth(lwd = 2, method = lm, formula = jitter(y, factor=.01*sqrt(max(y))) ~ poly(jitter(x, factor=.01),2), span=1e6, se=FALSE)
            #geom_smooth(lwd = 2, method = loess, formula = y ~ x, span = 1e7)  # Add linear regression line 
                             #  (by default includes 95% confidence region)
        all_batches <- unique(pmd_bench_table_combined$NumberClusters)
        temp_num_colors <- length(all_batches)
        temp_colors <- hue_pal()(temp_num_colors)
        for (temp_col_group_index in seq(1,temp_num_colors)){
            temp_batch <- all_batches[temp_col_group_index]
            temp_color <- temp_colors[temp_col_group_index]
            temp_subset <- pmd_bench_table_combined[pmd_bench_table_combined$NumberClusters==temp_batch,]
            #print(temp_batch)
            #print(head(temp_subset))
            #p6 <- p6 + ggplot(temp_subset, aes(x = percent_different_clusters_numeric, y = pmd, color = temp_color, linetype = dataset_sizes))
            if (length(unique(temp_subset$percent_different_clusters_numeric)) == 2){
                p8 <- p8 + stat_smooth(lwd = 2, data = temp_subset, method = lm, formula = y ~ x, span=1e6, se=TRUE)   # Add linear regression line 
            } else {
                p8 <- p8 + stat_smooth(lwd = 2, data = temp_subset, method = lm, formula = y ~ poly(x,2), span=1e6, se=TRUE)
            }
        }
        ############################
        print(p8)
        ############################
        p82 <- ggplot(pmd_bench_table_combined, aes(x = percent_different_clusters_numeric, y = cramers_v, color = NumberClusters, linetype = dataset_sizes)) +
            geom_point(shape = 1)   +
            xlab("Percent Different Clusters") +
            ylab("Bias Corrected Cramer's V")  + 
            theme(text = element_text(size = 18))   # Use hollow circles
            #stat_smooth(lwd = 2, method = lm, formula = jitter(y, factor=.01*sqrt(max(y))) ~ poly(jitter(x, factor=.01),2), span=1e6, se=FALSE)
            #geom_smooth(lwd = 2, method = loess, formula = y ~ x, span = 1e7)  # Add linear regression line 
                             #  (by default includes 95% confidence region)
        all_batches <- unique(pmd_bench_table_combined$NumberClusters)
        temp_num_colors <- length(all_batches)
        temp_colors <- hue_pal()(temp_num_colors)
        for (temp_col_group_index in seq(1,temp_num_colors)){
            temp_batch <- all_batches[temp_col_group_index]
            temp_color <- temp_colors[temp_col_group_index]
            temp_subset <- pmd_bench_table_combined[pmd_bench_table_combined$NumberClusters==temp_batch,]
            #print(temp_batch)
            #print(head(temp_subset))
            #p6 <- p6 + ggplot(temp_subset, aes(x = percent_different_clusters_numeric, y = pmd, color = temp_color, linetype = dataset_sizes))
            if (length(unique(temp_subset$percent_different_clusters_numeric)) == 2){
                p82 <- p82 + stat_smooth(lwd = 2, data = temp_subset, method = lm, formula = y ~ x, span=1e6, se=TRUE)   # Add linear regression line 
            } else {
                p82 <- p82 + stat_smooth(lwd = 2, data = temp_subset, method = lm, formula = y ~ poly(x,2), span=1e6, se=TRUE)
            }
        }
        ############################
        print(p82)
        ############################
        ############################
        ############################
        pInvSimp2 <- ggplot(pmd_bench_table_combined, aes(x = percent_different_clusters_numeric, y = inverse_simp, color = NumberClusters, linetype = dataset_sizes)) +
            geom_point(shape = 1)   +
            xlab("Percent Different Clusters") +
            ylab("Inverse Simpsons Index")  + 
            theme(text = element_text(size = 18))   # Use hollow circles
            #stat_smooth(lwd = 2, method = lm, formula = jitter(y, factor=.01*sqrt(max(y))) ~ poly(jitter(x, factor=.01),2), span=1e6, se=FALSE)
            #geom_smooth(lwd = 2, method = loess, formula = y ~ x, span = 1e7)  # Add linear regression line 
                             #  (by default includes 95% confidence region)
        all_batches <- unique(pmd_bench_table_combined$NumberClusters)
        temp_num_colors <- length(all_batches)
        temp_colors <- hue_pal()(temp_num_colors)
        for (temp_col_group_index in seq(1,temp_num_colors)){
            temp_batch <- all_batches[temp_col_group_index]
            temp_color <- temp_colors[temp_col_group_index]
            temp_subset <- pmd_bench_table_combined[pmd_bench_table_combined$NumberClusters==temp_batch,]
            #print(temp_batch)
            #print(head(temp_subset))
            #p6 <- p6 + ggplot(temp_subset, aes(x = percent_different_clusters_numeric, y = pmd, color = temp_color, linetype = dataset_sizes))
            if (length(unique(temp_subset$percent_different_clusters_numeric)) == 2){
                pInvSimp2 <- pInvSimp2 + stat_smooth(lwd = 2, data = temp_subset, method = lm, formula = y ~ x, span=1e6, se=TRUE)   # Add linear regression line 
            } else {
                pInvSimp2 <- pInvSimp2 + stat_smooth(lwd = 2, data = temp_subset, method = lm, formula = y ~ poly(x,2), span=1e6, se=TRUE)
            }
        }
        ############################
        print(pInvSimp2)
        ############################
        pShanEnt2 <- ggplot(pmd_bench_table_combined, aes(x = percent_different_clusters_numeric, y = shannon_entropy, color = NumberClusters, linetype = dataset_sizes)) +
            geom_point(shape = 1)   +
            xlab("Percent Different Clusters") +
            ylab("Shannon Entropy")  + 
            theme(text = element_text(size = 18))   # Use hollow circles
            #stat_smooth(lwd = 2, method = lm, formula = jitter(y, factor=.01*sqrt(max(y))) ~ poly(jitter(x, factor=.01),2), span=1e6, se=FALSE)
            #geom_smooth(lwd = 2, method = loess, formula = y ~ x, span = 1e7)  # Add linear regression line 
                             #  (by default includes 95% confidence region)
        all_batches <- unique(pmd_bench_table_combined$NumberClusters)
        temp_num_colors <- length(all_batches)
        temp_colors <- hue_pal()(temp_num_colors)
        for (temp_col_group_index in seq(1,temp_num_colors)){
            temp_batch <- all_batches[temp_col_group_index]
            temp_color <- temp_colors[temp_col_group_index]
            temp_subset <- pmd_bench_table_combined[pmd_bench_table_combined$NumberClusters==temp_batch,]
            #print(temp_batch)
            #print(head(temp_subset))
            #p6 <- p6 + ggplot(temp_subset, aes(x = percent_different_clusters_numeric, y = pmd, color = temp_color, linetype = dataset_sizes))
            if (length(unique(temp_subset$percent_different_clusters_numeric)) == 2){
                pShanEnt2 <- pShanEnt2 + stat_smooth(lwd = 2, data = temp_subset, method = lm, formula = y ~ x, span=1e6, se=TRUE)   # Add linear regression line 
            } else {
                pShanEnt2 <- pShanEnt2 + stat_smooth(lwd = 2, data = temp_subset, method = lm, formula = y ~ poly(x,2), span=1e6, se=TRUE)
            }
        }
        ############################
        print(pShanEnt2)
        #####################################################
        dev.off()
        outfile <- paste(directory, "pmd_characterization_merged%01d.png", sep = '/')
        png(outfile,width=18000, height=4500, res=600)
        print(plot_grid(p2, p3, p4, p42, pInvSimp1, pShanEnt1, p6, p7, p8, p82, pInvSimp2, pShanEnt2, nrow=2, ncol=6))
        dev.off()
    }
    return(pmd_bench_table_combined)
}


