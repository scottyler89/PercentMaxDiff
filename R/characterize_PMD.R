#' characterize_pmd
#' @description \code{characterize_pmd} Run the PMD characterization over a range of power conditions.
#' @param num_cells_b1 The number of cells in batch 1
#' @param num_cells_b2 The number of cells in batch 2
#' @param num_clust_shared the number of clusters that are shared across the two batches
#' @param iters the number of iterations for each condition
#' @return Null
#' @include PercentMaximumDifference.R
#' @examples
#'    characterize_pmd()
#' @name characterize_pmd
#' @export
characterize_pmd<-function(num_cells_b1=1000,
	                       num_cells_b2=1000,
	                       num_clust_shared=0:10,
	                       iters = 10){
	max_clust <- max(num_clust_shared)
	skew_vect <-c()
	iter_vect <- c()
	pmd_vect <- c()
	num_shared_vect <- c()
	for (iter in 1:iters){
		for (num_shared in num_clust_shared){
			clust_prob_1 <- rep(1,max_clust)/max_clust
			clust_prob_2 <- rep(1,max_clust)/max_clust
			clust_vect_1 <- get_random_sample_cluster(clust_prob_1, num_cells_b1)
			clust_vect_2 <- get_random_sample_cluster(clust_prob_2, num_cells_b2)+num_shared
			clusters<-c(clust_vect_1, clust_vect_2)
			batches <- c(rep(1,num_cells_b1),rep(2, num_cells_b2))
			temp_pmd <- get_percent_max_diff(batches, clusters)
			## now log it all
			iter_vect <- c(iter_vect, iter)
			pmd_vect <- c(pmd_vect, temp_pmd)
			num_shared_vect <- c(num_shared_vect, num_shared)
		}
	}
	out_df <- data.frame(iter = iter_vect,
					     num_non_shared_clusters = num_shared_vect, 
					     pmd = pmd_vect)
	return(out_df)
}



#' do_full_pmd_characterization
#' @description \code{do_full_pmd_characterization} Run the PMD characterization over a range of power conditions.
#' @param directory the directory for putting the output image
#' @include PercentMaximumDifference.R
#' @return Null
#' @importFrom ggplot2 ggplot geom_point geom_smooth aes
#' @importFrom grDevices dev.off png
#' @importFrom stats loess
#' @examples
#'    do_full_pmd_characterization()
#' @name do_full_pmd_characterization
#' @export
do_full_pmd_characterization<-function(directory=''){
	if (directory==''){
		directory=getwd()
	}
	pmd_bench_table_1k_vs_1k<-characterize_pmd()
	pmd_bench_table_10k_vs_1k<-characterize_pmd(num_cells_b1=10000)
	pmd_bench_table_10k_vs_100<-characterize_pmd(num_cells_b1=10000, num_cells_b2=500)
	pmd_bench_table_10k_vs_10k<-characterize_pmd(num_cells_b1=10000, num_cells_b2=10000)
	dataset_size_vect<-c(rep('1k_vs_1k',dim(pmd_bench_table_1k_vs_1k)[1]),
						rep('10k_vs_1k',dim(pmd_bench_table_10k_vs_1k)[1]),
						rep('10k_vs_500',dim(pmd_bench_table_10k_vs_100)[1]),
						rep('10k_vs_10k',dim(pmd_bench_table_10k_vs_10k)[1]))
	pmd_bench_table_combined<-rbind(pmd_bench_table_1k_vs_1k,
		                            pmd_bench_table_10k_vs_1k,
		                            pmd_bench_table_10k_vs_100,
		                            pmd_bench_table_10k_vs_10k)
	pmd_bench_table_combined<-cbind(dataset_size_vect,pmd_bench_table_combined)
	names(pmd_bench_table_combined)<-c('dataset_sizes',names(pmd_bench_table_1k_vs_1k))
	outfile <- paste(directory,"pmd_characterization.png",sep='/')
	print(outfile)
	png(outfile,width=3000,height=2100,res=600)
	p <- ggplot(pmd_bench_table_combined, aes(x=num_non_shared_clusters, y=pmd, color = dataset_sizes)) +
	            geom_point(shape=1) +    # Use hollow circles
	            geom_smooth(method=loess)   # Add linear regression line 
	                             #  (by default includes 95% confidence region)
	print(p)
	dev.off()
}

