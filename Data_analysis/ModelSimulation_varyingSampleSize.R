
simulated_data_over_time <- readRDS("simulated_data_over_time_filtered_100percent.rds")

# create sampled datasets
sample_dataset_varied <- function(data_all, sample_size_vec){
  sample_size_vec_long <- c(rep(sample_size_vec[1],4), rep(sample_size_vec[2],4), rep(sample_size_vec[3],4), rep(sample_size_vec[4],4), rep(sample_size_vec[5],4))
  new_data <- data.frame(matrix(NA, nrow = nrow(data_all), ncol = ncol(data_all)))
  for (j in 1:ncol(data_all)) {
    new_data[,j] <- tabulate(sample(length(data_all[,j]), size=sample_size_vec_long[j] * sum(data_all[,j]), replace=TRUE, prob = as.vector(data_all[,j])), length(data_all[,j]))
  }
  new_data
}

# define combinations of sampling frequencies (choose 5 lowest sampling frequencies we tested, do one increasing, one decreasing and two random ones)
sample_size_vec_1 <- c(0.001, 0.0025, 0.005, 0.01, 0.025)
sample_size_vec_2 <- rev(sample_size_vec_1)
sample_size_vec_3 <- c(0.0050, 0.0100, 0.0025, 0.0250, 0.0010) # sample(sample_size_vec_1, 5, replace = FALSE)
sample_size_vec_4 <- c(0.0100, 0.0025, 0.0010, 0.0250, 0.0050) # sample(sample_size_vec_1, 5, replace = FALSE)

simulated_data_over_time_filtered_varied_increasing <- sample_dataset_varied(simulated_data_over_time, sample_size_vec_1)
simulated_data_over_time_filtered_varied_decreasing <- sample_dataset_varied(simulated_data_over_time, sample_size_vec_2)
simulated_data_over_time_filtered_varied_rand1 <- sample_dataset_varied(simulated_data_over_time, sample_size_vec_3)
simulated_data_over_time_filtered_varied_rand2 <- sample_dataset_varied(simulated_data_over_time, sample_size_vec_4)

saveRDS(simulated_data_over_time_filtered_varied_increasing, "simulated_data_over_time_filtered_varied_increasing.rds")
saveRDS(simulated_data_over_time_filtered_varied_decreasing, "simulated_data_over_time_filtered_varied_decreasing.rds")
saveRDS(simulated_data_over_time_filtered_varied_rand1, "simulated_data_over_time_filtered_varied_rand1.rds")
saveRDS(simulated_data_over_time_filtered_varied_rand2, "simulated_data_over_time_filtered_varied_rand2.rds")
