library(odin.dust)
library(mcstate)
library(coda)

simulate_model_for_plot2 <- function(Nepal_mcmc_run, Nepal_params_loc){
  simulated_data_long <- list()
  rand_ind <- sample(x = 1:nrow(Nepal_mcmc_run), size = 200, replace = FALSE)
  WF_PPxSero <- odin.dust::odin_dust("NFDS_Model_PPxSero.R")
  
  empty_vec <- rep(0,Nepal_params_loc$species_no)
  cluster_samples_ParamVar <- array(rep(empty_vec,length(rand_ind)),dim = c(length(rand_ind),Nepal_params_loc$species_no))
  
  for (i in 1:length(rand_ind)) {
    pars <- Nepal_mcmc_run[rand_ind[i],4:7]
    Nepal_params_loc$sigma_f <- pars[1]
    Nepal_params_loc$prop_f <- pars[2]
    Nepal_params_loc$m <- pars[3]
    Nepal_params_loc$v <- pars[4]
    
    WFmodel_ppxSero <- WF_PPxSero$new(pars = Nepal_params_loc,
                                      time = 0,
                                      n_particles = 1L,
                                      n_threads = 4L,
                                      seed = 1L)
    cluster_samples_ParamVar[i,] <- (WFmodel_ppxSero$run(5*12)[(2:(Nepal_params_loc$species_no+1)),]) # time point 0 is 2014, then simulate 5 years (until 2019)
  }
  cluster_samples_ParamVar
}

simulate_model_for_plot2_null <- function(Nepal_mcmc_run, Nepal_params_loc){
  simulated_data_long <- list()
  rand_ind <- sample(x = 1:nrow(Nepal_mcmc_run), size = 200, replace = FALSE)
  WF_PPxSero <- odin.dust::odin_dust("NFDS_Model_PPxSero.R")
  
  empty_vec <- rep(0,Nepal_params_loc$species_no)
  cluster_samples_ParamVar <- array(rep(empty_vec,length(rand_ind)),dim = c(length(rand_ind),Nepal_params_loc$species_no))
  
  for (i in 1:length(rand_ind)) {
    pars <- Nepal_mcmc_run[rand_ind[i],4:5]
    Nepal_params_loc$m <- pars[1]
    Nepal_params_loc$v <- pars[2]
    
    WFmodel_ppxSero <- WF_PPxSero$new(pars = Nepal_params_loc,
                                      time = 0,
                                      n_particles = 1L,
                                      n_threads = 4L,
                                      seed = 1L)
    cluster_samples_ParamVar[i,] <- (WFmodel_ppxSero$run(5*12)[(2:(Nepal_params_loc$species_no+1)),])
  }
  cluster_samples_ParamVar
}



seq_clusters <- readRDS("Nepal_PP.rds")
intermed_gene_presence_absence_consensus <- readRDS(file = "Nepal_ggCaller_intermed_consensus.rds")
intermed_gene_presence_absence_consensus_matrix <- sapply(intermed_gene_presence_absence_consensus[-1,-1],as.double)
delta_ranking <- readRDS(file = "Nepal_delta_ranking.rds")
mass_clusters <- length(unique(seq_clusters$GPSC))
sero_no = length(unique(seq_clusters$Serotype))
model_start_pop <- readRDS(file = "Nepal_PPsero_startpop.rds")
mass_VT <- readRDS(file = "Nepal_SeroVT.rds")
mass_clusters <- length(unique(seq_clusters$GPSC))
avg_cluster_freq <- readRDS(file = "Nepal_PPsero_mig.rds")
dt <- 1/12
vacc_time <- 1
species_no <- mass_clusters
no_clusters <- mass_clusters
gene_no <- nrow(intermed_gene_presence_absence_consensus_matrix)
Pop_ini <- data.frame(model_start_pop)
Pop_eq <- rowSums(model_start_pop)
capacity <- sum(model_start_pop)
delta <- delta_ranking
vaccTypes <- mass_VT
migVec <- data.frame(avg_cluster_freq)

# 4-param model
Nepal_4param_mcmc2_probs <- readRDS("Nepal_4param_mcmc2_probs.rds")
#Nepal_param_mean4 <- apply(Nepal_4param_mcmc2_probs, 2, mean)
#Nepal_params_4 = list(species_no = species_no, Pop_ini = as.matrix(Pop_ini), Pop_eq = (Pop_eq), Genotypes = intermed_gene_presence_absence_consensus_matrix, capacity = capacity, delta = delta, vaccTypes = vaccTypes, gene_no = gene_no, vacc_time = vacc_time, dt = dt, migVec = as.matrix(migVec), sero_no = sero_no, sigma_f = (Nepal_param_mean4[4]), prop_f = (Nepal_param_mean4[5]), m = (Nepal_param_mean4[6]), v = (Nepal_param_mean4[7]))
#Nepal_params_4_model_data <- simulate_model_for_plot(Nepal_params_4)
Nepal_params_4_woFit <- list(species_no = species_no, Pop_ini = as.matrix(Pop_ini), Pop_eq = (Pop_eq), Genotypes = intermed_gene_presence_absence_consensus_matrix, capacity = capacity, delta = delta, vaccTypes = vaccTypes, gene_no = gene_no, vacc_time = vacc_time, dt = dt, migVec = as.matrix(migVec), sero_no = sero_no, sigma_f = NA, prop_f = NA, m = NA, v = NA)
Nepal_params_4_model_data <- simulate_model_for_plot2(Nepal_4param_mcmc2_probs, Nepal_params_4_woFit)

Nepal_params_4_model_data_rel <- Nepal_params_4_model_data/rowSums(Nepal_params_4_model_data)
saveRDS(Nepal_params_4_model_data_rel, "Nepal_params_4_model_data_rel.rds")

#mean_Nepal_params_4 <- apply(Nepal_params_4_model_data_rel,2,mean)
#Nepal_params_4_percentile95_up <- apply(Nepal_params_4_model_data_rel,2,quantile,.975)
#Nepal_params_4_percentile95_low <- apply(Nepal_params_4_model_data_rel,2,quantile,.025)

# Null model
Nepal_2param_mcmc2_probs <- readRDS("Nepal_2param_mcmc2_probs.rds")

Nepal_params_null_woFit = list(species_no = species_no, Pop_ini = as.matrix(Pop_ini), Pop_eq = (Pop_eq), Genotypes = intermed_gene_presence_absence_consensus_matrix, capacity = capacity, delta = delta, vaccTypes = vaccTypes, gene_no = gene_no, vacc_time = vacc_time, dt = dt, migVec = as.matrix(migVec), sero_no = sero_no, sigma_f = -1000, prop_f = 1, m = NA, v = NA)
Nepal_params_null_model_data <- simulate_model_for_plot2_null(Nepal_2param_mcmc2_probs, Nepal_params_null_woFit)

Nepal_params_null_model_data_rel <- Nepal_params_null_model_data/rowSums(Nepal_params_null_model_data)
saveRDS(Nepal_params_null_model_data_rel, "Nepal_params_null_model_data_rel.rds")
#mean_Nepal_params_null <- apply(Nepal_params_null_model_data_rel,2,mean)
#Nepal_params_null_percentile95_up <- apply(Nepal_params_null_model_data_rel,2,quantile,.975)
#Nepal_params_null_percentile95_low <- apply(Nepal_params_null_model_data_rel,2,quantile,.025)