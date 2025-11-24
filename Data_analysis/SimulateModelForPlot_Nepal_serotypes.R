# script can be run locally or on the cluster
# use commented commands for local version (it's pretty fast, so no need for the cluster really)

library(odin.dust)
library(mcstate)
library(coda)

## command line arguments
args <- commandArgs(trailingOnly = TRUE)

#stop script if no arguments
if(length(args)<2){
  print("Please let me know which version of the model you want to run!")
  print("As the first argument, 2 or 4, reflecting the parameters.")
  print("As the second argument, give the file of the model fit with parameter samples from MCMC.")
  stop("Requires command line argument.")
}

model_version <- NA
model_version <- args[1]
print(paste("Model version set to ", model_version, sep = ""))
fit_file <- args[2]
print(paste("Input file is ", fit_file, sep = ""))

simulate_model_for_plot2 <- function(mcmc_run, params_loc){
  simulated_data_long <- list()
  rand_ind <- sample(x = 1:nrow(mcmc_run), size = 200, replace = FALSE)
  WF_PPxSero <- odin.dust::odin_dust("NFDS_Model_PPxSero.R")
  
  empty_vec <- rep(0,params_loc$sero_no)
  cluster_samples_ParamVar <- array(rep(empty_vec,length(rand_ind)),dim = c(length(rand_ind),params_loc$sero_no))
  
  for (i in 1:length(rand_ind)) {
    pars <- mcmc_run[rand_ind[i],]
    params_loc$sigma_f <- pars[1]
    params_loc$prop_f <- pars[2]
    params_loc$m <- pars[3]
    params_loc$v <- pars[4]
    
    WFmodel_ppxSero <- WF_PPxSero$new(pars = params_loc,
                                      time = 0,
                                      n_particles = 1L,
                                      n_threads = 4L,
                                      seed = 1L)
    cluster_samples_ParamVar[i,] <- rowSums(matrix(WFmodel_ppxSero$run(5*12)[-(1:(params_loc$species_no+1)),], nrow = params_loc$sero_no, ncol = params_loc$species_no, byrow = TRUE)) # time point 0 is 2014, simulate until 2019
  }
  cluster_samples_ParamVar
}

simulate_model_for_plot2_null <- function(mcmc_run, params_loc){
  simulated_data_long <- list()
  rand_ind <- sample(x = 1:nrow(mcmc_run), size = 200, replace = FALSE)
  WF_PPxSero <- odin.dust::odin_dust("NFDS_Model_PPxSero.R")
  
  empty_vec <- rep(0,params_loc$sero_no)
  cluster_samples_ParamVar <- array(rep(empty_vec,length(rand_ind)),dim = c(length(rand_ind),params_loc$sero_no))
  
  for (i in 1:length(rand_ind)) {
    pars <- mcmc_run[rand_ind[i],]
    params_loc$m <- pars[1]
    params_loc$v <- pars[2]
    
    WFmodel_ppxSero <- WF_PPxSero$new(pars = params_loc,
                                      time = 0,
                                      n_particles = 1L,
                                      n_threads = 4L,
                                      seed = 1L)
    cluster_samples_ParamVar[i,] <- rowSums(matrix(WFmodel_ppxSero$run(5*12)[-(1:(params_loc$species_no+1)),], nrow = params_loc$sero_no, ncol = params_loc$species_no, byrow = TRUE))
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
mass_VT <- readRDS(file = "Nepal_SeroVT_6A.rds")
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


if(model_version == "4"){
  # 4-param model
  fit_4param_mcmc2_probs <- readRDS(fit_file) # something like "PPxSero_ggCaller_PopPUNK_4param_det_pmcmc_run2.rds"
  # fit_4param_mcmc2_probs <- readRDS("/Users/llorenz/Documents/PhD_Project/Code/1st_project/WF_plots_postTAC2/2025_11_05/Nepal_GPSC_VT_NewMetaData_6A_preFitParams_v3_BasedOnNepal_NewMetaData_6A_preFitParams/Nepal_PPxSero_ggCaller_PopPUNK_4param_det_pmcmc_run2.rds")
  fit_4param_mcmc2_probs_pars <- mcstate::pmcmc_thin(fit_4param_mcmc2_probs, burnin = 5000, thin = 1)$pars
  
  params_4_woFit <- list(species_no = species_no, Pop_ini = as.matrix(Pop_ini), Pop_eq = (Pop_eq), Genotypes = intermed_gene_presence_absence_consensus_matrix, capacity = capacity, delta = delta, vaccTypes = vaccTypes, gene_no = gene_no, vacc_time = vacc_time, dt = dt, migVec = as.matrix(migVec), sero_no = sero_no, sigma_f = NA, prop_f = NA, m = NA, v = NA)
  params_4_model_data <- simulate_model_for_plot2(fit_4param_mcmc2_probs_pars, params_4_woFit)
  
  params_4_model_data_rel <- params_4_model_data/rowSums(params_4_model_data)
  saveRDS(params_4_model_data_rel, "params_4_model_data_rel.rds")
  # saveRDS(params_4_model_data_rel, "/Users/llorenz/Documents/PhD_Project/Code/1st_project/WF_plots_postTAC2/2025_10_23/US_SimulateModelForPlots/params_4_model_data_rel.rds")
} else if(model_version == "2"){
  # Null model
  fit_2param_mcmc2_probs <- readRDS(fit_file) # something like "PPxSero_ggCaller_PopPUNK_2param_det_pmcmc_run2.rds"
  # fit_2param_mcmc2_probs <- readRDS("/Users/llorenz/Documents/PhD_Project/Code/1st_project/WF_plots_postTAC2/2025_08_04/US_newPrior/PPxSero_ggCaller_PopPUNK_Null_det_pmcmc_run2.rds")
  fit_2param_mcmc2_probs_pars <- mcstate::pmcmc_thin(fit_2param_mcmc2_probs, burnin = 5000, thin = 1)$pars
  
  params_null_woFit = list(species_no = species_no, Pop_ini = as.matrix(Pop_ini), Pop_eq = (Pop_eq), Genotypes = intermed_gene_presence_absence_consensus_matrix, capacity = capacity, delta = delta, vaccTypes = vaccTypes, gene_no = gene_no, vacc_time = vacc_time, dt = dt, migVec = as.matrix(migVec), sero_no = sero_no, sigma_f = -1000, prop_f = 1, m = NA, v = NA)
  params_null_model_data <- simulate_model_for_plot2_null(fit_2param_mcmc2_probs_pars, params_null_woFit)
  
  params_null_model_data_rel <- params_null_model_data/rowSums(params_null_model_data)
  saveRDS(params_null_model_data_rel, "params_null_model_data_rel.rds")
  # saveRDS(params_null_model_data_rel, "/Users/llorenz/Documents/PhD_Project/Code/1st_project/WF_plots_postTAC2/2025_10_23/US_SimulateModelForPlots/params_2_model_data_rel.rds")
}
