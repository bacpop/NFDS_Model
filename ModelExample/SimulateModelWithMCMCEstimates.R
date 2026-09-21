# script for simulating the model forward with parameter estimates from MCMC
# script can be run locally or on the cluster
# script has one command line - the file of the MCMC output (rds)
# 

library(odin.dust)
library(mcstate)
library(coda)

## command line arguments
args <- commandArgs(trailingOnly = TRUE)

#stop script if no arguments
if(length(args)<1){
  print("Please let me know which MCMC output you want me to use to simulate the model forward!")
  print("As the first and only argument, give the file of the model fit with parameter samples from MCMC.")
  stop("Requires command line argument.")
}

fit_file <- args[1] # MCMC trace / output file
print(paste("Input file is ", fit_file, sep = ""))

simulate_model_for_plot <- function(mcmc_run, params_loc){
  simulated_data_long <- list()
  rand_ind <- sample(x = 1:nrow(mcmc_run), size = 200, replace = FALSE)
  WF_PPxSero <- odin.dust::odin_dust("NFDS_Model.R")
  
  empty_vec <- rep(0,params_loc$species_no)
  cluster_samples_ParamVar <- array(rep(empty_vec,length(rand_ind)),dim = c(length(rand_ind),params_loc$species_no))
  cluster_samples_ParamVar_sero <- array(rep(empty_vec,length(rand_ind)),dim = c(length(rand_ind),params_loc$sero_no))
  
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
                                      seed = 1L, deterministic = FALSE)
    ModelSimResult <- WFmodel_ppxSero$run(5*12)
    cluster_samples_ParamVar[i,] <- (ModelSimResult[(2:(params_loc$species_no+1)),]) # time point 0 is 2014, simulate until 2019
    cluster_samples_ParamVar_sero[i,] <- rowSums(matrix(ModelSimResult[-(1:(params_loc$species_no+1)),], nrow = params_loc$sero_no, ncol = params_loc$species_no, byrow = TRUE)) # time point 0 is 2014, simulate until 2019
    
  }
  list(cluster_samples_ParamVar,cluster_samples_ParamVar_sero)
}

vacc_time <- 4 # time of vaccination, in years after start of dataset
PPsero_startpop <- readRDS("PPsero_startpop.rds") # start population for model (GPSC x serotypes)
# rows are GPSCs, columns are serotypes
# based on real data (counts) from first time point of dataset
# created by Poisson sampling from real data and scaling it up to model population size (15,000)
# GPSC-serotype combinations that appear in the dataset - but not in the first time point - are set to 1
no_GPSC <- nrow(PPsero_startpop) # number of GPSCs
no_sero <- ncol(PPsero_startpop) # number of serotypes
Genotypes_matrix <- readRDS("Genotypes_matrix.rds")
delta_test <- readRDS("gene_delta_ranking.rds") # delta statistic (computed as in Corander et al.)
gene_no_test <- length(delta_test) # number of intermediate-frequency genes
migMatr_test_mtx <- readRDS("PPsero_mig.rds") # immigration matrix (rows are GPSCs, columns are serotypes)
vaccTypes_test <- readRDS(file = "SeroVT.rds") # information on which types are affected by vaccine (1 = serotype affected by vaccine, 0 = serotypes not affected)
dt_test <- 1/12 # determines how many generations there are between two data points (here: 12 generations per year)

#intermed_gene_presence_absence_consensus_matrix <- sapply(intermed_gene_presence_absence_consensus[-1,-1],as.double)


# read in MCMC trace
fit_4param_mcmc2_probs <- readRDS(fit_file) # something like "det_mcmc2.rds"
# filter out burn-in (and if desired thin chains)
fit_4param_mcmc2_probs_pars <- fit_4param_mcmc2_probs[,-(1:3)]

# parameters for model, just fitted parameters (sigma_f, prop_f, v, m) are set to NA
params_4_woFit <- list(dt = dt_test, species_no = no_GPSC, sero_no = no_sero, gene_no = gene_no_test, Pop_ini = PPsero_startpop, Pop_eq = rowSums(PPsero_startpop), capacity = sum(PPsero_startpop), Genotypes = Genotypes_matrix, delta = (delta_test), migVec = migMatr_test_mtx, vaccTypes = vaccTypes_test, vacc_time = 4, sigma_f = NA, prop_f = NA, m = NA, v = NA)
# simulate model forward using parameters from MCMC
params_4_model_data_list <- simulate_model_for_plot(fit_4param_mcmc2_probs_pars, params_4_woFit)
params_4_model_data_GPSCS <- params_4_model_data_list[[1]] # simulated GPSC counts
params_4_model_data_sero <- params_4_model_data_list[[2]] # simulated serotype counts
  
params_4_model_data_rel_GPSC <- params_4_model_data_GPSCS/rowSums(params_4_model_data_GPSCS) # GPSC frequencies
saveRDS(params_4_model_data_rel_GPSC, "params_4_model_data_rel_GPSC.rds")
  
params_4_model_data_rel_sero <- params_4_model_data_sero/rowSums(params_4_model_data_sero) # sero frequencies
saveRDS(params_4_model_data_rel_sero, "params_4_model_data_rel_sero.rds")
