library(odin.dust)
#install.packages("mcstate")
#library(mcstate)
#install.packages("mcstate")
library(mcstate)
library(coda)


## command line arguments
args <- commandArgs(trailingOnly = TRUE)

#stop script if no arguments
if(length(args)==0){
  print("This is a script to fit models to simulated data.")
  print("As the first argument, specify location of simulated dataset.")
  print("As the second and third argument, location of the simulated data and the sampling frequency.")
  stop("Requires command line argument.")
}

if(args[1] == "US"){
  vacc_time <- 4
  #vacc_time <- 48
  PPsero_startpop6 <- readRDS("PPsero_startpop6.rds") 
  freq_clust <- which(rowSums(PPsero_startpop6)>(sum(PPsero_startpop6)/100))
  freq_sero <- which(colSums(PPsero_startpop6)>(sum(PPsero_startpop6)/100))
  species_no_test <- length(freq_clust) # 52 #20
  sero_no_test <- length(freq_sero) #32 #15
  Pop_ini_test_mtx <- matrix(PPsero_startpop6[freq_clust, freq_sero], nrow = length(freq_clust), ncol = length(freq_sero))
  
  Genotypes_test_matrix <- readRDS("Genotypes_test_matrix.rds")
  delta_test <- readRDS("delta_test.rds")
  gene_no_test <- length(delta_test)
  migMatr_test_mtx <- readRDS("migMatr_test_mtx.rds")
  vaccTypes_test_mtx1 <- readRDS(file = "SeroVT.rds")
  vaccTypes_test_mtx <- vaccTypes_test_mtx1[freq_sero]
  dt_test <- 1/12
}

if(length(args)>2){
  simulated_data <- readRDS(args[2])
  sampling_freq <- as.numeric(args[3])
  print(paste("Setting the sampling frequency to ", args[3]))
  years_avail <- 1:20
  years_for_sampling <- floor(max(years_avail)/sampling_freq)
  years_real <- years_avail[(1:years_for_sampling) * sampling_freq]
  
  fitting_test_data <- data.frame("year" = years_real, t(simulated_data[,years_real]))
  names(fitting_test_data) <- c("year", as.character(1:length(freq_clust)))
  fitting_test_data_2 <- mcstate::particle_filter_data(data = fitting_test_data,
                                                     time = "year",
                                                     rate = 1 / dt_test,
                                                     initial_time = 0)
  output_filename <- paste("US_based_simulation_with_freq_", as.character(sampling_freq),sep = "")
  
} else{
  print("Please give me the location of the simulated data and the sampling frequency.")
  stop("Requires command line argument.")
}

combined_compare <- function(state, observed, pars = NULL) {
  result <- 0
  #data_size <- sum(unlist(observed))
  data_size <- sum(unlist(observed[as.character(1:(length(unlist(observed))-4))]))
  model_size = sum(unlist(state[-1, , drop = TRUE]))
  exp_noise <- 1e6
  data_vals <- unlist(observed[as.character(1:(length(unlist(observed))-4))])
  #model_vals <- state[-1, , drop = TRUE]
  model_vals <- rep(0, length(unlist(observed))-4)
  data_missing <- FALSE
  for (i in 1:(length(unlist(observed))-4)){ 
    state_name <- paste("sum_clust", i, sep = "")
    model_vals[i] <- state[state_name, , drop = TRUE]
    if (is.na(observed[[as.character(i)]])) {
      #Creates vector of zeros in ll with same length, if no data
      #ll_obs <- numeric(length( state[state_name, , drop = TRUE]))
      data_missing <- TRUE
    } 
  }
  models_vals_err <- model_vals + rexp(n = length(model_vals), rate = exp_noise)
  if(data_missing){
    ll_obs <- 0
  }
  else{
    ll_obs <- dmultinom(x = (data_vals), prob = models_vals_err/model_size, log = TRUE)   
  }
  result <- ll_obs
  #for (i in 1:(length(unlist(observed))-4)){ 
  #  state_name <- paste("sum_clust", i, sep = "")
  #  if (is.na(observed[[as.character(i)]])) {
  #    #Creates vector of zeros in ll with same length, if no data
  #    ll_obs <- numeric(length( state[state_name, , drop = TRUE]))
  #  } else {
  #lambda <-  state[state_name, , drop = TRUE]/model_size * data_size + rexp(n = length( state[state_name, , drop = TRUE]/model_size * data_size), rate = exp_noise)
  #ll_obs <- dpois(x = observed[[as.character(i)]], lambda = lambda, log = TRUE)
  #    ll_obs <- dmultinom(x = (data_vals), prob = model_vals/model_size, log = TRUE)
  #  }
  
  #  result <- result + ll_obs
  #}
  result
}

Fit_model_to_sim_data <- function(fitting_data, fitting_params){
  #WF <- odin.dust::odin_dust("NFDS_Model_PPxSero.R")
  WF <- odin.dust::odin_dust("NFDS_Model_PPxSero_NewSigma.R")
  
  #print(fitting_mass_data)
  det_filter <- particle_deterministic$new(data = fitting_data,
                                           model = WF,
                                           compare = combined_compare)
  
  #complex_params = list(species_no = species_no, Pop_ini = Pop_ini, Pop_eq = Pop_eq, Genotypes = intermed_gene_presence_absence_consensus, capacity = capacity, delta = delta, vaccTypes = vaccTypes, gene_no = gene_no, vacc_time = vacc_time, dt = dt, sigma_w = pmcmc_sigma_w, migVec = (migVec), sero_no = sero_no)
  complex_params = fitting_params
  #complex_params <- c(Pop_ini, Pop_eq, Genotypes, capacity, delta, species_no, gene_no, vacc_time, dt, migVec,vT)
  
  make_transform <- function(m) {
    function(theta) {
      as_double_mtx <- function(x){
        sapply(x,as.double)
      }
      c(lapply(m, as_double_mtx), as.list(theta))
    }
  }
  
  
  proposal_matrix <- diag(0.1,4) # the proposal matrix defines the covariance-variance matrix for a mult normal dist
  #mcmc_pars$names()
  #mcmc_pars$model(mcmc_pars$initial())
  # read this: https://mrc-ide.github.io/mcstate/reference/pmcmc_parameters.html
  # it explains how to not fit all parameters but just the ones I want
  # non-scalar parameters have to be transformed for this.
  
  index <- function(info) {
    list(run = c(sum_clust = info$index$Pop_tot),
         state = c(Pop = info$index$Pop))
  }
  #mcmc_pars <- mcstate::pmcmc_parameters$new(list(mcstate::pmcmc_parameter("sigma_f", -0.597837, min = -1000, max = 0), mcstate::pmcmc_parameter("prop_f", 0.125, min = 0, max = 1), mcstate::pmcmc_parameter("m", -4, min = -1000, max = 0), mcstate::pmcmc_parameter("v", 0.05, min = 0, max = 1)), proposal_matrix, make_transform(complex_params))
  proposal_matrix <- diag(c(exp(1), 0.1, exp(1), 0.1))
  mcmc_pars <- mcstate::pmcmc_parameters$new(list(mcstate::pmcmc_parameter("sigma_f", runif(n=1, min=-10, max=0), min = -1000, max = 0), mcstate::pmcmc_parameter("prop_f", runif(n=1, min=0, max=1), min = 0, max = 1), mcstate::pmcmc_parameter("m", runif(n=1, min=-10, max=0), min = -1000, max = 0), mcstate::pmcmc_parameter("v", runif(n=1, min=0, max=1), min = 0, max = 1)), proposal_matrix, make_transform(complex_params))
  mcmc_pars <- mcstate::pmcmc_parameters$new(list(mcstate::pmcmc_parameter("sigma_f", runif(n=1, min=-10, max=0), min = -1000, max = 0, prior = function(a) 1/a), mcstate::pmcmc_parameter("prop_f", runif(n=1, min=0, max=1), min = 0, max = 1, prior = function(a) a), mcstate::pmcmc_parameter("m", runif(n=1, min=-10, max=0), min = -1000, max = 0, prior = function(a) 1/a), mcstate::pmcmc_parameter("v", runif(n=1, min=0, max=1), min = 0, max = 1, , prior = function(a) a)), proposal_matrix, make_transform(complex_params))
  
  #proposal_matrix <- diag(0.1,1)
  #mcmc_pars <- mcstate::pmcmc_parameters$new(list(mcstate::pmcmc_parameter("v", runif(n=1, min=0, max=1), min = 0, max = 1)), proposal_matrix, make_transform(complex_params))
  
  mcmc_pars$initial()
  #mcmc_pars$model(mcmc_pars$initial())
  
  #WF$public_methods$has_openmp()
  det_filter <- particle_deterministic$new(data = fitting_data,
                                           model = WF,
                                           index = index,
                                           compare = combined_compare)
  
  n_steps <- 5
  n_burnin <- 0
  
  
  control <- mcstate::pmcmc_control(
    n_steps,
    save_state = TRUE, 
    save_trajectories = TRUE,
    progress = TRUE,
    adaptive_proposal = TRUE,
    n_chains = 1)
  det_pmcmc_run <- mcstate::pmcmc(mcmc_pars, det_filter, control = control)
  
  n_steps <- 1000
  #n_steps <- 200
  n_burnin <- 0
  
  
  control <- mcstate::pmcmc_control(
    n_steps,
    save_state = TRUE, 
    save_trajectories = TRUE,
    progress = TRUE,
    adaptive_proposal = TRUE,
    n_chains =4, n_workers = 4,
    n_threads_total = 4)
  
  #n_chains = 8, n_workers = 8,
  #n_threads_total = 8
  
  det_pmcmc_run <- mcstate::pmcmc(mcmc_pars, det_filter, control = control)
  processed_chains <- mcstate::pmcmc_thin(det_pmcmc_run, burnin = 200, thin = 1)
  parameter_mean_hpd <- apply(processed_chains$pars, 2, mean)
  print(parameter_mean_hpd)
  
  det_mcmc1 <- coda::as.mcmc(cbind(det_pmcmc_run$probabilities, det_pmcmc_run$pars))
  #pdf(file = paste(output_filename,"4param_det_mcmc1.pdf",sep = "_"),   # The directory you want to save the file in
  #    width = 6, # The width of the plot in inches
  #    height = 12)
  #plot(det_mcmc1)
  #dev.off()
  print("det_mcmc_1 final log likelihood")
  processed_chains$probabilities[nrow(processed_chains$probabilities),2]
  print("det_mcmc_1 mean log likelihood")
  mean(processed_chains$probabilities[,2])
  det_proposal_matrix <- cov(processed_chains$pars)
  
  det_mcmc_pars <- mcstate::pmcmc_parameters$new(list(mcstate::pmcmc_parameter("sigma_f", parameter_mean_hpd[1], min = -1000, max = 0), mcstate::pmcmc_parameter("prop_f", parameter_mean_hpd[2], min = 0, max = 1),mcstate::pmcmc_parameter("m", parameter_mean_hpd[3], min = -1000, max = 0), mcstate::pmcmc_parameter("v", parameter_mean_hpd[4], min = 0, max = 1)), det_proposal_matrix, make_transform(complex_params))
  mcmc_pars <- mcstate::pmcmc_parameters$new(list(mcstate::pmcmc_parameter("sigma_f",  parameter_mean_hpd[1], min = -1000, max = 0, prior = function(a) 1/a), mcstate::pmcmc_parameter("prop_f",  parameter_mean_hpd[2], min = 0, max = 1, prior = function(a) a), mcstate::pmcmc_parameter("m",  parameter_mean_hpd[3], min = -1000, max = 0, prior = function(a) 1/a), mcstate::pmcmc_parameter("v",  parameter_mean_hpd[4], min = 0, max = 1, , prior = function(a) a)), proposal_matrix, make_transform(complex_params))
  
  det_filter <- particle_deterministic$new(data = fitting_data,
                                           model = WF,
                                           index = index,
                                           compare = combined_compare)
  
  n_steps <- 5
  n_burnin <- 0
  
  
  control <- mcstate::pmcmc_control(
    n_steps,
    save_state = TRUE, 
    save_trajectories = TRUE,
    progress = TRUE,
    adaptive_proposal = TRUE,
    n_chains = 1)
  det_pmcmc_run <- mcstate::pmcmc(det_mcmc_pars, det_filter, control = control)
  
  n_steps <- 20000
  n_steps <- 2000
  n_burnin <- 0
  
  
  control <- mcstate::pmcmc_control(
    n_steps,
    save_state = TRUE, 
    save_trajectories = TRUE,
    progress = TRUE,
    adaptive_proposal = TRUE,
    n_chains = 4, n_workers = 4, n_threads_total = 4)
  det_pmcmc_run2 <- mcstate::pmcmc(det_mcmc_pars, det_filter, control = control)
  processed_chains <- mcstate::pmcmc_thin(det_pmcmc_run2, burnin = 500, thin = 1)
  parameter_mean_hpd <- apply(processed_chains$pars, 2, mean)
  print(parameter_mean_hpd)
  par(mfrow = c(1,1))
  
  det_mcmc2 <- coda::as.mcmc(cbind(processed_chains$probabilities, processed_chains$pars))
  percentile95_low <- apply(det_mcmc2,2,quantile,.025)
  percentile95_up <- apply(det_mcmc2,2,quantile,.975)
  print("lower 95 percentile")
  print(percentile95_low)
  print("upper 95 percentile")
  print(percentile95_up)
  pdf(file = paste(output_filename,"_4param_det_mcmc2_simulation.pdf",sep = "_"),   # The directory you want to save the file in
      width = 6, # The width of the plot in inches
      height = 12)
  plot(det_mcmc2)
  dev.off()
  print("det_mcmc_2 final log likelihood")
  processed_chains$probabilities[nrow(processed_chains$probabilities),2]
  print("det_mcmc_2 mean log likelihood")
  mean(processed_chains$probabilities[,2])
  saveRDS(det_pmcmc_run2, paste(output_filename, "_4param_det_pmcmc_run2_simulation.rds", sep = ""))
  return(det_pmcmc_run2)
}

fitting_test_params <- list(species_no = species_no_test, Pop_ini = data.frame(Pop_ini_test_mtx), Pop_eq = rowSums(Pop_ini_test_mtx), Genotypes = as.data.frame(Genotypes_test_matrix), capacity = sum(Pop_ini_test_mtx), delta = delta_test, vaccTypes = vaccTypes_test_mtx, gene_no = gene_no_test, vacc_time = 4, dt = dt_test, pmcmc_sigma_w = -1000, migVec = data.frame(migMatr_test_mtx), sero_no = sero_no_test)

#fitting_test_params <- list(species_no = species_no_test, Pop_ini = data.frame(Pop_ini_test_mtx), Pop_eq = rowSums(Pop_ini_test_mtx), Genotypes = as.data.frame(Genotypes_test_matrix), capacity = sum(Pop_ini_test_mtx), delta = delta_test, vaccTypes = vaccTypes_test_mtx, gene_no = gene_no_test, vacc_time = 4, dt = dt_test, pmcmc_sigma_w = -1000, migVec = data.frame(migMatr_test_mtx), sero_no = sero_no_test, sigma_f = log(0.0345), prop_f = 0.3, m = log(0.0133))

MCMC_chain_sample1year_run <- Fit_model_to_sim_data(fitting_data = fitting_test_data_2, fitting_params = fitting_test_params)

processed_chains_1year <- mcstate::pmcmc_thin(MCMC_chain_sample1year_run, burnin = 1000, thin = 2)
MCMC_chain_sample1year <- coda::as.mcmc(cbind(processed_chains_1year$probabilities, processed_chains_1year$pars))

mean_sample1year <- apply(MCMC_chain_sample1year,2,mean)
percentile95_low_sample1year <- apply(MCMC_chain_sample1year,2,quantile,.025)
percentile95_up_sample1year <- apply(MCMC_chain_sample1year,2,quantile,.975)
