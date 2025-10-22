library(odin.dust)
library(mcstate)
library(coda)

## command line arguments
args <- commandArgs(trailingOnly = TRUE)

#stop script if no arguments
if(length(args)==0){
  print("This is a script to simulate the model based on fitting the model to simulated data. It should be run as a job array.")
  print("As the first argument, specify location of simulated dataset.")
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
  #dt_test <- 1
}

job_array_filename <- args[2] # job_array_filename <- "sim_fit_budget_fits_to_simulation.tsv"
row_index <- as.integer(args[3]) # row_index <- 1
print(job_array_filename)
print(row_index)
#output_folder <- args[4]

job_array_df = read.table(job_array_filename, header=FALSE)

# load the list of codons
budget <- job_array_df[row_index,2] # budget of sampling
sample_freq <- job_array_df[row_index,3] # budget of sampling
budget_dict <- c("01", "025", "05","1","2_5","5", "10", "25")
names(budget_dict) <- c("0.1", "0.25", "0.5","1","2.5","5", "10", "25")
rep_no <- job_array_df[row_index,4] # replicate number (1,2,3)
if(rep_no == 1){
  filepath <- "/nfs/research/jlees/leonie/WF_fitting_2025/September/ModelSimulation_newPriorSampleYearsAndSize/"
} else if(rep_no == 2){
  filepath <- "/nfs/research/jlees/leonie/WF_fitting_2025/September/ModelSimulation_newPriorSampleYearsAndSize/"
} else if(rep_no == 3){
  filepath <- "/nfs/research/jlees/leonie/WF_fitting_2025/September/ModelSimulation_newPriorSampleYearsAndSize/"
} else{
  print("replicate number invalid")
}
file_name <- paste("US_based_simulation_with_freq_", as.character(sample_freq), "_", budget_dict[as.character(budget)], "percent_4param_det_pmcmc_run2_simulation.rds", sep = "")
#m_fit <- log(job_array_df[row_index,4])
#v_fit <- job_array_df[row_index,5]
#sigmaf_fit <- log(job_array_df[row_index,6])
#propf_fit <- job_array_df[row_index,7]



simulate_model_for_plot2 <- function(ModSim1, params_fixed){
  simulated_data_long <- list()
  WF_PPxSero <- odin.dust::odin_dust("NFDS_Model_PPxSero.R")
  
  rand_ind <- sample(x = 1:nrow(ModSim1), size = 200, replace = FALSE)
  
  empty_vec <- rep(0,params_fixed$species_no)
  cluster_samples_ParamVar <- array(rep(empty_vec,20 * length(rand_ind)),dim = c(20*length(rand_ind),params_fixed$species_no))
  #cluster_samples_ParamVar <- array(rep(empty_vec,20),dim = c(20,params_fixed$species_no))
  
  for (i in 1:length(rand_ind)) {
    pars <- ModSim1[rand_ind[i],]
    params_fixed$sigma_f <- pars[1]
    params_fixed$prop_f <- pars[2]
    params_fixed$m <- pars[3]
    params_fixed$v <- pars[4]
    
    WFmodel_ppxSero <- WF_PPxSero$new(pars = params_fixed,
                                      time = 0,
                                      n_particles = 1L,
                                      n_threads = 4L,
                                      seed = 1L)
    for (j in 1:20) {
      cluster_samples_ParamVar[(i-1) * 20 + j,] <- (WFmodel_ppxSero$run(j*12)[(2:(params_fixed$species_no+1)),])
    }
  }
  cluster_samples_ParamVar
}


dt <- 1/12
vacc_time <- 1
species_no <- species_no_test
no_clusters <- species_no_test
gene_no <- nrow(Genotypes_test_matrix)
Pop_ini <- data.frame(Pop_ini_test_mtx)
Pop_eq <- rowSums(Pop_ini_test_mtx)
capacity <- sum(Pop_ini_test_mtx)
delta <- delta_test
vaccTypes <- vaccTypes_test_mtx
migVec <- data.frame(migMatr_test_mtx)


params_4_woFit <- list(species_no = species_no, Pop_ini = as.matrix(Pop_ini), Pop_eq = (Pop_eq), Genotypes = Genotypes_test_matrix, capacity = capacity, delta = delta, vaccTypes = vaccTypes, gene_no = gene_no, vacc_time = vacc_time, dt = dt, migVec = as.matrix(migVec), sero_no = sero_no_test, sigma_f = NA, prop_f = NA, m = NA, v = NA)
#params_fit <- c("m" = m_fit, "v" = v_fit, "propf" = propf_fit, "sigmaf" = sigmaf_fit)
ModSim_probs_1 <- readRDS(paste(filepath, file_name, sep = ""))
ModSim_probs_1_pars <- mcstate::pmcmc_thin(ModSim_probs_1, burnin = 5000, thin = 1)$pars
params_4_model_data <- simulate_model_for_plot2(ModSim_probs_1_pars, params_4_woFit)

params_4_model_data_rel <- params_4_model_data/rowSums(params_4_model_data)

# plot(params_4_model_data_rel[,1], ylim = c(0,.25))
# lines(params_4_model_data_rel[,2])

saveRDS(params_4_model_data_rel, paste(row_index, "_", "params_4_model_data_rel_paramdist.rds", sep = ""))
