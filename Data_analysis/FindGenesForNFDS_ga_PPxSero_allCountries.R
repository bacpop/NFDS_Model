### Loading packages
# install.packages("drat") # -- if you don't have drat installed
# drat:::add("ncov-ic")
# install.packages("odin.dust")
library(odin.dust)

## command line arguments
args <- commandArgs(trailingOnly = TRUE)
# import data

if(args[1] == "US"){
  seq_clusters <- readRDS("PopPUNK_clusters.rds")
  sero_no = length(unique(seq_clusters$Serotype))
  intermed_gene_presence_absence_consensus <- readRDS(file = "ggCPP_intermed_gene_presence_absence_consensus.rds")
  intermed_gene_presence_absence_consensus_matrix <- sapply(intermed_gene_presence_absence_consensus[-1,-1],as.double)
  #model_start_pop <- readRDS("PP_mass_cluster_freq_1_sero.rds")
  #model_start_pop <- model_start_pop / 133 * 15708
  model_start_pop <- readRDS("PPsero_startpop6.rds") 
  #model_start_pop <- readRDS("PP_mass_cluster_freq_1_sero.rds") # try using data directly
  #model_start_pop <- readRDS(file = "PPsero_startpop4.rds")
  #model_start_pop <- readRDS(file = "PPsero_startpop.rds")
  delta_ranking <- readRDS(file = "ggC_delta_ranking.rds")
  mass_cluster_freq_1 <- readRDS(file = "PP_mass_cluster_freq_1.rds")
  mass_cluster_freq_2 <- readRDS(file = "PP_mass_cluster_freq_2.rds")
  mass_cluster_freq_3 <- readRDS(file = "PP_mass_cluster_freq_3.rds")
  #mass_VT <- readRDS(file = "SeroVT.rds")
  mass_VT <- readRDS(file = "SeroVT.rds")
  mass_clusters <- length(unique(seq_clusters$Cluster))
  avg_cluster_freq <- readRDS(file = "PPsero_mig.rds")
  dt <- 1/36
  peripost_mass_cluster_freq <- data.frame("year" = c(1, 2), rbind(mass_cluster_freq_2, mass_cluster_freq_3))
  names(peripost_mass_cluster_freq) <- c("year", as.character(1:mass_clusters))
  
  data_list <- empty_list <- vector(mode = "list", length = 2)
  data_list[[1]] <- mass_cluster_freq_2
  data_list[[2]] <- mass_cluster_freq_3
  
  vacc_time <- 0
  output_filename <- "FindGenes_US"
  
  
  Pop_ini <- data.frame(model_start_pop)
  Pop_ini <- matrix(model_start_pop, nrow = mass_clusters, ncol = (sero_no))
  
  FindGenes_ggCPP_params <- list(dt = 1/36, species_no = mass_clusters,  gene_no = nrow(intermed_gene_presence_absence_consensus)-1, Pop_ini = Pop_ini, Pop_eq = rowSums(model_start_pop), capacity = sum(model_start_pop), Genotypes = intermed_gene_presence_absence_consensus_matrix, sigma_f = -3.6, sigma_w = 0, prop_f = 1, m = -4.2, migVec = avg_cluster_freq, vaccTypes = mass_VT, v = 0.087, vacc_time = 0, sero_no = sero_no)
  
  ### load model
  WF <- odin.dust::odin_dust("NFDS_Model_FindGenes_PPxSero.R")
  
} else if(args[1] == "Nepal"){
  seq_clusters <- readRDS("Nepal_PP.rds")
  intermed_gene_presence_absence_consensus <- readRDS(file = "Nepal_ggCaller_intermed_consensus.rds")
  intermed_gene_presence_absence_consensus_matrix <- sapply(intermed_gene_presence_absence_consensus[-1,-1],as.double)
  
  delta_ranking <- readRDS(file = "Nepal_delta_ranking.rds")
  #mass_cluster_freq_1 <- readRDS(file = "Nepal_cluster_freqs_1.rds")
  #mass_cluster_freq_2 <- readRDS(file = "Nepal_cluster_freqs_2.rds")
  #mass_cluster_freq_3 <- readRDS(file = "Nepal_cluster_freqs_3.rds")
  #mass_cluster_freq_4 <- readRDS(file = "Nepal_cluster_freqs_4.rds")
  #mass_cluster_freq_5 <- readRDS(file = "Nepal_cluster_freqs_5.rds")
  #mass_cluster_freq_6 <- readRDS(file = "Nepal_cluster_freqs_6.rds")
  #mass_cluster_freq_7 <- readRDS(file = "Nepal_cluster_freqs_7.rds")
  #mass_cluster_freq_8 <- readRDS(file = "Nepal_cluster_freqs_8.rds")
  #mass_cluster_freq_9 <- readRDS(file = "Nepal_cluster_freqs_9.rds")
  #mass_cluster_freq_10 <- readRDS(file = "Nepal_cluster_freqs_10.rds") # now ignoring the earlier because pre-vac is composed of first three pre-vac timepoints
  mass_cluster_freq_11 <- readRDS(file = "Nepal_cluster_freqs_11.rds")
  mass_cluster_freq_12 <- readRDS(file = "Nepal_cluster_freqs_12.rds")
  mass_cluster_freq_13 <- readRDS(file = "Nepal_cluster_freqs_13.rds")
  mass_cluster_freq_14 <- readRDS(file = "Nepal_cluster_freqs_14.rds")
  mass_cluster_freq_15 <- readRDS(file = "Nepal_cluster_freqs_15.rds")
  mass_cluster_freq_16 <- readRDS(file = "Nepal_cluster_freqs_16.rds")
  mass_cluster_freq_17 <- readRDS(file = "Nepal_cluster_freqs_17.rds")
  
  mass_clusters <- length(unique(seq_clusters$Cluster))
  sero_no = length(unique(seq_clusters$Serotype))
  
  model_start_pop <- readRDS(file = "Nepal_PPsero_startpop.rds")
  
  mass_VT <- readRDS(file = "Nepal_SeroVT.rds")
  mass_clusters <- length(unique(seq_clusters$GPSC))
  avg_cluster_freq <- readRDS(file = "Nepal_PPsero_mig.rds")
  dt <- 1/12
  #peripost_mass_cluster_freq <- data.frame("year" = 1:10, rbind(mass_cluster_freq_8,mass_cluster_freq_9, mass_cluster_freq_10, mass_cluster_freq_11,mass_cluster_freq_12, mass_cluster_freq_13, mass_cluster_freq_14,mass_cluster_freq_15,mass_cluster_freq_16,mass_cluster_freq_17))
  peripost_mass_cluster_freq <- data.frame("year" = 1:7, rbind(mass_cluster_freq_11,mass_cluster_freq_12, mass_cluster_freq_13, mass_cluster_freq_14,mass_cluster_freq_15,mass_cluster_freq_16,mass_cluster_freq_17))
  names(peripost_mass_cluster_freq) <- c("year", as.character(1:mass_clusters))
  
  data_list <- empty_list <- vector(mode = "list", length = 7)
  data_list[[1]] <- mass_cluster_freq_11
  data_list[[2]] <- mass_cluster_freq_12
  data_list[[3]] <- mass_cluster_freq_13
  data_list[[4]] <- mass_cluster_freq_14
  data_list[[5]] <- mass_cluster_freq_15
  data_list[[6]] <- mass_cluster_freq_16
  data_list[[7]] <- mass_cluster_freq_17
  
  #vacc_time <- 6
  vacc_time <- 1
  #cov_matrix <- readRDS("Nepal_cov_matrix.rds")
  output_filename <- "FindGenes_Nepal"
  
  Pop_ini <- data.frame(model_start_pop)
  Pop_ini <- matrix(model_start_pop, nrow = mass_clusters, ncol = (sero_no))
  
  # use params from Nepal fit from 09.07.2025
  FindGenes_ggCPP_params <- list(dt = 1/36, species_no = mass_clusters,  gene_no = nrow(intermed_gene_presence_absence_consensus)-1, Pop_ini = Pop_ini, Pop_eq = rowSums(model_start_pop), capacity = sum(model_start_pop), Genotypes = intermed_gene_presence_absence_consensus_matrix, sigma_f = -0.92227757, sigma_w = 0, prop_f = 1, m = -3.50741786, migVec = avg_cluster_freq, vaccTypes = mass_VT, v = 0.09500775, vacc_time = 0, sero_no = sero_no)
  ### load model
  WF <- odin.dust::odin_dust("NFDS_Model_FindGenes_PPxSero.R")
} else if(args[1] == "UK"){
  seq_clusters <- readRDS("UK_PP.rds")
  intermed_gene_presence_absence_consensus <- readRDS(file = "UK_ggCaller_intermed_consensus.rds")
  intermed_gene_presence_absence_consensus_matrix <- sapply(intermed_gene_presence_absence_consensus[-1,-1],as.double)
  mass_clusters <- length(unique(seq_clusters$Cluster))
  delta_ranking <- readRDS(file = "UK_delta_ranking.rds")
  mass_cluster_freq_1 <- readRDS(file = "UK_cluster_freqs_1.rds")
  mass_cluster_freq_2 <- readRDS(file = "UK_cluster_freqs_2.rds")
  mass_cluster_freq_3 <- readRDS(file = "UK_cluster_freqs_3.rds")
  mass_cluster_freq_4 <- readRDS(file = "UK_cluster_freqs_4.rds")
  mass_cluster_freq_5 <- readRDS(file = "UK_cluster_freqs_5.rds")
  mass_cluster_freq_6 <- readRDS(file = "UK_cluster_freqs_6.rds")
  mass_cluster_freq_7 <- readRDS(file = "UK_cluster_freqs_7.rds")
  #mass_VT1 <- readRDS(file = "UK_SeroVT.rds")
  
  avg_cluster_freq <- readRDS(file = "UK_PPsero_mig.rds")
  output_filename <- "UK_PPxSero_ggCaller_PopPUNK"
  
  mass_clusters <- length(unique(seq_clusters$Cluster))
  sero_no = length(unique(seq_clusters$Serotype))

  vaccTypes1 <- readRDS("UK_SeroVT.rds")
  vaccTypes2 <- readRDS("UK_SeroVT_PCV13.rds")
  
  dt <- 1/12
  peripost_mass_cluster_freq <- data.frame("year" = 1:6, rbind(mass_cluster_freq_2,mass_cluster_freq_3,mass_cluster_freq_4,mass_cluster_freq_5,mass_cluster_freq_6, mass_cluster_freq_7))
  names(peripost_mass_cluster_freq) <- c("year", as.character(1:mass_clusters))
  
  data_list <- empty_list <- vector(mode = "list", length = 6)
  data_list[[1]] <- mass_cluster_freq_2
  data_list[[2]] <- mass_cluster_freq_3
  data_list[[3]] <- mass_cluster_freq_4
  data_list[[4]] <- mass_cluster_freq_5
  data_list[[5]] <- mass_cluster_freq_6
  data_list[[6]] <- mass_cluster_freq_7
  
  vacc_time1 <- 0
  vacc_time2 <- 4
  sero_no = length(unique(seq_clusters$Serotype))
  model_start_pop <- readRDS(file = "UK_PPsero_startpop.rds")
  #model_start_pop <- readRDS(file = "UK_PPsero_startpop3.rds")
  #mass_VT2 <- readRDS(file = "UK_SeroVT_PCV13.rds")
  
  output_filename <- "FindGenes_UK"
  
  Pop_ini <- data.frame(model_start_pop)
  Pop_ini <- matrix(model_start_pop, nrow = mass_clusters, ncol = (sero_no))
  
  # use params from UK fit from 14.07.2025
  FindGenes_ggCPP_params <- list(dt = 1/36, species_no = mass_clusters,  gene_no = nrow(intermed_gene_presence_absence_consensus)-1, Pop_ini = Pop_ini, Pop_eq = rowSums(model_start_pop), capacity = sum(model_start_pop), Genotypes = intermed_gene_presence_absence_consensus_matrix, sigma_f = -2.7588411, sigma_w = 0, prop_f = 1, m = -4.6414890, migVec = avg_cluster_freq, vaccTypes1 = vaccTypes1, vaccTypes2 = vaccTypes2, v = 0.1213283, vacc_time1 = vacc_time1, vacc_time2 = vacc_time2, sero_no = sero_no)
  # read in model from file
  WF <- odin.dust::odin_dust("NFDS_Model_FindGenes_PPxSero2vacc.R")
}


# likelihood for fitting:
combined_compare_ga <- function(state, observed, pars = NULL) {
  result <- 0
  #data_size <- sum(unlist(observed))
  data_size <- sum(observed)
  model_size = sum(state)
  exp_noise <- 1e6
  data_vals <- unlist(observed)
  #model_vals <- state[-1, , drop = TRUE]
  model_vals <- rep(0, length(data_vals))
  data_missing <- FALSE
  #for (i in 1:(length(data_vals))){ 
    #state_name <- paste("sum_clust", i, sep = "")
  #  model_vals[i] <- observed[i]
  #  if (is.na(observed[i])) {
      #Creates vector of zeros in ll with same length, if no data
      #ll_obs <- numeric(length( state[state_name, , drop = TRUE]))
  #    data_missing <- TRUE
  #  } 
  #}
  model_vals <- state
  models_vals_err <- model_vals + rexp(n = length(model_vals), rate = exp_noise)
  if(data_missing){
    ll_obs <- 0
  }
  else{
    ll_obs <- dmultinom(x = (data_vals), prob = models_vals_err/model_size, log = TRUE)   
  }
  result <- ll_obs
  result
}


### try using genetic / evolutionary algorithms for finding best genes
#install.packages("GA", repos = 'https://cran.ma.imperial.ac.uk/')
library(GA)



decode2 <- function(x)
{ 
  x <- round(x)         
  return(x)
}

make_transform <- function(m) {
  function(theta) {
    as_double_mtx <- function(x){
      sapply(x,as.double)
    }
    c(lapply(m, as_double_mtx), as.list(theta))
  }
}


fitting_closure_max_decode <- function(all_other_params, data_list1){
  null_fit_dfoptim_fl <- function(fit_params){
    fit_params <- decode2(fit_params)
    rnd_vect_full <- fit_params
    
    all_other_params$delta_bool = rnd_vect_full
    #all_other_params$Pop_ini <- make_transform(Pop_ini)
    WFmodel_ggCPP <- WF$new(pars = all_other_params,
                            time = 1,
                            n_particles = 10L,
                            n_threads = 4L,
                            seed = 1L)
    #n_particles <- 10L
    #n_times <- 73
    #x <- array(NA, dim = c(WFmodel_ggCPP$info()$len, n_particles, n_times))
    
    #for (t in seq_len(n_times)) {
    #  x[ , , t] <- WFmodel_ggCPP$run(t)
    #}
    #time <- x[1, 1, ]
    #x <- x[-1, , ]
    combined_compare_res <- 0
    for (i in 1:length(data_list1)) {
      if(sum(data_list1[[i]])>0){
        simMeanggCPP2 <- rowMeans(WFmodel_ggCPP$run(i * 1/dt)[(2:(mass_clusters+1)),])
        combined_compare_res <- combined_compare_res + combined_compare_ga(simMeanggCPP2,data_list1[[i]])
      }
    }
    combined_compare_res
    
    #simMeanggCPP2 <- rowMeans(WFmodel_ggCPP$run(36)[(2:(mass_clusters+1)),])
    #simMeanggCPP3 <- rowMeans(WFmodel_ggCPP$run(72)[(2:(mass_clusters+1)),])
    #combined_compare_ga(simMeanggCPP2,data1) + combined_compare_ga(simMeanggCPP3,data2) 
    #- combined_compare(x[,1,37],data1) - combined_compare(x[,1,73],data2) 
  }
}

gene_number <-  nrow(intermed_gene_presence_absence_consensus)-1
library(doParallel)
library(foreach)
library(iterators)

monitor_fn <- function(obj) {
  cat(sprintf("Generation %d | Best Fitness: %f\n", obj@iter, max(obj@fitness)))
  flush.console()
}

ga_fit_FindGenes_ggCPP_dec <- fitting_closure_max_decode(FindGenes_ggCPP_params, data_list)
gann <- ga(type = "real-valued", fitness = ga_fit_FindGenes_ggCPP_dec, lower = rep(0, gene_number), upper = rep(1,gene_number), 
           seed = 123, elitism = 100, maxiter = 200, popSize = 1000, run = 10, parallel = 48, monitor = monitor_fn)
summary(gann)
as.vector(t(apply(gann@solution, 1, decode2)))
sum(as.vector(t(apply(gann@solution, 1, decode2))))/gene_number
saveRDS(gann,paste(output_filename, "_gann.rds", sep = ""))
