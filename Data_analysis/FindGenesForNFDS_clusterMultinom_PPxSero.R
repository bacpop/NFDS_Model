# idea randomly "mutate" 0-1-vector that I can map to real gene vector (based on common absence and presence, i.e. treat genes that always occur together as one gene)

#### run this part 

### Likelihood
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

  model_vals <- state
  models_vals_err <- model_vals + rexp(n = length(model_vals), rate = exp_noise)
  result <- dmultinom(x = (data_vals), prob = models_vals_err/model_size, log = TRUE) 
  result
}
#setwd("/nfs/research/jlees/leonie/WF_fitting_2024/run3")
### load model
WF <- odin.dust::odin_dust("NFDS_Model_FindGenes_PPxSero.R")

###########################################

### try binary ga

fitting_closure_max_binary <- function(all_other_params, data1, data2){
  null_fit_dfoptim_fl <- function(fit_params){
    rnd_vect_full <- fit_params
    
    all_other_params$delta_bool = rnd_vect_full
    WFmodel_ggCPP <- WF$new(pars = all_other_params,
                            time = 0,
                            n_particles = 1L,
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
    #simMeanggCPP2 <- rowMeans(WFmodel_ggCPP$run(36)[-1,])
    #simMeanggCPP3 <- rowMeans(WFmodel_ggCPP$run(72)[-1,])
    #combined_compare(simMeanggCPP2,data1) + combined_compare(simMeanggCPP3,data2) 
    sim1 <- WFmodel_ggCPP$run(36)[(2:(mass_clusters+1))]
    sim2 <- WFmodel_ggCPP$run(72)[(2:(mass_clusters+1))]
    combined_compare_ga(sim1,data1) + combined_compare_ga(sim2,data2)
    #- combined_compare(x[,1,37],data1) - combined_compare(x[,1,73],data2) 
  }
}


library(parallel)
#install.packages("doParallel",repos = "http://cran.us.r-project.org")
library(doParallel)
#install.packages("doSNOW",repos = "http://cran.us.r-project.org")
library(doSNOW)
#install.packages("GA",repos = "http://cran.us.r-project.org")
library(GA)

# import data
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
vacc_time <- 0
output_filename <- "PPxSero_ggCaller_PopPUNK"


Pop_ini <- data.frame(model_start_pop)
Pop_ini <- matrix(model_start_pop, nrow = mass_clusters, ncol = (sero_no))

FindGenes_ggCPP_params <- list(dt = 1/36, species_no = mass_clusters,  gene_no = nrow(intermed_gene_presence_absence_consensus)-1, Pop_ini = Pop_ini, Pop_eq = rowSums(model_start_pop), capacity = sum(model_start_pop), Genotypes = intermed_gene_presence_absence_consensus_matrix, sigma_f = -3.6, sigma_w = 0, prop_f = 1, m = -4.2, migVec = avg_cluster_freq, vaccTypes = mass_VT, v = 0.087, vacc_time = 0, sero_no = sero_no)


ga_fit_FindGenes_ggCPP_bin <- fitting_closure_max_binary(FindGenes_ggCPP_params, mass_cluster_freq_2, mass_cluster_freq_3)

start_pop <- function(){
  ga_results@population
}

monitor_fn <- function(obj) {
  cat(sprintf("Generation %d | Best Fitness: %f\n", obj@iter, max(obj@fitness)))
  flush.console()
}

gann <- ga(type = "binary", nBits = 1934, fitness = ga_fit_FindGenes_ggCPP_bin, lower = rep(0, 1934), upper = rep(1,1934), 
           elitism = 10000, maxiter = 100, popSize = 100000, run = 20, pcrossover = 0.8, pmutation = 0.5, crossover = gabin_spCrossover, mutation = gabin_raMutation, parallel = 48, monitor = monitor_fn)

#plot(gann)
print(summary(gann))
print(sum(gann@solution))
#[1] 905

saveRDS(gann, "gann.rds")


pdf(file = "gann_plot.pdf",   # The directory you want to save the file in
    width = 6, # The width of the plot in inches
    height = 12)
plot(gann)
dev.off()
