library(dplyr)

calc_post_distr <- function(observed_data, tolerance_schedule, prior_distr, distance_fct, perturb_kernels, data_generating_fct, N_particles) {
  ## Function calculates the posterior distribution of model and parameter
  #
  # observed_data: Observed data
  # tolerance_schedule: vector or tolerance_schedule
  # prior_distr: list -> 1: function of prior distribution for models #!# -> could be relaxed to selection of prior
  #                              2: functions of prior distributions for parameter of each model (length is model number)
  # distance_fct: function distance function 
  # perturb_kernels: list of functions as perturbation kernels
  #                 -> list(KM <- model perturbation kernel, KP <- parameter pertubation kernel)
  # N_particles: Integer of N particle samples
  # data_generating_fct: Data generating function
  ## 
  
  
  # Test that all inputs are properly defined 
  
  # 1. tolerance_schedule:
  stopifnot(is.unsorted(tolerance_schedule) | any(tolerance_schedule) <= 0)
  # 2. That model number is same as defined in prior_distr[[2]]:
  #
  # ...
  
  number_models <- length(prior_distr[[2]])
  number_tolerances <- length(tolerance_schedule)
  
  # Create dataframe of marginal model probabilities (length tolerance_schedule * number models)
  marginal_model_probs <- data.frame(matrix(NA, nrow = number_tolerances, ncol = number_models))
  # Create list of population_dataframe for each epsilon (length tolerance_schedule) -> particle (model, parameter), weight
  list_populations <- lapply(tolerance_schedule, FUN = function(tolerance) list())

  # Loop over tolerance_schedule (which define the populations)
  for (t in seq_along(tolerance_schedule)) {
    
    # Create population_dataframe for population t (N_particles * (model, parameter, weight))

    population_dataframe <- data.frame(matrix(NA, nrow = N_particles, ncol = 3))
    colnames(population_dataframe) <- c("model","parameter","weight")

    # Loop over number of particles per population
    for (i in seq_along(N_particles)) {
      
      population_dataframe[i] <- sim_i_data(list_populations, 
                                            observed_data, 
                                            tolerance_schedule, 
                                            prior_distr, 
                                            distance_fct, 
                                            perturb_kernels, 
                                            data_generating_fct, 
                                            marginal_model_probs, 
                                            epsilon = tolerance_schedule[t], population_index = t, number_models) 
      
    }
    list_populations[[t]] <- population_dataframe #store simulated population in list_populations
  }
  list_populations[[t]]$weight <- normalize_weights(list_populations[[t]]$weight) #Calculate normalized weights per model
  marginal_model_probs[t,c(-1)] <- calc_marginal_model_probs(list_populations[[t]]) #Calculate marginal model probabilities out of model parameter and weights 
  #!# adjust!: models doesn't have to survive - currently matching
}


sim_i_data <- function(list_populations, 
                       observed_data, 
                       tolerance_schedule, 
                       prior_distr, 
                       distance_fct, 
                       perturb_kernels, 
                       data_generating_fct, 
                       marginal_model_probs, epsilon, population_index, number_models) {
  ## Function simulate a candidate dataset and evaluates it with distance function
  ## If d(D0, D*) > epsilon: reject and draw again
  ## If d(D0, D*) <= epsilon: accept candidate and return
  ## Distinguish between t == 1 and t > 1: Use previous populations
  while(TRUE) {
    if (population_index == 1) { #!# Check here if this can be put outside of while loop (not as in paper) 
      sample_model <- prior_distr[[1]] #Get sample model m**
      sample_param <- prior_distr[[2]][[sample_model]] #Get sample param theta**
      sample_model_param <- c(sample_model, sample_param)
    }
    else {
      # m_star <- random draw from population t-1 with weights of previous population
      sample_model <- sample(number_models, 1, prob = marginal_model_probs[population_index - 1,])
      # model <- KM(m_star) -> Apply KM on any possible model and sample from permutated probabilities
      perturbated_model <- sample(number_models, 1, prob = sapply(1:number_models, FUN = perturb_kernels[[1]], sample_model)) 
      # param_star <- sample parameter with weights from previous population restricted on model drawn (model)
      param_probs_last_pop <- list_populations[[population_index - 1]][[2]] |>
        dplyr::filter(model == perturbated_model) |>
        dplyr::select(param, weight)
      
      sample_param <- sample(param_probs_last_pop$param, 1, param_probs_last_pop$weight)
      # parameter <- PM(param_star) 
      perturbated_param <- sample(param_probs_last_pop$param, 1, prob = sapply(param_probs_last_pop$param, perturb_kernels[[2]], sample_param)) #!# Check, not clear from paper -> still confused about this!!!
      sample_model_param <- c(perturbated_model, perturbated_param) 
    }
    
    prior_prob_param <- prior_distr[[1]](sample_model_param[1])*prior_dist[[2]][[sample_model_param[1]]](sample_model_param[2])
    if (prior_prob_param == 0) next 
    
    candidate_data <- data_generating_fct(sample_model_param)
    
    if (distance_fct(candidate_data, observed_data) > epsilon) next
    
    weight_particle <- calc_particle_weight(list_populations,
                                            candidate_data, 
                                            sample_model_param,
                                            population_index, 
                                            marginal_model_probs, 
                                            perturb_kernels,
                                            prior_prob_param,
                                            param_probs_last_pop) #calculate weight of particle if accepted
    return(c(sample_model_param, weight_particle))
  }
}
  
calc_particle_weight <- function(list_populations, candidate_data, sample_model_param, population_index, marginal_model_probs, perturb_kernels, prior_prob_param, param_probs_last_pop) {
  if (population_index == 1) {
    return(1) #!# Must include calculation of bt (didn't understand it -> number of replicae simulation runs)
  }
  else {
    bt <- 1 #!# Must include calculation of bt (didn't understand it -> number of replicae simulation runs)
    numerator <- prior_prob_param * bt
    
    # Sum over all possible models (j):
    # marginal_model_probs[population_index - 1, j] * perturb_kernels[[1]](sample_model_param(model, model j)) * SUM over parameter associated with sampled model in last population ((weight of parameter in last pop * KP) / prob of model in last population) #!# didn't understand properly
    
    rem_marginal_model_probs <- marginal_model_probs[population_index - 1, marginal_model_probs[population_index,] != 0][,-1] #Select all remaining models from previous population
    numbers_rem_model <- colnames(pos_marginal_model_probs)
    KM <- sapply(1:numbers_rem_model, FUN = perturb_kernels[[1]], sample_model_param[1]) #perturb model weight
    
    S_1 <- sum(rem_marginal_model_probs * KM)
    
    weights_rem_param <- param_probs_last_pop
    KP <- sapply(param_probs_last_pop$param, perturb_kernels[[2]], sample_model_param[2])
    
    S_2 <- sum(weights_rem_param * KP)/marginal_model_probs[population_index - 1, sample_model_param[1]]
    
    denominator <- S_1 * S_2
    
    weight_particle <- numerator / denominator # (prior_prob_sample * bt) / S
    return(weight_particle)
  }
  
  
}

normalize_weights <- function(weights) {
  # normalize as weights / sum(weights)
  norm_weights <- weights / sum(weights)
  return(norm_weights)

}

calc_marginal_model_probs <- function(population_results) {
  # sum weights for each model / sum of weights
  marginal_model_prob <- population |>
    dplyr::group_by(model) |>
    dplyr::summarize(model_marginal = sum(weight))
  
  # return marginal model probabilities
  return(marginal_model_prob)
}

