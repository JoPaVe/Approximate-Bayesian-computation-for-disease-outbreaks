
calc_post_distr <- function(observed_data, tolerance_schedule, prior_distr, distance_fct, perturb_kernels, data_generating_fct, N_particles) {
  ## Function calculates the posterior distribution of model and parameter
  ##
  # observed_data: Observed data
  # tolerance_schedule: vector or tolerance_schedule
  # prior_distr: function of prior distribution #!# -> could be relaxed to selection of prior
  # distance_fct: function distance function 
  # perturb_kernels: list of functions as perturbation kernels
  #                 -> list(KM <- model perturbation kernel, KP <- parameter pertubation kernel)
  # N_particles: Integer of N particle samples
  # data_generating_fct: Data generating function
  ## 
  
  
  # Test that all inputs are properly defined 
  # e.g. all epsilon in tolerance_schedule are positive and decreasing
  # ...
  
  # Create dataframe of marginal model probabilities (length tolerance_schedule * number models)
  marginal_model_probs
  # Create list of population_dataframe for each epsilon (length tolerance_schedule) -> particle (model, parameter), weight
  list_populations  

  # Loop over tolerance_schedule (which define the populations)
  for (t in seq_along(tolerance_schedule)) {
    
    # Create population_dataframe for population t (N_particles * number params)

    population_dataframe

    
    for (i in seq_along(N_particles)) {
      
      population_dataframe[i] <- sim_i_data(list_populations, 
                                            observed_data, 
                                            tolerance_schedule, 
                                            prior_distr, 
                                            distance_fct, 
                                            perturb_kernels, 
                                            data_generating_fct, 
                                            marginal_model_probs, epsilon = tolerance_schedule[t], population_index = t) 
      
    }
    list_popluations[[t]] <- population_dataframe
  }
  Normalized_weights <- normalize_weights(list_populations[[t]]$weights)
  list_populations[[t]] # add normalized_weights
  marginal_model_probs[t] <- calc_marginal_model_probs(list_populations[[t]]) 

}


sim_i_data <- function(list_populations, 
                       observed_data, 
                       tolerance_schedule, 
                       prior_distr, 
                       distance_fct, 
                       perturb_kernels, 
                       data_generating_fct, 
                       marginal_model_probs, epsilon, population_index) {
  ## Function simulate a candidate dataset and evaluates it with distance function
  ## If d(D0, D*) > epsilon: reject and draw again
  ## If d(D0, D*) <= epsilon: accept candidate and return
  ## Distinghish between t == 1
  while(TRUE) {
    if (population_index == 1) { #!# Check here if this can be put outside of while loop (not as in paper) 
      sample_model_param <- prior_distr()
    }
    else {
      # m_star <- random draw from population t-1 with weights of previous population
      # model <- KM(m_star)
      # param_star <- sample parameter with weights from previous population restricted on model drawn (model)
      # parameter <- PM(param_star) 
      sample_model_param # complete c(parameter, model) 
    }
    
    prior_prob_sample <- prior_distr(sample_model_param)
    if (prior_prob_sample == 0) next #!# Don't understand why we need it - is in paper p. 105
    
    candidate_data <- data_generating_fct(sample_model_param)
    
    if (distance_fct(candidate_data, observed_data) > epsilon) next
    
    weight_particle <- calc_particle_weight(candidate_data, 
                                   sample_model_param,
                                   population_index, 
                                   marginal_model_probs, 
                                   perturb_kernels,
                                   prior_prob_sample) #calculate weight of particle if accepted
    return(sample_model_param, weight_particle)
  }
}
  
calc_particle_weight <- function(candidate_data, sample_model_param, population_index, marginal_model_probs, perturb_kernels, prior_prob_sample) {
  if (population_index == 1) {
    return(1) #!# Must include calculation of bt (didn't understand it -> number of replicae simulation runs)
  }
  else {
    bt <- 1 #!# Must include calculation of bt (didn't understand it -> number of replicae simulation runs)
    S <- # Sum over all possible models (j):
          # marginal_model_probs[population_index - 1, j] * perturb_kernels[[1]](sample_model_param(model, model j)) * SUM over parameter associated with sampled model in last population ((weight of parameter in last pop * KP) / prob of model in last population) #!# didn't understand properly
    weight_particle <- x# (prior_prob_sample * bt) / S
    return(weight_particle)
  }
  
  
}

normalize_weights <- function(weights) {
  # normalize as weights / sum(weights)
}

calc_marginal_model_probs <- function(population) {
  # sum weights for each model / sum of weights
  # return marginal model probabilities
}
