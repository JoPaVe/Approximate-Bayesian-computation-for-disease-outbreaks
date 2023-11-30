
calc_post_distr_base <- function(observed_data, epsilon, prior_distr, model_number_params, distance_fct, data_generating_fct, N_particles) {
  ## Function calculates the posterior distribution of model and parameter
  #
  # observed_data: Observed data
  # prior_distr: list -> 1: function of prior distribution for models #!# -> could be relaxed to selection of prior
  #                      2: list of functions of prior distributions for parameter of each model (length is model number) -> first entry is number of parameters
  # model_number_params: vector or model parameter: c(modelparam1, modelparam2, ...)
  # distance_fct: function distance function (candidate_data, observed_data)
  # N_particles: Integer of N particle samples
  # data_generating_fct: Data generating function
  # tolerance: Integer defining epsilon.
  ## 
  
  
  # Test that all inputs are properly defined 
  
  
  
  number_models <- length(prior_distr[[2]])
  number_params <- sum(model_number_params)
  
  
  
  # Create dataframe of marginal model probabilities (length tolerance_schedule * number models)
  marginal_model_probs <- data.frame(matrix(NA, nrow = 1, ncol = number_models))
  # Create dataframe of accepted particles
  # list of data.frame with accepted particles for every model
  posterior_model_distributions <- list()
  for (model in 1:number_models) {
    posterior_model_distributions[[model]] <- data.frame(matrix(NA, nrow = 1, ncol = model_number_params[model]))
  }
  
  accepted <- 1
  attempted <- 1
  # Loop until accepted == N_particles 
  while(accepted <= N_particles) {
    
    # Draw m*
    model_draw <- prior_distr[[1]]()
    
    # Draw theta*
    param_draw <- prior_distr[[2]][[model_draw]]()
    
    # Create candidate dataset
    sample_data <- data_generating_fct(param_draw, observed_data)
    
    # compute distance
    distance_data <- distance_fct(observed_data, sample_data)
    
    
    # Store or next if distance > epsilon
    if (distance_data <= epsilon) {
      posterior_model_distributions[[model_draw]][accepted,] <- param_draw
      accepted <- accepted + 1
    } 
    attempted <- attempted + 1
  }
  
  # Compute marginal model probabilities
  marginal_model_probs <- sapply(1:number_models, FUN = function(model) {
    model_marginale_probability <- dim(posterior_model_distributions[[model]])[1] / N_particles
    return(model_marginale_probability)
  })
  
  return(list(marginal_model_probs, posterior_model_distributions, accepted / attempted))
}
