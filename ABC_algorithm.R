
CalculatePosteriorBase <- function(observed_data, kEpsilon, prior_distr, model_number_params, DistanceFct, DataGeneratingFct, kNparticles) {
  
  ## Function calculates the posterior distribution of parameter of a previously specified model and marginal model probabilities
  ## Implementation according to Toni T., Stumpf M.P.H. (2009) - Simulation-based model selection for dynamical systems in systems and population biology
  
  # Args:
  #   observed_data:       Observed data
  #   kEpsilon:           tolerance level that determines how close the simulated data should be to the observed data in integer
  #   prior_distr:         List -> 1: Function for prior distribution of models
  #                                2: List -> 1: Function of prior distributions of parameter for first model 
  #                                           2: Function of prior distributions of parameter for second model
  #                                           ... 
  #   model_number_params: Vector -> Number of parameters for each model
  #   DistanceFct:         Function that calculates the distance between simulated and observed data
  #   kNparticles:         Integer of Number of particle samples to draw
  #   DataGeneratingFct:   Function that generates data based on model parameters
  
  # Returns:
  #   List -> 1: marginal_model_probs: Data frame of marginal probabilities for each model
  #           2: List posterior_model_distributions -> 1: Data frame of accepted parameter for first model
  #                                                    2: Data frame of accepted parameter for second model
  #                                                    ...
  #           3: Ratio of accepted to attempted data samples
  # 
  
  
  ## Initialization
  
  # Number of models as number of different parameter priors
  number_models <- length(prior_distr[[2]])
  
  # Number of parameter
  number_params <- sum(model_number_params)
  
  # Data frame of marginal model probabilities -> nrows: 1, ncols: number models
  marginal_model_probs <- data.frame(matrix(NA, nrow = 1, ncol = number_models))
  
  # List of Data Frame of accepted parameter for every model -> nrows: 1, ncol: number of model parameter for every model
  posterior_model_distributions <- list()
  for (model_number in 1:number_models) {
    posterior_model_distributions[[model_number]] <- data.frame(matrix(NA, nrow = 1, ncol = model_number_params[model_number]))
  }
  
  # Number of attempted and accepted data samples
  accepted <- 1
  attempted <- 1
  
  
  ## Loop until kNparticles of parameters are accepted
  
  while(accepted <= kNparticles) {
    
    # Draw m*
    model_draw <- prior_distr[[1]]()
    
    # Draw theta*
    param_draw <- prior_distr[[2]][[model_draw]]()  # Draw theta* from model parameter
    
    # Generate candidate dataset
    sample_data <- DataGeneratingFct(param_draw, observed_data)
    
    # Compute distance
    distance_data <- DistanceFct(observed_data, sample_data)
    
    
    # Store if distance <= epsilon, next otherwise
    if (distance_data <= kEpsilon) {
      posterior_model_distributions[[model_draw]][accepted,] <- param_draw  # Store theta* for model m*
      accepted <- accepted + 1
      print(accepted)
    } 
    attempted <- attempted + 1
  }
  
  # Remove all NA rows
  posterior_model_distributions <- lapply(posterior_model_distributions, FUN = function(param_df) {
    return(param_df |> na.exclude()) 
  })
  
  # Compute marginal model probabilities
  marginal_model_probs <- sapply(1:number_models, FUN = function(model) {
    model_marginale_probability <- dim(posterior_model_distributions[[model]])[1] / kNparticles  # Number of accepted parameters for model m' / number of total accepted particles 
    return(model_marginale_probability)
  })
  
  return(list("marginal_model_probabilities" = marginal_model_probs, 
              "accepted_parameter" = posterior_model_distributions, 
              "ratio_accepted_parameter" = accepted / attempted))
}
