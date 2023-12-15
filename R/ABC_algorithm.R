CalculatePosteriorBase <- function(observed_data, kEpsilon, kMaxIterations, prior_distr, model_number_params, DistanceFct, DataGeneratingFct, kNparticles) {
  
  ## Function calculates the posterior distribution of parameter of a previously specified model and marginal model probabilities
  ## Implementation according to Toni T., Stumpf M.P.H. (2009) - Simulation-based model selection for dynamical systems in systems and population biology
  
  # Args:
  #   observed_data:       Observed data
  #   kEpsilon:           tolerance level that determines how close the simulated data should be to the observed data in integer
  #   kMaxIterations:      Integer defining the maximum number of iterations
  #   prior_distr:         List -> 1: Function for prior distribution of models
  #                                2: List -> 1: Function of prior distributions of parameter for first model 
  #                                           2: Function of prior distributions of parameter for second model
  #    
  #   model_number_params: Vector -> Number of parameters for each model
  #   DistanceFct:         Function that calculates the distance between simulated and observed data
  #   DataGeneratingFct:   Function that generates data based on model parameters 
  #   kNparticles:         Integer of Number of particle samples to draw
  
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
    posterior_model_distributions[[model_number]] <- data.frame(data.frame(
      replicate(number_params, vector("numeric", length = 0))
    ))
  }
  
  # Number of attempted and accepted data samples
  accepted <- 1
  attempted <- 1
  
  current_iteration <- 0
  
  ## Loop until kNparticles of parameters are accepted
  
  while(accepted <= kNparticles) {
    
    current_iteration <- current_iteration + 1

    # Break the loop if maximum iterations are reached
    if (current_iteration > kMaxIterations) {
      cat("Maximum iterations reached. Breaking out of the loop.\n")
      break
    }
    
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
      posterior_model_distributions[[model_draw]] <- rbind(posterior_model_distributions[[model_draw]], param_draw)  # Store theta* for model m*
      accepted <- accepted + 1
    } 
    attempted <- attempted + 1
  }
  

  # Compute marginal model probabilities
  marginal_model_probs <- sapply(1:number_models, FUN = function(model) {
    model_marginale_probability <- dim(posterior_model_distributions[[model]])[1] / (accepted - 1)  # Number of accepted parameters for model m' / number of total accepted particles 
    return(model_marginale_probability)
  })
  
  return(list("marginal_model_probabilities" = marginal_model_probs, 
              "accepted_parameter" = posterior_model_distributions, 
              "ratio_accepted_parameter" = accepted / attempted))
}
