##Algorithm.R

```{r}
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
  # tolerance: Integer defining epsilon
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
```

##Test.R

```{r}
library(matrixcalc)



######### Examples from Paper #############
library(GillespieSSA)
#1)
time <- seq(0,0.095,0.005)
Y <- c(3,8,15,18,23,27,32,35,37,37,
       38,39,40,41,41,42,42,42,42,42)
data <- data.frame(time = time, Y = Y)

tolerance_schedule <- c(3000,1400,600,140,40)
number_models <- 2

model_prior <- function(pdf = FALSE) {
  if (pdf == FALSE) {
    return(sample(number_models,1))
  }
  else {
    return(1/number_models)
  }
}

param_prior <- list(
  k1 = function(pdf = FALSE) {
    if (pdf == FALSE) {
      return(runif(1,0,100))
    }
    else {
      return(punif(pdf,0,100))
    }},
  k2 = function(pdf = FALSE) {
    if (pdf == FALSE) {
      return(runif(1,0,100))
    }
    else {
      return(punif(pdf,0,100))
    }}
)

prior_distr <- list(model_prior, param_prior)


distance_function <- function(candidate_data, observed_data) {
  return(mean((candidate_data - observed_data)^2))
}

KM <- function(m, m_base, list_populations, population_index) {
  if (m == m_base) {
    return(0.7)
  }
  else {
    return(0.3)
  }
}
  
PM <- function(theta, theta_base, list_populations, population_index) {
  delta <- 2 * (max(list_populations[[population_index - 1]]$param) - min(list_populations[population_index - 1]$param))
  return(runif(-delta, delta))
}


perturb_kernels <- list(KM, PM)


data_generating_fct <- function(parameter) {
  # define gillespie algorithm for dataset
  sample_data <- GillespieSSA::ssa(x0 = c(X = 40,Y = 3), a = c("c*X"), nu = matrix(c(-1,+1), nrow = 2), parms = c(c = parameter), tf = 0.2,method = ssa.etl(0.005), consoleInterval = 0.005)
  return(sample_data$data[1:20,3])
  
}

calc_post_distr(observed_data = data$Y, tolerance_schedule = tolerance_schedule,
                prior_distr = prior_distr, distance_fct = distance_function, perturb_kernels = perturb_kernels,
                data_generating_fct = data_generating_fct, N_particles = 100)
```


```{r}
##2)
H3N2_1977_78 <- c(66, 87, 25, 22,  4,
                  13, 14, 15,  9,  4,
                  NA,  4,  4,  9,  1, 
                  NA, NA,  4,  3,  1,
                  NA, NA, NA,  1,  1,
                  NA, NA, NA, NA,  0)
H3N2_1980_81 <- c(44, 62, 47, 38, 9,
                  10, 13, 8, 11, 5,
                  NA, 9, 2, 7, 3, 
                  NA, NA, 3, 5, 1, 
                  NA, NA, NA, 1, 0, 
                  NA, NA, NA, NA, 1)

InfB_1975_76 <- c( 9, 12, 18,  9,  4,
                   1,  6,  6,  4,  3,
                  NA,  2,  3,  4,  0, 
                  NA, NA,  1,  3,  2,
                  NA, NA, NA,  0,  0,
                  NA, NA, NA, NA,  0)
H1N1_1978_79 <- c(15, 12,  4,
                  11, 17,  4,
                  NA, 21,  5,
                  NA, NA, NA,
                  NA, NA, NA,
                  NA, NA, NA)


observed_data_Table2 <- list(matrix(H3N2_1977_78, nrow = 6, ncol = 5, byrow = TRUE),
                             matrix(H3N2_1980_81, nrow = 6, ncol = 5, byrow = TRUE))

observed_data_Table3 <- list(matrix(InfB_1975_76, nrow = 6, ncol = 5, byrow = TRUE),
                             matrix(H1N1_1978_79, nrow = 6, ncol = 3, byrow = TRUE))

prior_model <- function() {
  return(1)
}

prior_param_model1 <- function() {
  qc1 <- runif(1,0,1)
  qh1 <- runif(1,0,1)
  qc2 <- runif(1,0,1)
  qh2 <- runif(1,0,1)
  return(c(qc1, qh1, qc2, qh2))
}

prior_params <- list(prior_param_model1)
                     

prior_distr <- list(prior_model,
                    prior_params)

distance_fct <- function(observed_data, sample_data) {
  
  observed_data <- lapply(observed_data, FUN = function(dm) {
    df <- as.data.frame(dm) 
    df[is.na(df)] <- 0
    return(as.matrix(df))})
  
  sample_data <- lapply(sample_data, FUN = function(dm) {
    df <- as.data.frame(dm) 
    df[is.na(df)] <- 0
    return(as.matrix(df))})
  
  difference_sample1 <- observed_data[[1]]-sample_data[[1]]
  difference_sample2 <- observed_data[[2]]-sample_data[[2]]
  
  frobenius_sample1 <- sqrt(matrixcalc::matrix.trace(difference_sample1 %*% t(difference_sample1)))
  frobenius_sample2 <- sqrt(matrixcalc::matrix.trace(difference_sample2 %*% t(difference_sample2)))
  
  distance_observed_sample_data <- 1/2 * (frobenius_sample1 + frobenius_sample2)
  
  return(distance_observed_sample_data)
}

data_generating_fct <- function(params_draw, observed_data) {
  qc1 <- params_draw[1]
  qh1 <- params_draw[2]
  qc2 <- params_draw[3]
  qh2 <- params_draw[4]
  
  data_matrix_77_78 <- create_matrix(qc1, qh1, observed_data[[1]])
  data_matrix2_80_81 <- create_matrix(qc2, qh2, observed_data[[2]])
  
  return(list(data_matrix_77_78,
              data_matrix2_80_81))
}
  
create_matrix <- function(qc, qh, observed_data) {
  
  matrix_rows <- dim(observed_data)[1]
  matrix_columns <- dim(observed_data)[2]
  
  total_observations <- apply(observed_data, MARGIN = 2, FUN = sum, na.rm = TRUE)
  
  data_matrix <- matrix(NA, nrow = matrix_rows, ncol = matrix_columns)
  
  data_probs_matrix <- data_matrix
  data_probs_matrix[1,] <- qc^(1:matrix_columns)
  
  
  for (column in 1:matrix_columns) {
    for (row in 2:matrix_rows) {
      if (is.na(observed_data[row,column])) {
        data_probs_matrix[row, column] <- NA 
        next 
      }
      
      if ((row - 1) == column) data_probs_matrix[row, column] <- 1 - sum(data_probs_matrix[,column], na.rm = TRUE)
      else {
        data_probs_matrix[row, column] <- choose(column, row-1) * data_probs_matrix[row,row-1] * (qc * qh^(row-1))^(column-(row-1))
      }
      
    }
    
  }
  
  data_matrix <- round(sweep(data_probs_matrix, MARGIN = 2, total_observations, "*"))
  return(data_matrix)
}

N_particles <- 100

model_number_params <- c(4)

epsilon <- 20

results_Table2 <- calc_post_distr_base(observed_data = observed_data_Table2, model_number_params = model_number_params, epsilon = epsilon, prior_distr = prior_distr, distance_fct = distance_fct, data_generating_fct = data_generating_fct, N_particles = N_particles)

epsilon <- 8

results_Table3 <- calc_post_distr_base(observed_data = observed_data_Table3, model_number_params = model_number_params, epsilon = epsilon, prior_distr = prior_distr, distance_fct = distance_fct, data_generating_fct = data_generating_fct, N_particles = N_particles)


# Test

# Examine accepted parameters qith low qh1 values:
# Create grid and check which is smalles distance
grid <- expand.grid(qh1 = seq(0,1,0.1), qc1 = seq(0.75,1,0.05), gh2 = seq(0,1,0.1), qc2 = seq(0.75,1,0.05))

grid_results <- apply(grid, MARGIN = 1, FUN = function(row) {
  return(distance_fct(observed_data_Table2, data_generating_fct(row, observed_data_Table2)))
})

test_results <- cbind(grid, grid_results)
test_results[test_results$grid_results < 20,]

# Distance between itself should be 0
distance_fct(observed_data, observed_data)


test1 <- c(63, 87, 25, 22, 4,
                  13, 14, 15, 9, 4,
                  NA, 4, 4, 9, 1, 
                  NA, NA, 4, 3, 1,
                  NA, NA, NA, 1, 1,
                  NA, NA, NA, NA, 0)
test2 <- c(42, 62, 47, 38, 9,
                  10, 13, 8, 11, 5,
                  NA, 9, 2, 7, 3, 
                  NA, NA, 3, 5, 1, 
                  NA, NA, NA, 1, 0, 
                  NA, NA, NA, NA, 1)

test_matrix <- list(matrix(test1, nrow = 6, ncol = 5, byrow = TRUE),
                      matrix(test2, nrow = 6, ncol = 5, byrow = TRUE))

distance_fct(observed_data, test_matrix)


## Probs should add up to 1
replicate(10, round(create_matrix(runif(1,0,1), runif(1,0,1), observed_data[[1]])))

test_data <- lapply(observed_data, FUN = function(dm) {
  df <- as.data.frame(dm) 
  df[is.na(df)] <- 0
  return(as.matrix(df))})
observed_data_test <- observed_data[[1]][is.na(observed_data[[1]])] <- 0
```

##Plots.R

```{r echo=FALSE}
library(ggplot2)

colnames(results_Table2[[2]][[1]]) <- c("qc1", "qh1", "qc2", "qh2")

ggplot(data = results_Table2[[2]][[1]]) +
  geom_point(aes(qh1, qc1, colour = "red")) +
  geom_point(aes(qh2, qc2, colour = "blue")) +
  labs(x = expression(q[c]), y = expression(q[h])) + 
  xlim(c(0, 1)) + ylim(c(0,1)) +
  ggtitle(label = "Figure 3(a)" , 
          subtitle = "ABC posterior distributions for parameters modeling 
  Supplementary Table 2")

colnames(results_Table3[[2]][[1]]) <- c("qc1", "qh1", "qc2", "qh2")

ggplot(data = results_Table3[[2]][[1]]) +
  geom_point(aes(qh1, qc1, colour = "red")) +
  geom_point(aes(qh2, qc2, colour = "blue")) +
  labs(x = expression(q[c]), y = expression(q[h])) + 
  xlim(c(0, 1)) + ylim(c(0,1)) +
  ggtitle(label = "Figure 3(b)" , 
          subtitle = "ABC posterior distributions for parameters modeling
  Supplementary Table 3")
```