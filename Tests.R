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

#2)
H3N2_1977_78 <- c(66, 87, 25, 22, 4,
                  13, 14, 15, 9, 4,
                  NA, 4, 4, 9, 1, 
                  NA, NA, 4, 3, 1,
                  NA, NA, NA, 1, 1,
                  NA, NA, NA, NA, 0)
H3N2_1980_81 <- c(44, 62, 47, 38, 9,
                  10, 13, 8, 11, 5,
                  NA, 9, 2, 7, 3, 
                  NA, NA, 3, 5, 1, 
                  NA, NA, NA, 1, 0, 
                  NA, NA, NA, NA, 1)

observed_data <- list(matrix(H3N2_1977_78, nrow = 6, ncol = 5, byrow = TRUE),
                      matrix(H3N2_1980_81, nrow = 6, ncol = 5, byrow = TRUE))
epsilon <- 25

prior_model <- function() {
  return(1)
}

prior_param_model1 <- function() {
  qc1 <- runif(1,0,1)
  qh1 <- runif(1,0,1)
  qc2 <- runif(1,0,1)
  qh2 <- runif(1,0,1)
  return(c(qh1, qc1, qh2, qc2))
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
  
  data_matrix <- sweep(data_probs_matrix, MARGIN = 2, total_observations, "*")
  return(data_matrix)
}

N_particles <- 100

model_number_params <- c(4)

undebug(calc_post_distr_base)
results <- calc_post_distr_base(observed_data = observed_data, model_number_params = model_number_params, epsilon = epsilon, prior_distr = prior_distr, distance_fct = distance_fct, data_generating_fct = data_generating_fct, N_particles = N_particles)
plot(results[[2]][[1]]$X1, results[[2]][[1]]$X2,
     xlim = c(0,1),
     ylim = c(0,1))
points(results[[2]][[1]]$X3, results[[2]][[1]]$X4, pch = 16)

# Test

# Distance between itself should be 0
distance_fct(observed_data, observed_data)


test1 <- c(60, 80, 10, 22, 4,
                  13, 14, 15, 9, 4,
                  NA, 4, 4, 9, 1, 
                  NA, NA, 4, 3, 1,
                  NA, NA, NA, 1, 1,
                  NA, NA, NA, NA, 0)
test2 <- c(40, 60, 45, 38, 9,
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
