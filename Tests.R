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


#### Example 3.3 Influenza infection outbreaks of Toni T., Stumpf M.P.H. (2009) - Simulation-based model selection for dynamical systems in systems and population biology

# Observed data according to Supplementary Table 2
# Influenza A (H3N2) infection in 1977-78 (middle column) and 1980-81 (right column) epidemics, Tecumseh, Michigan (Addy C, Jr IL and Haber M. A generalized stochastic model for the analysis of infectious disease nal size data. Biometrics, 961-974, 1991.)
H3N2_1977_78 <- c(66, 87, 25, 22,  4,
                  13, 14, 15,  9,  4,
                  NA,  4,  4,  9,  1, 
                  NA, NA,  4,  3,  1,
                  NA, NA, NA,  1,  1,
                  NA, NA, NA, NA,  0)
H3N2_1980_81 <- c(44, 62, 47, 38,  9,
                  10, 13,  8, 11,  5,
                  NA,  9,  2,  7,  3, 
                  NA, NA,  3,  5,  1, 
                  NA, NA, NA,  1,  0, 
                  NA, NA, NA, NA,  1)

# Observed data according to Supplementary Table 3
# Inuenza B infection in 1975-76 epidemic (middle column) and influenza A (H1N1) infection in 1978-79 epidemic (right column), Seattle, Washington (Jr IL and Koopman J. Household and community transmission parameters from nal distribu-tions of infections in households. Biometrics, 115-126, 1982.)
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

# Observed data stored in List -> 1: First data matrix 
#                                 2: Second data matrix
observed_data_Table2 <- list(matrix(H3N2_1977_78, nrow = 6, ncol = 5, byrow = TRUE),
                             matrix(H3N2_1980_81, nrow = 6, ncol = 5, byrow = TRUE))

observed_data_Table3 <- list(matrix(InfB_1975_76, nrow = 6, ncol = 5, byrow = TRUE),
                             matrix(H1N1_1978_79, nrow = 6, ncol = 3, byrow = TRUE))

# Prior for model
PriorModel <- function() {
  return(1)
}

# Prior for parameter of model 1 (only model in example)
PriorParameterModel1 <- function() {
  qc1 <- runif(1,0,1)
  qh1 <- runif(1,0,1)
  qc2 <- runif(1,0,1)
  qh2 <- runif(1,0,1)
  return(c(qc1, qh1, qc2, qh2))
}

# Parameter prior stored in List -> 1: Prior of parameter for model 1
prior_params <- list(PriorParameterModel1)
                     
# Prior of model and parameter in List -> 1: Prior of model 1
#                                         2: Prior of parameter for model 1
prior_distr <- list(PriorModel,
                    prior_params)


# Distance function as described in Example 3.3 
DistanceFct <- function(observed_data, sample_data) {
  
  # Replace NAs in data matrices with 0 (to compute Frobenius norm)
  observed_data <- lapply(observed_data, FUN = function(dm) {
    df <- as.data.frame(dm) 
    df[is.na(df)] <- 0
    return(as.matrix(df))})
  
  sample_data <- lapply(sample_data, FUN = function(dm) {
    df <- as.data.frame(dm) 
    df[is.na(df)] <- 0
    return(as.matrix(df))})
  
  # Compute distances between observed and sampled data matrices
  difference_sample1 <- observed_data[[1]]-sample_data[[1]] 
  difference_sample2 <- observed_data[[2]]-sample_data[[2]]
  
  # Compute frobenius norms for distances 
  # Attention: Frobenius norm (F(A)) calculated as F(A) = sqrt(trace(t(A)%*%A))
  frobenius_sample1 <- sqrt(matrixcalc::matrix.trace(difference_sample1 %*% t(difference_sample1)))
  frobenius_sample2 <- sqrt(matrixcalc::matrix.trace(difference_sample2 %*% t(difference_sample2)))
  
  # Return distance as 1/2 * sum of frobenius norms
  distance_observed_sample_data <- 1/2 * (frobenius_sample1 + frobenius_sample2)
  
  return(distance_observed_sample_data)
}


# Data generating process as described in Example 3.3 
DataGeneratingFct <- function(params_draw, observed_data) {
  qc1 <- params_draw[1]
  qh1 <- params_draw[2]
  qc2 <- params_draw[3]
  qh2 <- params_draw[4]
  
  data_matrix_77_78 <- CreateMatrix(qc1, qh1, observed_data[[1]])
  data_matrix2_80_81 <- CreateMatrix(qc2, qh2, observed_data[[2]])
  
  return(list(data_matrix_77_78,
              data_matrix2_80_81))
}

# Helper function to create sample matrix
CreateMatrix <- function(qc, qh, observed_data) {
  
  # Define matrix rows and columns as in observed data
  matrix_rows <- dim(observed_data)[1] 
  matrix_columns <- dim(observed_data)[2]
  
  # Calculate total observations per column
  total_observations <- apply(observed_data, MARGIN = 2, FUN = sum, na.rm = TRUE)
  
  # Define data matrix with same dimensions as observed data 
  data_matrix <- matrix(NA, nrow = matrix_rows, ncol = matrix_columns)
  
  # Define matrix of probabilities
  data_probs_matrix <- data_matrix
  data_probs_matrix[1,] <- qc^(1:matrix_columns)  # First row of probability matrix (w_0s) is defined as qc^s 
  
  # Fill matrix of probabilities row wise
  for (column in 1:matrix_columns) {  # Start at first column
    for (row in 2:matrix_rows) {  # Start at second row
      if (is.na(observed_data[row,column])) {
        data_probs_matrix[row, column] <- NA  # Skip weight calculation if entry is NA in observed data matrix (i.e. more infected persons as susceptible)
        next 
      }
      
      if ((row - 1) == column) data_probs_matrix[row, column] <- 1 - sum(data_probs_matrix[,column], na.rm = TRUE)  # Calculate entries w_jj as 1 - sum_i^j-1(entries ij)
      else {
        data_probs_matrix[row, column] <- choose(column, row-1) * data_probs_matrix[row,row-1] * (qc * qh^(row-1))^(column-(row-1))  # Calculate entries w_js according to formula given in example
      }
      
    }
    
  }
  
  # Calculate sampled data matrix as Total observations in column * probabilities in column
  data_matrix <- round(sweep(data_probs_matrix, MARGIN = 2, total_observations, "*"))  # Round solution to reduce distances due to decimal points
  return(data_matrix)
}

# Number of accepted particles required
kNparticles <- 100

# Number of maximum iterations
kMaxIterations <- 10000

# Number of parameters per model
model_number_params <- c(4)

kEpsilon <- 20

results_Table2 <- CalculatePosteriorBase(observed_data = observed_data_Table2, model_number_params = model_number_params, kEpsilon = kEpsilon, kMaxIterations = kMaxIterations, prior_distr = prior_distr, DistanceFct = DistanceFct, DataGeneratingFct = DataGeneratingFct, kNparticles = kNparticles)

kEpsilon <- 8

results_Table3 <- CalculatePosteriorBase(observed_data = observed_data_Table3, model_number_params = model_number_params, kEpsilon = kEpsilon, kMaxIterations = kMaxIterations, prior_distr = prior_distr, DistanceFct = DistanceFct, DataGeneratingFct = DataGeneratingFct, kNparticles = kNparticles)


## Tests for arguments of calc_post_distr_base(...)
# Strategy: Test each Function in arguments and use debugging to test if calc_post_distr_base works properly


# Examine accepted parameters qh low qh1 values:
# Create grid and check which is smallest distance
grid <- expand.grid(qh1 = seq(0,1,0.1), qc1 = seq(0.75,1,0.05), gh2 = seq(0,1,0.1), qc2 = seq(0.75,1,0.05))

grid_results <- apply(grid, MARGIN = 1, FUN = function(row) {
  return(DistanceFct(observed_data_Table2, data_generating_fct(row, observed_data_Table2)))
})

test_results <- cbind(grid, grid_results)
test_results[test_results$grid_results < 20,]

## DistanceFct
# Distance between itself should be 0
DistanceFct(observed_data_Table2, observed_data_Table2)

# Manually create small differences and assess differences
test1 <- c(64, 87, 25, 22,  4,
           13, 14, 15,  9,  4,
           NA,  4,  4,  9,  1, 
           NA, NA,  4,  3,  1,
           NA, NA, NA,  1,  1,
           NA, NA, NA, NA,  0)

test2 <- c(44, 62, 47, 38,  9,
           10, 13,  8, 11,  5,
           NA,  9,  2,  7,  3, 
           NA, NA,  3,  5,  1, 
           NA, NA, NA,  1,  0, 
           NA, NA, NA, NA,  1)

test_matrix <- list(matrix(test1, nrow = 6, ncol = 5, byrow = TRUE),
                    matrix(test2, nrow = 6, ncol = 5, byrow = TRUE))


## Data generating process function

nrows <- dim(observed_data_Table2[[1]])[1]
ncols <- dim(observed_data_Table2[[1]])[2]

# Sums per column in data matrix should add up number of observed data matrix when matrix is created

for (i in 1:10) {
  random_matrix <- round(matrix(replicate(nrows * ncols,
                                    runif(1,0,5)), 
                          nrow = nrows, 
                          ncol = ncols) * observed_data_Table2[[1]])
  
  generated_matrix <- round(CreateMatrix(runif(1,0,1), 
                     runif(1,0,1), 
                     random_matrix
  ))
  
  print(all(colSums(random_matrix, na.rm = TRUE), colSums(generated_matrix, na.rm = TRUE)))
  
}          


CreateMatrixProbTest <- function(qc, qh, observed_data) {
  
  # Define matrix rows and columns as in observed data
  matrix_rows <- dim(observed_data)[1] 
  matrix_columns <- dim(observed_data)[2]
  
  # Calculate total observations per column
  total_observations <- apply(observed_data, MARGIN = 2, FUN = sum, na.rm = TRUE)
  
  # Define data matrix with same dimensions as observed data 
  data_matrix <- matrix(NA, nrow = matrix_rows, ncol = matrix_columns)
  
  # Define matrix of probabilities
  data_probs_matrix <- data_matrix
  data_probs_matrix[1,] <- qc^(1:matrix_columns)  # First row of probability matrix (w_0s) is defined as qc^s 
  
  # Fill matrix of probabilities row wise
  for (column in 1:matrix_columns) {  # Start at first column
    for (row in 2:matrix_rows) {  # Start at second row
      if (is.na(observed_data[row,column])) {
        data_probs_matrix[row, column] <- NA  # Skip weight calculation if entry is NA in observed data matrix (i.e. more infected persons as susceptibles)
        next 
      }
      
      if ((row - 1) == column) data_probs_matrix[row, column] <- 1 - sum(data_probs_matrix[,column], na.rm = TRUE)  # Calculate entries w_jj as 1 - sum_i^j-1(entries ij)
      else {
        data_probs_matrix[row, column] <- choose(column, row-1) * data_probs_matrix[row,row-1] * (qc * qh^(row-1))^(column-(row-1))  # Calculate entries w_js according to formula given in example
      }
      
    }
    
  }
  
  # Calculate sampled data matrix as Total observations in column * probabilities in column
  data_matrix <- round(sweep(data_probs_matrix, MARGIN = 2, total_observations, "*"))  # Round solution to reduce distances due to decimal points
  return(list(data_matrix, data_probs_matrix))
}

# Sums per column in probability matrix should add up to 1 when matrix is created
for (i in 1:10) {
  random_matrix <- round(matrix(replicate(nrows * ncols,
                                          runif(1,0,5)), 
                                nrow = nrows, 
                                ncol = ncols) * observed_data_Table2[[1]])
  
  generated_matrix <- CreateMatrixProbTest(runif(1,0,1), 
                                         runif(1,0,1), 
                                         random_matrix)
                            
  
  print(colSums(generated_matrix[[2]], na.rm = TRUE))
  
}


## Test with two or more models (Fig. 3 b)
# Prior for model
PriorModel <- function() {
  return(sample(1:2,1))
}

# Prior for parameter of model 2 (4 parameter)
PriorParameterModel2 <- function() {
  qc1 <- runif(1,0,1)
  qh1 <- runif(1,0,1)
  qc2 <- runif(1,0,1)
  qh2 <- runif(1,0,1)
  return(c(qc1, qh1, qc2, qh2))
}

# Prior for parameter of model 1 (2 parameter)
PriorParameterModel1 <- function() {
  qc1 <- runif(1,0,1)
  qh1 <- runif(1,0,1)
  qc2 <- qc1
  qh2 <- qh1
  return(c(qc1, qh1, qc2, qh2))
}

# Parameter prior stored in List -> 1: Prior of parameter for model 1
#                                   2: Prior of parameter for model 2
prior_params <- list(PriorParameterModel1, PriorParameterModel2)

# Prior of model and parameter in List -> 1: Prior of model 1
#                                         2: Prior of parameter for model 1 and model 2
prior_distr <- list(PriorModel,
                    prior_params)
  
# Number of parameters per model
model_number_params <- c(4, 4)

kEpsilon <- 15

kNparticles <- 100

results_Table2_test <- CalculatePosteriorBase(observed_data = observed_data_Table2, model_number_params = model_number_params, kEpsilon = kEpsilon, kMaxIterations = kMaxIterations, prior_distr = prior_distr, DistanceFct = DistanceFct, DataGeneratingFct = DataGeneratingFct, kNparticles = kNparticles)
results_Table2_test$marginal_model_probabilities


## Test with two or more models (Fig. 3 d)
# Prior for parameter of model 1 (3 parameter)
PriorParameterModel1 <- function() {
  qc1 <- runif(1,0,1)
  qh1 <- runif(1,0,1)
  qc2 <- runif(1,0,1)
  qh2 <- qh1
  return(c(qc1, qh1, qc2, qh2))
}

# Prior for parameter of model 1 (2 parameter)
PriorParameterModel2 <- function() {
  qc1 <- runif(1,0,1)
  qh1 <- runif(1,0,1)
  qc2 <- qc1
  qh2 <- qh1
  return(c(qc1, qh1, qc2, qh2))
}

# Parameter prior stored in List -> 1: Prior of parameter for model 1
#                                   2: Prior of parameter for model 2
prior_params <- list(PriorParameterModel1, PriorParameterModel2)

# Prior of model and parameter in List -> 1: Prior of model 1
#                                         2: Prior of parameter for model 1 and model 2
prior_distr <- list(PriorModel,
                    prior_params)


results_Table3_test <- CalculatePosteriorBase(observed_data = observed_data_Table3, model_number_params = model_number_params, kEpsilon = kEpsilon, kMaxIterations = kMaxIterations, prior_distr = prior_distr, DistanceFct = DistanceFct, DataGeneratingFct = DataGeneratingFct, kNparticles = kNparticles)
results_Table3_test$marginal_model_probabilities
