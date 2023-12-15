rm(list=ls())
source("R/ABC_algorithm.R")
library(matrixcalc)


######### Examples from Paper #############
library(GillespieSSA)

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
# Influenza B infection in 1975-76 epidemic (middle column) and influenza A (H1N1) infection in 1978-79 epidemic (right column), Seattle, Washington (Jr IL and Koopman J. Household and community transmission parameters from nal distribu-tions of infections in households. Biometrics, 115-126, 1982.)

InfB_1975_76 <- c( 9, 12, 18,  9,  4,
                   1,  6,  6,  4,  3,
                   NA,  2,  3,  4,  0, 
                   NA, NA,  1,  3,  2,
                   NA, NA, NA,  0,  0,
                   NA, NA, NA, NA,  0)
H1N1_1978_79 <- c(15, 12,  4,
                  11, 17,  4,
                  NA, 21,  4,
                  NA, NA,  5,
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
  
  data_matrix_77_78 <- CreateInfectionMatrix(qc1, qh1, observed_data[[1]])
  data_matrix2_80_81 <- CreateInfectionMatrix(qc2, qh2, observed_data[[2]])
  
  return(list(data_matrix_77_78,
              data_matrix2_80_81))
}

# Helper function to create sample matrix
CreateInfectionMatrix <- function(qc, qh, observed_data) {
  
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
kMaxIterations <- 500000

# Number of parameters per model
model_number_params <- c(4)

kEpsilon <- 20

results_Table2 <- CalculatePosteriorBase(observed_data = observed_data_Table2, model_number_params = model_number_params, kEpsilon = kEpsilon, kMaxIterations = kMaxIterations, prior_distr = prior_distr, DistanceFct = DistanceFct, DataGeneratingFct = DataGeneratingFct, kNparticles = kNparticles)

kEpsilon <- 7

results_Table3 <- CalculatePosteriorBase(observed_data = observed_data_Table3, model_number_params = model_number_params, kEpsilon = kEpsilon, kMaxIterations = kMaxIterations, prior_distr = prior_distr, DistanceFct = DistanceFct, DataGeneratingFct = DataGeneratingFct, kNparticles = kNparticles)

source("R/Plotting.R")
#source("R/Comparison.R")
