







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
observed_data <- matrix()
tolerance <- 
prior_distr <- 
distance_fct <- 
data_generating_fct <- 
N_particles <- 200
