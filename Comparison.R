rm(list = ls())
#setwd("C:/Users/ssara/Documents/IU Courses/Stats Degree/Courses/FA23/STAT-S 610 Intro to Statistical Computing/Project/Project/Approximate-Bayesian-computation-for-disease-outbreaks")
source("ABC_algorithm.R")
source("Tests.R")

# epsilon values
epsilon_values <- c(25, 20, 15, 10, 7, 4)

# Initialize variables to store execution times and results
execution_times_3a <- numeric(length(epsilon_values))
execution_times_3b <- numeric(length(epsilon_values))
results_list_3a <- list()
results_list_3b <- list()

# Loop over epsilon values for 3a
for (i in seq_along(epsilon_values)) {
  start_time <- Sys.time()
  results <- CalculatePosteriorBase(observed_data = observed_data_Table2, model_number_params = model_number_params, kEpsilon = epsilon_values[i], prior_distr = prior_distr, DistanceFct = DistanceFct, DataGeneratingFct = DataGeneratingFct, kNparticles = kNparticles)
  end_time <- Sys.time()
  execution_times_3a[i] <- difftime(end_time, start_time, units = "secs")
  results_list_3a[[i]] <- results  # Store results for plotting
  
  cat("Completed epsilon value:", epsilon_values[i], "in", execution_times_3a[i], "seconds\n")
}

# Loop over epsilon values for 3b
for (i in seq_along(epsilon_values)) {
  start_time <- Sys.time()
  results <- CalculatePosteriorBase(observed_data = observed_data_Table3, model_number_params = model_number_params, kEpsilon = epsilon_values[i], prior_distr = prior_distr, DistanceFct = DistanceFct, DataGeneratingFct = DataGeneratingFct, kNparticles = kNparticles)
  end_time <- Sys.time()
  execution_times_3b[i] <- difftime(end_time, start_time, units = "secs")
  results_list_3b[[i]] <- results  # Store results for plotting
  
  cat("Completed epsilon value:", epsilon_values[i], "in", execution_times_3b[i], "seconds\n")
}

epsilon_table_3a <- data.frame(Epsilon = epsilon_values, ExecutionTime = execution_times_3a)
epsilon_table_3b <- data.frame(Epsilon = epsilon_values, ExecutionTime = execution_times_3b)

plots_list_3a <- list()
plots_list_3b <- list()


for (i in seq_along(results_list_3a)) {
  
  colnames(results_list_3a[[i]][[2]][[1]]) <- c("qc1", "qh1", "qc2", "qh2")
  
  plots_list_3a[[i]] <- ggplot(data = results_list_3a[[i]][[2]][[1]]) + 
    geom_point(aes(qh1, qc1), color = "red", shape = 1, size = 2.2) +
    geom_point(aes(qh2, qc2), color = "blue", shape = 1, size = 2.2) +
    labs(x = "qh", y = "qc") + 
    scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0,1)) +
    scale_x_continuous(breaks = seq(0, 1, 0.2), limits = c(0,1)) +
    theme(
      panel.grid.major = element_blank(),
      panel.background = element_rect(fill = "white", linetype = "solid", colour = "black"),
      plot.title = element_text(face = "bold")
    ) +
    ggtitle(paste("Epsilon =", epsilon_values[i]))
}

for (i in seq_along(results_list_3b)) {
  
  colnames(results_list_3b[[i]][[2]][[1]]) <- c("qc1", "qh1", "qc2", "qh2")
  
  plots_list_3b[[i]] <- ggplot(data = results_list_3b[[i]][[2]][[1]]) + 
    geom_point(aes(qh1, qc1), color = "red", shape = 1, size = 2.2) +
    geom_point(aes(qh2, qc2), color = "blue", shape = 1, size = 2.2) +
    labs(x = "qh", y = "qc") + 
    scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0,1)) +
    scale_x_continuous(breaks = seq(0, 1, 0.2), limits = c(0,1)) +
    theme(
      panel.grid.major = element_blank(),
      panel.background = element_rect(fill = "white", linetype = "solid", colour = "black"),
      plot.title = element_text(face = "bold")
    ) +
    ggtitle(paste("Epsilon =", epsilon_values[i]))
}

print(epsilon_table_3a)
print(epsilon_table_3b)


do.call(grid.arrange, c(plots_list_3a, ncol = length(epsilon_values)))

do.call(grid.arrange, c(plots_list_3b, ncol = length(epsilon_values)))
