# Approximate-Bayesian-computation-for-disease-outbreaks

## R folder
Contains the following files
  1. ABC_algorithm.R -> Used for code of ABC for model selection on the joint space
  2.  Plotting.R -> Used to create Plots
  3.  Comparison.R -> Used to compare plots of different epsilon values.
  4.  Main.R -> Running of all code and description (Calling Comparison.R has been commented as it takes a lot of time to compute at half a million iterations)

## tests folder
  1. Added tests folder to the repository.
  2. Folder contains test-functions.R file, which contain all tests we conducted.
  3. Ideally, using devtools::test_dir(’.’), the tests should run and your code should pass the tests.
