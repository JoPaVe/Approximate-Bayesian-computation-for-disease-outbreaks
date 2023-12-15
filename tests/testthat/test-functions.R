## Tests for arguments of calc_post_distr_base(...)
# Strategy: Test each Function in arguments and use debugging to test if calc_post_distr_base works properly

# Test that distance between identical data sets is 0
test_that("Distance between identical datasets is 0", {
  expect_equal(DistanceFct(observed_data_Table2, observed_data_Table2), 0)
})

# Test that manually created small differences yield expected distance
test_that("Distance between manually created test matrices is 1", {
  test1 <- c(66, 87, 25, 22,  4,
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
  
  distance <- DistanceFct(test_matrix[[1]], test_matrix[[2]])
  expect_equal(distance, 1)
})

test_that("Column sums of generated matrix match column add upto 1", {
  nrows <- dim(observed_data_Table2[[1]])[1]
  ncols <- dim(observed_data_Table2[[1]])[2]
  
  for (i in 1:10) {
    random_matrix <- round(matrix(replicate(nrows * ncols, runif(1,0,5)), 
                                  nrow = nrows, 
                                  ncol = ncols) * observed_data_Table2[[1]])
    
    generated_matrix <- round(CreateMatrix(runif(1,0,1), 
                                           runif(1,0,1), 
                                           random_matrix
    ))
    
    expect_true(all(colSums(random_matrix, na.rm = TRUE) == colSums(generated_matrix, na.rm = TRUE)))
  }
})
