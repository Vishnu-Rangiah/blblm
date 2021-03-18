test_that("blblm package functions tests cases", {

  #Test for Blblm AND time

  bnch = bench::mark(

    fit <- blblm(mpg ~ wt * hp, data = mtcars, m = 3, B = 1000, parallel = TRUE, nthreads = 2), #limited to only 2 cores
    fit2 <- blblm(mpg ~ wt * hp, data = mtcars, m = 3, B = 1000),
    check = FALSE

    )
  #View(bnch)
  #expect_lte(bnch$total_time[1],bnch$total_time[2])

  #length of coeffs
  fit_len = length(coef(fit))
  expect_identical(fit_len,4)
  fit2_len = length(coef(fit2))
  expect_identical(fit2_len,4)

  #class of coeffs
  fit_coef = class(coef(fit))
  fit2_coef = class(coef(fit2))
  expect_identical(fit_coef,fit2_coef)

  #a represents (fit) a single process and b is the (fit2) parallel process
  #we should expect the parallel process to have a quicker total run time for larger datasets
  #however since we have a smaller data-set the initiation time for parallel computing may cause fit to have a longer time than fit2

  #We will similarly test the length of the other exported functions in our package to find any potential errors.


  #length of confint
  cfl = length(confint(fit, c("wt", "hp")))
  expect_identical(cfl,4)
  cfl2 = length(confint(fit2, c("wt", "hp"))) #4
  expect_identical(cfl2,4)

  #length of sigma
  sl = length(sigma(fit))
  expect_identical(sl,1)
  sl2 = length(sigma(fit2)) #1
  expect_identical(sl2,1)

  #length of sigma with confidence
  scl = length(sigma(fit, confidence = TRUE)) #3
  expect_identical(scl,3)
  scl2 = length(sigma(fit2, confidence = TRUE))
  expect_identical(scl2,3)

  #length of sigma with confidence
  pl = length(predict(fit, data.frame(wt = c(2.5, 3), hp = c(150, 170)))) #2
  expect_identical(pl,2)
  pl2 = length(predict(fit, data.frame(wt = c(2.5, 3), hp = c(150, 170)))) #2
  expect_identical(pl2,2)

  #length of sigma with confidence
  pcl = length(predict(fit, data.frame(wt = c(2.5, 3), hp = c(150, 170)), confidence = TRUE))#6
  expect_identical(pcl,6)
  pcl2 = length(predict(fit, data.frame(wt = c(2.5, 3), hp = c(150, 170)), confidence = TRUE))#6
  expect_identical(pcl2,6)



  })