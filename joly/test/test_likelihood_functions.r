library(testthat)


setwd("/home/rasmus-emil/github/thesis/joly/")
source("likelihood_functions.R")


test_that("case_1_likelihood_calculates_correctly_for_linear_splines", {
  V_0 <- 1
  V_healthy <- 2
  T_obs <- 4
  
  a01 <- function(x) 1

  a02 <- function(x) 1/2
  
  a12 <- function(x) 1/4
  
  A01 <- function(x) x
  
  A02 <- function(x) 1/2*x
  
  A12 <- function(x) 1/4*x
  
  calculated <- case_1_likelihood(V_0, V_healthy, T_obs, a01, a02, a12, A01, A02, A12)
  
  expected <- 1/(exp(-1*1-1/2*1)) *(exp(-1*4-1/2*4)+4*(exp(5/2)-1)/(5*exp(6)) )
  
  expect_equal(calculated, expected, tolerance = 10^(-6))
})

test_that("case_2_likelihood_calculates_correctly_for_linear_splines", {
  V_0 <- 1
  V_healthy <- 2
  T_obs <- 4

  a01 <- function(x) 1

  a02 <- function(x) 1/2

  a12 <- function(x) 1/4

  A01 <- function(x) x

  A02 <- function(x) 1/2*x

  A12 <- function(x) 1/4*x
  
  expected <- 1/(exp(-1*1-1/2*1)) *(exp(-1*4-1/2*4)*1/2+4*1/4*(exp(5/2)-1)/(5*exp(6)) )
  
  calculated <- case_2_likelihood(V_0, V_healthy, T_obs, a01, a02, a12, A01, A02, A12)
  
  expect_equal(calculated, expected, tolerance = 10^(-6))
  
})
 
test_that("case_3_likelihood_calculates_correctly_for_linear_splines", {
  V_0 <- 1
  V_healthy <- 2
  V_ill <- 3
  T_obs <- 4
  
  a01 <- function(x) 1

  a02 <- function(x) 1/2

  a12 <- function(x) 1/4

  A01 <- function(x) x

  A02 <- function(x) 1/2*x

  A12 <- function(x) 1/4*x

  
  expected <- 1/(exp(-1*1-1/2*1)) *4*(exp(5/4)-1)/(5*exp(19/4)) 
  
  calculated <- case_3_likelihood(V_0, V_healthy, V_ill, T_obs, a01, a02, a12, A01, A02, A12)
  
  expect_equal(calculated, expected, tolerance = 10^(-6))
})

test_that("case_4_likelihood_calculates_correctly_for_linear_splines", {
  V_0 <- 1
  V_healthy <- 2
   V_ill <- 3
  T_obs <- 4

  a01 <- function(x) 1

  a02 <- function(x) 1/2

  a12 <- function(x) 1/4

  A01 <- function(x) x

  A02 <- function(x) 1/2*x

  A12 <- function(x) 1/4*x

  expected <- 1/(exp(-1*1-1/2*1)) *(exp(5/4)-1)/(5*exp(19/4)) 
  
  calculated <- case_4_likelihood(V_0, V_healthy, V_ill, T_obs, a01, a02, a12, A01, A02, A12)
  
  expect_equal(calculated, expected, tolerance = 10^(-6))
})