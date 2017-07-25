library(ape)
library(sampler)

context("Nee_May_1997")

set.seed(50)
tree <- rcoal(50)
n <- 10


test_that("test run_sampler_geo", {
  selection <- Nee_May_1997(x = tree, n = n)
  expect_equal(class(selection), "phylo")
  expect_equal(Ntip(selection), n)

  expect_true(sum(selection$edge.length) > sum(sample(tree$edge.length, 10)))

  expect_error(Nee_May_1997(x = 1, n = n))
  expect_error(Nee_May_1997(x = tree, n = 1))
  expect_error(Nee_May_1997(x = tree, n = c(1, 2)))
  expect_error(Nee_May_1997(x = tree, n = "a"))
  expect_error(Nee_May_1997(x = tree, n = 2.3))
  expect_error(Nee_May_1997(x = tree, n = 50000))
})

