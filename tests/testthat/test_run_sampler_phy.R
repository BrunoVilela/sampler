library(ape)
library(sampler)
context("run_sampler_phy")

tree <- rcoal(100)
## Calculate the distance
## Highly overdispersed 50% resample design (alpha = 50)
n <- 10

test_that("test run_sampler_phy", {
  selection <- run_sampler_phy(x = tree, n = n, alpha = 100, starting = "t10")
  expect_equal(class(selection), "phylo")
  expect_equal(Ntip(selection), n)
  expect_true(any(selection$tip.label %in%  "t10"))
  expect_error(run_sampler_phy(x = 1, n = n, alpha = 100, starting = "t10"))
  expect_error(run_sampler_phy(x = tree, n = n, alpha = 100, starting = "t10",
                               dist.func = 1))
  expect_error(run_sampler_phy(x = tree, n = n, alpha = 100, starting = "t10",
                               dist.func = function(x){x}))
  if (require(adephylo)) {
    selection <- run_sampler_phy(x = tree, n = n, alpha = 100, starting = "t10",
                                 dist.func = function(x){as.matrix(distTips(x))})
    expect_equal(class(selection), "phylo")
    expect_equal(Ntip(selection), n)
    expect_true(any(selection$tip.label %in%  "t10"))
  }

})

