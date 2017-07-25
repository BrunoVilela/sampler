library(sp)
library(sampler)
library(raster)
library(maptools)
context("run_sampler_geo")

data(wrld_simpl)  # World map
Brazil <- wrld_simpl[wrld_simpl$NAME == "Brazil", ]  # Brazil (polygon)
coords <- slot(spsample(Brazil, 100, "regular"), "coords")
rownames(coords) <- paste0("t", 1:nrow(coords))

## Calculate the distance
## Highly overdispersed 50% resample design (alpha = 50)
n <- 10

test_that("test run_sampler_geo", {
  selection <- run_sampler_geo(x = coords, n = n, alpha = 100, starting = "t10")
  expect_equal(class(selection), "matrix")
  expect_equal(nrow(selection), n)
  expect_equal(ncol(selection), 2)
  expect_true(any(rownames(selection) %in%  "t10"))
  expect_error(run_sampler_geo(x = 1, n = n, alpha = 100, starting = "t10"))
  expect_error(run_sampler_geo(x = coords, n = n, alpha = 100, starting = "t10",
                               dist.func = 1))
  expect_error(run_sampler_geo(x = coords, n = n, alpha = 100, starting = "t10",
                               dist.func = function(x){x}))

    selection <- run_sampler_geo(x = coords, n = n, alpha = 100, starting = "t10",
                                 dist.func = function(x){as.matrix(dist(x))})
    expect_equal(class(selection), "matrix")
    expect_equal(nrow(selection), n)
    expect_equal(ncol(selection), 2)
    expect_true(any(rownames(selection) %in%  "t10"))

  selection <- run_sampler_geo(x = coords, n = n, alpha = 100, starting = "t10",
                               return_start = TRUE)
  expect_equal(class(selection[[1]]), "matrix")
  expect_equal(nrow(selection[[1]]), n)
  expect_true(any(rownames(selection[[1]]) %in%  "t10"))
  expect_true(all(selection[[2]] ==  "t10"))
})

