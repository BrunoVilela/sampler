library(ape)
library(sampler)
context("run_sampler")

tree <- rcoal(100)
## Calculate the distance
dist2 <- cophenetic(tree)
## Highly overdispersed 50% resample design (alpha = 50)
n <- 10

test_that("test run_sampler", {
  selection <- run_sampler(x = dist2, n = n, alpha = 100, starting = "t10")
  expect_equal(class(selection), "character")
  expect_equal(length(selection), n)
  expect_true(any(selection %in%  "t10"))

  selection <- run_sampler(x = dist2, n = n, alpha = 100)
  expect_equal(class(selection), "character")
  expect_equal(length(selection), n)

  selection <- run_sampler(x = dist2, n = n, alpha = 100,
                           starting = c("t10", "t3"), n_start = 2)
  expect_equal(class(selection), "character")
  expect_equal(length(selection), n)
  expect_true(all(c("t10", "t3") %in%  selection))

  selection <- run_sampler(x = dist2, n = n, alpha = 100,
                           starting = c("t10", "t3"),
                           n_start = 2, return_start = TRUE)
  expect_equal(class(selection[[1]]), "character")
  expect_equal(length(selection[[1]]), n)
  expect_true(all(c("t10", "t3") %in%  selection[[1]]))
  expect_true(all(selection[[2]] ==  c("t10", "t3")))
})

test_that("test run_sampler alpha and aggregation versus overdispersed", {
  selection <- run_sampler(x = dist2, n = n, alpha = 50, starting = "t10")
  selection2 <- run_sampler(x = dist2, n = n, alpha = -50, starting = "t10")
  pd1 <- sum(drop.tip(tree, tree$tip.label[!tree$tip.label %in% selection])$edge.length)
  pd2 <- sum(drop.tip(tree, tree$tip.label[!tree$tip.label %in% selection2])$edge.length)
  expect_true(pd1 > pd2)
})

test_that("test run_sampler errors", {
  expect_error(run_sampler(x = dist2, n = n, alpha = 50, starting = "t102"))
  expect_error(run_sampler(x = dist2, n = n, alpha = 10000, starting = "t99"))
  expect_error(run_sampler(x = dist2, n = nrow(dist) + 1, alpha = 0))
  expect_error(run_sampler(x = dist2[, -1], n = nrow(dist) + 1, alpha = 0))
  expect_error(run_sampler(x = dist2, n = 1, alpha = 0))
  expect_error(run_sampler(x = dist2, n = 2.2, alpha = 0))
  dist3 <- dist2
  colnames(dist3) <- NULL
  dist3 <- dist2
  colnames(dist3)[1:2] <- rownames(dist3)[1:2] <- c("t1", "t1")
  expect_error(run_sampler(x = dist3, n = n, alpha = 0))
  expect_error(run_sampler(x = dist2, n = "a", alpha = 0))
  expect_error(run_sampler(x = dist2, n = c(3, 4), alpha = 0))
  expect_error(run_sampler(x = dist2, n = n, alpha = c(2, 3)))
  expect_error(run_sampler(x = dist2, n = n, alpha = 0, n_start = 1002))
  expect_error(run_sampler(x = dist2, n = n, alpha = 0, n_start = 2.2))
  dist3 <- dist2
  dist3[1, 1] <- NA
  expect_error(run_sampler(x = dist3, n = n, alpha = 0, n_start = 2))
  dist3[1, 1] <- Inf
  expect_error(run_sampler(x = dist3, n = n, alpha = 0, n_start = 2))
  dist3[1, 1] <- -1
  expect_error(run_sampler(x = dist3, n = n, alpha = 0, n_start = 2))
})

