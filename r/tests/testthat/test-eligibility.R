brute <- function(value, k) {
  best <- -Inf; arg <- NULL
  for (mask in seq_len(2^k - 1)) {
    inc <- bitwAnd(mask, 2^(0:(k - 1))) > 0
    val <- value(inc)
    if (is.finite(val) && val > best) { best <- val; arg <- inc }
  }
  list(value = best, included = arg)
}

test_that("power objective equals exhaustive search", {
  set.seed(11)
  for (rep in 1:25) {
    k <- 10
    g <- rnorm(k, 0.3, 0.6); v <- runif(k, 0.3, 3); w <- runif(k)
    if (!any(g > 0)) next
    fit <- eligibility(g, v, w)
    val <- function(inc) sum((w * g)[inc]) / sqrt(sum((w * v)[inc]))
    expect_equal(val(fit$included), brute(val, k)$value)
    expect_gte(fit$threshold, 0.5)
  }
})

test_that("cost and benefit objectives equal exhaustive search", {
  set.seed(12)
  for (rep in 1:25) {
    k <- 10
    g <- rnorm(k, 0.3, 0.6); v <- runif(k, 0.3, 3); w <- runif(k); w <- w / sum(w)
    h <- runif(k, 0, 2)
    if (!any(g > 0)) next
    cost <- function(inc) {
      G <- sum((w * g)[inc]); if (G <= 0) return(-Inf)
      -(4 * sum((w * v)[inc]) / G^2) * (1 + 2 * sum(w[inc]))
    }
    fit <- eligibility(g, v, w, objective = "cost", screen_cost = 1, patient_cost = 2)
    expect_equal(cost(fit$included), brute(cost, k)$value)
    ben <- function(inc) {
      G <- sum((w * g)[inc]); if (G <= 0) return(-Inf)
      pnorm(sqrt(50) * G / sqrt(4 * sum((w * v)[inc])) - qnorm(0.975)) * sum((w * h)[inc])
    }
    fit <- eligibility(g, v, w, objective = "benefit", benefit = h, n_screened = 50)
    expect_equal(ben(fit$included), brute(ben, k)$value)
  }
})

test_that("unadjusted analysis equals exhaustive search and penalises atypical prognosis", {
  set.seed(13)
  for (rep in 1:25) {
    k <- 10
    g <- abs(rnorm(k, 0.3, 0.3)); v <- runif(k, 0.5, 2); w <- runif(k); mu <- rnorm(k)
    val <- function(inc) {
      f <- sum(w[inc]); m0 <- sum((w * mu)[inc]) / f; m1 <- sum((w * (mu + g))[inc]) / f
      noise <- v + 0.5 * (mu - m0)^2 + 0.5 * (mu + g - m1)^2
      sum((w * g)[inc]) / sqrt(sum((w * noise)[inc]))
    }
    fit <- eligibility(g, v, w, prognosis = mu)
    expect_equal(val(fit$included), brute(val, k)$value, tolerance = 1e-8)
  }
  x <- seq(-3, 3, length.out = 601)
  adj <- eligibility(effect = 1, variance = 1, weight = dnorm(x))
  unadj <- eligibility(effect = 1, variance = 1, weight = dnorm(x), prognosis = 2 * x)
  expect_true(all(adj$included))
  expect_lt(unadj$enrolled, 1)
})

test_that("two-subgroup results agree with the published efficiency table", {
  # Number screened, untargeted relative to targeted design, when the effect in
  # marker-negative patients is half that in marker-positive patients.
  for (case in list(c(0.75, 0.98), c(0.5, 0.89), c(0.25, 0.64))) {
    prev <- case[1]
    ratio <- relative_duration(c(1, 0.5), weight = c(prev, 1 - prev)) /
      relative_duration(c(1, 0.5), weight = c(prev, 1 - prev), included = c(TRUE, FALSE))
    expect_equal(round(ratio, 2), case[2])
  }
  fit <- eligibility(c(1, 0.34), weight = c(0.25, 0.75))
  expect_true(all(fit$included))
  fit <- eligibility(c(1, 0.32), weight = c(0.25, 0.75))
  expect_equal(fit$included, c(TRUE, FALSE))
})

test_that("enrolling everyone never costs more than the bound", {
  set.seed(14)
  for (rep in 1:20) {
    g <- rexp(200); v <- runif(200, 0.5, 2); w <- runif(200)
    expect_lte(relative_duration(g, v, w), inclusion_bound(g, v, w) + 1e-12)
  }
  t <- seq(0.0005, 0.9995, by = 0.001)
  expect_equal(inclusion_bound(1 - t), 4 / 3, tolerance = 1e-4)
  expect_equal(relative_duration(1 - t), (4 / 9)^2 / (2 / 3) / 0.25, tolerance = 1e-3)
})

test_that("inputs are validated", {
  expect_error(eligibility(c(-1, -2)))
  expect_error(eligibility(c(1, NA)))
  expect_error(eligibility(c(1, 2), variance = c(1, -1)))
  expect_error(eligibility(c(1, 2), weight = c(1, 2, 3)))
  expect_error(eligibility(c(1, 2), objective = "benefit"))
  expect_error(relative_duration(c(1, 2), included = c(TRUE, FALSE, TRUE)))
  expect_equal(relative_duration(1, weight = c(1, 1), included = c(TRUE, FALSE)), 2)
})

test_that("design_power matches the textbook two-sample formula", {
  out <- design_power(effect = 0.3, variance = 1, n_screened = 400)
  expect_equal(out$ncp, 0.3 * sqrt(400) / 2)
  expect_equal(out$power, pnorm(3 - qnorm(0.975)))
  score <- design_power(effect = 0.5, information = 0.25 * 0.2 * 0.8, n_screened = 1000)
  expect_equal(score$ncp, 0.5 * sqrt(1000 * 0.25 * 0.16))
  expect_error(design_power(effect = c(1, 2), included = c(FALSE, FALSE), n_screened = 10))
})
