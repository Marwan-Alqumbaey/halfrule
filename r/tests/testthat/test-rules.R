test_that("thresholds match their closed forms", {
  expect_equal(rule_threshold("power"), 0.5)
  expect_equal(rule_threshold("cost", cost_share = 0.5), 0.75)
  e <- power_elasticity(0.8)
  expect_equal(e, (qnorm(0.975) + qnorm(0.8)) * dnorm(qnorm(0.8)) / 0.8)
  expect_equal(rule_threshold("benefit", power = 0.8), e / (2 * (1 + e)))
  expect_equal(round(rule_threshold("benefit", power = 0.8), 3), 0.248)
  expect_lt(rule_threshold("coverage", power = 0.8), 0)
  expect_error(rule_threshold("cost", cost_share = 2))
  expect_error(power_elasticity(1))
})

test_that("subgroup threshold tends to one half and matches the two-group formula", {
  expect_equal(subgroup_threshold(0), 0.5)
  prev <- c(0.25, 0.5, 0.75)
  expect_equal(subgroup_threshold((1 - prev) / prev), sqrt(prev) / (1 + sqrt(prev)))
  expect_error(subgroup_threshold(-1))
  expect_error(subgroup_threshold(1, reference = "best"))
})

test_that("group thresholds agree with a direct comparison of the two sets", {
  set.seed(1)
  kappa <- runif(2000, -0.5, 1.5); rho <- exp(runif(2000, log(0.01), log(20)))
  better <- (1 + kappa * rho) / sqrt(1 + rho) > 1
  expect_equal(better, kappa > subgroup_threshold(rho))
  kappa_all <- kappa * (1 + rho) / (1 + kappa * rho)
  keep <- 1 + kappa * rho > 0
  expect_equal(better[keep], kappa_all[keep] > subgroup_threshold(rho[keep], reference = "all"))
  expect_true(all(subgroup_threshold(rho) <= 0.5))
  expect_true(all(subgroup_threshold(rho, reference = "all") >= 0.5))
})

test_that("linear decay with even arrival closes at two thirds of the zero-effect time", {
  expect_equal(time_window(function(t) 1 - t / 6.3, upper = 8)$close, 4.2, tolerance = 2e-3)
  expect_equal(time_window(function(t) 1 - t / 9, upper = 9)$close, 6, tolerance = 2e-3)
})

test_that("exponential decay with even arrival closes where exp(u) = 1 + 2u", {
  u <- uniroot(function(u) exp(u) - 1 - 2 * u, c(0.5, 3))$root
  expect_equal(time_window(function(t) exp(-t / 2), upper = 20)$close, 2 * u, tolerance = 5e-3)
})

test_that("benefit objective reproduces the quarter rule", {
  t <- seq(0.0005, 0.9995, by = 0.001)
  fit <- eligibility(effect = 1 - t, objective = "benefit", n_screened = 60)
  e <- fit$elasticity
  expect_equal(fit$threshold, e / (2 * (1 + e)), tolerance = 5e-3)
  expect_gt(fit$enrolled, eligibility(effect = 1 - t)$enrolled)
})

test_that("cost objective reproduces the cost rule", {
  t <- seq(0.0005, 0.9995, by = 0.001)
  fit <- eligibility(effect = 1 - t, objective = "cost", screen_cost = 1, patient_cost = 3)
  expect_equal(fit$threshold, 1 - fit$cost_share / 2, tolerance = 5e-3)
  free <- eligibility(effect = 1 - t, objective = "cost", screen_cost = 1, patient_cost = 0)
  expect_equal(free$included, eligibility(effect = 1 - t)$included)
})

test_that("with a constant effect, enrolment follows information", {
  set.seed(1)
  info <- rgamma(4000, 2, 4)
  fit <- eligibility(effect = 0.3, information = info, objective = "cost",
                     screen_cost = 1, patient_cost = 5)
  cut <- (1 - fit$cost_share) * mean(info[fit$included])
  expect_equal(min(info[fit$included]), cut, tolerance = 0.02)
  expect_true(all(eligibility(effect = 0.3, information = info)$included))
})
