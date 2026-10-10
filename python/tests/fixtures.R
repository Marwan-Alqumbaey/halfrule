# Writes fixtures.json: results of the R package for random problems, used to check that the
# Python package gives the same answers. Run from this folder with the R package installed.
library(halfrule)
set.seed(2026)
num <- function(x) unname(as.numeric(x))
cases <- lapply(1:12, function(i) {
  k <- 40
  g <- rnorm(k, 0.3, 0.5); v <- runif(k, 0.3, 3); w <- runif(k)
  prog <- rnorm(k); h <- runif(k, 0, 2); info <- runif(k, 0.1, 1.5)
  idx <- function(fit) which(fit$included) - 1L
  a <- eligibility(g, v, w)
  b <- eligibility(g, v, w, objective = "cost", screen_cost = 1, patient_cost = 2.5)
  cc <- eligibility(g, v, w, objective = "benefit", benefit = h, n_screened = 80)
  d <- eligibility(abs(g), v, w, prognosis = prog)
  e <- eligibility(g, information = info, weight = w, objective = "cost", screen_cost = 1, patient_cost = 4)
  list(effect = g, variance = v, weight = w, prognosis = prog, benefit = h, information = info,
       power = list(included = I(idx(a)), ncp = a$ncp, threshold = a$threshold),
       cost = list(included = I(idx(b)), cost_share = b$cost_share, relative_duration = b$relative_duration),
       benefit_fit = list(included = I(idx(cc)), power = cc$power, elasticity = cc$elasticity),
       unadjusted = list(included = I(idx(d)), ncp = d$ncp),
       score_cost = list(included = I(idx(e)), ncp = e$ncp),
       all_duration = relative_duration(g, v, w),
       bound = inclusion_bound(g, v, w),
       design = design_power(g, v, w, included = g > 0.2, n_screened = 300)$power)
})
writeLines(jsonlite::toJSON(cases, digits = NA, auto_unbox = TRUE), "fixtures.json")
