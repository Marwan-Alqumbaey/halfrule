# Simulated trials. For nested eligibility sets, compare simulated power with the power predicted from the
# non-centrality parameter, and the set with the highest simulated power with the set given by the rule.
suppressMessages(library(halfrule))
set.seed(20261010)
dir.create("output", showWarnings = FALSE)
reps <- 20000
z_crit <- qnorm(0.975)

# Welch z statistics for nested sets: patients are ordered by `key`, and set k holds the first cuts[k] patients.
nested_z <- function(key, y, arm, cuts) {
  o <- order(key); y <- y[o]; a <- arm[o]
  n1 <- cumsum(a)[cuts]; n0 <- cuts - n1
  s1 <- cumsum(y * a)[cuts]; s0 <- cumsum(y * (1 - a))[cuts]
  q1 <- cumsum(y^2 * a)[cuts]; q0 <- cumsum(y^2 * (1 - a))[cuts]
  v1 <- (q1 - s1^2 / n1) / (n1 - 1); v0 <- (q0 - s0^2 / n0) / (n0 - 1)
  (s1 / n1 - s0 / n0) / sqrt(v1 / n1 + v0 / n0)
}

# z statistics for the treatment coefficient in a linear regression of y on arm and x, for nested sets.
nested_z_adjusted <- function(key, y, arm, x, cuts) {
  o <- order(key); X <- cbind(1, arm[o], x[o]); y <- y[o]
  cx <- apply(X[, c(1, 1, 1, 2, 2, 3)] * X[, c(1, 2, 3, 2, 3, 3)], 2, cumsum)
  cy <- apply(X * y, 2, cumsum); cyy <- cumsum(y^2)
  sapply(cuts, function(m) {
    xtx <- matrix(cx[m, c(1, 2, 3, 2, 4, 5, 3, 5, 6)], 3)
    inv <- solve(xtx); b <- inv %*% cy[m, ]
    s2 <- (cyy[m] - sum(b * cy[m, ])) / (m - 3)
    b[2] / sqrt(s2 * inv[2, 2])
  })
}

# ---- Scenarios A and B: time window, continuous outcome, two arrival patterns --------------------------------
run_window <- function(label, rdelay, ddelay, m, d0 = 0.30, zero = 6.3, upper = 8) {
  windows <- seq(1, upper, by = 0.25)
  eff <- function(t) d0 * (1 - t / zero)
  reject <- matrix(0, reps, length(windows))
  for (r in seq_len(reps)) {
    t <- rdelay(m); arm <- rbinom(m, 1, 0.5)
    y <- rnorm(m) + arm * eff(t)
    cuts <- findInterval(windows, sort(t))
    reject[r, ] <- nested_z(t, y, arm, cuts) > z_crit
  }
  grid <- (seq_len(4000) - 0.5) * upper / 4000
  pred <- sapply(windows, function(w) design_power(eff(grid), 1, ddelay(grid), prognosis = 0,
                                                    included = grid <= w, n_screened = m)$power)
  gain <- sapply(windows, function(w) sum((eff(grid) * ddelay(grid))[grid <= w]) / sum(ddelay(grid)))
  rule_power <- time_window(eff, ddelay, upper = upper, prognosis = 0)$close
  rule_benefit <- time_window(eff, ddelay, upper = upper, prognosis = 0, objective = "benefit", n_screened = m)$close
  data.frame(scenario = label, set = windows, simulated = colMeans(reject), predicted = pred,
             benefit_simulated = colMeans(reject) * gain, benefit_predicted = pred * gain,
             rule_power = rule_power, rule_benefit = rule_benefit)
}
sA <- run_window("A: even arrival", function(n) runif(n, 0, 8), function(t) rep(1, length(t)), m = 1500)
sB <- run_window("B: later arrival more common", function(n) 8 * sqrt(runif(n)), function(t) t, m = 4000)

# ---- Scenario C: binary outcome, odds ratio falls with delay, prognosis varies, unadjusted analysis -----------
run_binary <- function(m = 1300, b0 = 0.75, zero = 6.3, upper = 8) {
  windows <- seq(1, upper, by = 0.25)
  reject <- matrix(0, reps, length(windows))
  for (r in seq_len(reps)) {
    t <- runif(m, 0, upper); s <- rnorm(m); arm <- rbinom(m, 1, 0.5)
    y <- rbinom(m, 1, plogis(-0.6 + s + arm * b0 * (1 - t / zero)))
    cuts <- findInterval(windows, sort(t))
    reject[r, ] <- nested_z(t, y, arm, cuts) > z_crit
  }
  # Patient types are delays; the prognostic score s is averaged over within each delay.
  nt <- 400; tg <- (seq_len(nt) - 0.5) * upper / nt; sg <- qnorm((seq_len(400) - 0.5) / 400)
  p0 <- outer(tg, sg, function(t, s) plogis(-0.6 + s))
  p1 <- outer(tg, sg, function(t, s) plogis(-0.6 + s + b0 * (1 - t / zero)))
  eff <- rowMeans(p1 - p0); mu0 <- rowMeans(p0)
  # Within-delay variance, averaged over arms: binomial variance plus the spread of risk across s.
  var <- rowMeans((p0 * (1 - p0) + p1 * (1 - p1)) / 2) + (apply(p0, 1, var) + apply(p1, 1, var)) / 2
  pred <- sapply(windows, function(w) design_power(eff, var, prognosis = mu0, included = tg <= w, n_screened = m)$power)
  fit <- eligibility(eff, var, prognosis = mu0)
  data.frame(scenario = "C: binary outcome", set = windows, simulated = colMeans(reject), predicted = pred,
             benefit_simulated = NA, benefit_predicted = NA,
             rule_power = max(tg[fit$included]) + upper / nt / 2, rule_benefit = NA)
}
sC <- run_binary()

# ---- Scenario D: constant effect, strong prognostic covariate, unadjusted or adjusted analysis --------------
run_band <- function(m = 1500, d = 0.2, slope = 1.5) {
  bands <- seq(0.4, 3, by = 0.2)
  rej_un <- rej_adj <- matrix(0, reps, length(bands))
  for (r in seq_len(reps)) {
    x <- rnorm(m); arm <- rbinom(m, 1, 0.5); y <- slope * x + d * arm + rnorm(m)
    cuts <- findInterval(bands, sort(abs(x)))
    rej_un[r, ] <- nested_z(abs(x), y, arm, cuts) > z_crit
    rej_adj[r, ] <- nested_z_adjusted(abs(x), y, arm, x, cuts) > z_crit
  }
  x <- qnorm((seq_len(4000) - 0.5) / 4000)
  pred_un <- sapply(bands, function(b) design_power(d, 1, prognosis = slope * x, included = abs(x) <= b, n_screened = m)$power)
  pred_adj <- sapply(bands, function(b) design_power(d, 1, included = abs(x) <= b, n_screened = m)$power)
  fit <- eligibility(d, 1, prognosis = slope * x)
  rbind(data.frame(scenario = "D: unadjusted analysis", set = bands, simulated = colMeans(rej_un), predicted = pred_un,
                   benefit_simulated = NA, benefit_predicted = NA, rule_power = max(abs(x[fit$included])), rule_benefit = NA),
        data.frame(scenario = "D: adjusted analysis", set = bands, simulated = colMeans(rej_adj), predicted = pred_adj,
                   benefit_simulated = NA, benefit_predicted = NA, rule_power = Inf, rule_benefit = NA))
}
sD <- run_band()

sim <- rbind(sA, sB, sC, sD)
write.csv(sim, "output/simulation_curves.csv", row.names = FALSE)

summary_table <- do.call(rbind, lapply(split(sim, sim$scenario), function(d) {
  at_rule <- which.min(abs(d$set - min(d$rule_power[1], max(d$set))))
  data.frame(scenario = d$scenario[1],
             rule = d$rule_power[1],
             best_simulated = d$set[which.max(d$simulated)],
             power_best_simulated = max(d$simulated),
             power_nearest_rule = d$simulated[at_rule],
             largest_error = max(abs(d$simulated - d$predicted)),
             rule_benefit = d$rule_benefit[1],
             best_benefit_simulated = if (all(is.na(d$benefit_simulated))) NA else d$set[which.max(d$benefit_simulated)])
}))
write.csv(summary_table, "output/simulation_summary.csv", row.names = FALSE)
print(summary_table, row.names = FALSE, digits = 3)
