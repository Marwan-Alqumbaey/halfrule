# Checks of the software and of the rules by exhaustive search.
#  1. For random problems with 12 patient types, compare the set returned by eligibility()
#     with the best of all 4095 possible sets, for each objective.
#  2. Count how often more than one set satisfies the half rule (it is a condition
#     that the best set must meet, not a rule that picks it out).
#  3. Check the two-group rules on random pairs of groups.
suppressMessages(library(halfrule))
set.seed(20261010)
dir.create("output", showWarnings = FALSE)

# ---- 1 and 2. Random problems ------------------------------------------------------------------------------
k <- 12; n_problems <- 1000
member <- sapply(0:(k - 1), function(j) bitwAnd(seq_len(2^k - 1), 2^j) > 0) * 1   # 4095 x 12
z <- qnorm(0.975)
agree <- c(power = 0, cost = 0, benefit = 0, unadjusted = 0)
done <- 0; worst_gap <- 0; prefix_ok <- 0; several_stationary <- 0
for (i in seq_len(n_problems)) {
  g <- rnorm(k, 0.3, 0.6); v <- runif(k, 0.3, 3); w <- runif(k); w <- w / sum(w)
  h <- runif(k, 0, 2); mu <- rnorm(k)
  if (!any(g > 0)) next
  done <- done + 1
  G <- drop(member %*% (w * g)); V <- drop(member %*% (w * 4 * v)); F <- drop(member %*% w); H <- drop(member %*% (w * h))
  val <- list(
    power = ifelse(G > 0, G / sqrt(V), -Inf),
    cost = ifelse(G > 0, -(V / G^2) * (1 + 2 * F), -Inf),
    benefit = ifelse(G > 0, pnorm(sqrt(50) * G / sqrt(V) - z) * H, -Inf))
  fits <- list(power = eligibility(g, v, w),
               cost = eligibility(g, v, w, objective = "cost", screen_cost = 1, patient_cost = 2),
               benefit = eligibility(g, v, w, objective = "benefit", benefit = h, n_screened = 50))
  for (o in names(fits)) {
    row <- sum(2^(which(fits[[o]]$included) - 1))
    agree[o] <- agree[o] + (val[[o]][row] >= max(val[[o]]) - 1e-10 * abs(max(val[[o]])))
  }
  # The best set for power is a top slice of the types ranked by effect-to-variance ratio.
  r <- g / (4 * v); a <- fits$power$included
  prefix_ok <- prefix_ok + (min(r[a]) > max(c(-Inf, r[!a])))
  # Number of top slices that meet the half rule (every enrolled type at least half the pooled ratio,
  # every excluded type at most half).
  o <- order(r, decreasing = TRUE)
  stationary <- sapply(seq_len(k), function(m) {
    s <- o[seq_len(m)]; pooled <- sum((w * g)[s]) / sum((w * 4 * v)[s])
    pooled > 0 && all(r[s] >= pooled / 2) && all(r[-s] <= pooled / 2)
  })
  several_stationary <- several_stationary + (sum(stationary) > 1)
  # Unadjusted analysis (positive effects).
  ga <- abs(g)
  Gu <- drop(member %*% (w * ga))
  m0 <- drop(member %*% (w * mu)) / F; m1 <- drop(member %*% (w * (mu + ga))) / F
  Vu <- drop(member %*% (w * v)) + 0.5 * (drop(member %*% (w * mu^2)) - F * m0^2) +
    0.5 * (drop(member %*% (w * (mu + ga)^2)) - F * m1^2)
  vu <- Gu / sqrt(Vu)
  fu <- eligibility(ga, v, w, prognosis = mu)
  got <- vu[sum(2^(which(fu$included) - 1))]
  agree["unadjusted"] <- agree["unadjusted"] + (got >= max(vu) * (1 - 1e-10))
  worst_gap <- max(worst_gap, 1 - got / max(vu))
}
checks <- data.frame(
  check = c("power: best set found", "cost: best set found", "benefit: best set found",
            "unadjusted analysis: best set found", "power: best set is a top slice of the ranking",
            "more than one top slice meets the half rule"),
  problems = done,
  count = c(agree, prefix_ok, several_stationary),
  largest_shortfall = c(0, 0, 0, worst_gap, NA, NA))
write.csv(checks, "output/exhaustive_checks.csv", row.names = FALSE)
print(checks, row.names = FALSE)

# ---- 3. Two-group rules ---------------------------------------------------------------------------------------
# Group A is enrolled; group B is a candidate. kappa = ratio of B / pooled ratio of A, rho = noise of B / noise of A.
n_pairs <- 100000
kappa <- runif(n_pairs, -0.5, 1.5); rho <- exp(runif(n_pairs, log(0.01), log(20)))
gain <- (1 + kappa * rho) / sqrt(1 + rho) > 1                  # adding B raises G / sqrt(V)
exact <- kappa > subgroup_threshold(rho)
pooled_with_b <- (1 + kappa * rho) / (1 + rho)                  # pooled ratio of A and B together, relative to A
kappa_c <- kappa / pooled_with_b                                 # ratio of B relative to the pooled ratio with B
pairs <- data.frame(
  check = c("exact threshold agrees with direct comparison",
            "B at least half of A: adding B shortens the trial",
            "B below half of A and B together: dropping B shortens the trial"),
  cases = c(n_pairs, sum(kappa >= 0.5), sum(pooled_with_b > 0 & kappa_c < 0.5)),
  correct = c(sum(gain == exact), sum(gain[kappa >= 0.5]), sum(!gain[pooled_with_b > 0 & kappa_c < 0.5])))
write.csv(pairs, "output/two_group_checks.csv", row.names = FALSE)
print(pairs, row.names = FALSE)
