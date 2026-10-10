# Checks against published results, and closed-form thresholds.
suppressMessages(library(halfrule))
dir.create("output", showWarnings = FALSE)
res <- list()
add <- function(example, setting, value, published, note = "") {
  res[[length(res) + 1]] <<- data.frame(example = example, setting = setting, value = value,
                                        published = published, note = note)
}

# 1. Simon and Maitournam (2004), equation (E) with f = 1: patients randomised in the untargeted design
#    divided by patients screened in the targeted design, when marker-negative patients get half the effect.
#    The paper prints the value for 50% prevalence (16/18).
for (prev in c(0.75, 0.5, 0.25)) {
  ratio <- relative_duration(c(1, 0.5), weight = c(prev, 1 - prev)) /
    relative_duration(c(1, 0.5), weight = c(prev, 1 - prev), included = c(TRUE, FALSE))
  add("Simon and Maitournam 2004, eq. (E)", sprintf("prevalence %.2f", prev), ratio,
      if (prev == 0.5) 16 / 18 else NA,
      sprintf("threshold for marker-negative effect %.3f", subgroup_threshold((1 - prev) / prev)))
}

# 2. Rudser et al. (2014), main example: latent groups with prevalence 0.4, biomarker N(2, 1) in the group
#    that benefits fully and N(0, 1) in the other group, effect 0.28 SD in the first group, unadjusted analysis.
#    Published: enrolling everyone screens the fewest patients when the second group keeps at least 0.295 of
#    the effect.
rudser_threshold <- function(prev, mu, effect, n_grid = 4000) {
  q <- (seq_len(n_grid) - 0.5) / n_grid
  x <- sapply(q, function(p) uniroot(function(z) prev * pnorm(z - mu) + (1 - prev) * pnorm(z) - p, c(-10, 15))$root)
  ppv <- prev * dnorm(x - mu) / (prev * dnorm(x - mu) + (1 - prev) * dnorm(x))
  all_best <- function(a) {
    eff <- effect * (ppv + (1 - ppv) * a)
    # Within a biomarker value, treated patients are a mixture of the two groups.
    var <- 1 + 0.5 * ppv * (1 - ppv) * (1 - a)^2 * effect^2
    eligibility(eff, var, prognosis = 0)$enrolled > 0.9999
  }
  uniroot(function(a) all_best(a) - 0.5, c(0, 1), tol = 1e-4)$root
}
add("Rudser et al. 2014", "main example", rudser_threshold(0.4, 2, 0.28), 0.295,
    "largest retained fraction at which restriction still helps")
# Complete separation of the groups, prevalence 0.95: the exact two-group threshold with an adjusted
# analysis, and with an unadjusted analysis for two standardised effects.
add("Rudser et al. 2014", "separation, prevalence 0.95, adjusted", subgroup_threshold(0.05 / 0.95), NA,
    "two-group threshold")
for (d in c(0.28, 1)) {
  f <- function(a) (eligibility(c(d, a * d), 1, c(0.95, 0.05), prognosis = 0)$enrolled > 0.999) - 0.5
  add("Rudser et al. 2014", sprintf("separation, prevalence 0.95, unadjusted, effect %.2f SD", d),
      uniroot(f, c(0, 1), tol = 1e-6)$root, NA,
      "published 0.53 is the largest value across their settings")
}

# 3. Prognostic enrichment example of Kerr et al. (2017): event rate 20%, 30% relative reduction,
#    marker AUC 0.72, 90% power, screening cost 1000, cost per enrolled patient 10000.
#    Published: lowest cost near 85% of patients screened out.
mu <- sqrt(2) * qnorm(0.72)
q <- (seq_len(20000) - 0.5) / 20000
x <- sapply(q, function(p) uniroot(function(z) 0.2 * pnorm(z - mu) + 0.8 * pnorm(z) - p, c(-10, 12))$root)
risk <- 0.2 * dnorm(x - mu) / (0.2 * dnorm(x - mu) + 0.8 * dnorm(x))
p0 <- risk; p1 <- 0.7 * risk
fit <- eligibility(effect = p0 - p1, variance = (p0 * (1 - p0) + p1 * (1 - p1)) / 2, prognosis = p0,
                   objective = "cost", screen_cost = 1000, patient_cost = 10000)
add("Kerr et al. 2017", "share screened out, cost rule", 1 - fit$enrolled, 0.85,
    sprintf("cost share w = %.3f", fit$cost_share))
if (requireNamespace("BioPET", quietly = TRUE)) {
  u <- seq(0, 0.99, by = 0.0025)
  cost <- sapply(u, function(v) {
    ss <- BioPET:::sample_size(mean(risk[q > v]), 0.3, 0.025, 0.9, "one.sided")
    10000 * ss + 1000 * ss / (1 - v)
  })
  add("Kerr et al. 2017", "share screened out, their cost formula", u[which.min(cost)], 0.85,
      "direct minimisation on a grid of 0.25%")
}

# 4. Time windows with even arrival.
for (zero in c(5.0, 6.3, 13.8)) {
  add("Linear decay, even arrival", sprintf("zero effect at %.1f h", zero),
      time_window(function(t) 1 - t / zero, upper = 24, n_grid = 48000)$close, 2 * zero / 3, "closing time, hours")
}
u <- uniroot(function(u) exp(u) - 1 - 2 * u, c(0.5, 3), tol = 1e-10)$root
add("Exponential decay, even arrival", "time constant 1", time_window(function(t) exp(-t), upper = 20, n_grid = 40000)$close,
    u, "closing time; exp(u) = 1 + 2u")
t <- (seq_len(100000) - 0.5) / 100000
add("Linear decay, arrivals up to zero-effect time", "everyone versus best window",
    relative_duration(1 - t), 32 / 27, "relative duration")
add("Linear decay, arrivals up to zero-effect time", "bound 1 + CV^2", inclusion_bound(1 - t), 4 / 3, "")

out <- do.call(rbind, res)
write.csv(out, "output/published_examples.csv", row.names = FALSE)
print(out, row.names = FALSE, digits = 4)

# 5. Thresholds by planned power (equal variances).
pw <- c(0.5, 0.7, 0.8, 0.9, 0.95)
tab <- data.frame(power = pw, elasticity = power_elasticity(pw),
                  benefit = rule_threshold("benefit", pw), coverage = rule_threshold("coverage", pw))
write.csv(tab, "output/thresholds_by_power.csv", row.names = FALSE)
print(tab, row.names = FALSE, digits = 3)
