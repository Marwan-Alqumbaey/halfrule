# Illustration on the Third International Stroke Trial (IST-3).
# The data are not distributed with this repository. Download them from https://doi.org/10.7488/ds/1350
# and set IST3_PATH to the file datashare_aug2015.sas7bdat.
suppressMessages({library(halfrule); library(haven)})
set.seed(20261010)
path <- Sys.getenv("IST3_PATH", "data/datashare_aug2015.sas7bdat")
if (!file.exists(path)) stop("IST-3 data not found. Set IST3_PATH.")
dir.create("output", showWarnings = FALSE)
raw <- as.data.frame(read_sas(path))
ist <- data.frame(treat = as.integer(raw$itt_treat == 0),          # source coding: 0 = alteplase, 1 = control
                  good = as.integer(raw$alivefav6 == 1),            # alive with Oxford Handicap Score 0-1 at six months
                  indep = as.integer(raw$aliveind6 == 1),           # alive with score 0-2 at six months
                  imputed = as.integer(raw$imputed6 == 1),          # six-month score imputed by the trial team
                  delay = raw$randdelay, age = raw$age, nihss = raw$nihss)
flow <- data.frame(randomised = nrow(ist), missing_data = sum(!complete.cases(ist)),
                   delay_over_6h = sum(ist$delay > 6, na.rm = TRUE))
ist <- ist[complete.cases(ist) & ist$delay <= 6, ]
flow$analysed <- nrow(ist); flow$outcome_imputed <- sum(ist$imputed); flow$imputed_good_outcome <- sum(ist$imputed & ist$good)
flow$share_over_80 <- mean(ist$age > 80)
write.csv(flow, "output/ist3_flow.csv", row.names = FALSE)
print(flow, row.names = FALSE)

# Information per patient, p(1 - p)/4, from a prognostic model without treatment terms.
info_of <- function(d, outcome) {
  m <- glm(as.formula(paste(outcome, "~ delay + age + nihss + I(nihss^2)")), binomial, d)
  p <- fitted(m)
  p * (1 - p) / 4
}

# ---- A. Time window with the effect profile of the pooled analysis -----------------------------------------
# The log odds ratio falls linearly with delay and reaches zero at `zero` hours. Only the shape is used.
# For expected benefit, the health gain of a patient is the absolute risk difference, about 4 * beta * i.
window_analysis <- function(d, outcome, zero) {
  info <- info_of(d, outcome)
  shape <- zero - d$delay
  best <- eligibility(shape, information = info)
  run <- function(b) eligibility(b * shape, information = info, objective = "benefit",
                                 benefit = b * shape * 4 * info, n_screened = nrow(d))
  scale <- sapply(c(0.8, 0.9), function(target)
    uniroot(function(b) run(b)$power - target, c(1e-3, 2), extendInt = "upX")$root)
  benefit_close <- sapply(scale, function(b) max(d$delay[run(b)$included]))
  # Expected health gain of the 6-hour window relative to the best window, at 80% power.
  gain <- function(b, inc) design_power(b * shape, information = info, included = inc,
                                        n_screened = nrow(d))$power * sum((b * shape * 4 * info)[inc])
  gain_6h <- gain(scale[1], d$delay <= 6) / gain(scale[1], run(scale[1])$included)
  dur <- sapply(c(3, 4.5, 6), function(w) relative_duration(shape, information = info, included = d$delay <= w))
  data.frame(outcome = outcome, zero_effect_h = zero, power_window_h = max(d$delay[best$included]),
             enrolled_share = best$enrolled, benefit_window_80 = benefit_close[1], benefit_window_90 = benefit_close[2],
             duration_3h = dur[1], duration_4.5h = dur[2], duration_6h = dur[3], benefit_6h_80 = gain_6h,
             bound_6h = inclusion_bound(shape, information = info))
}
tabA <- do.call(rbind, lapply(c("good", "indep"), function(o)
  do.call(rbind, lapply(c(5.0, 6.3, 13.8), function(z) window_analysis(ist, o, z)))))
write.csv(tabA, "output/ist3_window_external.csv", row.names = FALSE)
cat("\nA. Time window with the external profile\n"); print(tabA, row.names = FALSE, digits = 3)

# ---- B. Time window with the profile estimated from IST-3 alone ------------------------------------------
internal_window <- function(d) {
  m <- glm(good ~ treat * delay + age + nihss + I(nihss^2), binomial, d)
  b <- coef(m)
  eff <- b["treat"] + b["treat:delay"] * d$delay
  out <- c(slope = unname(b["treat:delay"]), close = NA, all = 0)
  if (!any(eff > 0)) return(out)
  fit <- eligibility(eff, information = info_of(d, "good"))
  out["all"] <- as.numeric(fit$enrolled > 0.9999)
  if (b["treat:delay"] < 0) out["close"] <- max(d$delay[fit$included])
  out
}
est <- internal_window(ist)
boot <- t(replicate(2000, internal_window(ist[sample(nrow(ist), replace = TRUE), ])))
falls <- boot[, "slope"] < 0
tabB <- data.frame(interaction_per_hour = est["slope"], window_h = est["close"],
                   boot_share_effect_rises_with_delay = mean(!falls),
                   boot_share_enrol_all = mean(boot[, "all"] == 1),
                   boot_window_q10 = quantile(boot[falls, "close"], 0.10, na.rm = TRUE),
                   boot_window_q50 = quantile(boot[falls, "close"], 0.50, na.rm = TRUE),
                   boot_window_q90 = quantile(boot[falls, "close"], 0.90, na.rm = TRUE))
write.csv(tabB, "output/ist3_window_internal.csv", row.names = FALSE)
cat("\nB. Time window with the profile estimated from IST-3\n"); print(tabB, row.names = FALSE, digits = 3)

# ---- C. Upper age limit under a constant odds ratio ------------------------------------------------------
info <- info_of(ist, "good")
tabC <- data.frame(share_over_80 = mean(ist$age > 80),
                   mean_info_80_or_younger = mean(info[ist$age <= 80]), mean_info_over_80 = mean(info[ist$age > 80]),
                   duration_if_80_or_younger = relative_duration(1, information = info, included = ist$age <= 80))
write.csv(tabC, "output/ist3_age.csv", row.names = FALSE)
cat("\nC. Age limit\n"); print(tabC, row.names = FALSE, digits = 3)

# ---- Data for Figure 3 -------------------------------------------------------------------------------------
br <- seq(0, 6, by = 0.5)
write.csv(data.frame(delay_mid = head(br, -1) + 0.25, patients = as.vector(table(cut(ist$delay, br, include.lowest = TRUE)))),
          "output/ist3_arrival.csv", row.names = FALSE)
grid <- seq(1, 6, by = 0.05)
curve <- do.call(rbind, lapply(c(5.0, 6.3, 13.8), function(z)
  data.frame(zero_effect_h = z, window = grid,
             relative_duration = sapply(grid, function(w) relative_duration(z - ist$delay, information = info,
                                                                              included = ist$delay <= w)))))
write.csv(curve, "output/ist3_duration_curve.csv", row.names = FALSE)
