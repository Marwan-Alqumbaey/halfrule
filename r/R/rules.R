#' Threshold of the half rule and its relatives
#'
#' Gives the effect of the marginal patient type at the best eligibility set,
#' as a fraction of the average effect among enrolled patients. Types above the
#' threshold are enrolled and types below it are not. All types are assumed to
#' be equally noisy and small; for a whole group of patients use
#' [subgroup_threshold()].
#'
#' @details
#' * `"power"` (or `"cost"`): `1 - cost_share / 2`. With `cost_share = 1`
#'   (all variable cost scales with screening or time) this is one half.
#' * `"benefit"`: `e / (2 * (1 + e))`, where `e` is the elasticity of power at
#'   the planned power. About one quarter at 80% power.
#' * `"coverage"`: `1/2 - 1/e`. This is the threshold when the pay-off is the
#'   number of patients covered by a positive trial, whatever their benefit.
#'   It is negative at usual power, so patients with no effect are enrolled.
#'
#' @param objective One of `"power"`, `"cost"`, `"benefit"` or `"coverage"`.
#' @param power Planned power (benefit and coverage only).
#' @param alpha One-sided significance level.
#' @param cost_share Share of variable trial cost that scales with the number
#'   of screened patients or with time.
#' @return The threshold as a fraction of the average effect among enrolled patients.
#' @examples
#' rule_threshold("power")
#' rule_threshold("benefit", power = 0.8)
#' rule_threshold("cost", cost_share = 0.5)
#' @export
rule_threshold <- function(objective = c("power", "cost", "benefit", "coverage"),
                           power = 0.8, alpha = 0.025, cost_share = 1) {
  objective <- match.arg(objective)
  if (objective %in% c("power", "cost")) {
    if (any(cost_share < 0 | cost_share > 1)) stop("`cost_share` must lie between 0 and 1.", call. = FALSE)
    return(1 - cost_share / 2)
  }
  e <- power_elasticity(power, alpha)
  if (objective == "benefit") e / (2 * (1 + e)) else 0.5 - 1 / e
}

#' Elasticity of power
#'
#' Percentage change in power for a one percent change in the non-centrality
#' parameter, for a one-sided z-test.
#'
#' @param power Power of the test.
#' @param alpha One-sided significance level.
#' @return The elasticity.
#' @examples
#' power_elasticity(0.8)
#' @export
power_elasticity <- function(power = 0.8, alpha = 0.025) {
  if (any(power <= 0 | power >= 1)) stop("`power` must lie strictly between 0 and 1.", call. = FALSE)
  zb <- stats::qnorm(power)
  (stats::qnorm(1 - alpha) + zb) * stats::dnorm(zb) / power
}

#' Threshold for adding or dropping a whole group
#'
#' Exact condition for whether a whole group of patients shortens the trial.
#' Let `rho` be the information the group carries, divided by the information
#' of the other enrolled patients (for equally noisy patients, the ratio of the
#' two group sizes).
#'
#' * `reference = "others"`: adding the group shortens the trial if and only if
#'   its effect-to-noise ratio, divided by the pooled ratio of the other
#'   enrolled patients, exceeds `1 / (1 + sqrt(1 + rho))`. This is never more
#'   than one half, so a group whose effect is at least half the average effect
#'   of the others always shortens the trial.
#' * `reference = "all"`: keeping an enrolled group lengthens the trial if and
#'   only if its ratio, divided by the pooled ratio of all enrolled patients
#'   including the group, is below `sqrt(1 + rho) / (1 + sqrt(1 + rho))`. This is
#'   never less than one half, so a group whose effect is below half the average
#'   effect of all enrolled patients should always be dropped.
#'
#' Both values tend to one half as the group becomes small.
#'
#' @param rho Information of the group divided by the information of the other
#'   enrolled patients.
#' @param reference `"others"` or `"all"`, the patients whose average the
#'   group is compared with.
#' @return The threshold.
#' @examples
#' # Marker-positive patients are 25% of those screened. Marker-negative
#' # patients are worth enrolling if their effect exceeds one third of the
#' # effect in marker-positive patients.
#' subgroup_threshold(0.75 / 0.25)
#' subgroup_threshold(0.75 / 0.25, reference = "all")
#' @export
subgroup_threshold <- function(rho, reference = c("others", "all")) {
  reference <- match.arg(reference)
  if (any(rho < 0)) stop("`rho` must be non-negative.", call. = FALSE)
  s <- sqrt(1 + rho)
  if (reference == "others") 1 / (1 + s) else s / (1 + s)
}

#' Optimal time window
#'
#' Finds the best treatment time window when the effect depends on the delay
#' from onset to treatment.
#'
#' @param effect Function giving the expected effect at a delay.
#' @param density Function giving the relative number of patients who present
#'   at a delay. Defaults to a constant.
#' @param upper Longest delay considered.
#' @param lower Shortest delay considered.
#' @param variance Function giving the outcome variance at a delay, or a number.
#' @param n_grid Number of grid points.
#' @param ... Further arguments passed to [eligibility()], for example
#'   `objective`, `screen_cost`, `patient_cost`, `n_screened`.
#' @return A list with the delays on the grid (`delay`), the fit from
#'   [eligibility()] (`fit`), and the latest enrolled delay (`close`).
#' @examples
#' # Linear decay to zero at 6.3 hours, patients spread evenly over time.
#' time_window(function(t) 1 - t / 6.3, upper = 8)$close   # 4.2 hours
#' @export
time_window <- function(effect, density = function(t) rep(1, length(t)), upper,
                        lower = 0, variance = 1, n_grid = 4000, ...) {
  if (!is.function(effect) || !is.function(density)) stop("`effect` and `density` must be functions.", call. = FALSE)
  step <- (upper - lower) / n_grid
  t <- lower + (seq_len(n_grid) - 0.5) * step
  v <- if (is.function(variance)) variance(t) else variance
  fit <- eligibility(effect = effect(t), variance = v, weight = density(t), ...)
  list(delay = t, fit = fit, close = max(t[fit$included]) + step / 2)
}
