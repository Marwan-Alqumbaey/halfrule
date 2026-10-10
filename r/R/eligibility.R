#' Optimal eligibility criteria
#'
#' Finds which patient types a two-arm randomised trial should enrol.
#' Enrolling more types recruits faster but mixes in patients with smaller
#' effects. The function returns the set that is best for the chosen objective.
#'
#' @details
#' Each element of `effect` describes one patient type (or one patient).
#' Three objectives are available.
#'
#' * `"power"`: highest power for a fixed number of screened patients, which is
#'   the same as the shortest trial for a fixed power. The best set consists of
#'   the types with the largest effect-to-noise ratios. In it, every enrolled
#'   type has a ratio of at least half the pooled ratio of the enrolled patients
#'   and every excluded type at most half (the half rule). The function ranks
#'   the types and compares all candidate sets, so it returns the best set even
#'   when several sets meet the half rule.
#' * `"cost"`: lowest cost for a fixed power when each screened patient costs
#'   `screen_cost` (this includes the cost of time) and each enrolled patient
#'   costs `patient_cost`.
#' * `"benefit"`: highest expected benefit, defined as the probability of a
#'   positive trial times the total `benefit` of the enrolled types. With
#'   `benefit = 1` this is the expected number of patients covered by a
#'   positive result.
#'
#' By default the outcome is compared between arms as a difference in means
#' with 1:1 allocation and `variance` is the outcome variance within an arm.
#' If `information` is given, `effect` is read as a model coefficient (for
#' example a log odds ratio or log hazard ratio) and `information` as the
#' Fisher information that one enrolled patient carries about it.
#'
#' If `prognosis` is given, the analysis is taken to be unadjusted. Types whose
#' expected outcome is far from the average of the enrolled patients then add
#' noise and need a larger effect to be worth enrolling. The set is then found
#' by a search that is not guaranteed to be exact.
#'
#' @param effect Expected treatment effect of each patient type.
#' @param variance Outcome variance within an arm for each type.
#' @param weight Frequency of each type among screened patients.
#' @param information Optional Fisher information per enrolled patient (score form).
#' @param objective One of `"power"`, `"cost"` or `"benefit"`.
#' @param screen_cost Cost per screened patient, including the cost of time.
#' @param patient_cost Cost per enrolled patient.
#' @param benefit Health gain per future patient of each type. Defaults to `effect`.
#' @param n_screened Number of patients screened during the trial. Needed for `"benefit"`.
#' @param alpha One-sided significance level.
#' @param prognosis Optional expected outcome under control for each type.
#'
#' @return An object of class `halfrule`: a list with
#' \item{included}{logical vector, `TRUE` for enrolled types}
#' \item{ratio}{effect divided by noise (4 times the variance, or the Fisher information in score form) for each type}
#' \item{threshold}{lowest ratio among enrolled types divided by the pooled ratio}
#' \item{enrolled}{fraction of screened patients who are enrolled}
#' \item{mean_effect}{average effect among enrolled patients}
#' \item{ncp}{non-centrality parameter per square root of screened patients}
#' \item{relative_duration}{trial duration relative to the power-optimal set, for equal power}
#' \item{cost_share}{share of variable cost that scales with screening (cost objective)}
#' \item{power, elasticity}{power and its elasticity (benefit objective)}
#'
#' @examples
#' # Effect falls linearly with treatment delay and vanishes at 6 hours.
#' delay <- seq(0.005, 5.995, by = 0.01)
#' fit <- eligibility(effect = 1 - delay / 6)
#' max(delay[fit$included])   # 4 hours: two thirds of 6
#' @export
eligibility <- function(effect, variance = 1, weight = 1, information = NULL,
                        objective = c("power", "cost", "benefit"),
                        screen_cost = 1, patient_cost = 0, benefit = NULL,
                        n_screened = NULL, alpha = 0.025, prognosis = NULL) {
  objective <- match.arg(objective)
  chk <- .check_inputs(effect, variance, weight, information, prognosis)
  effect <- chk$effect
  w <- chk$weight
  if (!any(effect[w > 0] > 0)) stop("No patient type has a positive expected effect.", call. = FALSE)
  if (objective == "cost" && (screen_cost < 0 || patient_cost < 0 || screen_cost + patient_cost <= 0))
    stop("Costs must be non-negative and not both zero.", call. = FALSE)
  h <- NULL
  if (objective == "benefit") {
    if (is.null(n_screened) || n_screened <= 0) stop("`n_screened` is needed for the benefit objective.", call. = FALSE)
    h <- if (is.null(benefit)) effect else rep_len(benefit, length(effect))
  }
  z <- stats::qnorm(1 - alpha)
  profile <- !is.null(chk$prognosis) && is.null(chk$information)
  solve <- function(obj, included = NULL, centre = NULL) {
    sn <- .signal_noise(effect, chk$variance, chk$information, chk$prognosis, w, included, centre)
    .solve_fixed(sn$g, sn$v, w, h, obj, screen_cost, patient_cost, n_screened, z)
  }
  value_of <- function(included, obj) {
    if (!any(included)) return(-Inf)
    sn <- .signal_noise(effect, chk$variance, chk$information, chk$prognosis, w, included)
    .set_value(included, sn$g, sn$v, w, h, obj, screen_cost, patient_cost, n_screened, z)
  }
  find <- function(obj) {
    pos <- w > 0
    if (!profile) return(solve(obj, included = effect > 0 & pos))
    # Unadjusted analysis: the noise of a type depends on the mean outcome of
    # the enrolled patients. Search over that mean, then refine.
    prog <- chk$prognosis[pos]
    centres <- if (sum(pos) <= 100) {
      u <- sort(unique(prog)); sort(unique(c(u, (u[-1] + u[-length(u)]) / 2)))
    } else {
      stats::quantile(prog, seq(0.01, 0.99, length.out = 65), names = FALSE)
    }
    shifts <- unique(stats::quantile(effect[effect > 0 & pos], c(0.1, 0.5, 0.9), names = FALSE))
    starts <- list(effect > 0 & pos)
    for (m in centres) for (s in shifts) starts[[length(starts) + 1L]] <- solve(obj, centre = c(m, m + s))
    best <- NULL; best_val <- -Inf
    for (inc in unique(starts)) {
      val <- value_of(inc, obj)
      repeat {
        new <- solve(obj, included = inc)
        nv <- value_of(new, obj)
        if (!(nv > val + 1e-12)) break
        inc <- new; val <- nv
      }
      if (val > best_val) { best <- inc; best_val <- val }
    }
    # Local search: single flips, and swaps when the number of types is small.
    idx <- which(pos)
    repeat {
      improved <- FALSE
      for (k in idx) {
        trial <- best; trial[k] <- !trial[k]
        tv <- value_of(trial, obj)
        if (tv > best_val + 1e-12) { best <- trial; best_val <- tv; improved <- TRUE }
      }
      if (length(idx) <= 60) {
        for (a in idx[best[idx]]) for (b in idx[!best[idx]]) {
          if (!best[a] || best[b]) next
          trial <- best; trial[a] <- FALSE; trial[b] <- TRUE
          tv <- value_of(trial, obj)
          if (tv > best_val + 1e-12) { best <- trial; best_val <- tv; improved <- TRUE }
        }
      }
      if (!improved) break
    }
    best
  }
  inc <- find(objective)
  ref <- if (objective == "power") inc else find("power")
  sn <- .signal_noise(effect, chk$variance, chk$information, chk$prognosis, w, inc)
  t <- .totals(inc, sn$g, sn$v, w, h)
  snr <- .signal_noise(effect, chk$variance, chk$information, chk$prognosis, w, ref)
  tr <- .totals(ref, snr$g, snr$v, w)
  ratio <- sn$g / sn$v
  out <- list(included = inc, ratio = ratio,
              threshold = min(ratio[inc]) / (t$G / t$V),
              enrolled = t$F,
              mean_effect = sum((w * effect)[inc]) / t$F,
              ncp = t$G / sqrt(t$V),
              relative_duration = (tr$G^2 / tr$V) / (t$G^2 / t$V),
              objective = objective)
  if (objective == "cost") out$cost_share <- screen_cost / (screen_cost + patient_cost * t$F)
  if (objective == "benefit") {
    nu <- sqrt(n_screened) * out$ncp
    out$power <- stats::pnorm(nu - z)
    out$elasticity <- nu * stats::dnorm(nu - z) / out$power
  }
  class(out) <- "halfrule"
  out
}

#' @export
print.halfrule <- function(x, digits = 3, ...) {
  cat("Optimal eligibility (objective: ", x$objective, ")\n", sep = "")
  cat("  enrolled share of screened patients: ", format(x$enrolled, digits = digits), "\n", sep = "")
  cat("  mean effect among enrolled:          ", format(x$mean_effect, digits = digits), "\n", sep = "")
  cat("  lowest ratio / pooled ratio:         ", format(x$threshold, digits = digits), "\n", sep = "")
  cat("  duration relative to power-optimal:  ", format(x$relative_duration, digits = digits), "\n", sep = "")
  if (!is.null(x$cost_share)) cat("  share of cost that scales with screening: ", format(x$cost_share, digits = digits), "\n", sep = "")
  if (!is.null(x$power)) cat("  power: ", format(x$power, digits = digits), "\n", sep = "")
  invisible(x)
}

#' Duration of a trial with a given eligibility set
#'
#' Compares any eligibility set with the power-optimal set. The result is the
#' factor by which the trial must be longer to reach the same power.
#'
#' @inheritParams eligibility
#' @param included Logical vector marking the enrolled types. Defaults to all types.
#' @return Relative duration (1 means as fast as the optimal set).
#' @examples
#' delay <- seq(0.005, 5.995, by = 0.01)
#' relative_duration(effect = 1 - delay / 6)                   # enrol up to 6 hours
#' relative_duration(effect = 1 - delay / 6, included = delay <= 4.5)
#' @export
relative_duration <- function(effect, variance = 1, weight = 1, information = NULL,
                              prognosis = NULL, included = NULL) {
  chk <- .check_inputs(effect, variance, weight, information, prognosis, length(included))
  effect <- chk$effect
  w <- chk$weight
  if (is.null(included)) included <- rep(TRUE, length(effect))
  if (length(included) != length(effect)) stop("`included` must have one value per patient type.", call. = FALSE)
  included <- as.logical(included) & w > 0
  best <- eligibility(effect, chk$variance, w, chk$information, prognosis = chk$prognosis)
  sn <- .signal_noise(effect, chk$variance, chk$information, chk$prognosis, w, included)
  t <- .totals(included, sn$g, sn$v, w)
  if (t$G <= 0) return(Inf)
  best$ncp^2 / (t$G^2 / t$V)
}

#' Upper bound on the cost of enrolling everyone
#'
#' A trial that enrols every patient type needs at most `1 + CV^2` times the
#' duration of the best restricted trial, where `CV` is the coefficient of
#' variation of the effect-to-noise ratio.
#'
#' @inheritParams eligibility
#' @return The bound `1 + CV^2`.
#' @examples
#' delay <- seq(0.005, 5.995, by = 0.01)
#' inclusion_bound(effect = 1 - delay / 6)   # 4/3
#' @export
inclusion_bound <- function(effect, variance = 1, weight = 1, information = NULL) {
  chk <- .check_inputs(effect, variance, weight, information, NULL)
  sn <- .signal_noise(chk$effect, chk$variance, chk$information, NULL, chk$weight, TRUE)
  m <- chk$weight * sn$v
  r <- sn$g / sn$v
  m1 <- sum(m * r) / sum(m)
  if (m1 <= 0) return(Inf)
  sum(m * r^2) / sum(m) / m1^2
}

#' Power of a trial with a given eligibility set
#'
#' @inheritParams eligibility
#' @param included Logical vector marking the enrolled types. Defaults to all types.
#' @return A list with the enrolled share (`enrolled`), the mean effect among
#'   enrolled patients (`mean_effect`), the non-centrality parameter (`ncp`)
#'   and the power of a one-sided test (`power`).
#' @examples
#' delay <- seq(0.005, 5.995, by = 0.01)
#' design_power(effect = 0.3 * (1 - delay / 6), n_screened = 1500, included = delay <= 4)
#' @export
design_power <- function(effect, variance = 1, weight = 1, information = NULL,
                         prognosis = NULL, included = NULL, n_screened, alpha = 0.025) {
  chk <- .check_inputs(effect, variance, weight, information, prognosis, length(included))
  effect <- chk$effect
  w <- chk$weight
  if (is.null(included)) included <- rep(TRUE, length(effect))
  if (length(included) != length(effect)) stop("`included` must have one value per patient type.", call. = FALSE)
  included <- as.logical(included) & w > 0
  if (!any(included)) stop("No patient type is enrolled.", call. = FALSE)
  sn <- .signal_noise(effect, chk$variance, chk$information, chk$prognosis, w, included)
  t <- .totals(included, sn$g, sn$v, w)
  ncp <- sqrt(n_screened) * t$G / sqrt(t$V)
  list(enrolled = t$F, mean_effect = sum((w * effect)[included]) / t$F,
       ncp = ncp, power = stats::pnorm(ncp - stats::qnorm(1 - alpha)))
}
