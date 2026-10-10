# Internal helpers. All sums are over patient types weighted by their
# frequency among screened patients.

.check_inputs <- function(effect, variance, weight, information, prognosis, n_min = 1L) {
  n <- max(length(effect), length(variance), length(weight), length(information), length(prognosis), n_min)
  rec <- function(x, name) {
    if (is.null(x)) return(NULL)
    if (!is.numeric(x) || anyNA(x)) stop("`", name, "` must be numeric without missing values.", call. = FALSE)
    if (length(x) == 1L) x <- rep(x, n)
    if (length(x) != n) stop("`", name, "` must have length 1 or the common length ", n, ".", call. = FALSE)
    x
  }
  effect <- rec(effect, "effect")
  variance <- rec(variance, "variance")
  weight <- rec(weight, "weight")
  information <- rec(information, "information")
  prognosis <- rec(prognosis, "prognosis")
  if (any(variance <= 0)) stop("`variance` must be positive.", call. = FALSE)
  if (any(weight < 0) || sum(weight) <= 0) stop("`weight` must be non-negative with a positive sum.", call. = FALSE)
  if (!is.null(information) && any(information <= 0)) stop("`information` must be positive.", call. = FALSE)
  list(effect = effect, variance = variance, weight = weight / sum(weight),
       information = information, prognosis = prognosis)
}

# Signal and noise of each type. Difference in means with 1:1 allocation:
# signal = effect, noise = 4 * variance, so that the non-centrality parameter
# is sqrt(n_screened) * G / sqrt(V). Score form: signal = effect * information.
.signal_noise <- function(effect, variance, information, prognosis, weight, included, centre = NULL) {
  if (!is.null(information)) {
    return(list(g = effect * information, v = information))
  }
  v <- 4 * variance
  if (!is.null(prognosis) && (any(included) || !is.null(centre))) {
    mu1 <- prognosis + effect
    if (is.null(centre)) {
      f <- sum(weight[included])
      centre <- c(sum((weight * prognosis)[included]) / f, sum((weight * mu1)[included]) / f)
    }
    v <- 4 * (variance + 0.5 * (prognosis - centre[1])^2 + 0.5 * (mu1 - centre[2])^2)
  }
  list(g = effect, v = v)
}

.totals <- function(included, g, v, w, h = NULL) {
  list(G = sum((w * g)[included]), V = sum((w * v)[included]), F = sum(w[included]),
       H = if (is.null(h)) NA_real_ else sum((w * h)[included]))
}

# Value to maximise for a set with totals G, V, F, H.
.value <- function(G, V, F, H, objective, screen_cost, patient_cost, n_screened, z) {
  ok <- G > 0 & V > 0
  out <- rep(-Inf, length(G))
  if (objective == "power") {
    out[ok] <- log(G[ok]) - 0.5 * log(V[ok])
  } else if (objective == "cost") {
    out[ok] <- -(log(V[ok]) - 2 * log(G[ok]) + log(screen_cost + patient_cost * F[ok]))
  } else {
    ok <- ok & H > 0
    nu <- sqrt(n_screened) * G[ok] / sqrt(V[ok])
    out[ok] <- stats::pnorm(nu - z, log.p = TRUE) + log(H[ok])
  }
  out
}

# Best prefix when types are ranked by `index` (largest first).
.best_prefix <- function(index, g, v, w, h, objective, screen_cost, patient_cost, n_screened, z) {
  o <- order(index, decreasing = TRUE)
  G <- cumsum((w * g)[o]); V <- cumsum((w * v)[o]); F <- cumsum(w[o])
  H <- if (is.null(h)) rep(NA_real_, length(o)) else cumsum((w * h)[o])
  val <- .value(G, V, F, H, objective, screen_cost, patient_cost, n_screened, z)
  k <- which.max(val)
  inc <- logical(length(o)); inc[o[seq_len(k)]] <- TRUE
  list(included = inc, value = val[k])
}

# Marginal value of adding a little of each type to the current set.
.marginal <- function(included, g, v, w, h, objective, screen_cost, patient_cost, n_screened, z) {
  t <- .totals(included, g, v, w, h)
  if (objective == "power") return(g / t$G - v / (2 * t$V))
  if (objective == "cost") return(2 * g / t$G - v / t$V - patient_cost / (screen_cost + patient_cost * t$F))
  nu <- sqrt(n_screened) * t$G / sqrt(t$V)
  eps <- nu * exp(stats::dnorm(nu - z, log = TRUE) - stats::pnorm(nu - z, log.p = TRUE))
  eps * (g / t$G - v / (2 * t$V)) + h / t$H
}

.set_value <- function(included, g, v, w, h, objective, screen_cost, patient_cost, n_screened, z) {
  if (!any(included)) return(-Inf)
  t <- .totals(included, g, v, w, h)
  .value(t$G, t$V, t$F, t$H, objective, screen_cost, patient_cost, n_screened, z)
}

# Solve for fixed signal and noise.
.solve_fixed <- function(g, v, w, h, objective, screen_cost, patient_cost, n_screened, z) {
  args <- list(g = g, v = v, w = w, h = h, objective = objective, screen_cost = screen_cost,
               patient_cost = patient_cost, n_screened = n_screened, z = z)
  pos <- w > 0
  if (objective == "power" || (objective == "cost" && patient_cost == 0)) {
    idx <- ifelse(pos, g / v, -Inf)
    args$objective <- "power"
    return(do.call(.best_prefix, c(list(index = idx), args))$included & pos)
  }
  # Candidate sets from the two-parameter family implied by the first-order condition.
  cands <- list()
  if (objective == "cost") {
    for (a in seq(0, max(g), length.out = 60)) cands[[length(cands) + 1L]] <- ifelse(pos, (g - a) / v, -Inf)
  } else {
    for (th in seq(0, pi / 2, length.out = 60)) cands[[length(cands) + 1L]] <- ifelse(pos, (cos(th) * g + sin(th) * h) / v, -Inf)
  }
  best <- NULL
  for (idx in cands) {
    b <- do.call(.best_prefix, c(list(index = idx), args))
    if (is.null(best) || b$value > best$value) best <- b
  }
  inc <- best$included & pos; val <- best$value
  # Polish: apply the rule until the set stops changing, then try single flips.
  for (it in seq_len(200)) {
    new <- (do.call(.marginal, c(list(included = inc), args)) >= 0) & pos
    nv <- do.call(.set_value, c(list(included = new), args))
    if (identical(new, inc) || !(nv > val)) break
    inc <- new; val <- nv
  }
  repeat {
    t <- .totals(inc, g, v, w, h)
    s <- ifelse(inc, -1, 1) * w
    Hk <- if (is.null(h)) rep(NA_real_, length(g)) else t$H + s * h
    vals <- .value(t$G + s * g, t$V + s * v, t$F + s, Hk, objective, screen_cost,
                   patient_cost, n_screened, z)
    vals[!pos] <- -Inf
    k <- which.max(vals)
    if (!(vals[k] > val + 1e-12)) break
    inc[k] <- !inc[k]; val <- vals[k]
  }
  inc
}
