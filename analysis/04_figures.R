# Figures. Run after scripts 01 to 03. Figures are 170 mm wide.
suppressMessages(library(halfrule))
dir.create("figures", showWarnings = FALSE)
blue <- "#2a78d6"; orange <- "#eb6834"; aqua <- "#1baf7a"; ink <- "#0b0b0b"; muted <- "#52514e"; grid_col <- "#e4e3df"
width <- 170 / 25.4

draw <- function(name, height, plot_fun) {
  cairo_pdf(file.path("figures", paste0(name, ".pdf")), width = width, height = height, pointsize = 9, family = "Helvetica")
  plot_fun(); invisible(dev.off())
  png(file.path("figures", paste0(name, ".png")), width = width, height = height, units = "in", res = 600, pointsize = 9)
  plot_fun(); invisible(dev.off())
}
base_par <- function(...) par(mar = c(4, 4.3, 2, 0.8), mgp = c(2.5, 0.6, 0), tcl = -0.3, las = 1,
                              col.axis = muted, col.lab = ink, fg = muted, bty = "l", ...)
frame_grid <- function() grid(col = grid_col, lty = 1, lwd = 0.6)
letter <- function(l) mtext(l, side = 3, adj = -0.2, line = 0.6, font = 2, col = ink, cex = 1.1)

# ---- Figure 1: thresholds for whole groups, and time windows ---------------------------------------------
fig1 <- function() {
  layout(matrix(1:2, 1)); base_par()
  # A. Threshold for adding or dropping a group, by its size.
  rho <- exp(seq(log(0.01), log(100), length.out = 400))
  s <- sqrt(1 + rho)
  plot(NA, xlim = c(0.01, 100), ylim = c(0, 1), log = "x", xaxt = "n", xaxs = "i", yaxs = "i",
       xlab = "Size of the group relative to those already enrolled",
       ylab = "Effect of the group (fraction of the average)")
  axis(1, at = c(0.01, 0.1, 1, 10, 100), labels = c("0.01", "0.1", "1", "10", "100"))
  frame_grid(); abline(h = 0.5, col = muted, lty = 3)
  lines(rho, 1 / (1 + s), col = blue, lwd = 2)
  lines(rho, s / (1 + s), col = orange, lwd = 2, lty = 2)
  text(0.013, 0.30, "Add a group if its effect is above\nthis fraction of the average effect\nof the patients already enrolled",
       adj = c(0, 1), col = ink, cex = 0.85)
  text(0.013, 0.97, "Drop a group if its effect is below\nthis fraction of the average effect\nof all enrolled patients,\nthe group included",
       adj = c(0, 1), col = ink, cex = 0.85)
  letter("A")
  # B. Duration against window when the effect falls linearly to zero at time 1.
  w <- seq(0.15, 1.2, by = 0.005); t <- (seq_len(6000) - 0.5) * 1.2 / 6000
  patterns <- list("even arrival" = function(t) rep(1, length(t)), "more late arrivals" = function(t) t,
                   "more early arrivals" = function(t) 1.2 - t)
  cols <- c(blue, orange, aqua); ltys <- c(1, 2, 4)
  plot(NA, xlim = c(0.15, 1.2), ylim = c(1, 2), xaxs = "i", yaxs = "i",
       xlab = "Window closes at (fraction of the zero-effect time)", ylab = "Duration relative to the best window")
  frame_grid(); abline(v = 1, col = muted, lty = 3)
  for (k in seq_along(patterns)) {
    dens <- patterns[[k]](t)
    dur <- sapply(w, function(x) relative_duration(1 - t, weight = dens, included = t <= x))
    lines(w, dur, col = cols[k], lwd = 2, lty = ltys[k])
    points(w[which.min(dur)], 1, pch = 21, bg = cols[k], col = "white", cex = 1.4, lwd = 1.5, xpd = NA)
  }
  legend(0.39, 2.0, names(patterns), col = cols, lty = ltys, lwd = 2, bty = "n", text.col = ink, seg.len = 3)
  letter("B")
}
draw("figure1", 3.3, fig1)

# ---- Figure 2: simulated and predicted power -------------------------------------------------------------
sim <- read.csv("output/simulation_curves.csv")
fig2 <- function() {
  layout(matrix(1:4, 2, byrow = TRUE)); base_par()
  panel <- function(d, title, xlab, l, extra = NULL) {
    plot(NA, xlim = range(d$set), ylim = c(0, 1), xlab = xlab, ylab = "Power", yaxs = "i")
    frame_grid()
    if (is.finite(d$rule_power[1])) abline(v = d$rule_power[1], col = muted, lty = 3)
    lines(d$set, d$predicted, col = blue, lwd = 2)
    points(d$set, d$simulated, pch = 16, cex = 0.6, col = ink)
    if (!is.null(extra)) {
      lines(extra$set, extra$predicted, col = orange, lwd = 2, lty = 2)
      points(extra$set, extra$simulated, pch = 1, cex = 0.7, col = ink)
    }
    mtext(title, side = 3, line = 0.5, col = ink); letter(l)
  }
  panel(sim[sim$scenario == "A: even arrival", ], "Continuous outcome, even arrival", "Window closes at (hours)", "A")
  legend("bottomright", c("predicted", "simulated", "rule"), col = c(blue, ink, muted), lty = c(1, NA, 3),
         pch = c(NA, 16, NA), lwd = c(2, NA, 1), bty = "n", text.col = ink)
  panel(sim[sim$scenario == "B: later arrival more common", ], "Continuous outcome, later arrival more common",
        "Window closes at (hours)", "B")
  panel(sim[sim$scenario == "C: binary outcome", ], "Binary outcome, unadjusted analysis", "Window closes at (hours)", "C")
  panel(sim[sim$scenario == "D: unadjusted analysis", ], "Constant effect, strong prognostic covariate",
        "Enrol patients within this many SD of the mean", "D", extra = sim[sim$scenario == "D: adjusted analysis", ])
  text(2.95, 0.90, "adjusted analysis", adj = 1, col = ink); text(2.95, 0.48, "unadjusted analysis", adj = 1, col = ink)
}
draw("figure2", 5.6, fig2)

# ---- Figure 3: IST-3 ---------------------------------------------------------------------------------------
if (file.exists("output/ist3_arrival.csv")) {
  arr <- read.csv("output/ist3_arrival.csv"); cur <- read.csv("output/ist3_duration_curve.csv")
  fig3 <- function() {
    layout(matrix(1:2, 1)); base_par()
    plot(NA, xlim = c(0, 6), ylim = c(0, 500), xlab = "Delay from stroke to randomisation (hours)",
         ylab = "Patients", xaxs = "i", yaxs = "i")
    frame_grid()
    rect(arr$delay_mid - 0.23, 0, arr$delay_mid + 0.23, arr$patients, col = blue, border = NA)
    letter("A")
    plot(NA, xlim = c(1, 6), ylim = c(1, 2.5), xlab = "Window closes at (hours)",
         ylab = "Duration relative to the best window", xaxs = "i", yaxs = "i")
    frame_grid(); abline(v = c(3, 4.5), col = muted, lty = 3)
    zs <- c(5, 6.3, 13.8); cols <- c(orange, blue, aqua); ltys <- c(2, 1, 4)
    for (k in 1:3) {
      d <- cur[cur$zero_effect_h == zs[k], ]
      lines(d$window, pmin(d$relative_duration, 9), col = cols[k], lwd = 2, lty = ltys[k])
    }
    legend(4.0, 2.5, c("5.0 h", "6.3 h", "13.8 h"), title = "Benefit ends at", title.adj = 0, col = cols,
           lty = ltys, lwd = 2, bty = "n", text.col = ink, title.col = ink, seg.len = 3)
    letter("B")
  }
  draw("figure3", 3.3, fig3)
}

# Record the software versions used by the analysis scripts.
for (p in c("haven", "BioPET")) if (requireNamespace(p, quietly = TRUE)) loadNamespace(p)
writeLines(capture.output(sessionInfo()), "output/session_info.txt")
