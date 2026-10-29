## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")
spdbl_version <- tryCatch(utils::packageVersion("spDBL"),
                          error = function(e) NULL)
if (is.null(spdbl_version) || spdbl_version < "1.1.1") {
  stop(
    "This vignette requires spDBL >= 1.1.1; the active installation is ",
    if (is.null(spdbl_version)) "missing" else as.character(spdbl_version), ". ",
    "From the spDBL source directory, run ",
    "'bash vignettes/renderEmulation.sh quarto calibration', ",
    "or install the updated package and restart R.",
    call. = FALSE
  )
}

## -----------------------------------------------------------------------------
library(spDBL)
data("dt_calibration", package = "spDBL")
d <- dt_calibration
dim(d$training$input)
dim(d$training$responses[[1]])
range(d$time)

## -----------------------------------------------------------------------------
n_emulator <- 500L
em <- emulator_learn(
  d$training$input, d$training$responses, d$Nx, d$Ny,
  covariance = "sigma2R", episode_window = 2,
  nsam = n_emulator, seed = 1234
)
cal <- calibration_prepare(em, d$lower, d$upper)
stopifnot(identical(cal$sampling_method, "marginal"),
          isFALSE(cal$posterior$joint_trajectories))
c(times = cal$T, locations = cal$S, episodes = cal$E, draws = cal$M)

## ----eval=FALSE---------------------------------------------------------------
# cal_more <- calibration_prepare(em, d$lower, d$upper,
#                                 n_draws = 1000, seed = 1235)

## -----------------------------------------------------------------------------
predicted <- emulator_predict(cal, d$validation$input)
dim(predicted[[1]])
validation <- calibration_validate(
  cal, d$validation$input, d$validation$responses
)
validation$field

## -----------------------------------------------------------------------------
sites <- seq_len(cal$S)
z <- d$field$z[, sites, drop = FALSE]
dim(z)
range(diag(d$field$U))

## -----------------------------------------------------------------------------
recorded <- readRDS(system.file(
  "extdata", "calibration-example.rds", package = "spDBL", mustWork = TRUE
))
samples <- recorded$samples
stopifnot(recorded$provenance$profile == "oneiter500",
          identical(recorded$locations, sites),
          identical(recorded$provenance$emulator_sampling_method, "marginal"),
          isFALSE(recorded$provenance$joint_trajectories))
data.frame(
  emulator_draws = length(unique(samples$emulator_draw)),
  chains_per_draw = recorded$settings$chains,
  iterations_per_chain = recorded$settings$iter,
  warmup = recorded$settings$warmup,
  retained_samples = nrow(samples$mcmc_eta),
  observed_sites = length(recorded$locations),
  emulator_sampling = recorded$provenance$emulator_sampling_method
)

## ----eta-histograms, fig.width=11, fig.height=4.2, message=FALSE, warning=FALSE----
library(ggplot2)
active <- which(d$upper > d$lower)
eta_quantiles <- apply(samples$mcmc_eta, 2, quantile, c(0.025, 0.5, 0.975))
hist_colors <- c("#D9A450", "#78B0A1", "#8CA4BF")
line_colors <- c("True value" = "#C44D42", "Median" = "#26384A",
                 "95% interval" = "#637381")
line_types <- c("True value" = "solid", "Median" = "solid",
                "95% interval" = "dashed")
eta_plots <- lapply(seq_along(active), function(k) {
  j <- active[k]
  values <- samples$mcmc_eta[, j]
  q <- eta_quantiles[, j]
  span <- range(c(values, d$truth$eta[j]))
  padding <- max(diff(span) * 0.1, (d$upper[j] - d$lower[j]) * 0.005)
  limits <- c(max(d$lower[j], span[1] - padding),
              min(d$upper[j], span[2] + padding))
  references <- data.frame(
    value = c(q[1], q[3], q[2], d$truth$eta[j]),
    reference = factor(c("95% interval", "95% interval", "Median", "True value"),
                       levels = names(line_colors))
  )
  ggplot(data.frame(value = values), aes(x = value)) +
    geom_histogram(bins = 28, fill = hist_colors[k],
                   color = "white", linewidth = 0.25) +
    geom_vline(data = references,
               aes(xintercept = value, color = reference, linetype = reference),
               linewidth = 0.75, show.legend = TRUE) +
    scale_color_manual(values = line_colors, name = NULL, drop = FALSE) +
    scale_linetype_manual(values = line_types, name = NULL, drop = FALSE) +
    scale_x_continuous(breaks = scales::breaks_extended(n = 5),
                       labels = scales::label_number()) +
    scale_y_continuous(breaks = scales::breaks_pretty(n = 4),
                       expand = expansion(mult = c(0, 0.1))) +
    coord_cartesian(xlim = limits) +
    labs(title = bquote(eta[.(j)]),
         subtitle = sprintf("Median %.3f  |  95%% [%.3f, %.3f]", q[2], q[1], q[3]),
         x = "Parameter value", y = if (k == 1) "Count" else NULL) +
    theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(size = 16, color = "#26384A"),
      plot.subtitle = element_text(size = 9, color = "#5D6975",
                                   margin = margin(b = 10)),
      axis.title = element_text(color = "#44515D"),
      axis.text = element_text(color = "#44515D"),
      axis.line = element_line(color = "#AEB8C0", linewidth = 0.3),
      axis.ticks = element_line(color = "#AEB8C0", linewidth = 0.3),
      panel.grid.major.y = element_line(color = "#E8ECEF", linewidth = 0.35),
      legend.position = "bottom",
      legend.text = element_text(size = 10, color = "#44515D"),
      legend.key.width = grid::unit(1.4, "lines"),
      plot.margin = margin(8, 10, 8, 8)
    )
})
eta_panel <- ggpubr::ggarrange(
  plotlist = eta_plots, ncol = 3, align = "hv",
  common.legend = TRUE, legend = "bottom"
)
ggpubr::annotate_figure(
  eta_panel,
  top = ggpubr::text_grob("Mechanistic parameter samples",
                         face = "bold", size = 14, color = "#26384A"),
  bottom = ggpubr::text_grob(
    sprintf("%d emulator draws  |  %d transition per draw  |  %d warmup",
            length(unique(samples$emulator_draw)),
            recorded$settings$iter, recorded$settings$warmup),
    size = 9.5, color = "#5D6975"
  )
)

## -----------------------------------------------------------------------------
eta_table <- data.frame(
  parameter = colnames(samples$mcmc_eta)[active],
  truth = unname(d$truth$eta[active]),
  median = eta_quantiles[2, active],
  lower = eta_quantiles[1, active],
  upper = eta_quantiles[3, active],
  row.names = NULL
)
knitr::kable(eta_table, digits = 5)

## ----tau2-u-z-panel, fig.width=12, fig.height=4-------------------------------
tau_q <- apply(samples$mcmc_tau2, 2, quantile, c(0.025, 0.5, 0.975))
tau_frame <- data.frame(
  time = d$time, truth = d$truth$tau2,
  lower = tau_q[1, ], median = tau_q[2, ], upper = tau_q[3, ]
)
plot_tau2 <- ggplot(tau_frame, aes(x = time)) +
  geom_ribbon(aes(ymin = lower, ymax = upper),
              fill = "#A6CEE3", alpha = 0.6) +
  geom_line(aes(y = median), color = "#1F78B4",
            linewidth = 0.7, linetype = "dotdash") +
  geom_line(aes(y = truth), color = "#E31A1C", linewidth = 0.7) +
  labs(title = expression(Variance~(tau^2)), x = "Time", y = "Value") +
  theme_bw(base_size = 11)

comparison_panel <- function(q, reference, title, x_label) {
  frame <- data.frame(
    reference = as.vector(reference),
    lower = as.vector(q[1, , ]), median = as.vector(q[2, , ]),
    upper = as.vector(q[3, , ])
  )
  frame <- frame[order(frame$reference), ]
  ggplot(frame, aes(x = reference)) +
    geom_ribbon(aes(ymin = lower, ymax = upper),
                fill = "#A6CEE3", alpha = 0.6) +
    geom_point(aes(y = median), color = "#1F78B4", size = 0.6) +
    geom_abline(slope = 1, intercept = 0, color = "#E31A1C",
                linewidth = 0.7) +
    labs(title = title, x = x_label, y = "Recovered value") +
    theme_bw(base_size = 11)
}
plot_u <- comparison_panel(
  recorded$field_quantiles$u, d$truth$u[, sites, drop = FALSE],
  expression(Spatial~bias~(u)), "True spatial bias"
)
plot_z <- comparison_panel(
  recorded$field_quantiles$z_rep, z,
  expression(Field~value~(z)), "Observed field value"
)
ggpubr::ggarrange(plot_tau2, plot_u, plot_z, ncol = 3, align = "hv")

## -----------------------------------------------------------------------------
data.frame(
  divergences = sum(recorded$diagnostics$divergences),
  treedepth_hits = sum(recorded$diagnostics$max_treedepth_hits)
)

## -----------------------------------------------------------------------------
system.file("examples", "run_calibration.R", package = "spDBL")

