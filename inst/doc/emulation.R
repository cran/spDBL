## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
require(spDBL)
require(magrittr)
require(ggplot2)
require(ggpubr)
seed <- 1234
set.seed(seed)

# Load example data
data("dt_emulation")
pde_para_train = dt_emulation$pde_para_train
pde_para_test  = dt_emulation$pde_para_test
dt_pde_train   = dt_emulation$dt_pde_train
dt_pde_test    = dt_emulation$dt_pde_test

Nx <- 10
Ny <- 10

## -----------------------------------------------------------------------------
## Fit the emulator and predict at held-out inputs ----
emulator <- emulator_learn(pde_para_train = pde_para_train,
                         # pde_para_test = pde_para_test,
                         dt_pde_train = dt_pde_train,
                         # dt_pde_test = dt_pde_test, 
                         Nx = Nx,
                         Ny = Ny)

res_pre_exact <- emulator_predict(emulator = emulator,
                 input_new = pde_para_test,
                 dt_pde_test = dt_pde_test)

res_pre_MC <- emulator_predict(emulator = emulator,
                 input_new = pde_para_test,
                 dt_pde_test = dt_pde_test,
                 MC = TRUE)

## ----plot 1 prepare-----------------------------------------------------------
para_ffbs <- emulator$para_ffbs
res_ffbs <- emulator$res_ffbs
nT_ori <- emulator$setup$nT_ori
N_sp <- emulator$setup$N_sp
N_people <- emulator$setup$N_people
# Plot settings ----
# Set ggplot theme
col_epa <- c("#00e400", "#ffff00", "#ff7e00", "#ff0000", "#99004c", "#7e0023")
col_bgr <- c("#d5edfc", "#a5d9f6", "#7eb4e0", "#588dc8", "#579f8b", "#5bb349",
             "#5bb349", "#f3e35a", "#eda742", "#e36726", "#d64729", "#c52429",
             "#a62021", "#871b1c")

## Plot PDE results ----
### Heat map ----

input_num <- 1
tstamp <- as.integer(seq(1, nT_ori, length.out = 9))
dat <- dt_pde_test
max_y <- max(as.vector(unlist(dat))) # set max limit for all plots

{
  plot_ls <- list()
  ind_sp <- data.frame(row = rep(1:Ny, times = Nx), col = rep(1:Nx, each = Ny))
  ind_plot <- 1

  for (i in tstamp) {
    temp <- dat[[i]][input_num,]
    rownames(temp) <- NULL
    colnames(temp) <- NULL
    dt <- data.frame(row = ind_sp$row, col = ind_sp$col, sol = temp)%>%
      as.data.frame()

    p <- ggplot(dt, aes(x = col, y = row, fill = sol)) +
      geom_raster() +
      scale_fill_gradientn(colours = col_bgr,
                           limits = c(0, max_y),
                           oob = scales::squish) +
      labs(x = "x", y = "y", fill = "Value") 

    plot_ls[[ind_plot]] <- p
    ind_plot <- ind_plot + 1
  }
}

pde_heat <- plot_ls

## ----plot 1, fig.width=12, fig.height=10, out.width="100%", warning=FALSE, message=FALSE----
labels <- paste0("PDE: t = ", tstamp[1:9] - 1)

ggarrange(
  plotlist = pde_heat[1:9],
  ncol = 3,
  nrow = 3,
  labels = labels,
  font.label = list(size = 14, face = "bold"),
  hjust = -0.1,
  vjust = 1.2,
  align = "hv",
  common.legend = TRUE,
  legend = "right"
)

## ----plot 2 prepare-----------------------------------------------------------
dat <- res_pre_exact 
{
  plot_ls <- list()
  ind_sp <- data.frame(row = rep(1:Ny, times = Nx), col = rep(1:Nx, each = Ny))
  ind_plot <- 1

  for (i in tstamp) {
    temp <- dat[[i]][input_num,]
    rownames(temp) <- NULL
    colnames(temp) <- NULL
    dt <- data.frame(row = ind_sp$row, col = ind_sp$col, sol = temp)%>%
      as.data.frame()

    p <- ggplot(dt, aes(x = col, y = row, fill = sol)) +
      geom_raster() +
      scale_fill_gradientn(colours = col_bgr,
                           limits = c(0, max_y),
                           oob = scales::squish) +
      labs(x = "x", y = "y", fill = "Value")

    plot_ls[[ind_plot]] <- p
    ind_plot <- ind_plot + 1
  }
}

ffbs_heat <- plot_ls

## ----plot 2, fig.width=12, fig.height=7, out.width="100%", warning=FALSE, message=FALSE----
labels <- c(
  paste0("PDE solution: t = ", tstamp[c(4, 6, 8)] - 1),
  paste0("FFBS emulation: t = ", tstamp[c(4, 6, 8)] - 1)
)

ggarrange(
  plotlist = list(
    pde_heat[[4]], pde_heat[[6]], pde_heat[[8]],
    ffbs_heat[[4]], ffbs_heat[[6]], ffbs_heat[[8]]
  ),
  ncol = 3,
  nrow = 2,
  labels = labels,
  font.label = list(size = 16, face = "bold"),
  hjust = -0.1,
  vjust = 1.2,
  align = "hv",
  common.legend = TRUE,
  legend = "right"
)

## ----plot 3 prepare-----------------------------------------------------------
{
  ### Error plot ----
  #### Specific time, one spatial location, all inputs ----
  {
    alpha_error <- 0.5
    res_pre <- res_pre_MC
    
    time_num <- 15
    sp_num <- 2
    y_true <- dt_pde_test[[time_num]][,sp_num]
    y_pre <- res_pre[[time_num]][,sp_num,]
    y_pre_stat <- data.frame(y_true = y_true,
                             med = apply(X = y_pre, MARGIN = 1, FUN = median),
                             lower = apply(X = y_pre, MARGIN = 1, FUN = quantile, prob = 0.025),
                             upper = apply(X = y_pre, MARGIN = 1, FUN = quantile, prob = 0.975))
    error_width <- (max(y_true) - min(y_true)) / 30 # denominator is for aesthetic 
    y_pre_stat %>% ggplot(aes(x = y_true, y = med)) + 
      geom_pointrange(aes(ymin = lower, ymax = upper), size =.2)+
      geom_errorbar(aes(ymin = lower, ymax = upper), width = error_width) + 
      geom_abline(col = "red")
  }
  
  #### Specific time, one input, all spatial locations ----
  {
    time_num <- 20
    sp_num <- c(1:N_sp)
    y_true <- dt_pde_test[[time_num]][input_num,sp_num]
    y_pre <- res_pre[[time_num]][input_num,sp_num,]
    y_pre_stat <- data.frame(y_true = y_true,
                             med = apply(X = y_pre, MARGIN = 1, FUN = median),
                             lower = apply(X = y_pre, MARGIN = 1, FUN = quantile, prob = 0.025),
                             upper = apply(X = y_pre, MARGIN = 1, FUN = quantile, prob = 0.975))
    y_pre_stat <- y_pre_stat / N_people
    error_width <- (max(y_pre_stat["y_true"]) - min(y_pre_stat["y_true"])) / 30
    plot_error_1 <- y_pre_stat %>% ggplot(aes(x = y_true, y = med)) + 
      geom_pointrange(aes(ymin = lower, ymax = upper), size =.2, alpha = alpha_error)+
      geom_errorbar(aes(ymin = lower, ymax = upper), width = error_width, alpha = alpha_error) + 
      geom_abline(col = "red") + 
      labs(x = "PDE solution", y = "FFBS prediction")
  }
  
  # panel
  time_p <- tstamp[c(4, 6, 8)]
  sp_num <- c(1:N_sp)
  plot_error_comp_ls <- list()
  plot_band_ls <- list()
  y_pre_stat_error_comp_ls <- list()
  for (t in 1:length(time_p)) {
    time_num <- time_p[t]
    y_true <- dt_pde_test[[time_num]][input_num,sp_num]
    y_pre <- res_pre[[time_num]][input_num,sp_num,]
    y_pre_stat <- data.frame(y_true = y_true,
                             med = apply(X = y_pre, MARGIN = 1, FUN = median),
                             lower = apply(X = y_pre, MARGIN = 1, FUN = quantile, prob = 0.025),
                             upper = apply(X = y_pre, MARGIN = 1, FUN = quantile, prob = 0.975))
    y_pre_stat <- y_pre_stat / N_people
    y_pre_stat_error_comp_ls[[t]] <- y_pre_stat
  }
  
  for (t in 1:length(time_p)) {
    # scatter plot
    error_width <- (max(y_pre_stat_error_comp_ls[[t]]["y_true"]) - min(y_pre_stat_error_comp_ls[[t]]["y_true"])) / 30
    plot_error_1 <- y_pre_stat_error_comp_ls[[t]] %>% ggplot(aes(x = y_true, y = med)) + 
      geom_pointrange(aes(ymin = lower, ymax = upper), size =.2, alpha = alpha_error / 3)+
      geom_errorbar(aes(ymin = lower, ymax = upper), width = error_width, alpha = alpha_error / 3) + 
      geom_abline(col = "red") + 
      labs(x = "PDE solution", y = "FFBS emulation")
    plot_error_comp_ls[[t]] <- plot_error_1
    
    # error band plot
    plot_band_predict <- y_pre_stat_error_comp_ls[[t]] %>%
      ggplot(aes(x = y_true, y = med)) + 
      geom_ribbon(aes(ymin = lower, ymax = upper, x = y_true), fill = "#A6CEE3", alpha = 1) +
      geom_point(aes(y = med), color = "#1F78B4", size = 0.7, alpha = 1) +
      geom_abline(color = "#E31A1C", linewidth = 1) + 
      labs(x = "PDE solution", y = "FFBS emulation") +
      theme_minimal()
    plot_band_ls[[t]] <- plot_band_predict
  }
}

## ----plot 3, fig.width=13, fig.height=4.5, out.width="100%", warning=FALSE, message=FALSE----
labels <- paste0("t = ", tstamp[c(4, 6, 8)] - 1)

ggarrange(
  plotlist = plot_error_comp_ls[1:3],
  ncol = 3,
  nrow = 1,
  labels = labels,
  font.label = list(size = 16, face = "bold"),
  hjust = -0.1,
  vjust = 1.2,
  align = "hv",
  common.legend = TRUE,
  legend = "right"
)

