#Freshwater ecological responses to warming conditions and increased
#competition in an anadromous keystone species
#Christman, Schindler & Scheuerell 

#1. Load data          (single file: Christman_et_al_BiolLett_data.xlsx)
#2. August lake temperature (Results)
#3. Growth data cleaning
#4. Growth summary statistics + Welch's t-tests (Results)
#5. Figure 1  Environmental context
#6. Figure 2  Freshwater growth anomalies
#7. MARSS growth models, AICc, model averaging (Figure 3 and Appendix 2)


###### 0. data setup #####
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(scales)
library(cowplot)
library(ggplotify)
library(ggbreak)
library(readxl)
library(MARSS)

data_file <- "Christman_et_al_BiolLett_data.xlsx"

historical_end_year <- 1950        # historical period: FW1_Year <= 1950
modern_start_year   <- 1960        # modern period:     FW1_Year >= 1960
temp_years          <- 1964:2022   # years with lake-temperature record
outlier_sd          <- 2.5         # drop fish > 2.5 SD from River x FW1_Year x age mean
z90                 <- qnorm(0.95) # 90% CIs for MARSS covariate effects

river_names <- c(N = "Nushagak", W = "Wood")

dir.create("figures", showWarnings = FALSE)
dir.create("tables",  showWarnings = FALSE)
save_fig <- function(plot, name, width, height) {
  ggsave(file.path("figures", paste0(name, ".png")), plot,
         width = width, height = height, units = "in", dpi = 300, bg = "white")
}
save_tab <- function(df, name) {
  write.csv(df, file.path("tables", paste0(name, ".csv")), row.names = FALSE)
}

##### 1. Load in data #####

read_sheet <- function(sheet) read_excel(data_file, sheet = sheet, guess_max = 1e6)

growth_raw       <- read_sheet("growth")
lake_temp        <- read_sheet("lake_temp")
escapement_model <- read_sheet("escapement_model")
esc_plot         <- read_sheet("escapement_fig1") %>% rename(EscpYear = Year)
snap             <- read_sheet("snap_air_temp")
density          <- read_sheet("spawner_density")

# Escapement covariate: log10 spawners, lagged so that spawners in return
# year t are matched to juvenile lake growth (FW1) in year t + 1
log10_escapement <- function(watershed) {
  escapement_model %>%
    filter(Watershed == watershed) %>%
    transmute(FW1_Year = Return_Year + 1, log10_Esc = log10(Escapement)) %>%
    filter(FW1_Year <= max(temp_years))
}
wbrood <- log10_escapement("Wood")
nbrood <- log10_escapement("Nushagak")


##### 2. Lake Temperatures #####
# August = day of year 209-225, Annual mean of the sampled values

aug_obs <- tibble(FW1_Year = temp_years) %>%
  left_join(lake_temp %>%
              filter(DayOfYear >= 209, DayOfYear <= 225, YearSampled %in% temp_years) %>%
              group_by(FW1_Year = YearSampled) %>%
              summarise(Aug_Temp_obs = mean(AvgOfTemperature, na.rm = TRUE), .groups = "drop"),
            by = "FW1_Year")

# Univariate biased random walk state-space model:
#   x_t = x_{t-1} + u + w_t,  w_t ~ N(0, q)   (August temperature)
#   y_t = x_t + v_t,          v_t ~ N(0, r)   (sampled annual mean)
aug_brw_model <- list(B = matrix(1), U = matrix("u"), Q = matrix("q"),
                      Z = matrix(1), A = matrix(0),   R = matrix("r"),
                      x0 = matrix("x0"), tinitx = 1)
y_aug <- matrix(aug_obs$Aug_Temp_obs, nrow = 1)

fit_aug_brw <- MARSS(y_aug, model = aug_brw_model,
                     control = list(maxit = 5000, allow.degen = TRUE))

# q is estimated at 0, so the final model fixes q = 0 (identical likelihood)
fit_aug <- MARSS(y_aug, model = modifyList(aug_brw_model, list(Q = matrix(0))),
                 control = list(maxit = 5000))

aug_q_check <- tibble(model  = c("q estimated", "q fixed at 0"),
                      q      = c(coef(fit_aug_brw)$Q[1, 1], 0),
                      u      = c(coef(fit_aug_brw)$U[1, 1], coef(fit_aug)$U[1, 1]),
                      logLik = c(fit_aug_brw$logLik, fit_aug$logLik),
                      AICc   = c(fit_aug_brw$AICc, fit_aug$AICc))

# Trend u (*C/year) with 95% CI; change 1964-2022
ci_aug  <- MARSSparamCIs(fit_aug, method = "hessian", alpha = 0.05)
n_steps <- diff(range(temp_years))
aug_trend <- tibble(u_per_year   = ci_aug$par$U[1, 1],
                    u_lo95       = ci_aug$par.lowCI$U[1, 1],
                    u_hi95       = ci_aug$par.upCI$U[1, 1]) %>%
  mutate(per_decade   = u_per_year * 10,
         total_change = u_per_year * n_steps,
         total_lo95   = u_lo95 * n_steps,
         total_hi95   = u_hi95 * n_steps,
         n_years_obs  = sum(!is.na(aug_obs$Aug_Temp_obs)))
print(aug_trend)

# Kalman-smoothing fills missing years
aug_states <- tsSmooth(fit_aug, type = "xtT", interval = "confidence", level = 0.95) %>%
  as_tibble() %>%
  transmute(FW1_Year = temp_years[t], Aug_Temp_xT = .estimate,
            Aug_Temp_se = .se, lo95 = .conf.low, hi95 = .conf.up) %>%
  left_join(aug_obs, by = "FW1_Year")

aug_filled_years <- aug_states %>% filter(is.na(Aug_Temp_obs)) %>%
  dplyr::select(FW1_Year, Aug_Temp_xT, Aug_Temp_se, lo95, hi95)
print(aug_filled_years)

# Covariate series: observed where sampled, smoothed state where not
lt_interp <- aug_states %>%
  transmute(FW1_Year,
            Aug_Temp          = coalesce(Aug_Temp_obs, Aug_Temp_xT),
            Aug_Temp_observed = !is.na(Aug_Temp_obs))


##### 3. Freshwater growth data #####
# One row per fish; FW1_g and FW2_g in mm. Outliers (> 2.5 SD from the
# River x FW1_Year x age mean) and fish without FW1 growth are removed

growth <- growth_raw %>%
  mutate(age_group = if_else(age %in% c(12, 13), "12/13", "22/23")) %>%
  group_by(River, FW1_Year, age) %>%
  mutate(z = as.numeric(scale(FW1_g))) %>%
  filter(abs(z) < outlier_sd) %>%
  ungroup() %>%
  dplyr::select(-z)

cat("Removed", nrow(growth_raw) - nrow(growth), "fish (outliers + missing FW1)\n")


##### 4. Growth summary stats #####
# Annual mean FW1 growth (µm) by river, age class and modern vs historical
annual_means <- growth %>%
  group_by(River, age_group, FW1_Year) %>%
  summarise(n_fish         = sum(!is.na(FW1_g)),
            annual_mean_um = mean(FW1_g, na.rm = TRUE) * 1000, .groups = "drop") %>%
  mutate(period = case_when(FW1_Year <= historical_end_year ~ "Historical",
                            FW1_Year >= modern_start_year   ~ "Modern")) %>%
  filter(!is.na(period))

# Summary by period: number of fish, mean and SD of annual means (µm).
# The Results text reports the Modern rows.
growth_summary <- annual_means %>%
  group_by(River = river_names[River], age_group, period) %>%
  summarise(n_fish  = sum(n_fish),
            n_years = n(),
            mean_um = mean(annual_mean_um),
            sd_um   = sd(annual_mean_um),
            .groups = "drop")
print(growth_summary)

growth_period_tests <- annual_means %>%
  group_by(River, age_group) %>%
  group_modify(function(d, key) {
    h  <- d$annual_mean_um[d$period == "Historical"]
    m  <- d$annual_mean_um[d$period == "Modern"]
    tt <- t.test(h, m, var.equal = FALSE)
    tibble(hist_n_years = length(h), hist_mean_um = mean(h),
           modern_n_years = length(m), modern_mean_um = mean(m),
           t_stat = unname(tt$statistic), df = unname(tt$parameter), p_value = tt$p.value)
  }) %>%
  ungroup() %>%
  mutate(River = river_names[River])
print(growth_period_tests)


##### 4. Figure 1 (four panels) #####
#Note, map is added afterward, not in this code

snap_5yr <- snap %>%
  mutate(Lake = factor(Lake, levels = c("Aleknagik", "Nuyakuk")),
         Year_5 = floor(Year / 5) * 5) %>%
  group_by(Year_5, Lake) %>%
  summarise(Temp = mean(Temp, na.rm = TRUE), .groups = "drop")

aug_mean <- mean(lt_interp$Aug_Temp, na.rm = TRUE)
lt_anom  <- lt_interp %>% mutate(AnomSign = if_else(Aug_Temp >= aug_mean, "pos", "neg"))

make_x <- function(years, by) {
  lims <- range(years, na.rm = TRUE)
  list(scale_x_continuous(breaks = seq(floor(lims[1] / by) * by,
                                       ceiling(lims[2] / by) * by, by = by)),
       coord_cartesian(xlim = lims))
}
x_left  <- make_x(c(snap_5yr$Year_5, density$Year), by = 30)
x_right <- make_x(c(lt_anom$FW1_Year, esc_plot$EscpYear), by = 20)

base_theme <- theme_minimal(base_size = 18) +
  theme(panel.grid        = element_blank(),
        axis.line         = element_line(color = "black"),
        axis.ticks        = element_line(color = "black"),
        legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5))
no_x <- theme(axis.text.x = element_blank(), axis.title.x = element_blank())
legend_topleft <- theme(legend.position = "inside", legend.position.inside = c(0.02, 0.98),
                        legend.justification = c(0, 1))
watershed_colors <- c("Nushagak" = "grey60", "Wood" = "black")

p_snap <- ggplot(snap_5yr, aes(Year_5, Temp, color = Lake)) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = c("Aleknagik" = "black", "Nuyakuk" = "grey60")) +
  x_left +
  labs(y = "SNAP air\ntemperature (°C)", color = "Lake") +
  base_theme + no_x + legend_topleft

p_lake <- ggplot(lt_anom, aes(x = FW1_Year, color = AnomSign)) +
  geom_hline(yintercept = aug_mean, linetype = "dashed", color = "gray60") +
  geom_segment(aes(xend = FW1_Year, y = aug_mean, yend = Aug_Temp, alpha = Aug_Temp_observed),
               linewidth = 2) +
  scale_color_manual(values = c("pos" = "#d73027", "neg" = "#2c728e")) +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.35), guide = "none") +  # pale = filled year
  scale_y_continuous(breaks = pretty_breaks(n = 5)) +
  x_right +
  labs(y = "August lake\ntemperature (°C)") +
  base_theme + no_x + theme(legend.position = "none")

p_dens <- ggplot(density, aes(Year, Spawnerperkm2, color = Watershed)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = watershed_colors) +
  scale_y_continuous(breaks = seq(0, 20000, by = 5000)) +
  x_left +
  labs(x = "Year", y = "Spawners\nper km²", color = "Watershed") +
  base_theme + legend_topleft

p_esc <- ggplot(esc_plot, aes(EscpYear, Escapement, color = Watershed)) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = watershed_colors) +
  scale_y_log10(breaks = c(1e5, 3e5, 1e6, 3e6, 1e7), labels = function(x) x / 1e6) +
  x_right +
  labs(x = "Year", y = "Spawner escapement\n(millions, log scale)") +
  base_theme + theme(legend.position = "none")

fig1 <- wrap_plots(p_snap, p_lake, p_dens, p_esc, ncol = 2)


##### 5. Figure 2 #####
#Rows = river; columns = age 1.X (FW1) | age 2.X (FW1 + FW2 stacked).
#Baseline (0) = mean of the river's annual means (historical + modern).
#Historical point = mean of historical annual anomalies; whiskers = range.
#Historical annual means are grouped by sampling Year; modern by FW1_Year.


growth_named <- growth %>% mutate(RiverName = unname(river_names[River]))

shared_theme <- theme_minimal(base_size = 21) +
  theme(legend.position = "bottom", legend.title = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "grey80", linewidth = 0.3),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background  = element_rect(fill = "white", color = NA),
        axis.line = element_line(color = "black"),
        axis.line.x.top = element_blank(), axis.text.x.top = element_blank(),
        axis.ticks.x.top = element_blank(),
        axis.ticks = element_line(color = "black"), axis.ticks.length = unit(4, "pt"),
        plot.title = element_text(hjust = 0))

hist_x_breaks <- c(1932, seq(1960, 2020, by = 10))
hist_x_labels <- c("1912–1939", seq(1960, 2020, by = 10))
sign_colors   <- c("Above" = "#950404", "Below" = "#0F542F")

# Annual means (µm) for one age group, historical and modern
yearly_means <- function(age_grp, vars) {
  summarise_vars <- function(d) {
    d %>% summarise(across(all_of(vars), ~ mean(.x, na.rm = TRUE) * 1000, .names = "{.col}_mean"),
                    .groups = "drop")
  }
  d <- filter(growth_named, age_group == age_grp)
  list(hist   = d %>% filter(FW1_Year <= historical_end_year) %>%
                  group_by(RiverName, Year) %>% summarise_vars(),
       modern = d %>% filter(FW1_Year >= modern_start_year) %>%
                  group_by(RiverName, Year = FW1_Year) %>% summarise_vars())
}

##### 6, panels a-d. Freshwater growth anomalies #####
ym1 <- yearly_means("12/13", "FW1_g")
baseline_age1 <- ym1$modern %>%
  group_by(RiverName) %>% summarise(baseline = mean(FW1_g_mean, na.rm = TRUE), .groups = "drop")

modern_age1 <- ym1$modern %>% left_join(baseline_age1, by = "RiverName") %>%
  mutate(FW1_anomaly = FW1_g_mean - baseline)
hist_age1 <- ym1$hist %>% left_join(baseline_age1, by = "RiverName") %>%
  mutate(FW1_anomaly = FW1_g_mean - baseline) %>%
  group_by(RiverName) %>%
  summarise(point_anomaly = mean(FW1_anomaly, na.rm = TRUE), ymin = min(FW1_anomaly, na.rm = TRUE),
            ymax = max(FW1_anomaly, na.rm = TRUE), .groups = "drop")

make_panel_age1 <- function(river, title_label, show_x_label = FALSE) {
  mod_dat  <- modern_age1 %>% filter(RiverName == river) %>%
    mutate(Sign = if_else(FW1_anomaly >= 0, "Above", "Below"))
  hist_dat <- hist_age1 %>% filter(RiverName == river) %>%
    mutate(Sign = if_else(point_anomaly >= 0, "Above", "Below"))
  ggplot() +
    geom_hline(yintercept = 0, color = "grey40", linewidth = 0.4, linetype = "dashed") +
    geom_errorbar(data = hist_dat, aes(x = 1932, ymin = ymin, ymax = ymax, color = Sign),
                  width = 2, linewidth = 0.8) +
    geom_point(data = hist_dat, aes(x = 1932, y = point_anomaly, color = Sign), size = 3) +
    geom_col(data = mod_dat, aes(x = Year, y = FW1_anomaly, fill = Sign), width = 1) +
    scale_fill_manual(values = sign_colors, guide = "none") +
    scale_color_manual(values = sign_colors, guide = "none") +
    scale_x_continuous(breaks = hist_x_breaks, labels = hist_x_labels) +
    scale_y_continuous(limits = c(-100, 100), breaks = seq(-100, 100, by = 20)) +
    scale_x_break(c(1937, 1959), space = 0.05) +
    labs(x = if (show_x_label) "Year" else "", y = "Growth anomaly (μm)", title = title_label) +
    shared_theme + theme(legend.position = "none")
}


ym2 <- yearly_means("22/23", c("FW1_g", "FW2_g"))
baseline_age2 <- ym2$modern %>%
  group_by(RiverName) %>%
  summarise(baseline_FW1 = mean(FW1_g_mean, na.rm = TRUE),
            baseline_FW2 = mean(FW2_g_mean, na.rm = TRUE), .groups = "drop")

modern_age2 <- ym2$modern %>% left_join(baseline_age2, by = "RiverName") %>%
  transmute(RiverName, Year,
            FW1_anomaly = FW1_g_mean - baseline_FW1,
            FW2_anomaly = FW2_g_mean - baseline_FW2) %>%
  pivot_longer(c(FW1_anomaly, FW2_anomaly), names_to = "Zone", values_to = "Anomaly") %>%
  mutate(Zone = factor(Zone, levels = c("FW2_anomaly", "FW1_anomaly")))   # FW2 stacked on top

hist_age2 <- ym2$hist %>% left_join(baseline_age2, by = "RiverName") %>%
  transmute(RiverName, FW1 = FW1_g_mean - baseline_FW1, FW2 = FW2_g_mean - baseline_FW2) %>%
  pivot_longer(c(FW2, FW1), names_to = "Metric", values_to = "anom") %>%
  group_by(RiverName, Metric) %>%
  summarise(anomaly = mean(anom, na.rm = TRUE), ymin = min(anom, na.rm = TRUE),
            ymax = max(anom, na.rm = TRUE), .groups = "drop")

stack_totals <- modern_age2 %>% group_by(RiverName, Year) %>%
  summarise(total = sum(Anomaly), .groups = "drop")
y_lim_age2 <- c(floor(min(c(stack_totals$total, hist_age2$ymin), na.rm = TRUE) / 50) * 50,
                ceiling(max(c(stack_totals$total, hist_age2$ymax), na.rm = TRUE) / 50) * 50)

zone_sign_fill    <- c("FW1 above" = "#950404", "FW1 below" = "#0F542F",
                       "FW2 above" = "white",   "FW2 below" = "white")
zone_sign_outline <- c("FW1 above" = "#950404", "FW1 below" = "#0F542F",
                       "FW2 above" = "#950404", "FW2 below" = "#0F542F")

make_panel_age2 <- function(river, title_label, show_x_label = TRUE) {
  mod_dat <- modern_age2 %>% filter(RiverName == river) %>%
    mutate(Zone_sign = paste(if_else(Zone == "FW1_anomaly", "FW1", "FW2"),
                             if_else(Anomaly >= 0, "above", "below")))
  hist_dat <- hist_age2 %>% filter(RiverName == river) %>%
    mutate(Zone_sign = paste(Metric, if_else(anomaly >= 0, "above", "below")),
           x_pos     = if_else(Metric == "FW2", 1930, 1934))
  ggplot() +
    geom_hline(yintercept = 0, color = "grey40", linewidth = 0.4, linetype = "dashed") +
    geom_col(data = mod_dat,
             aes(x = Year, y = Anomaly, fill = Zone_sign, color = Zone_sign, group = Zone),
             position = "stack", width = 1, linewidth = 0.35) +
    geom_errorbar(data = hist_dat, aes(x = x_pos, ymin = ymin, ymax = ymax, color = Zone_sign),
                  width = 1.5, linewidth = 0.8, show.legend = FALSE) +
    geom_point(data = hist_dat, aes(x = x_pos, y = anomaly, fill = Zone_sign, color = Zone_sign),
               size = 3, shape = 21, stroke = 0.8, show.legend = FALSE) +
    scale_fill_manual(values = zone_sign_fill, name = NULL) +
    scale_color_manual(values = zone_sign_outline, name = NULL) +
    scale_x_continuous(breaks = hist_x_breaks, labels = hist_x_labels) +
    scale_y_continuous(limits = y_lim_age2, breaks = seq(y_lim_age2[1], y_lim_age2[2], by = 50)) +
    scale_x_break(c(1937, 1959)) +
    labs(x = if (show_x_label) "Year" else "", y = "Growth anomaly (μm)", title = title_label) +
    guides(fill = guide_legend(ncol = 1), color = "none") +
    shared_theme + theme(legend.position = "right")
}


panel_nush_age1 <- make_panel_age1("Nushagak", "a) Nushagak 1.X")
panel_nush_age2 <- make_panel_age2("Nushagak", "b) Nushagak 2.X", show_x_label = FALSE)
panel_wood_age1 <- make_panel_age1("Wood",     "c) Wood 1.X", show_x_label = TRUE)
panel_wood_age2 <- make_panel_age2("Wood",     "d) Wood 2.X")

no_legend    <- theme(legend.position = "none")
render_break <- function(p) ggplotify::as.ggplot(print(p))   # ggbreak plot -> arrangeable grob

fig2 <- cowplot::plot_grid(
  cowplot::plot_grid(render_break(panel_nush_age1 + no_legend),
                     render_break(panel_nush_age2 + no_legend),
                     render_break(panel_wood_age1 + no_legend),
                     render_break(panel_wood_age2 + no_legend), ncol = 2),
  cowplot::get_legend(panel_nush_age2),
  ncol = 2, rel_widths = c(1, 0.12)
)

##### 7. Marss modelling of FW1 growth #####

# Two-river state-space model (state rows: W, N):
#y_t = x_t + v_t,                 v_t ~ MVN(0, R_t); R_t diagonal
#x_t = B x_{t-1} + C c_t + w_t,   w_t ~ MVN(0, Q); Q diagonal and equal
#B diagonal: AR(1) estimated per river, or random walk (B = 1) for comparison
#Response and covariates are z-scored. Covariates: August lake temperature
#(both rivers) and log10 brood-year escapement (each river its own).
#16 covariate combinations x {AR(1), RW} = 32 models. Results are averaged
#across the 16 AR(1) models with AICc weights:
#theta_bar  = sum_i w_i theta_i   (theta_i = 0 when absent from model i)
#se_unc     = sum_i w_i sqrt(var_i + (theta_i - theta_bar)^2)


zscore <- function(x) {
  m <- mean(x, na.rm = TRUE); s <- sd(x, na.rm = TRUE)
  list(z = (x - m) / s, sd = s)
}

# Years kept = both rivers have growth, escapement and temperature, so every
# model is fit to identical data
annual_fw1 <- growth %>%
  filter(age_group == "12/13") %>%
  group_by(River, FW1_Year) %>%
  summarise(FW1_mean = mean(FW1_g, na.rm = TRUE),
            FW1_sd   = sd(FW1_g, na.rm = TRUE),
            n        = sum(!is.na(FW1_g)), .groups = "drop")

river_series <- function(rv) {
  annual_fw1 %>%
    filter(River == rv) %>%
    left_join(if (rv == "W") wbrood else nbrood, by = "FW1_Year") %>%
    left_join(dplyr::select(lt_interp, FW1_Year, Aug_Temp), by = "FW1_Year") %>%
    filter(is.finite(log10_Esc), !is.na(FW1_mean)) %>%
    arrange(FW1_Year)
}
sW <- river_series("W"); sN <- river_series("N")
stopifnot(!anyDuplicated(sW$FW1_Year), !anyDuplicated(sN$FW1_Year))

resp_mu <- c(W = mean(sW$FW1_mean), N = mean(sN$FW1_mean))
resp_sd <- c(W = sd(sW$FW1_mean),   N = sd(sN$FW1_mean))

years <- sort(union(sW$FW1_Year, sN$FW1_Year))
fW <- tibble(FW1_Year = years) %>% left_join(sW, by = "FW1_Year")
fN <- tibble(FW1_Year = years) %>% left_join(sN, by = "FW1_Year")
clim <- coalesce(fW$Aug_Temp, fN$Aug_Temp)
keep <- !is.na(clim) & is.finite(fW$log10_Esc) & is.finite(fN$log10_Esc)
fW <- fW[keep, ]; fN <- fN[keep, ]; clim <- clim[keep]
model_years <- years[keep]
cat("MARSS years:", length(model_years), "(", min(model_years), "-", max(model_years), ")\n")

z_clim <- zscore(clim); z_eW <- zscore(fW$log10_Esc); z_eN <- zscore(fN$log10_Esc)

Y <- rbind(W = (fW$FW1_mean - resp_mu[["W"]]) / resp_sd[["W"]],
           N = (fN$FW1_mean - resp_mu[["N"]]) / resp_sd[["N"]])
covariates <- rbind(Climate = z_clim$z, Esc_W = z_eW$z, Esc_N = z_eN$z)

# Known observation variances (squared SE of annual mean, z-scale);
# zero or missing values replaced with that river's median
fill_med <- function(v) ifelse(is.na(v) | v == 0, median(v, na.rm = TRUE), v)
R_array <- array(0, dim = c(2, 2, length(model_years)))
R_array[1, 1, ] <- fill_med((fW$FW1_sd^2 / fW$n) / resp_sd[["W"]]^2)
R_array[2, 2, ] <- fill_med((fN$FW1_sd^2 / fN$n) / resp_sd[["N"]]^2)

#  (effect = theta * SD_response / SD_covariate)
param_info <- tibble(
  parameter    = c("CW_climate", "CW_esc", "CN_climate", "CN_esc"),
  RiverName    = c("Wood", "Wood", "Nushagak", "Nushagak"),
  short        = rep(c("August\ntemperature", "Spawner\nescapement"), 2),
  label        = paste0(RiverName, ": ", rep(c("August temperature", "Spawner escapement"), 2)),
  sd_response  = resp_sd[c("W", "W", "N", "N")] * 1000,                    # µm
  sd_covariate = c(z_clim$sd, z_eW$sd, z_clim$sd, z_eN$sd),
  unit         = rep(c("µm per 1 °C", "µm per 10-fold increase in spawners"), 2)
)
params <- param_info$parameter

# C matrix template (columns: Climate, Esc_W, Esc_N; 0)
C_template <- matrix(c("CW_climate", "CW_esc", "0",
                       "CN_climate", "0",      "CN_esc"), nrow = 2, byrow = TRUE)

specs <- tidyr::expand_grid(dynamics = c("AR1", "RW"),
                            CW_climate = c(FALSE, TRUE), CW_esc = c(FALSE, TRUE),
                            CN_climate = c(FALSE, TRUE), CN_esc = c(FALSE, TRUE))
cov_label <- function(clim_on, esc_on) {
  case_when(clim_on & esc_on ~ "Temperature + escapement",
            clim_on          ~ "Temperature",
            esc_on           ~ "Escapement",
            TRUE             ~ "None")
}
specs <- specs %>%
  mutate(Wood_covariates     = cov_label(CW_climate, CW_esc),
         Nushagak_covariates = cov_label(CN_climate, CN_esc),
         model = paste0(dynamics, " | W: ", Wood_covariates, " / N: ", Nushagak_covariates))

build_C <- function(spec) {
  C <- matrix(list(0), 2, 3)
  for (i in 1:2) for (j in 1:3) {
    p <- C_template[i, j]
    if (p != "0" && isTRUE(spec[[p]])) C[[i, j]] <- p
  }
  C
}
build_B <- function(dynamics) {
  if (dynamics == "AR1") matrix(list("BW", 0, 0, "BN"), 2, 2) else diag(2)
}

fits <- setNames(vector("list", nrow(specs)), specs$model)
for (i in seq_len(nrow(specs))) {
  fits[[i]] <- MARSS(Y, model = list(B = build_B(specs$dynamics[i]), U = matrix(0, 2, 1),
                                     Q = "diagonal and equal", Z = "identity", A = "scaling",
                                     R = R_array, C = build_C(specs[i, ]), c = covariates),
                     control = list(maxit = 2000, allow.degen = TRUE), silent = TRUE)
}
conv <- sapply(fits, function(f) f$convergence)
if (any(conv != 0)) warning("Non-zero convergence: ", paste(names(conv)[conv != 0], collapse = "; "))

aicc_weights <- function(tab) {
  tab %>% mutate(delta_AICc = AICc - min(AICc),
                 weight     = exp(-delta_AICc / 2) / sum(exp(-delta_AICc / 2))) %>%
    arrange(AICc)
}

aicc_all <- specs %>%
  dplyr::select(model, dynamics, Wood_covariates, Nushagak_covariates) %>%
  mutate(K      = sapply(fits[model], function(f) length(unlist(f$par))),
         logLik = sapply(fits[model], function(f) f$logLik),
         AICc   = sapply(fits[model], function(f) f$AICc)) %>%
  aicc_weights()


dynamics_support <- aicc_all %>%
  group_by(dynamics) %>%
  summarise(total_weight = sum(weight), best_delta_AICc = min(delta_AICc), .groups = "drop")
print(dynamics_support)


aicc_ar1 <- aicc_all %>% filter(dynamics == "AR1") %>% aicc_weights()

appendix2 <- aicc_ar1 %>%
  transmute(Rank = row_number(), `Wood covariates` = Wood_covariates,
            `Nushagak covariates` = Nushagak_covariates, K,
            delta_AICc = round(delta_AICc, 2), weight = signif(weight, 3))
print(appendix2, n = Inf)


coefs <- bind_rows(lapply(aicc_ar1$model, function(m) {
  fit <- fits[[m]]
  U   <- fit$par$U
  est <- se <- setNames(rep(0, length(params)), params)
  present <- intersect(rownames(U), params)
  if (length(present)) {
    ci  <- MARSSparamCIs(fit, method = "hessian")
    pos <- match(present, rownames(U))
    est[present] <- U[pos, 1]
    se[present]  <- ci$par.se$U[pos, 1]
  }
  tibble(model = m, parameter = params, estimate = est, se = se, in_model = params %in% present)
}))
stopifnot(all(is.finite(coefs$se[coefs$in_model]) & coefs$se[coefs$in_model] > 0))

avg <- coefs %>%
  left_join(dplyr::select(aicc_ar1, model, weight), by = "model") %>%
  group_by(parameter) %>%
  summarise(theta_bar  = sum(weight * estimate),
            se_unc     = sum(weight * sqrt(se^2 + (estimate - theta_bar)^2)),
            importance = sum(weight[in_model]), .groups = "drop") %>%
  mutate(lo90 = theta_bar - z90 * se_unc, hi90 = theta_bar + z90 * se_unc) %>%
  left_join(param_info, by = "parameter") %>%
  mutate(effect      = theta_bar * sd_response / sd_covariate,
         effect_lo90 = lo90      * sd_response / sd_covariate,
         effect_hi90 = hi90      * sd_response / sd_covariate) %>%
  arrange(match(parameter, params))

model_avg_table <- avg %>%
  transmute(Covariate = label,
            theta = round(theta_bar, 3), SE = round(se_unc, 3),
            `90% CI` = sprintf("%.3f to %.3f", lo90, hi90),
            Importance = round(importance, 2),
            `Effect (90% CI)` = sprintf("%.1f (%.1f to %.1f)", effect, effect_lo90, effect_hi90),
            Unit = unit)
print(model_avg_table)


plot_river <- function(river, show_y) {
  d <- avg %>% filter(RiverName == river) %>%
    mutate(short = factor(short, levels = rev(unique(avg$short))))
  p <- ggplot(d, aes(theta_bar, short)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_errorbar(aes(xmin = lo90, xmax = hi90), width = 0.15, linewidth = 0.9,
                  orientation = "y") +
    geom_point(size = 4, shape = 21, fill = "steelblue", color = "black") +
    coord_cartesian(xlim = c(-0.4, 0.4)) +
    labs(title = river, x = "Effect size (standardized)", y = NULL) +
    theme_bw(base_size = 16) +
    theme(panel.grid.minor = element_blank())
  if (!show_y) p <- p + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  p
}
fig3 <- (plot_river("Nushagak", TRUE) + plot_river("Wood", FALSE)) +
  plot_annotation(tag_levels = "a", tag_suffix = ")")



