#####################################################################################################################################
#
#                                 Load libraries - functions & simulation results
#
######################################################################################################################################

# Libraries
library(ggplot2)
library(tidyverse)
library(scales)
library(ggh4x)

# source functions
source("MA Models & Help Functions.R")

# simulation results
load("simulation_results\\simulation_results.RData")

########################################################################################################################################
#
#                                                     Transform data in long format
#
########################################################################################################################################

# Make a scenario id column for each simulation scenario
scenarios$id <- paste0(sprintf("S%02d", 1:dim(scenarios)[1]))

# index column to merge results with simulated scenarios
scenarios$scenario <- paste0("rand=", scenarios$rand, "_rate=", scenarios$n_rct_rate, "_ROB=", scenarios$rob_structure, "_Bias=", scenarios$bias_pattern)

# Long format
sim_long <- bind_rows(
  results %>%
    transmute(scenario, method = "Bias-adjusted SSVS", bias = bias, coverage = as.numeric(coverage), ci_width = CIw),
  results %>%
    transmute(scenario, method = "Naive: all studies", bias = bias_naive, coverage = as.numeric(coverage_naive), ci_width = CIw_naive),
  results %>%
    transmute(scenario, method = "Naive: RCT only", bias = bias_rct_est, coverage = as.numeric(coverage_rct), ci_width = CIw_rct),
  results %>%
    transmute(scenario, method = "Naive: OBS only", bias = bias_obs_est, coverage = as.numeric(coverage_obs), ci_width = CIw_obs)
) %>%
  left_join(scenarios, by = "scenario") %>%
  mutate(
    method = factor(method, levels = c("Bias-adjusted SSVS", "Naive: all studies", "Naive: RCT only", "Naive: OBS only")),
    bias_pattern = factor(bias_pattern, levels = c("No Bias", "NRE moderate", "NRE strong", "Rct-NRE Moderate", "Rct20_NRE20-35"),  labels = c("No Bias", "NRE moderate", "NRE strong", "RCT/NRE Moderate", "RCT 20%, NRE 20-35%")),
    id = factor(id, levels = scenarios$id)
  )

# Summary table with Performance metrics
summary_tbl <- sim_long %>%
  group_by(id, scenario, rand, n_rct_rate, rob_structure, bias_pattern, method) %>%
  summarise(
    n_sim = n(),
    mean_bias = mean(bias, na.rm = TRUE),
    abs_bias = mean(abs(bias), na.rm = TRUE),
    coverage = mean(coverage, na.rm = TRUE),
    mcse_coverage = sqrt(coverage * (1 - coverage) / n_sim),
    ci_width = mean(ci_width, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    coverage_lcl = pmax(0, coverage - 1.96 * mcse_coverage),
    coverage_ucl = pmin(1, coverage + 1.96 * mcse_coverage)
  )

# Manuscript's Table S3 
Table_s3 <- sim_long %>%
  group_by(id, scenario, rand, n_rct_rate, rob_structure, bias_pattern, method) %>%
  summarise(
    n_sim = n(),
    mean_bias = mean(bias, na.rm = TRUE),
    mean_bias_mcse = sqrt((1/(n_sim*(n_sim-1)))*sum((bias - mean(bias))^2)),
    abs_bias = mean(abs(bias), na.rm = TRUE),
    coverage = mean(coverage, na.rm = TRUE),
    mcse_coverage = sqrt(coverage * (1 - coverage) / n_sim),
    ci_wid = mean(ci_width, na.rm = TRUE),
    ci_widht_mcse = sqrt((1/(n_sim*(n_sim-1)))*sum((ci_width - mean(ci_width))^2)),
    .groups = "drop"
  )

# writexl::write_xlsx(Table_s3, "TS3.xlsx")

############################################################################################################################
#
#                                     Visualization of performance metrics
#
#############################################################################################################################

summary_tbl$index <- paste0("Rate: ", summary_tbl$n_rct_rate, " ROB: ", summary_tbl$rob_structure) # index column

# Helper to avoid repeating the percentage formatting
coverage_scale <- function(breaks) {scale_y_continuous(breaks = breaks,labels = scales::label_percent(accuracy = 1))}

# Manuscript's Figure S8
Fig_S8 <- ggplot(summary_tbl, aes(x= index, y=coverage, colour = rand)) +
  geom_point(size = 1.8, position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = coverage_lcl, ymax = coverage_ucl), position = position_dodge(width = 0.9),linewidth = 0.4, width=0.4, alpha=0.9, size=1.3) +  
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red", linewidth = 0.7) +
  facet_grid(method~bias_pattern, scales = "free")+
  coord_flip() +
  facetted_pos_scales(
    y = list(
      COL == 1 ~ coverage_scale(c(0.94, 0.96, 0.98, 1.00)),
      COL == 2 ~ coverage_scale(c(0.80, 0.85, 0.9, 0.95, 1.00)),
      COL == 3 ~ coverage_scale(c(0.60, 0.70, 0.80, 0.9, 1)),
      COL == 4 ~ coverage_scale(c(0.60, 0.70, 0.80, 0.9, 1)),
      COL == 5 ~ coverage_scale(c(0.25, 0.50, 0.75, 0.9))
    )) +
  xlab("") + ylab("Coverage") + labs(colour = "Randomization") +
  theme(legend.position = "bottom", text = element_text(size = 18), axis.text.x = element_text(size = 16) )

tiff("Figure S8.tiff", width = 12000, height = 8280, res = 600, compression = "lzw")
Fig_S8
dev.off()

sim_long$index <- paste0("Rate: ", sim_long$n_rct_rate, " ROB: ", sim_long$rob_structure) # index column 

# Figure S5 - Absolute Bias for each simulated scenario
sim_long$abs_bias <- abs(sim_long$bias)
Fig_S5 <- performance_boxplot(sim_long, "abs_bias", "Distribution of absolute bias across simulation replicates", "Absolute Bias", hline = 0)
tiff("Figure S5.tiff", width = 14000, height = 8280, res = 600, compression = "lzw")
Fig_S5
dev.off()


# Figure S6 - Bias for each simulated scenario
Fig_S6 <- performance_boxplot(sim_long, "bias", "Distribution of bias across simulation replicates", "Bias on OR scale", hline = 0)
tiff("Figure S6.tiff", width = 14000, height = 8280, res = 600, compression = "lzw")
Fig_S6
dev.off()

# Figure S7 - Credible interval width for each simulated scenario
Fig_S7 <- performance_boxplot(sim_long, "ci_width", "Distribution of credible interval width across simulation replicates", "Credible interval width")
tiff("Figure S7.tiff", width = 14000, height = 8280, res = 600, compression = "lzw")
Fig_S7
dev.off()

################################################################
# Summarize results by method used and Bias scenario
################################################################

# Table 3
Table_3 <- sim_long %>%
  group_by(bias_pattern, method) %>%
  summarise(
    n = n(),
    
    mean_abs_bias = mean(abs_bias, na.rm = TRUE),
    se_abs_bias = sd(abs_bias, na.rm = TRUE) / sqrt(n),
    abs_bias_lcl = mean_abs_bias - 1.96 * se_abs_bias,
    abs_bias_ucl = mean_abs_bias + 1.96 * se_abs_bias,
    
    mean_CIw = mean(ci_width, na.rm = TRUE),
    se_CIw = sd(ci_width, na.rm = TRUE) / sqrt(n),
    CIw_lcl = mean_CIw - 1.96 * se_CIw,
    CIw_ucl = mean_CIw + 1.96 * se_CIw,
    
    coverage = mean(coverage, na.rm = TRUE),
    coverage_lcl = coverage - 1.96 * sqrt(coverage * (1 - coverage) / n),
    coverage_ucl = coverage + 1.96 * sqrt(coverage * (1 - coverage) / n),
    
    .groups = "drop"
  ) %>%
  mutate(coverage_lcl = pmax(0, coverage_lcl), coverage_ucl = pmin(1, coverage_ucl),
    abs_bias_text = sprintf("%.2f (%.2f, %.2f)", mean_abs_bias, abs_bias_lcl, abs_bias_ucl),
    CIw_text = sprintf("%.2f (%.2f, %.2f)", mean_CIw, CIw_lcl, CIw_ucl),
    coverage_text = sprintf("%.1f%% (%.1f%%, %.1f%%)", 100 * coverage, 100 * coverage_lcl, 100 * coverage_ucl)
  ) %>%
  select(bias_pattern, method, n, `Absolute bias, mean (95% CI)` = abs_bias_text, `CI width, mean (95% CI)` = CIw_text, `Coverage, % (95% CI)` = coverage_text)

Table_3$method <- factor(Table_3$method, levels = c("Bias-adjusted SSVS",  "Naive: OBS only", "Naive: RCT only", "Naive: all studies"))
Table_3 <- Table_3[order(Table_3$bias_pattern, Table_3$method), ]
writexl::write_xlsx(Table_3, "Table 3.xlsx") 

# Figure 5 - Absolute bias summarized by bias assign scenario
Fig_5 <- ggplot(sim_long, aes(x = bias_pattern, y = abs_bias, fill = method)) +
  scale_y_continuous(limits = c(0, 1.2), breaks = seq(0,1.2,0.1))+
  geom_boxplot(position = position_dodge(width = 0.80), width = 0.65, outlier.alpha = 0.25, outlier.size = 0.7) +
  scale_fill_manual(values = method_cols, drop = FALSE) +
  labs(x = "Bias scenario", y = "Absolute bias", fill = NULL) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom",
        panel.grid.major.x = element_blank(),
        axis.text.x = element_text(size=11, face = "bold"),
        axis.title = element_text(size=14, face = "bold"),
        plot.title = element_text(face = "bold"),
        legend.text = element_text(size=12, face = "bold"),
        panel.grid.minor.y = element_blank())

tiff("Figure 5.tiff", width = 6000, height = 4080, res = 600, compression = "lzw")
Fig_5
dev.off()

# Figure 7 - Credible interval width summarized by bias assign scenario
Fig_7 <- ggplot(sim_long, aes(x = bias_pattern, y = ci_width, fill = method)) +
  scale_y_continuous(breaks = seq(0, 3,0.5))+
  geom_boxplot(position = position_dodge(width = 0.80), width = 0.65, outlier.alpha = 0.25, outlier.size = 0.7) +
  scale_fill_manual(values = method_cols, drop = FALSE) +
  labs(x = "Bias scenario", y = "Credible interval width", fill = NULL) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(size=11, face = "bold"),
        axis.title = element_text(size=14, face = "bold"),
        plot.title = element_text(face = "bold"),
        legend.text = element_text(size=12, face = "bold")
  )

tiff("Figure 7.tiff", width = 6000, height = 4080, res = 600, compression = "lzw")
Fig_7
dev.off()

# Figure 6 - Coverage summarized by bias assign scenario
coverage_df <- sim_long %>%
  group_by(bias_pattern, method) %>%
  summarise(
    n = sum(!is.na(coverage)),
    coverage = mean(coverage, na.rm = TRUE),
    coverage_se = sqrt(coverage * (1 - coverage) / n),
    coverage_lcl = pmax(0, coverage - 1.96 * coverage_se),
    coverage_ucl = pmin(1, coverage + 1.96 * coverage_se),
    .groups = "drop"
  )


Fig_6 <- ggplot(coverage_df, aes(x = bias_pattern, y = coverage, fill = method, group = method)) +
  geom_col(position = position_dodge(width = 0.85), width = 0.75, color = "white", linewidth = 0.25) +
  geom_errorbar(aes(ymin = coverage_lcl, ymax = coverage_ucl), position = position_dodge(width = 0.85),
    width = 0.20, linewidth = 0.75, color = "black") +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red", linewidth = 0.7) +
  scale_fill_manual(values = method_cols, drop = FALSE) +
  scale_y_continuous(limits = c(0, 1.02), breaks = seq(0, 1, 0.05), labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Bias scenario", y = "Coverage", fill = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text.x = element_text(size=11, face = "bold"),
    axis.title = element_text(size=14, face = "bold"),
    plot.title = element_text(face = "bold"),
    legend.text = element_text(size=12, face = "bold")
  )

tiff("Figure 6.tiff", width = 6000, height = 4080, res = 600, compression = "lzw")
Fig_6
dev.off()


