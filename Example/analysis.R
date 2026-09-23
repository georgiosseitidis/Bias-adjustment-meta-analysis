###############################
# Load libraries
###############################
library(R2jags)
library(ggplot2)
library(readxl)

# Load Meta-analysis models and Data
source("MA model and Help functions.R") # Help function and Meta-analysis models

# Read data
d <- read_xlsx("Data.xlsx") 

###############################################################################################################
#
#                                                   Meta-analysis
#
# Apply 4 Meta-analysis models
#     1) Bias adjusted MA model
#     2) Naive meta-analysis ignoring the study design
#     3) Naive meta-analysis using only the RCTs
#     4) Naive meta-analysis using only the non-RCTs
#
###############################################################################################################

###########################
# 1) Bias adjusted model
###########################

# Design specific RoB category 
d$cat_id <- NA
for(i in 1:dim(d)[1]){d$cat_id[i] <- get_cat_id(d$design[i], d$Risk[i])}

# MCMC model data
data_MCMC_bias <- list(
  y = d$TE,
  prec_s = 1/d$seTE^2,
  RoB = d$cat_id,
  RCT = c(1, 1, 1, 0 ,0 , 0), 
  ROB_Mod = c(0, 1, 0, 0, 1, 0), 
  ROB_High = c(0, 0, 1, 0, 0, 1), 
  n = dim(d)[1],
  n_des = 6,
  design = as.factor(d$RCT)
)

set.seed(2026)
fit <- jags(data = data_MCMC_bias,
            parameters.to.save = c("mu", "b", "theta", "tau", "gamma", "w","w_star", "w_NRE", "D", 'D_NRE_bias'),
            model.file = textConnection(model_ma_bias),
            n.chains = 3, n.iter = 150000, n.burnin = 50000)

w_ssvs <- 100*fit$BUGSoutput$median$w_star # Relative weights

########################################
# 2) Naive model using all studies
########################################

# MCMC model data
data_MCMC <- list(y = d$TE, prec_s = 1/d$seTE^2, n = dim(d)[1])

set.seed(2026)
fit_ma_all <- jags(data = data_MCMC,
                   parameters.to.save = c("mu", "theta", "tau", "w", "w_star", "w_NRE"),
                   model.file = textConnection(model_ma),
                   n.chains = 3, n.iter = 150000, n.burnin = 50000)

w_naive <- 100*fit_ma_all$BUGSoutput$median$w_star # Relative weights

####################################
# 3) Naive model using only RCTs
#####################################

# MCMC model data
data_MCMC_RCT <- list(y = d$TE[which(d$RCT == 1)], prec_s = 1/d$seTE[which(d$RCT == 1)], n = length(which(d$RCT == 1)))
set.seed(2026)
fit_ma_rct <- jags(data = data_MCMC_RCT,
                   parameters.to.save = c("mu", "theta", "tau"),
                   model.file = textConnection(model_ma_2),
                   n.chains = 3, n.iter = 150000, n.burnin = 50000)

#######################################
# 4) Naive model using only non-RCTs
########################################

# MCMC model data
data_MCMC_OBS <- list(y = d$TE[-which(d$RCT == 1)], prec_s = 1/d$seTE[-which(d$RCT == 1)], n = dim(d)[1] - length(which(d$RCT == 1)))
set.seed(2026)
fit_ma_os <- jags(data = data_MCMC_OBS,
                  parameters.to.save = c("mu", "theta", "tau"),
                  model.file = textConnection(model_ma_2),
                  n.chains = 3, n.iter = 150000, n.burnin = 50000)


######################################
#              Results
######################################

# Figure 3 - Forest plot 
forest_fit_all(fit_bias = fit, fit_naive = fit_ma_all, fit_ma_os = fit_ma_os, fit_ma_rct = fit_ma_rct, d) + scale_x_log10(breaks = c(0.1, 0.2, 0.3, 0.4, 0.5,0.6, 0.8, 1, 1.2, 1.5, 2, 2.5))

# Data-frame with study weights
df_w <- data.frame(element = rep(seq_along(w_ssvs), 2), vector = rep(c("Bias Adjusted", "Naive RE meta-analysis"), each = length(w_ssvs)), value = c(fit$BUGSoutput$median$w, fit_ma_all$BUGSoutput$median$w))

# Figure 4 - Study Weights
ggplot(df_w, aes(x = factor(element), y = value, fill = vector)) +
  geom_col(position = "dodge") +
  labs(x = "Study", y = "Weights", fill = "",) +
  theme_minimal()+
  theme(legend.position = "bottom",
        axis.text.x = element_text(size = 13, color = "black"),
        axis.text.y = element_text(size = 13, color = "black"),
        legend.text = element_text(size = 13, color = "black"),
        axis.title = element_text(size = 14, color = "black"),
        axis.title.x = element_text(vjust = 0.5, size = 14, color = "black"),
        legend.title = element_text(size = 14, color = "black"),
        panel.grid.minor = element_blank())

# Table S2
Table_S2 <- data.frame("study" = d$name, "design" = d$design, "Rob" = d$Risk, "w_ba" = fit$BUGSoutput$median$w, 
                       "omega_ba" = fit$BUGSoutput$median$w_star, "Di_bias" = fit$BUGSoutput$median$D, 
                       "Di_total" = fit$BUGSoutput$median$w/fit_ma_all$BUGSoutput$median$w)

# Figure S2 - Posterior Inclusion Probabilities
lolipo_pip(fit)

# Figure S3 - Spike & Slab
spike_slab(fit)

# Data-frame with relative study weights
df_w <- data.frame(element = rep(seq_along(w_ssvs), 2), vector = rep(c("Bias Adjusted", "Naive RE meta-analysis"), each = length(w_ssvs)), value = c(w_ssvs, w_naive))

# Figure S4 - Relative weights
ggplot(df_w, aes(x = factor(element), y = value, fill = vector)) +
  geom_col(position = "dodge") +
  labs(x = "Study", y = "Relative Weights", fill = "",) +
  theme_minimal()+
  theme(legend.position = "bottom",
        axis.text.x = element_text(size = 13, color = "black"),
        axis.text.y = element_text(size = 13, color = "black"),
        legend.text = element_text(size = 13, color = "black"),
        axis.title = element_text(size = 14, color = "black"),
        legend.title = element_text(size = 14, color = "black"),
        panel.grid.minor = element_blank())

