##############################################################################################################################
#
#                                          Meta-analysis models
#
##############################################################################################################################

# Proposed Bias adjusted Meta-analysis model
model_ma_bias_adj <- "
model {

  # Likelihood
  for (i in 1:n) {
    y[i] ~ dnorm(theta[i], prec_s[i])
    theta[i] ~ dnorm(mu + b[RoB[i]], prec_tau[design[i]])
  }
    
    # SSVS prior on b[i]
    for(j in 1:n_des){
    
      b[j] ~ dnorm(0, prec_b[j])
    
      # Define precision based on gamma
      prec_b[j] <- equals(gamma[j], 0) * inv_psi2[j] + equals(gamma[j], 1) * inv_c2psi2[j]
    
      # Prior inclusion probability
      gamma[j] <- equals(RCT[j], 0)*1 + equals(RCT[j], 1)*g[j] # NON RCTs heve always bias correction
    
      # Tuning parameters
      inv_psi2[j] <- pow(psi[j], -2) # spike
      inv_c2psi2[j] <- pow(c[j]*c[j] * psi[j]*psi[j], -1) #slab
    
      # Set ci based on the ROB
      c[j] <- equals(ROB_Mod[j], 0)*equals(ROB_High[j], 0)*1 + 
            equals(ROB_Mod[j], 1)*equals(ROB_High[j], 0)*3.57 + 
            equals(ROB_Mod[j], 0)*equals(ROB_High[j], 1)*7
      psi[j] <- 0.1
    
    
    # RCTs have inclusion probabilities
    g[j] ~ dbern(p[j])
    logit(p[j]) <- l_0 + l_1*ROB_Mod[j] + l_2*ROB_High[j]
    
    }

  # Priors
  mu ~ dnorm(0, 0.001)
  
  #heterogeneity 1 RCT, 2 OBS
  for(i in 1:2){
    tau[i] ~ dunif(0.001, 4)
    prec_tau[i] <- pow(tau[i], -2)
  }
  
  l_0 ~ dnorm(-2.928, pow(0.514, -2)) 
  l_1 ~ dnorm(2.233, pow(0.286, -2))
  l_2 ~ dnorm(4.8, pow(0.486, -2))
}

"

# Standard Bayesian random-effects meta-analysis model
model_ma <- "
model {

  for (i in 1:n) {
    y[i] ~ dnorm(theta_i[i], prec_s[i])
    theta_i[i] ~ dnorm(mu, prec_tau)
  }

  mu ~ dnorm(0, 0.001)

  tau ~ dunif(0.001, 4)
  prec_tau <- pow(tau, -2)
}
"

##################################################################################################################################
#
#                                                            Help functions
#
##################################################################################################################################

# Number of studies based on RoB allocation
get_rob_studies <- function(n, type) {
  if (type == "LH") {return(c(low = n/2, unclear = 0, high = n/2))}
  
  if (type == "L_UH") {return(c(low = n/2, unclear = n/4, high = n/4))}
}

# Get the category id based on study-design and RoB
get_cat_id <- function(design, rob) {
  ifelse(design == "RCT" & rob == "low", 1,
         ifelse(design == "RCT" & rob == "unclear", 2,
                ifelse(design == "RCT" & rob == "high", 3,
                       ifelse(design == "OBS" & rob == "low", 4,
                              ifelse(design == "OBS" & rob == "unclear", 5, 6)))))
}

# Get designs names
get_cat_name <- function(design, rob) {
  cat_names <- c("low_rct", "unclear_rct", "high_rct", "low_obs", "unclear_obs", "high_obs")
  cat_names[get_cat_id(design, rob)]
}

# Generate study level data using 1-1 randomization
gen_studies <- function(n.stud, design, rob, mean_bias, sd_bias, theta = 0, tau_rct = 0.20, tau_obs = 0.30) {
  
  n.stud <- as.integer(n.stud)
  
  if (n.stud == 0) {
    return(NULL)
  }
  
  cat_id <- get_cat_id(design, rob)
  cat_name <- get_cat_name(design, rob)
  
  if (design == "RCT") {
    mc <- if (rob == "high") {
      sample(30:50, n.stud, replace = TRUE)
    } else {
      sample(50:80, n.stud, replace = TRUE)
    }
    mt <- mc
    tau_design <- tau_rct
  }
  
  if (design == "OBS") {
    mc <- if (rob == "high") {
      sample(300:700, n.stud, replace = TRUE)
    } else {
      sample(500:1000, n.stud, replace = TRUE)
    }
    mt <- mc
    tau_design <- tau_obs
  }
  
  true_theta <- rnorm(n.stud, theta, tau_design)
  
  true_bias <- rnorm(n.stud, mean = mean_bias[cat_name], sd = sd_bias[cat_name])
  
  observed_theta <- true_theta + true_bias
  
  pic <- runif(n.stud, 0.20, 0.45)
  pit <- pic * (exp(observed_theta) / (1 - pic + pic * exp(observed_theta)))
  
  r.c <- rbinom(n.stud, mc, pic)
  r.e <- rbinom(n.stud, mt, pit)
  
  data.frame(
    n.e = mt,
    n.c = mc,
    r.e = r.e,
    r.c = r.c,
    true_theta = true_theta,
    true_bias = true_bias,
    observed_theta = observed_theta,
    design = rep(design, n.stud),
    Risk = rep(rob, n.stud),
    cat_id = rep(cat_id, n.stud),
    cat_name = rep(cat_name, n.stud)
  )
}

# Generate study level data using 1-2 randomization for observational studies (#control = 2*#Treatment)
gen_studies_12obs <- function(n.stud, design, rob, mean_bias, sd_bias, theta = 0, tau_rct = 0.20, tau_obs = 0.30) {
  
  n.stud <- as.integer(n.stud)
  
  if (n.stud == 0) {
    return(NULL)
  }
  
  cat_id <- get_cat_id(design, rob)
  cat_name <- get_cat_name(design, rob)
  
  if (design == "RCT") {
    mt <- if (rob == "high") {
      sample(30:50, n.stud, replace = TRUE)
    } else {
      sample(50:80, n.stud, replace = TRUE)
    }
    mc <- mt
    tau_design <- tau_rct
  }
  
  if (design == "OBS") {
    mt <- if (rob == "high") {
      sample(300:700, n.stud, replace = TRUE)
    } else {
      sample(500:1000, n.stud, replace = TRUE)
    }
    mc <- 2*mt
    tau_design <- tau_obs
  }
  
  true_theta <- rnorm(n.stud, theta, tau_design)
  
  true_bias <- rnorm(n.stud, mean = mean_bias[cat_name], sd = sd_bias[cat_name])
  
  observed_theta <- true_theta + true_bias
  
  pic <- runif(n.stud, 0.20, 0.45)
  pit <- pic * (exp(observed_theta) / (1 - pic + pic * exp(observed_theta)))
  
  r.c <- rbinom(n.stud, mc, pic)
  r.e <- rbinom(n.stud, mt, pit)
  
  data.frame(
    n.e = mt,
    n.c = mc,
    r.e = r.e,
    r.c = r.c,
    true_theta = true_theta,
    true_bias = true_bias,
    observed_theta = observed_theta,
    design = rep(design, n.stud),
    Risk = rep(rob, n.stud),
    cat_id = rep(cat_id, n.stud),
    cat_name = rep(cat_name, n.stud)
  )
}

# Get scenarios labels
make_scenario_label <- function(scen) {paste0("rand=", scen$rand, "_rate=", scen$n_rct_rate, "_ROB=", scen$rob_structure, "_Bias=", scen$bias_pattern)}

# Get MCMC performance metrics for parameter μ. The function returns: Bias, Coverage, CI length
get_mu_metrics <- function(fit, theta_true, alpha = 0.05) {
  s <- fit$BUGSoutput$summary
  
  mu_med <- s["mu", "50%"] # logOR Median
  mu_l <- s["mu", "2.5%"] # logOR LB
  mu_u <- s["mu", "97.5%"] #logOR UB
  
  target_or <- exp(theta_true) # true OR
  l_or <- exp(mu_l) # LB of OR
  u_or <- exp(mu_u) # UB of OR
  
  data.frame(bias_or = target_or - exp(mu_med), coverage = theta_true >= mu_l & theta_true <= mu_u, CIw = u_or - l_or, Rhat = s["mu", "Rhat"])
}

# Boxplots with performance metrics of the simulation's results
performance_boxplot <- function(data, metric, title, ylab, hline = NULL) {
  
  p <- ggplot(data, aes(x = method, y = .data[[metric]], fill = method, alpha = factor(rand), group = interaction(method, rand))) +
    geom_boxplot( position = position_dodge2(width = 0.8, preserve = "single"), width = 0.70, outlier.alpha = 0.20, outlier.size = 0.5) +
    facet_wrap(bias_pattern ~ index, ncol = 6, scales = "free_y") +
    scale_fill_manual(values = method_cols, drop = FALSE) +
    scale_alpha_manual(values = c(0.4, 1), name = "Randomization", labels = function(x) gsub("-", ":", x)) +
    guides(fill = "none", alpha = guide_legend(override.aes = list(fill = "grey40"))) +
    labs(title = title, x = NULL, y = ylab) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom", text = element_text(size = 16), 
          axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
          strip.text = element_text(face = "bold"), plot.title = element_text(face = "bold"))
  
  if (!is.null(hline)) {p <- p + geom_hline(yintercept = hline, linetype = "dashed", colour = "red")}
  
  p
}

# Colors for the model used
method_cols <- c("Bias-adjusted SSVS" = "#1B9E77", "Naive: all studies" = "#D95F02", "Naive: RCT only" = "#7570B3", "Naive: OBS only" = "#E7298A")
