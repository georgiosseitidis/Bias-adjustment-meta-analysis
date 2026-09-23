# Random-Effects meta-analysis model
model_ma <- "
model {
  # Likelihood
  for (i in 1:n) {
    y[i] ~ dnorm(theta[i], prec_s[i])
    theta[i] ~ dnorm(mu, prec_tau)
    w[i] <- pow(pow(prec_s[i], -1) + pow(prec_tau, -1), -1) 
  }

  # Priors
  mu ~ dnorm(0, 0.001)
  tau ~ dunif(0, 5)
  prec_tau <- pow(tau, -2)
  
  for(i in 1:n){w_star[i] <- w[i]*pow(sum(w), -1)}
  w_NRE <- sum(w_star[1:16])
}
"

# Random-Effects meta-analysis model (no weights)
model_ma_2 <- "
model {
  # Likelihood
  for (i in 1:n) {
    y[i] ~ dnorm(theta[i], prec_s[i])
    theta[i] ~ dnorm(mu, prec_tau)
  }

  # Priors
  mu ~ dnorm(0, 0.001)
  tau ~ dunif(0, 5)
  prec_tau <- pow(tau, -2)
}
"

# Bias adjusted meta-analysis model
model_ma_bias <- "
model {
  # Likelihood
  for (i in 1:n) {
    y[i] ~ dnorm(theta[i], prec_s[i])
    theta[i] ~ dnorm(mu + b[RoB[i]], prec_tau[design[i]])
    
    # weights
    w[i] <- pow(pow(prec_s[i], -1) + pow(prec_tau[design[i]], -1) + (1-g[RoB[i]])*psi[RoB[i]]*psi[RoB[i]] + g[RoB[i]]*c[RoB[i]]*c[RoB[i]]*psi[RoB[i]]*psi[RoB[i]], -1)
    
    w_ds[i] <- pow(pow(prec_s[i], -1) + pow(prec_tau[design[i]], -1), -1)  # design-specific weights

    D_ar[i] <- (1-g[RoB[i]])*psi[RoB[i]]*psi[RoB[i]] + g[RoB[i]]*(c[RoB[i]]*c[RoB[i]]*psi[RoB[i]]*psi[RoB[i]])
    D[i] <- 1- D_ar[i]*w[i]
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
  
  # Heterogeneity 1 RCT, 2 OBS
  for(i in 1:2){
    tau[i] ~ dunif(0.001, 4)
    prec_tau[i] <- pow(tau[i], -2)
  }
  
  l_0 ~ dnorm(-2.928, pow(0.514, -2)) 
  l_1 ~ dnorm(2.233, pow(0.286, -2))
  l_2 ~ dnorm(4.8, pow(0.486, -2))
  
  # Calculation of Weights D_NRE
  for(i in 1:n){w_star[i] <- w[i]*pow(sum(w), -1)}
  w_NRE <- sum(w_star[1:16])
  
  D_NRE_bias <- sum(w[1:16])*pow(sum(w_ds[1:16]),-1)
}

"
############################################################################################################################################
#
#                                                               Help functions
#
############################################################################################################################################

get_cat_id <- function(design, rob) {
  ifelse(design == "RCT" & rob == "Low", 1,
         ifelse(design == "RCT" & rob == "Unclear", 2,
                ifelse(design == "RCT" & rob == "High", 3,
                       ifelse(design == "OBS" & rob == "Low", 4,
                              ifelse(design == "OBS" & rob == "Unclear", 5, 6)))))
}

# Forest plot comparing the estimates from the MCMC run of the employed MA models 
# d: data.frame (RCTs must be at the bottom)
forest_fit_all <- function(fit_bias, fit_naive, fit_ma_rct, fit_ma_os, d){
  
  n_studies = dim(d)[1] #Total number of studies
  
  ###################################################
  # posterior means and intervals for the studies
  ###################################################
  
  ############################
  # study specific estimates
  ############################
  
  # Bias model 
  fit_summary <- fit$BUGSoutput$summary
  adjusted_mean <- fit_summary[grep("theta", row.names(fit_summary)), "50%"]
  adjusted_lower <- fit_summary[grep("theta", row.names(fit_summary)), "2.5%"]
  adjusted_upper <- fit_summary[grep("theta", row.names(fit_summary)), "97.5%"]
  
  # naive MA model
  fit_summary_naive <- fit_naive$BUGSoutput$summary
  ma_mean_naive <- fit_summary_naive[grep("theta", row.names(fit_summary_naive)), "50%"]
  ma_lower_naive <- fit_summary_naive[grep("theta", row.names(fit_summary_naive)), "2.5%"]
  ma_upper_naive <- fit_summary_naive[grep("theta", row.names(fit_summary_naive)), "97.5%"]
  
  # MA model for RCT
  fit_summary_ma_RCT <- fit_ma_rct$BUGSoutput$summary
  ma_mean_rct <- fit_summary_ma_RCT[grep("theta", row.names(fit_summary_ma_RCT)), "50%"]
  ma_lower_rct <- fit_summary_ma_RCT[grep("theta", row.names(fit_summary_ma_RCT)), "2.5%"]
  ma_upper_rct <- fit_summary_ma_RCT[grep("theta", row.names(fit_summary_ma_RCT)), "97.5%"]
  
  
  # MA model for Observational
  fit_summary_ma_os <- fit_ma_os$BUGSoutput$summary
  ma_mean_os <- fit_summary_ma_os[grep("theta", row.names(fit_summary_ma_os)), "50%"]
  ma_lower_os <- fit_summary_ma_os[grep("theta", row.names(fit_summary_ma_os)), "2.5%"]
  ma_upper_os <- fit_summary_ma_os[grep("theta", row.names(fit_summary_ma_os)), "97.5%"]
  
  
  # Study labels
  study_labels <- as.character(d$name)
  
  ma_mean <- c(ma_mean_os, ma_mean_rct)
  ma_lower <- c(ma_lower_os, ma_lower_rct)
  ma_upper <- c(ma_upper_os, ma_upper_rct)
  
  #########################################
  # Bias adjusted Meta-analysis
  ##########################################
  mu_summary <- c(mean =  fit_summary[grep("mu", row.names(fit_summary)), "50%"], 
                  `2.5%` = fit_summary[grep("mu", row.names(fit_summary)), "2.5%"], 
                  `97.5%` = fit_summary[grep("mu", row.names(fit_summary)), "97.5%"])
  
  # Raw naive meta-analysis
  mu_summary_ma_naive <- c(mean =  fit_summary_naive[grep("mu", row.names(fit_summary_naive)), "50%"], 
                           `2.5%` = fit_summary_naive[grep("mu", row.names(fit_summary_naive)), "2.5%"], 
                           `97.5%` = fit_summary_naive[grep("mu", row.names(fit_summary_naive)), "97.5%"])
  
  # Raw meta-analysis using RCTs
  mu_summary_ma_rct <- c(mean =  fit_summary_ma_RCT[grep("mu", row.names(fit_summary_ma_RCT)), "50%"], 
                         `2.5%` = fit_summary_ma_RCT[grep("mu", row.names(fit_summary_ma_RCT)), "2.5%"], 
                         `97.5%` = fit_summary_ma_RCT[grep("mu", row.names(fit_summary_ma_RCT)), "97.5%"])
  
  # Raw meta-analysis using Observations
  mu_summary_ma_os <- c(mean =  fit_summary_ma_os[grep("mu", row.names(fit_summary_ma_os)), "50%"], 
                        `2.5%` = fit_summary_ma_os[grep("mu", row.names(fit_summary_ma_os)), "2.5%"], 
                        `97.5%` = fit_summary_ma_os[grep("mu", row.names(fit_summary_ma_os)), "97.5%"])
  
  # -------------------------------
  # Create plotting dataframe
  # -------------------------------
  
  # Main study rows (raw and adjusted stacked)
  plot_df <- data.frame(
    Study = c(study_labels, study_labels),
    Type = c(rep("Adjusted", times = n_studies), rep("Only non-RCT", times = n_studies- sum(d$RCT) ), rep("Only RCT", times = sum(d$RCT) )),
    y_pos = c(1:n_studies, 1:n_studies),
    Effect = c(adjusted_mean, ma_mean),
    Lower = c(adjusted_lower, ma_lower),
    Upper = c(adjusted_upper, ma_upper),
    Design = plyr::mapvalues(c(d$RCT, d$RCT), c(1, 0), c("RCT", "Non-RCT")),
    RoB = c(d$Risk, d$Risk)
  )
  
  plot_df_naive <- data.frame(
    Study = c(study_labels),
    Type = c(rep("Naive MA", times = n_studies) ),
    y_pos = c(1:n_studies),
    Effect = c(ma_mean_naive),
    Lower = c(ma_lower_naive),
    Upper = c(ma_upper_naive),
    Design = plyr::mapvalues(c(d$RCT), c(1, 0), c("RCT", "Non-RCT")),
    RoB = c(d$Risk)
  )  
  
  plot_df <- rbind(plot_df, plot_df_naive) 
  plot_df <- plot_df[order(plot_df$y_pos), ]
  
  # Meta-analysis rows
  meta_df <- data.frame(
    Study = rep("Meta-analysis", 4),
    Type = c("Adjusted", "Naive MA",  "Only non-RCT", "Only RCT"),
    y_pos = rep(dim(d)[1] + 1, 4),
    Effect = c(mu_summary["mean"], mu_summary_ma_naive["mean"],  mu_summary_ma_os["mean"], mu_summary_ma_rct["mean"]),
    Lower = c(mu_summary["2.5%"], mu_summary_ma_naive["2.5%"], mu_summary_ma_os["2.5%"], mu_summary_ma_rct["2.5%"]),
    Upper = c(mu_summary["97.5%"], mu_summary_ma_naive["97.5%"], mu_summary_ma_os["97.5%"], mu_summary_ma_rct["97.5%"]),
    Design = c("-", "-", "-", "-"),
    RoB = c("-", "-", "-", "-")
  )
  
  # Combine everything
  library(dplyr)
  
  plot_df_full <- bind_rows(plot_df, meta_df)
  
  # Fix order of factor levels for plotting
  plot_df_full$Study <- factor(plot_df_full$Study, levels = rev(c(as.character(study_labels), "Meta-analysis")))
  
  # Reassign y-positions based on correct ordering
  plot_df_full$y_pos <- as.numeric(plot_df_full$Study)
  
  # -------------------------------
  # 4. Plot: full forest plot
  # -------------------------------
  
  plot_df_full$Effect <- exp(plot_df_full$Effect)
  plot_df_full$Lower <- exp(plot_df_full$Lower)
  plot_df_full$Upper <- exp(plot_df_full$Upper)
  
  ggplot(plot_df_full, aes(x = Effect, y = Study , color = Type)) +
    geom_point(size = 3, position = position_dodge(width = 0.5)) +
    geom_errorbarh(aes(xmin = Lower, xmax = Upper), position = position_dodge(width = 0.5), height = 0.2) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
    scale_color_manual(values = c("Only non-RCT" = "red", "Naive MA"="purple", "Only RCT" = "blue", "Adjusted" = "darkgreen")) +
    scale_x_log10()+
    labs(
      title = "",
      x = "Odds Ratio",
      y = NULL
    ) +
    theme_minimal()+
    theme(
      legend.position = "bottom",
      axis.text.x = element_text(size = 13, color = "black"),
      axis.text.y = element_text(size = 13, color = "black"),
      legend.text = element_text(size = 13, color = "black"),
      axis.title = element_text(size = 14, color = "black"),
      legend.title = element_text(size = 14, color = "black"),
      panel.grid.minor = element_blank()
    )
  
}

# Lollipop plot of the posterior inclusion probabilities
lolipo_pip <- function(fit){
  
  # Extract posterior matrix
  samples <- as.mcmc(fit)
  post_mat <- do.call(rbind, samples)
  gammas <- post_mat[, grep("gamma", colnames(post_mat))]
  
  PIP <- apply(gammas, 2, mean)
  s_PIP <- as.integer(gsub("gamma\\[|\\]", "", labels(PIP)))
  
  df_PIP <- data.frame(PIP, s_PIP)
  df_PIP <- df_PIP[order(df_PIP$s_PIP), ]
  
  df_PIP$study <- c("RCT-Low", "RCT-Unclear", "RCT-High", "OBS-Low", "OBS-Unclear", "OBS-High")
  df_PIP$RCT <- c(rep("RCT", 3), rep("OBS", 3))
  df_PIP$ROB <- c("Low", "Unclear", "High", "Low", "Unclear", "High")
  
  
  # Lolipo plot
  ggplot(df_PIP[-which(df_PIP$ROB == "Unclear" | (df_PIP$ROB != "High" & df_PIP$RCT == "OBS") ), ], aes(x = study, y = PIP, color = ROB)) +
    geom_segment(aes(xend = study, y = 0, yend = PIP), color = "gray60") +
    geom_point(size = 5) +
    facet_wrap(~ RCT, ncol = 1, scales = "free_y") +
    coord_flip() +
    labs(
      title = "",
      x = NULL,
      y = "Posterior Inclusion Probability of Bias term"
    ) +
    theme_minimal(base_size = 14)+
    theme(legend.position = "bottom")
}

# Spike and Slab for bias term
spike_slab <- function(fit){
  
  # Extract posterior samples from fit
  samples <- as.mcmc(fit)
  post_mat <- do.call(rbind, samples)
  
  # Detect how many b[i] terms there are
  b_names <- grep("^b\\[", colnames(post_mat), value = TRUE)
  g_names <- gsub("b", "gamma", b_names)  # match gamma[i]
  
  # Extract numeric indices for ordering
  study_ids <- as.integer(gsub("[^0-9]", "", b_names))
  
  # Initialize empty list to collect data
  density_list <- list()
  
  # Loop safely through studies
  for (i in seq_along(b_names)) {
    b_name <- b_names[i]
    g_name <- g_names[i]
    study_label <- paste0("Study ", study_ids[i])
    
    b_vals <- post_mat[, b_name]
    g_vals <- post_mat[, g_name]
    
    # Initialize empty df for this study
    df_study <- data.frame()
    
    if (sum(g_vals == 0) > 10) {
      d0 <- density(b_vals[g_vals == 0])
      df0 <- data.frame(x = d0$x, y = d0$y, Gamma = "γ = 0", Study = study_label)
      df_study <- bind_rows(df_study, df0)
    }
    
    if (sum(g_vals == 1) > 10) {
      d1 <- density(b_vals[g_vals == 1])
      df1 <- data.frame(x = d1$x, y = d1$y, Gamma = "γ = 1", Study = study_label)
      df_study <- bind_rows(df_study, df1)
    }
    
    density_list[[i]] <- df_study
  }
  
  # Combine all studies into one dataframe
  density_df <- bind_rows(density_list)
  density_df$Study <- factor(density_df$Study, levels = unique(density_df$Study))
  
  study_info <- data.frame(
    StudyID = paste0("Study ", 1:6),
    Design = c(rep("RCT", 3), rep("OBS", 3)),
    RoB = c("Low", "Unclear", "High", "Low", "Unclear", "High")
  )
  
  # Combine info for new label
  study_info <- study_info %>%
    mutate(Label = paste0(StudyID, " (", Design, ", ", RoB, ")"))
  
  density_df <- density_df %>%
    left_join(study_info, by = c("Study" = "StudyID")) %>%
    mutate(FacetLabel = factor(Label, 
                               levels = study_info$Label, 
                               label = c("RCT - Low", "RCT - Unclear", "RCT - High",  "OBS - Low", "OBS - Unclear", "OBS - High" )))
  
  # here we do not have RCT of unclear ROB and OBS of low-unclear (exclude them)
  ggplot(density_df[-which(density_df$Study %in% paste0("Study ", c(2,4,5))), ], aes(x = x, y = y, fill = Gamma, color = Gamma)) +
    geom_area(alpha = 0.2, position = "identity") +
    geom_line(size = 1) +
    facet_wrap(~ FacetLabel, scales = "free") +
    scale_fill_manual(values = c("γ = 0" = "#377eb8", "γ = 1" = "#4daf4a")) +
    scale_color_manual(values = c("γ = 0" = "#377eb8", "γ = 1" = "#4daf4a")) +
    labs(
      title = expression("Posterior Densities of " * b[i] * " Conditional on " * gamma[i]),
      x = expression(b[i]),
      y = "Density",
      fill = "Condition",
      color = "Condition"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      legend.position = "top",
      strip.text = element_text(face = "bold"),
      axis.title = element_text(face = "bold")
    )
  
}
