###############################################################################################################################
#
#                                                       SIMULATION SETUP
#
################################################################################################################################

# Load libraries
library(R2jags)
library(doParallel)
library(foreach)

# Load Meta-analysis models and help functions
source("MA Models & Help Functions.R")

# Simulation setup
nsim <- 1000 # number of iterations
theta <- 0 # theta value
N <- 24 # Number of studies

n.iter <- 100000 # MCMC iterations
n.burnin <- 20000 # MCMC burn-in
n.chains <- 2 # number of chains
n.thin <- 2 # thinning

################################################################################################################################
#
#                                                                 Set simulations scenarios
#
################################################################################################################################


# Bias scenarios
bias_patterns <- data.frame( bias_pattern = c("No Bias", "NRE moderate", "NRE strong", "Rct-NRE Moderate", "Rct20_NRE20-35"),
                             low_rct = c(0, 0, 0, 0, log(1.20)),
                             unclear_rct = c(0, log(1.05), log(1.05), log(1.10), log(1.20)),
                             high_rct = c(0, log(1.10), log(1.10), log(1.15), log(1.20)),
                             low_obs = c(0, log(1.05), log(1.05), log(1.05), log(1.20)),
                             unclear_obs = c(0, log(1.10), log(1.20), log(1.20), log(1.20)),
                             high_obs = c(0, log(1.15), log(1.35), log(1.35), log(1.35)),
                             stringsAsFactors = FALSE)
# Simulation scenarios
scenarios <- expand.grid(n_rct_rate = c(1.5, 2, 3),
                         rob_structure = c("LH", "L_UH"),
                         bias_pattern = bias_patterns$bias_pattern,
                         stringsAsFactors = FALSE)
# Add bias values
scenarios <- merge(scenarios, bias_patterns, by = "bias_pattern", all.x = TRUE)

scenarios$rob_structure <- factor(scenarios$rob_structure, levels = c("LH", "L_UH"))
scenarios$bias_pattern <- factor(scenarios$bias_pattern, levels = bias_patterns$bias_pattern)

# Order by RCT rate, ROB & Bias pattern
scenarios <- scenarios[order(scenarios$n_rct_rate, scenarios$rob_structure, scenarios$bias_pattern), ]
rownames(scenarios) <- NULL

# Add randomization patterns (1-1, 1-2 for OBS)
scenarios <- rbind(scenarios, scenarios)
scenarios <- data.frame("rand" = c(rep("1-1", 30), rep("1-2", 30)), scenarios)  

# SD used for the bias terms
sd_bias_default <- c(low_rct = 0.005, unclear_rct = 0.03, high_rct = 0.05, low_obs = 0.03, unclear_obs = 0.06, high_obs = 0.08)

################################################################################################################################
#
#                                                       Parallel simulation
#
################################################################################################################################

n_cores <- parallel::detectCores()  # number of parallel cores

cl <- makeCluster(n_cores)
registerDoParallel(cl)

out_dir <- "simulation_results" # folder for saving simulations results
dir.create(out_dir, showWarnings = FALSE)

# Run simulation
results <- foreach(s = seq_len(nrow(scenarios)), .combine = rbind, .packages = c("R2jags", "meta", "netmeta")) %dopar% {
  
  # Simulated Scenario 
  scen <- scenarios[s, ]
  
  # Scenario file
  scenario_file <- file.path(out_dir, paste0("scenario_", s, ".csv"))
  
  scenario_results <- vector("list", nsim)
  
  # Mean bias for the simulated scenario
  mean_bias_scenario <- c(low_rct = scen$low_rct, unclear_rct = scen$unclear_rct, high_rct = scen$high_rct, low_obs = scen$low_obs, unclear_obs = scen$unclear_obs, high_obs = scen$high_obs)
  
  # Calculate the Number of RCT and OBS studies
  n_rct <- N / scen$n_rct_rate # RCTs
  n_obs <- N - n_rct # OBS studies
  
  # RoB allocation to the studies
  rct_stud <- get_rob_studies(n_rct, scen$rob_structure) #RCTs
  obs_stud <- get_rob_studies(n_obs, scen$rob_structure) #OBSs
  
  for (j in seq_len(nsim)) {
    
    set.seed(nsim + s * 10000 + j)
    
    # Generate RCTs
    rcts <- lapply(names(rct_stud), function(rob) {gen_studies(n.stud = rct_stud[rob], design = "RCT", rob = rob, mean_bias = mean_bias_scenario, sd_bias = sd_bias_default, theta = theta, tau_rct = 0.20, tau_obs = 0.30)})
    
    # Generate OBSs
    if(scen$rand == "1-1"){
      obs <- lapply(names(obs_stud), function(rob) {gen_studies(n.stud = obs_stud[rob], design = "OBS", rob = rob, mean_bias = mean_bias_scenario, sd_bias = sd_bias_default, theta = theta, tau_rct = 0.20, tau_obs = 0.30)})
    }else{
      obs <- lapply(names(obs_stud), function(rob) {gen_studies_12obs(n.stud = obs_stud[rob], design = "OBS", rob = rob, mean_bias = mean_bias_scenario, sd_bias = sd_bias_default, theta = theta, tau_rct = 0.20, tau_obs = 0.30)})
    }
    
    # Combine the studies
    df <- do.call(rbind, c(rcts, obs))
    df$id <- seq_len(nrow(df))
    
    # Calculate Treatment Effects
    z <- pairwise(treat = list(t1 = rep("Exp", nrow(df)), t2 = rep("Con", nrow(df))), event = list(r.e, r.c), n = list(n.e, n.c), studlab = id, sm = "OR", data = df)
    
    # Data frame containing study-level data
    d <- data.frame(TE = z$TE, seTE = z$seTE, n.v = z$n1, n.c = z$n2, design = df$design, Risk = df$Risk, 
                    true_theta = df$true_theta, true_bias = df$true_bias, observed_theta = df$observed_theta, 
                    true_bias = df$true_bias, cat_id = df$cat_id)
    
    d$RCT <- ifelse(d$design == "RCT", 1, 0)
    d$design_id <- ifelse(d$design == "RCT", 1, 2)
    d$ROB_Mod <- ifelse(d$Risk == "unclear", 1, 0)
    d$ROB_High <- ifelse(d$Risk == "high", 1, 0)
    
    #########################################
    #           Bias-adusted model
    ##########################################
    
    # MCMC data for the Bias adjusted model
    data_jags_ssvs <- list(y = d$TE, 
                           prec_s = 1/d$seTE^2,
                           n = nrow(d),
                           n_des = 6,
                           RoB = d$cat_id,
                           RCT = c(1, 1, 1, 0 ,0 , 0),
                           design = plyr::mapvalues(d$RCT, from = c(0, 1), to = c(2, 1)),
                           ROB_Mod = c(0, 1, 0, 0, 1, 0),
                           ROB_High = c(0, 0, 1, 0, 0, 1))
    
    
    # Fit the Bias-adusted model
    fit_ssvs <- jags(data = data_jags_ssvs, parameters.to.save = "mu", model.file = textConnection(model_ma_bias_adj),
                     n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin, n.thin =  n.thin, DIC = FALSE)
    
    #########################################
    #           Naive meta-analysis models
    ##########################################
    
    # MCMC data for the naive meta-analysis model
    data_jags_all <- list(y = d$TE, prec_s = 1 / d$seTE^2, n = nrow(d))
    
    # Fit naive meta-analysis
    fit_ma_all <- jags(data = data_jags_all, parameters.to.save = "mu", model.file = textConnection(model_ma), 
                       n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin, n.thin =  n.thin, DIC = FALSE)
    
    # MCMC data for the sensitivity RCT meta-analysis model
    data_jags_rct <- list(y = d$TE[d$RCT == 1], prec_s = 1 / d$seTE[d$RCT == 1]^2, n = sum(d$RCT == 1))
    
    # Fit naive RCT sensitivity
    fit_ma_rct <- jags(data = data_jags_rct, parameters.to.save = "mu", model.file = textConnection(model_ma),
                       n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin, n.thin =  n.thin, DIC = FALSE)
    
    # MCMC data for the sensitivity OBS meta-analysis model
    data_jags_obs <- list(y = d$TE[d$RCT == 0], prec_s = 1 / d$seTE[d$RCT == 0]^2, n = sum(d$RCT == 0))
    
    # Fit naive OBS sensitivity
    fit_ma_obs <- jags(data = data_jags_obs, parameters.to.save = "mu", model.file = textConnection(model_ma),
                       n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin, n.thin =  n.thin, DIC = FALSE)
    
    # Get models' estimates
    perf_ssvs <- get_mu_metrics(fit_ssvs, theta) # bias-adjusted
    perf_naive <- get_mu_metrics(fit_ma_all, theta) # naive MA
    perf_rct <- get_mu_metrics(fit_ma_rct, theta) # RCT sensitivity
    perf_obs <- get_mu_metrics(fit_ma_obs, theta) # OBS scensitivity
    
    
    scenario_label <- make_scenario_label(scen)
    
    # Export results
    iter_export <- data.frame(sim = j, # Simulation id
                              scenario = scenario_label, # id label
                              
                              # Simulation Argument set-up
                              Rand = scen$rand, # Randomisation setup
                              Nrct_rate = scen$n_rct_rate, # RCT rate
                              ROB_srt = scen$rob_structure, # ROB allocation
                              bias_pattern = as.character(scen$bias_pattern), # Bias scenario

                              # bias arguments set up
                              low_rct = scen$low_rct, # RCT of low RoB
                              unclear_rct = scen$unclear_rct, # RCT of unclear RoB
                              high_rct = scen$high_rct, # RCT of High RoB
                              low_obs = scen$low_obs, # OBS of low RoB
                              unclear_obs = scen$unclear_obs, # OBS of unclear RoB
                              high_obs = scen$high_obs, # OBS of High RoB
                              
                              mean_true_bias = mean(d$true_bias), # Average bias assigned to studies
                              mean_abs_true_bias = mean(abs(d$true_bias)), # Average absolute bias assigned to studies
                              
                              # Models' Perfomance metrics
                              
                              # Bias adjusted
                              bias = perf_ssvs$bias_or, # bias estimating true OR
                              coverage = perf_ssvs$coverage, #coverage
                              CIw = perf_ssvs$CIw, # CI width
                              Rhat = perf_ssvs$Rhat, # Rhat
                              
                              # Naive meta-analysis
                              bias_naive = perf_naive$bias_or, # bias estimating true OR
                              coverage_naive = perf_naive$coverage, # coverage
                              CIw_naive = perf_naive$CIw, #CI width
                              R_naive = perf_naive$Rhat, # Rhat
                              
                              # RCT sensitivity
                              bias_rct_est = perf_rct$bias_or, # bias estimating true OR
                              coverage_rct = perf_rct$coverage, # coverage
                              CIw_rct = perf_rct$CIw, #CI width
                              R_rct = perf_rct$Rhat, # Rhat
                              
                              # OBS sensitivity 
                              bias_obs_est = perf_obs$bias_or, # bias estimating true OR
                              coverage_obs = perf_obs$coverage, # coverage
                              CIw_obs = perf_obs$CIw, #CI width,
                              R_rct = perf_obs$Rhat # Rhat
                              )
      
    
    first_write <- !file.exists(scenario_file)
    
    # Export results
    write.table(iter_export, file = scenario_file, sep = ",", row.names = FALSE, col.names = first_write, append = !first_write)
    
    scenario_results[[j]] <- iter_export
  }
  
  do.call(rbind, scenario_results)
}

stopCluster(cl)

# Export all the results
write.csv(results, file = file.path(out_dir, "simulation_results.csv"), row.names = FALSE)

# Save R-session
save.image(file.path(out_dir, "simulation_results.RData"))
