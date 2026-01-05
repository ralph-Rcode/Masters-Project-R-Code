
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)



# Generates the shape of a single invasion trajectory
get_flicker_traj <- function(b, c, I, lambda, R_GL, R_LT_eq, p0=0.001) {
  p <- numeric(2500)
  p[1] <- p0
  for(t in 1:2499) {
    # Hierarchy: R_total = R_GL + (1 - R_GL) * R_LT(t)
    R_LT <- R_LT_eq * (1 - exp(-lambda * (t - 1)))
    R_tot <- R_GL + (1 - R_GL) * R_LT
    
    p_now <- p[t]
    p_star_C <- p_now + (1 - p_now) * R_tot
    p_star_A <- p_now * (1 - R_tot)
    
    # Fitness logic (Additive)
    wC <- 1 + b * I * (1 - p_star_C)
    wA <- 1 - c * I + b * I * (1 - p_star_A)
    w_bar <- p_now * wC + (1 - p_now) * wA
    
    p[t+1] <- max(0, min(1, p_now * wC / w_bar))
    if(p[t+1] < 1e-9) break
  }
  return(p[p > 0])
}

# Poisson Superposition for Mean Load, CV, and State Occupancy
#Gemini AI was used to debug the rpois function here
run_poisson_analysis <- function(mu, b, c, I, lambda, R_GL, R_LT_eq, T_total = 1000000) {
  traj <- get_flicker_traj(b, c, I, lambda, R_GL, R_LT_eq)
  len_traj <- length(traj)
  
  # Ensure mu is numeric to avoid rpois error
  hits <- rpois(T_total, lambda = as.numeric(mu))
  hit_times <- which(hits > 0)
  
  total_p <- numeric(T_total + len_traj)
  lineage_counts <- integer(T_total + len_traj)
  
  for (i in seq_along(hit_times)) {
    t_start <- hit_times[i]
    num_mutations <- hits[t_start]
    t_end <- t_start + len_traj - 1
    total_p[t_start:t_end] <- total_p[t_start:t_end] + (traj * num_mutations)
    lineage_counts[t_start:t_end] <- lineage_counts[t_start:t_end] + num_mutations
  }
  
  final_p <- total_p[1:T_total]
  final_counts <- lineage_counts[1:T_total]
  
  # 1. Summary Stats
  m_p <- mean(final_p)
  cv_p <- if(m_p > 0) sd(final_p) / m_p else 0
  
  # 2. State Binning (0 to 6+)
  state_table <- table(factor(final_counts, levels = 0:20))
  props <- as.data.frame(state_table / T_total)
  colnames(props) <- c("Category", "Prop")
  props <- props %>%
    mutate(Category = as.character(Category)) %>%
    mutate(Category = ifelse(as.numeric(Category) >= 6, "6+", Category)) %>%
    group_by(Category) %>%
    summarise(Prop = sum(Prop)) %>%
    ungroup()
  
  return(list(stats = data.frame(Mean_p = m_p, CV = cv_p), occupancy = props))
}

#plotting

plot_sensitivity_stats <- function(df, param_name, param_expr) {
  p1 <- ggplot(df, aes(x = .data[[param_name]], y = Mean_p)) +
    geom_line(color = "steelblue", size = 1) + geom_point(color = "steelblue", size = 3) +
    labs(title = "Effect on Mean Cheater Frequency", x = param_expr, y = expression(paste("Mean frequency ", hat(p)))) +
    theme_minimal() + theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))
  
  p2 <- ggplot(df, aes(x = .data[[param_name]], y = CV)) +
    geom_line(color = "firebrick3", size = 1) + geom_point(color = "firebrick3", size = 3) +
    labs(title = "Effect on Coefficient of Variation", x = param_expr, y = "CV") +
    theme_minimal() + theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))
  return(p1 / p2)
}

plot_state_transition <- function(df, x_var, x_label) {
  df$Category <- factor(df$Category, levels = c("6+", "5", "4", "3", "2", "1", "0"))
  ggplot(df, aes(x = .data[[x_var]], y = Prop, fill = Category)) +
    geom_area(alpha = 0.9, color = "white", size = 0.1) +
    scale_fill_manual(values = c("0" = "#f7fbff", "1" = "#deebf7", "2" = "#c6dbef", 
                                 "3" = "#9ecae1", "4" = "#6baed6", "5" = "#3182bd", "6+" = "#08519c")) +
    labs(title = "Number of Segregating Cheater Alleles Across the Metapopulation", x = x_label, y = "Proportion of Time", fill = "Active Lineages") +
    theme_minimal() + theme(panel.border = element_rect(colour = "black", fill = NA, size = 1), panel.grid = element_blank())
}

#baseline parameters
BASE <- list(mu=0.006, lambda=0.006, I=0.5, R_LT_eq=0.5, R_GL=0.15, b=4, c=1)

# Function to run both analyses for a given parameter range
run_full_sweep <- function(param_name, vals, param_expr) {
  res_stats <- data.frame()
  res_states <- data.frame()
  
  for(v in vals) {
    current <- BASE
    current[[param_name]] <- v
    out <- run_poisson_analysis(mu=current$mu, b=current$b, c=current$c, I=current$I, 
                                lambda=current$lambda, R_GL=current$R_GL, R_LT_eq=current$R_LT_eq)
    # Collect Stats
    row_stats <- out$stats
    row_stats[[param_name]] <- v
    res_stats <- rbind(res_stats, row_stats)
    # Collect States
    row_states <- out$occupancy
    row_states[[param_name]] <- v
    res_states <- rbind(res_states, row_states)
  }
  
  p_stats <- plot_sensitivity_stats(res_stats, param_name, param_expr)
  p_states <- plot_state_transition(res_states, param_name, param_expr)
  return(list(p_stats = p_stats, p_states = p_states))
}

#run all sensitivity profiles

# A: Mutation Hit Rate (mu)
sweep_mu <- run_full_sweep("mu", seq(0.001, 0.02, length.out = 20), 
                           expression(paste("Mutation Rate (", mu, ")")))

# B: Lag Rate (lambda)
sweep_lambda <- run_full_sweep("lambda", seq(0.003, 0.015, length.out = 20), 
                               expression(paste("Rate of Structural Accumulation (", lambda, ")")))

# C: Investment Level (I)
sweep_I <- run_full_sweep("I", seq(0.1, 1.0, length.out = 20), 
                          "Investment Level (I)")

# D: Equilibrium Structure (R_LT_eq)
sweep_RLT <- run_full_sweep("R_LT_eq", seq(0.3, 0.9, length.out = 14), 
                            expression(R[LT]^eq))

# E: Local Relatedness (R_GL)
sweep_RGL <- run_full_sweep("R_GL", seq(0.05, 0.25, length.out = 20), 
                            expression(R[GL]))

# F: Benefit (b)
sweep_b <- run_full_sweep("b", seq(2, 6, length.out = 9), 
                          "Benefit (b)")

# G: Cost (c)
sweep_c <- run_full_sweep("c", seq(0.5, 2.0, length.out = 10), 
                          "Cost (c)")

#run results



# lag rate sweep
sweep_lambda$p_stats
sweep_lambda$p_states

# mutation rate
sweep_mu$p_stats
sweep_mu$p_states



#investment
sweep_I$p_stats
sweep_I$p_states

#R_GL,b,c 
sweep_RGL$p_states
sweep_b$p_states
sweep_c$p_states

#R_LT
sweep_RLT$p_states
