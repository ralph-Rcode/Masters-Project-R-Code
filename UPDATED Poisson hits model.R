library(ggplot2)
library(dplyr)


# this function returns the vector of frequencies for a single mutation event from its birth until it is purged (extinction).

get_additive_flicker_trajectory <- function(
    b = 4,              # Benefit
    c = 1,              # Cost
    I = 0.5,            # Investment Level
    lambda = 0.006,      # Lag rate
    R_GL = 0.15,         # Local Relatedness
    R_LT_eq = 0.5,      # Equilibrium Structure
    p0 = 0.001,         # Initial mutation size
    threshold = 1e-9    # Extinction threshold where the allele is deemed extinct
) {
  
  # We simulate enough generations to ensure extinction
  
  max_gen_guess <- 2000 
  p <- numeric(max_gen_guess)
  p[1] <- p0
  
  for(t in 1:(max_gen_guess-1)) {
    
    # 1. Calculate Time-Dependent Structure (The Age of the Lineage)
    # The 't' here represents the age of this specific mutation
    R_LT <- R_LT_eq * (1 - exp(-lambda * (t - 1)))
    R_total <- R_GL + (1 - R_GL) * R_LT
    
    # 2. Subjective Frequencies (Additive Model)
    p_now <- p[t]
    
    # Cheater sees: p + (1-p)R
    p_star_C <- p_now + (1 - p_now) * R_total
    # Altruist sees: p(1-R)
    p_star_A <- p_now * (1 - R_total)
    
    # 3. Fitness
    w_C <- 1 + b * I * (1 - p_star_C)
    w_A <- 1 - c * I + b * I * (1 - p_star_A)
    w_bar <- p_now * w_C + (1 - p_now) * w_A
    
    # 4. Update
    p_next <- (p_now * w_C) / w_bar
    
    # Check for extinction
    if(p_next < threshold) {
      p[t+1] <- 0
      break
    } else {
      p[t+1] <- p_next
    }
  }
  
  
  return(p[p > 0])
}

#Poisson superposition
# This function places the "flickers" onto a timeline based on mutation rate mu.

run_poisson_hits <- function(
    mu,                # Mutation Rate (Events per generation)
    T_total = 5000,    # Total Time Steps to simulate
    flicker_traj       # The shape of one flicker 
) {
  
  # 1. Generate Mutation Events
  # Draw number of mutations at each time step
  # (Most steps will be 0, some will be 1 or more)
  n_mutations_per_step <- rpois(T_total, lambda = mu)
  
  # Identify the times where mutations occurred
  mutation_times <- which(n_mutations_per_step > 0)
  mutation_counts <- n_mutations_per_step[mutation_times]
  
  # 2. Superposition (Add them up)
  # Initialize the global frequency vector
  total_p <- numeric(T_total + length(flicker_traj)) 
  
  for (i in seq_along(mutation_times)) {
    start_t <- mutation_times[i]
    count <- mutation_counts[i] # In case 2 mutations happen at once
    
    # Add the flicker trajectory to the timeline at the correct start time
    # We multiply by 'count' if multiple mutations happened at the exact same timestep
    end_t <- start_t + length(flicker_traj) - 1
    
    total_p[start_t:end_t] <- total_p[start_t:end_t] + (flicker_traj * count)
  }
  
  # Trim back to T_total (remove the "run-off" after simulation ends)
  total_p <- total_p[1:T_total]
  
  return(data.frame(
    Time = 1:T_total,
    Frequency = total_p,
    MutationRate = paste0("Mutation Rate: ", mu)
  ))
}

# ==============================================================================
# 3. RUN SIMULATIONS
# ==============================================================================

# Step A: Pre-calculate the shape of the flicker 
# (We assume all mutations have the same underlying dynamics)
base_flicker <- get_additive_flicker_trajectory(
  b = 4, c = 1, I = 0.5, lambda = 0.004, R_GL = 0.15, R_LT_eq = 0.5, p0 = 0.001
)

# Step B: Run for 3 different mutation rates
# Note: mu needs to be high enough to see events, but low enough for the 'rare' assumption
sim_low <- run_poisson_hits(mu = 0.002, T_total = 2000, flicker_traj = base_flicker)
sim_med <- run_poisson_hits(mu = 0.006,  T_total = 2000, flicker_traj = base_flicker)
sim_high <- run_poisson_hits(mu = 0.05,  T_total = 2000, flicker_traj = base_flicker)

# Combine Data
all_results <- rbind(sim_low, sim_med, sim_high)

# Factor ordering to ensure plots appear Low -> High
all_results$MutationRate <- factor(all_results$MutationRate, levels = unique(all_results$MutationRate))

# ==============================================================================
# 4. PLOTTING
# ==============================================================================

p_final <- ggplot(all_results, aes(x = Time, y = Frequency)) +
  # Use a filled area or line. Area looks nice for "Load".
  geom_area(alpha = 0.7) +
  geom_line(color = "black", size = 0.3) +
  
  # Facet into 3 panels vertically
  facet_wrap(~MutationRate, ncol = 1, scales = "free_y") +
  
  labs(
    title = "Poisson Hits Model",
    
    y = "Total cheater frequency",
    x = "Time (Generations)"
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    strip.text = element_text(face = "bold", size = 12),
    panel.grid.minor = element_blank()
  )

print(p_final)