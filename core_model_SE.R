run_simulation = function(
                          group_size = 7,                        # Average group size, base = 20
  generations ,                           
  ngroups = 20,                           # Number of groups, base = 20
  b = 2,
  c = 1,
  e = 0.001,
  m = 0.2,
  k = 0.25,
  inmS = 0.1,
  inmT = 0.1,
  base = 10
){

   results_a = results_n = results_s = results_t = meta_results = rep(NA,times=generations)

  n <- group_size * ngroups                # Total number of agents
  
  # Initialize agents dataframe
  agents <- data.frame(
    group = rep(1:ngroups, each = group_size),
    allele = 'N',
    payoff = base
  )
  
  
  # Initialize group_stats dataframe
  group_stats <- data.frame(s = rep(0, ngroups), t = rep(0, ngroups))
  
  
  # Simulation loop
  for (time in 1:generations) {
    #  print(paste("steps ",time,"pop total=",nrow(agents)))
    # Restore each generation to base fitness
    agents$payoff <- base
  if(time%%100==0)print(time)
    # Interaction Phase
    for (g in 1:ngroups) {
      group_agents <- agents[agents$group == g, ]
      group_s <- group_stats$s[g]
      group_t <- group_stats$t[g]
      

      
      #### PAIRING ####
      
      # positions inside group
      g_pos <- seq_len(nrow(group_agents))
      
      # segmented set ( s * group_size)
      seg_size <- round(group_s * length(g_pos))
      seg_pos  <- if (seg_size > 0) sample(g_pos, seg_size, replace = FALSE) else integer(0)
      
      # --- SEGMENTED: pair only A's among the segmented ---
      if (length(seg_pos) > 0) {
        seg_A_pos <- seg_pos[group_agents$allele[seg_pos] == "A"]
        if (length(seg_A_pos) > 1) {
          if (length(seg_A_pos) %% 2 == 1) {
            seg_A_pos <- seg_A_pos[-sample(length(seg_A_pos), 1)]
          }
          seg_pairs <- split(sample(seg_A_pos), rep(1:(length(seg_A_pos)/2), each = 2))
          for (pr in seg_pairs) {
            i <- pr[1]; j <- pr[2]
            # both are A 
            group_agents$payoff[i] <- group_agents$payoff[i] - c + b
            group_agents$payoff[j] <- group_agents$payoff[j] - c + b
          }
        }
      }
      
      # --- NON-SEGMENTED: pair everyone else at random ---
      nonseg_pos <- setdiff(g_pos, seg_pos)
      if (length(nonseg_pos) > 1) {
        if (length(nonseg_pos) %% 2 == 1) {
          nonseg_pos <- nonseg_pos[-sample(length(nonseg_pos), 1)]
        }
        ns_pairs <- split(sample(nonseg_pos), rep(1:(length(nonseg_pos)/2), each = 2))
        for (pr in ns_pairs) {
          i <- pr[1]; j <- pr[2]
          p_i <- group_agents$allele[i]
          p_j <- group_agents$allele[j]
          
          # i's payoff
          if (!is.na(p_i) && p_i == "A") group_agents$payoff[i] <- group_agents$payoff[i] - c
          if (!is.na(p_j) && p_j == "A") group_agents$payoff[i] <- group_agents$payoff[i] + b
          
          # j's payoff
          if (!is.na(p_j) && p_j == "A") group_agents$payoff[j] <- group_agents$payoff[j] - c
          if (!is.na(p_i) && p_i == "A") group_agents$payoff[j] <- group_agents$payoff[j] + b
        }
      }
      
      
      #      if (n_nonseg %% 2 == 1) {
      #        agents$payoff[as.numeric(rownames(non_segmented[excluded_N, ]))] <- 
      #          sum(non_segmented$payoff) / nrow(non_segmented)
      #      }
      
      #### TAXATION ####
      
      total_tax <- sum(group_agents$payoff * group_t)
      group_agents$payoff <- (group_agents$payoff * (1 - group_t)) + (total_tax / nrow(group_agents))
      
      
      
      #### REPRODUCTION ####
      
      avg_payoff <- mean(group_agents$payoff)
      offspring <- sample(1:nrow(group_agents), nrow(group_agents), replace = TRUE, 
                          prob = group_agents$payoff / avg_payoff)
      new_agents <- group_agents[offspring, ]
      
      
      
      #### MUTATION ####
      
      mutate <- runif(nrow(new_agents)) < e
      if (any(mutate)) {
        new_agents$allele[mutate] <- sample(c("N", "A"), sum(mutate), replace = TRUE)
      }
      
      
      # Update agents dataframe
      agents[agents$group == g, ] <- new_agents
      
    }
    
    ### Institution Mutation ###
    group_stats$s = group_stats$s + sample(c(-0.1, 0, 0.1),1, prob = c( (inmS/2), (1-inmS), (inmS/2)))
    group_stats$t = group_stats$t + sample(c(-0.1, 0, 0.1),1, prob = c( (inmT/2), (1-inmT), (inmT/2)))
    group_stats$s[which(group_stats$s>0.5)] = 0.5
    group_stats$s[which(group_stats$s<0)] = 0
    group_stats$t[which(group_stats$t>1)] = 1
    group_stats$t[which(group_stats$t<0)] = 0    
    
    
    #### MIGRATION ####
    migrants <- sample(1:nrow(agents), size = round(m * nrow(agents)))
    if (length(migrants) > 0) {
      for (migrant in migrants) {
        if (sum(agents$group == agents$group[migrant]) > 4) {
          possible_groups <- setdiff(1:ngroups, agents$group[migrant])
          agents$group[migrant] <- sample(possible_groups, 1)
        }
      }
    }
    
    
    #### CONFLICT  ####
    
    # Which groups go to war (prob k)
    tempwar       <- rbinom(ngroups, 1, k) == 1
    wargoers      <- which(tempwar)
    non_wargoers  <- setdiff(seq_len(ngroups), wargoers)
    
    # If odd, pick ONE MORE group to join the competition
    if (length(wargoers) %% 2 == 1 && length(non_wargoers) > 0) {
      wargoers <- c(wargoers, sample(non_wargoers, 1))
    }
    
    if (length(wargoers) >= 2) {
      # Randomly pair the wargoers
      wargoers <- sample(wargoers)
      l_groups <- wargoers[1:(length(wargoers)/2)]
      m_groups <- wargoers[(length(wargoers)/2 + 1):length(wargoers)]
      
      for (war in seq_along(l_groups)) {
        gl <- l_groups[war]
        gm <- m_groups[war]
        
        # Winner by TOTAL group payoff 
        payoff_l <- sum(agents$payoff[agents$group == gl])  ### ADD institution costs ###
        payoff_m <- sum(agents$payoff[agents$group == gm])  ### ADD institution costs ###
        
        winner <- if (payoff_l > payoff_m) gl else gm
        loser  <- if (payoff_l > payoff_m) gm else gl
        
        # --- Winner repopulates loser; winner becomes temporarily enlarged ---
        winner_rows <- which(agents$group == winner)
        loser_rows  <- which(agents$group == loser)
        gw          <- length(winner_rows)
        glos        <- length(loser_rows)
        
        # Apply winner's population frequency to the ENLARGED group (no quadratic bias):
        pA      <- mean(agents$allele[winner_rows] == "A")  # winner's A-share
        n_total <- gw + glos
        nA      <- round(pA * n_total)
        nN      <- n_total - nA
        
        # Build the enlarged group by copying random winner rows for other fields,
        # then set alleles to exactly match the desired frequencies.
        take_A_from <- sample(winner_rows, nA, replace = TRUE)
        take_N_from <- sample(winner_rows, nN, replace = TRUE)
        
        combined <- rbind(agents[take_A_from, ], agents[take_N_from, ])
        combined$allele <- c(rep("A", nA), rep("N", nN))
        
        # Randomly split the MERGED group back into two groups of original sizes
        assign_to_winner <- c(rep(TRUE, gw), rep(FALSE, glos))
        assign_to_winner <- sample(assign_to_winner)
        
        combined$group[assign_to_winner]  <- winner
        combined$group[!assign_to_winner] <- loser
        
        # Write back both groups at once
        agents[c(winner_rows, loser_rows), ] <- combined
        
        # Copy institutions from winner to loser (has no effect when inm=0 and s=t=0)
        group_stats$s[loser] <- group_stats$s[winner]
        group_stats$t[loser] <- group_stats$t[winner]
      }
    
      
      
    }
    
    results_a[time] <- nrow(agents[agents$allele=="A",])/nrow(agents)
#    results_n[time] <- nrow(agents[agents$allele=="N",])/nrow(agents)
    results_s[time] <- mean(group_stats$s)
    results_t[time] <- mean(group_stats$t)
  }
	
  return(rbind(results_a,results_s,results_t))
  
}
