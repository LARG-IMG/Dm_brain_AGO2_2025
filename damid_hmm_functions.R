run_hmm_on_granges <- function(granges_obj, n_states = 2, by_chr = TRUE, score_col = "score", 
                               chromosomes = c("2L", "2R", "3L", "3R", "X", "Y", "4"), 
                               collapse_enriched = TRUE) {
  # Load necessary libraries
  library(HiddenMarkov)
  library(GenomicRanges)
  library(dplyr)
  library(cluster)
  
  # Filter to only include specified chromosomes
  granges_obj <- keepSeqlevels(granges_obj, chromosomes,
                               pruning.mode = 'coarse')
  # Process by chromosome if specified
  if (by_chr) {
    # Split the GRanges object by chromosome
    gr_list <- split(granges_obj, seqnames(granges_obj))
    # Apply HMM to each chromosome
    results_list <- lapply(gr_list, function(chrom_data) {
      # Run HMM on this chromosome's data
      # print(head(chrom_data))
      process_granges_hmm(chrom_data, n_states = n_states, score_col = score_col)
    })
    
    # Combine the results back into a single GRanges object
    granges_with_states <- unlist(GRangesList(results_list))
    
  } else {
    # Run HMM on the entire GRanges object
    granges_with_states <- process_granges_hmm(granges_obj, n_states = n_states, score_col = score_col)
  }
  print("2")
  
  if (collapse_enriched) {
    # Identify enriched regions and collapse them
    enriched_granges <- granges_with_states[granges_with_states$state == "Enriched"]
    # Reduce to merge overlapping or adjacent enriched regions
    enriched_regions <- GenomicRanges::reduce(enriched_granges)
    return(enriched_regions)
  } else {
    return(granges_with_states)
  }
}

process_granges_hmm <- function(gr_obj, n_states = 2, score_col = "score") {
  # Extract the score values
  # print(head(gr_obj))
  observations <- mcols(gr_obj)[[score_col]]
  # print(head(observations))
  # Remove NA values
  valid_indices <- which(!is.na(observations))
  observations <- observations[valid_indices]
  gr_obj <- gr_obj[valid_indices]
  
  # Ensure observations are numeric
  observations <- as.numeric(observations)
  # Initialize transition probabilities
  trans_prob <- matrix(1 / n_states, nrow = n_states, ncol = n_states)
  # Initialize emission parameters
  # For the normal distribution, we need means and standard deviations for each state
  # We'll initialize the means based on quantiles of the data
  print(length(observations) >= n_states)
  print(length(observations))
  print(n_states)
  if (length(observations) >= n_states) {
    # init_means <- quantile(observations,
    #                        probs = seq(0, 1, length.out = n_states + 1))[-c(1, n_states + 1)]
    clar <- clara(observations, n_states)
    init_means <- clar$medoids
    init_sds <- tapply(observations, clar$clustering, sd)*0.5
    print(init_means)
  } else {
    # If not enough observations, use the overall mean
    init_means <- rep(mean(observations), n_states)
    init_sds <- rep(sd(observations), n_states)
  }
  
  
  # Define HMM
  hmm_model <- dthmm(
    observations,
    Pi = trans_prob,
    delta = rep(1 / n_states, n_states),
    distn = "norm",
    pm = list(mean = init_means, sd = init_sds),
    discrete = FALSE
  )
  
  # Train HMM
  trained_hmm <- BaumWelch(hmm_model)
  
  # Decode states
  decoded_states <- Viterbi(trained_hmm)
  
  # Assign state labels
  state_labels <- assign_state_labels(observations, decoded_states, n_states)
  
  # Add the state labels to the GRanges object
  mcols(gr_obj)$state <- state_labels
  
  return(gr_obj)
}

assign_state_labels <- function(observations, decoded_states, n_states) {
  # Calculate the mean observation for each state
  print("5")
  state_means <- tapply(observations, decoded_states, mean)
  
  # Rank the states based on their means (from lowest to highest)
  state_order <- order(state_means)
  
  # Assign labels based on the rank
  labels <- rep(NA, length(decoded_states))
  
  # For 2 states
  if (n_states == 2) {
    labels[decoded_states == state_order[1]] <- "Background"
    labels[decoded_states == state_order[2]] <- "Enriched"
  } else if (n_states == 3) {
    labels[decoded_states == state_order[1]] <- "Depleted"
    labels[decoded_states == state_order[2]] <- "Background"
    labels[decoded_states == state_order[3]] <- "Enriched"
  } else {
    # For more than 3 states
    label_names <- paste0("State", seq_len(n_states))
    labels <- label_names[match(decoded_states, state_order)]
  }
  
  return(labels)
}

