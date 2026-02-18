#--------------------------------------------------------
#----------------- General Functions --------------------
#--------------------------------------------------------
#
kendall_tau_matrix_to_pearson_correlation_matrix <- function(kendall_tau_matrix) {
  sin((pi/2)*kendall_tau_matrix)
}


positive_semi_definite_maker <- function(A) {
  # making positive definite
  minimum_eigen_value <- min(eigen(A)$values)
  d <- nrow(A)
  if(minimum_eigen_value < 0) {
    A <- A - diag(minimum_eigen_value, nrow = d, ncol = d)
  }
  A
}


# Estimate pearson correlation using the mapping of the empirical kendall's tau.
## The random_set_ratio arguments contorols the ratio of samples incorporated in estimation.
## It is used in validation step for parameter tuning.
weighted_rank_based_pearson_correlation_estimator <- function(datasets, random_set_ratio=1) {
  number_of_nodes <- ncol(datasets[[1]])
  empirical_kendall <- list()
  w <- NULL
  for(p in 1:length(datasets)) {
    number_of_samples <- nrow(datasets[[p]])
    w <- c(w, number_of_samples)
    # Create random sub set of samples
    indices <- sample(1:number_of_samples, size = random_set_ratio*number_of_samples, replace = F)
    # Estimate empirical kendall's tau using weighted mean of empirical kandall's tau of datasets
    empirical_kendall[[p]] <- cor(datasets[[p]][indices,],
                                  method = "kendall",
                                  use = "pairwise.complete.obs");
  }
  mat_array <- array(unlist(empirical_kendall), dim = c(number_of_nodes,
                                                        number_of_nodes, 
                                                        length(empirical_kendall)))
  weighted_mean_empirical_kendall <- apply(mat_array, c(1,2), 
                                           function(x) weighted.mean(x, w, na.rm = TRUE))
  empirical_pearson = kendall_tau_matrix_to_pearson_correlation_matrix(weighted_mean_empirical_kendall)
  # make it positive semi definite
  empirical_pearson <- positive_semi_definite_maker(empirical_pearson)
  empirical_pearson
}


#--------------------------------------------------------
#-------- Functions in Real Data Experiments ------------
#--------------------------------------------------------
#
instability_evaluator_of_solution_pathes <- function(solution_pathes, number_of_nodes){
  number_of_random_set <- length(solution_pathes)
  
  instability_results <- data.frame()
  set_index <- rep(1, number_of_random_set)
  while(T){
    max_index = NA;
    max_lambda_value = 0
    edge_frequency = rep(0, number_of_nodes^2)
    number_of_results <- 0;
    for(i in 1:number_of_random_set){
      edge_existence = rep(0, number_of_nodes^2)
      if(set_index[i] <= length(solution_pathes[[i]])) {
        if(set_index[i] > 1) {
          edge_existence[(solution_pathes[[i]][[set_index[i]-1]]$active_set+1)] <- 1
        }
        edge_frequency <- edge_frequency + edge_existence
        
        number_of_results = number_of_results+1
        
        if(max_lambda_value < solution_pathes[[i]][[set_index[i]]]$knots_lambdas){
          max_index = i
          max_lambda_value = solution_pathes[[i]][[set_index[i]]]$knots_lambdas
        }
      }
    }
    if(number_of_results < number_of_random_set){
      break;
    }
    
    if(max_lambda_value == 0){
      break;
    }
    
    Adj <- matrix(edge_frequency/number_of_random_set, nrow = number_of_nodes)
    edge_existence_probability <- Adj[upper.tri(Adj)]
    instability <- sum(2*(1-edge_existence_probability)*edge_existence_probability)/choose(number_of_nodes,2)
    
    instability_results <- rbind(instability_results,
                                 data.frame(lambda=max_lambda_value,
                                            instability=instability))
    
    set_index[max_index] = set_index[max_index]+1;
  }
  instability_results
}


library(stats)
fisher_test <- function(all_genes, hub_genes, functional_genes) {
  functional_index <- all_genes %in% functional_genes
  hub_index <- all_genes %in% hub_genes
  
  pDNA_test <-
    matrix(c(sum(functional_index&hub_index), sum(!functional_index&hub_index),
             sum(functional_index&!hub_index), sum(!functional_index&!hub_index)),
           nrow = 2,
           dimnames = list(inDR = c("T", "F"),
                           isHub = c("T", "F")))
  fisher.test(pDNA_test, alternative = "greater")$p.value
}



#--------------------------------------------------------
#------ Functions in Synthetic Data Experiments ---------
#--------------------------------------------------------
#
library("GGMselect")
library("igraph")
# create two complete precision matrix
generate_reference_models <- function(number_of_nodes, number_of_samples,
                                      type="ScaleFree", density_of_graph = 0.2,
                                      power = 1, number_of_changes, mult=1,
                                      change_type="Random")
{
  #######################################
  ########### Generate Model ############
  #######################################
  
  # zero matrix generation
  A <- matrix(rep(0,number_of_nodes*number_of_nodes), nrow=number_of_nodes)
  B <- matrix(rep(0,number_of_nodes*number_of_nodes), nrow=number_of_nodes)
  total_possible_edges <- choose(number_of_nodes, 2)
  
  # create common model
  random_weights <- rnorm(total_possible_edges)
  random_weights <- random_weights + mult*sign(random_weights)
  
  
  # define base of precision matrix
  A[lower.tri(A)] <- random_weights
  A <- t(A)
  A[lower.tri(A)] <- random_weights
  
  B <- A
  
  
  change_mask <- matrix(data = 0, nrow = number_of_nodes,
                        ncol = number_of_nodes)
  if(change_type == "Random") {
    indices <- 1:choose(number_of_nodes,2)
    change_mask[lower.tri(change_mask)] <- indices
    change_mask <- t(change_mask)
    change_mask[lower.tri(change_mask)] <- indices
    mask <- sample(indices,number_of_changes)
    change_mask[change_mask %in% mask] <- -1
    change_mask[change_mask != -1] <- 0
  } else if(change_type == "Hub") {
    index_hub <- sample(1:number_of_nodes,1)
    index_leaf <- sample(setdiff(1:number_of_nodes, index_hub), number_of_changes)
    change_mask[index_hub, index_leaf] <- -1
    change_mask[index_leaf, index_hub] <- -1
  } else {
    errorCondition(message = "The change_type value is not standard!")
  }
  # create structure
  if(type == "ScaleFree") {
    g <- barabasi.game(number_of_nodes, power, directed = F);
    graph_matrix <- as_adjacency_matrix(g, sparse = F);
    A <- A*graph_matrix;
    B <- B*(graph_matrix+change_mask);
  } else {
    errorCondition(message = "The type value is not standard!")
  }
  
  
  # making positive definite
  minimum_eigen_value <- min(c(eigen(A)$values, eigen(B)$values))
  
  if(minimum_eigen_value < 1)
  {
    A <- A + diag(1-minimum_eigen_value, nrow = number_of_nodes,
                  ncol = number_of_nodes)
    B <- B + diag(1-minimum_eigen_value, nrow = number_of_nodes,
                  ncol = number_of_nodes)
  }
  
  
  #########################################
  ########### Generate Samples ############
  #########################################
  covariance_matrix_A <- solve(A)
  samples_A <- rmvnorm(number_of_samples, mean=rep(0,number_of_nodes), sigma=covariance_matrix_A)
  covariance_matrix_B <- solve(B)
  samples_B <- rmvnorm(number_of_samples, mean=rep(0,number_of_nodes), sigma=covariance_matrix_B)
  
  list(samples_A = samples_A, precision_matrix_A = A, covariance_matrix_A = covariance_matrix_A,
       samples_B = samples_B, precision_matrix_B = B, covariance_matrix_B = covariance_matrix_B)
}


# Estimated differential structure evaluator
result_evaluator <- function(real_differences, estimated_differences){
  real_differences_support <- sign(real_differences)
  estimated_differences_support  <- sign(estimated_differences)
  
  real_differences_support <- real_differences_support[upper.tri(real_differences_support)]
  estimated_differences_support <- estimated_differences_support[upper.tri(estimated_differences_support)]
  
  # sign evaluator
  ST <- sum(real_differences_support == estimated_differences_support)
  STN <- sum(real_differences_support == estimated_differences_support & real_differences_support == 0)
  STP <- ST - STN
  
  SF <- sum(real_differences_support != estimated_differences_support)
  SFP <- sum(real_differences_support != estimated_differences_support & real_differences_support == 0)
  SFN <- SF - SFP
  
  STPR <- STP/(STP+SFN)
  SFPR <- SFP/(STN+SFP)
  
  SACC <- ST/(ST+SF)
  
  SPrecision <- STP/(STP+SFP)
  SRecall <- STP/(STP+SFN)
  
  # support evaluator
  NT <- sum(abs(real_differences_support) == abs(estimated_differences_support))
  NTN <- sum(abs(real_differences_support) == abs(estimated_differences_support) & real_differences_support == 0)
  NTP <- NT - NTN
  
  NF <- sum(abs(real_differences_support) != abs(estimated_differences_support))
  NFP <- sum(abs(real_differences_support) != abs(estimated_differences_support) & real_differences_support == 0)
  NFN <- NF - NFP
  
  NTPR <- NTP/(NTP+NFN)
  NFPR <- NFP/(NTN+NFP)
  
  NACC <- NT/(NT+NF)
  
  NPrecision <- NTP/(NTP+NFP)
  NRecall <- NTP/(NTP+NFN)
  
  data.frame(STP=STP, STN=STN, SFP=SFP, SFN=SFN, STPR=STPR, SFPR=SFPR, SACC=SACC, SPrecision=SPrecision, SRecall=SRecall,
             NTP=NTP, NTN=NTN, NFP=NFP, NFN=NFN, NTPR=NTPR, NFPR=NFPR, NACC=NACC, NPrecision=NPrecision, NRecall=NRecall)
}


solution_path_performance_evaluator <- function(solution_path, differential_structure) {
  performances <- NULL
  for(i in 1:length(solution_path$solution_path)) {
    if(i == 1) {
      estimated_differences <- rep(0, length = nrow(differential_structure) * ncol(differential_structure))
    } else {
      estimated_differences <- rep(0, length = nrow(differential_structure) * ncol(differential_structure))
      estimated_differences[1+solution_path$solution_path[[i-1]]$active_set] <- solution_path$solution_path[[i]]$active_set_values
    }
    estimated_differences <- matrix(estimated_differences, nrow=nrow(differential_structure))
    performances <- rbind(performances, data.frame(lambda = solution_path$solution_path[[i]]$knots_lambdas,
                                                   result_evaluator(real_differences = differential_structure,
                                                                   estimated_differences = estimated_differences)))
  }
  performances
}


aggregate_dtrace_solution_path_resuts <- function(number_of_repetition,
                                              dtrace_solution_path_performances){
  # mean of dtrace solution path results
  set_index <- rep(1, number_of_repetition)
  mean_dtrace_solution_path_performances <- NULL;
  remaining_states <- sum(sapply(dtrace_solution_path_performances, nrow))
  
  number_of_results = number_of_repetition
  while(remaining_states > 0 && number_of_results == number_of_repetition){
    max_index = NA;
    max_lambda_value = 0
    temp_performance_record <- NULL;
    number_of_results <- 0;
    for(i in 1:number_of_repetition){
      if(set_index[i] <= nrow(dtrace_solution_path_performances[[i]])){
        temp_performance_record <- rbind(temp_performance_record,
                                         dtrace_solution_path_performances[[i]][set_index[i],])
        number_of_results = number_of_results+1
        
        if(max_lambda_value < dtrace_solution_path_performances[[i]][set_index[i],1]){
          max_index = i
          max_lambda_value = dtrace_solution_path_performances[[i]][set_index[i],1]
        }
      }
    }
    
    if(max_lambda_value == 0){
      break;
    }
    
    temp_performance_record <- apply(temp_performance_record,2,mean, na.rm=T);
    temp_performance_record[1] <- max_lambda_value
    mean_dtrace_solution_path_performances <- rbind(mean_dtrace_solution_path_performances,
                                                  temp_performance_record)
    set_index[max_index] = set_index[max_index]+1;
    remaining_states <- remaining_states - 1;
  }
  colnames(mean_dtrace_solution_path_performances) <- colnames(dtrace_solution_path_performances[[1]])
  mean_dtrace_solution_path_performances <- head(mean_dtrace_solution_path_performances,-1)
  mean_dtrace_solution_path_performances
}

compute_gamma_submatrix <- function(pre_active_set_indices, CA, CB, nVar){
  # Pre-allocate the symmetric matrix g
  active_set_size = length(pre_active_set_indices)
  g <- matrix(0, nrow = active_set_size, ncol = active_set_size)
  
  for (i in 1:active_set_size) {
    # Get the value for the current row from the list
    # Using double brackets [[ ]] to access list elements by index
    current_index_val <- pre_active_set_indices[i]
    
    # Calculate the 2D coordinates for the current element (for the row)
    columnA_current <- floor((current_index_val) / nVar) + 1 # R indices start at 1
    columnB_current <- ((current_index_val) %% nVar) + 1      # R indices start at 1
    
    # Iterate for the columns (like the inner iterator 'it')
    for (j in 1:active_set_size) {
      # Get the value for the current column from the list
      inner_index_val <- pre_active_set_indices[j]
      
      # Calculate the 2D coordinates for the inner element (for the column)
      rowA_inner <- floor((inner_index_val) / nVar) + 1
      rowB_inner <- ((inner_index_val) %% nVar) + 1
      
      # Perform the calculation
      # Note: R uses [row, column] indexing for matrices
      val <- CA[rowA_inner, columnA_current] * CB[rowB_inner, columnB_current] + 
        CB[rowA_inner, columnA_current] * CA[rowB_inner, columnB_current]
      
      # Assign the value to the matrix
      g[i, j] <- val
    }
  }
  
  return(g)
}


find_max_index <- function(sp, lambda) {
  low <- 1
  high <- length(sp)
  
  # Edge cases
  if (high == 0)  # Empty list
    stop("Error: List is empty!")
  if (sp[[1]]$knots_lambdas <= lambda) return(0)  # No element satisfies condition
  if (sp[[high]]$knots_lambdas > lambda) return(NULL)  # All elements satisfy condition

  # Binary search
  while (low <= high) {
    mid <- floor((low + high) / 2)
    
    if (sp[[mid]]$knots_lambdas > lambda) {
      # Check if this is the last element that satisfies the condition
      if (mid == high || sp[[mid + 1]]$knots_lambdas <= lambda) {
        return(mid)
      } else {
        low <- mid + 1
      }
    } else {
      high <- mid - 1
    }
  }
  
  return(0)  # Not found
}

library("Matrix")
get_knot_matrix <- function(i, sp, dim, eps=1e-12) {
  row_indices <- NULL
  column_indices <- NULL
  values <- NULL
  m <- NULL
  
  if (i == 1) {
    m = Matrix(0, dim, dim, sparse = TRUE)
  } else {
    active_set_values = sp[[i]]$active_set_values
    active_set_values[abs(active_set_values) < eps] = 0
    active_set = sp[[i-1]]$active_set
    active_set_size = length(active_set_values)
    
    for (j in 1:active_set_size) {
      # Get the value for the current column from the list
      index_val <- active_set[j]
      val <- active_set_values[j]
      
      # Calculate the 2D coordinates for the inner element (for the column)
      row_idx <- floor((index_val) / dim) + 1
      column_idx <- ((index_val) %% dim) + 1
      
      row_indices <- c(row_indices, row_idx)
      column_indices <- c(column_indices, column_idx)
      values <- c(values, val)
    }
    m = sparseMatrix(i = row_indices, j = column_indices, x = values, dims = c(dim, dim))
  }
  m
}

get_estimated_delta <- function(sp, lambda, dim) {
  i = find_max_index(sp, lambda)
  
  if (is.null(i)) return(NULL)
  if (i == 0) {
    Matrix(0, dim, dim, sparse = TRUE)
  } else {
    delta_knot_pre = get_knot_matrix(i, sp, dim)
    lambda_knot_pre = sp[[i]]$knots_lambdas
    delta_knot_past = get_knot_matrix(i+1, sp, dim)
    lambda_knot_past = sp[[i+1]]$knots_lambdas
    
    knot_pre_weight = (lambda - lambda_knot_past)/(lambda_knot_pre - lambda_knot_past)
    estimated_delta = (1-knot_pre_weight)*delta_knot_past + knot_pre_weight*delta_knot_pre
    estimated_delta
  }
}

get_obj_value <- function(delta, covA, covB, lambda) {
  if (is.null(delta)) return(NA)
  ( sum(diag(covA%*%delta%*%covB%*%delta)) + sum(diag(covB%*%delta%*%covA%*%delta))) / 4 -
    sum(diag(delta%*%(covA-covB))) +
    lambda*sum(abs(delta))
}

get_norm1_residuals <- function(delta, covA, covB, lambda) {
  if (is.null(delta)) return(NA)
  delta = as.matrix(delta)
  derivatives = (covA%*%delta%*%covB + covB%*%delta%*%covA)/2 + (covB-covA)
  residuals_matrix = matrix(NA, nrow = nrow(derivatives), ncol = ncol(derivatives))
  
  zero_ind = delta == 0
  residuals_matrix[!zero_ind] = abs(derivatives[!zero_ind] + 
                                         lambda*sign(delta[!zero_ind]))
  residuals_matrix[zero_ind] = pmax(abs(derivatives[zero_ind])-lambda, 0)
  sum(abs(residuals_matrix))
}
