source('./SPDtrace/Libraray.R')
source('./Simulation/editedDtrace.R')
Rcpp::sourceCpp('./SPDtrace/SolutionPath.cpp')
Rcpp::sourceCpp('./Simulation/CrossFDTL.cpp')

library(dplyr)
library(ggplot2)
library(plotly)
library(stringr)
library(ggpubr)
library(Matrix)
library(scales)
library(peakRAM)

#--------------------------------------------------------
#------------------ Initialization ----------------------
#--------------------------------------------------------
#
set.seed(1);
number_of_nodes = c(50, 100);
number_of_changes = c(20);
number_of_tuning_parameters = c(50);
maximum_legal_sparsity_in_solution_path = c(100);
number_of_samples = c(5)
number_of_sources = c(1)
maximum_iteration = c(50);
number_of_repetition = 100;

settings <- expand.grid(number_of_nodes = number_of_nodes,
                        number_of_changes = number_of_changes,
                        number_of_tuning_parameters = number_of_tuning_parameters,
                        maximum_legal_sparsity_in_solution_path = maximum_legal_sparsity_in_solution_path,
                        number_of_samples = number_of_samples,
                        number_of_sources = number_of_sources,
                        maximum_iteration = maximum_iteration)

setting_results <- list();



#--------------------------------------------------------
#----- Repetition of evaluation on synthetic data ------
#--------------------------------------------------------
#
for(setting_index in 1:nrow(settings)) {
  print(paste("setting index: ", setting_index, " of ", nrow(settings)))
  
  # Experiment Setting
  number_of_nodes = settings$number_of_nodes[setting_index]
  number_of_changes = settings$number_of_changes[setting_index]
  number_of_tuning_parameters = settings$number_of_tuning_parameters[setting_index]
  maximum_legal_sparsity_in_solution_path = settings$maximum_legal_sparsity_in_solution_path[setting_index]
  number_of_samples = settings$number_of_samples[setting_index]
  number_of_sources = settings$number_of_sources[setting_index]
  maximum_iteration = settings$maximum_iteration[setting_index]
  
  dtrace_solution_path_performances <- list();
  dtrace_solution_path_res <- list();
  dtrace_solution_path_times <- vector();
  dtrace_performances <- list();
  dtrace_res <- list();
  
  dtrace_solution_path_res <- list();
  dtrace_res <- list();

  
  # Generate tuning parameters
  tuning_parameters <- 0.001+0.399*(0:number_of_tuning_parameters)/(number_of_tuning_parameters)
  
  for(k in 1:number_of_repetition) {
    # Generate Models
    reference <- generate_reference_models(
      number_of_nodes = number_of_nodes,
      number_of_samples = number_of_samples*number_of_sources,
      number_of_changes = number_of_changes
    )
    print(paste("Repetition number: (", setting_index, " : ", k, ")"))
    
    model_precision_A = reference$precision_matrix_A
    model_precision_B = reference$precision_matrix_B
    
    differential_structure <- sign(model_precision_B - model_precision_A)
    diag(differential_structure) <- 0
    
    sample_covariances <- list(reference$covariance_matrix_A, reference$covariance_matrix_B)
    
    # Differential Network Learning
    # Inference Using Solution Path Method
    dtrace_solution_path_mem <- peakRAM({
      dtrace_solution_path_start_time <- Sys.time()
      dtrace_solution_path_results <- inference_Dtrace_solution_path(sample_covariances[[1]], sample_covariances[[2]],
                                                                     maximum_legal_sparsity_in_solution_path)
      dtrace_solution_path_stop_time <- Sys.time()
    })
    dtrace_solution_path_mem_total_mb <- dtrace_solution_path_mem$Peak_RAM_Used_MiB
    
    dtrace_solution_path_time <- dtrace_solution_path_stop_time - dtrace_solution_path_start_time
    dtrace_solution_path_times[k] <- as.numeric(dtrace_solution_path_time);
    dtrace_solution_path_performances[[k]] <- solution_path_performance_evaluator(solution_path = dtrace_solution_path_results,
                                                                                  differential_structure = differential_structure)
    
    # Inference Using D-trace and CrossFDTL Methods
    temp_dtrace_performances <- NULL;
    temp_dtrace_time <- vector();
    
    temp_dtrace_solution_path_res <- data.frame();
    temp_dtrace_res <- data.frame();

    for(iterator in 1:length(tuning_parameters)) {
      lambda <- tuning_parameters[iterator]
      dtrace_mem <- peakRAM({
        dtrace_results <- edited_Dtrace(X = sample_covariances, lambda = lambda, maxiter = maximum_iteration)
      })

      # Record the estimation of solution path corresponding to each tuning parameters
      sp = dtrace_solution_path_results$solution_path
      estimated_delta=get_estimated_delta(sp, lambda, number_of_nodes)
      obj_solution_path = get_obj_value(estimated_delta,
                                        sample_covariances[[1]],
                                        sample_covariances[[2]],
                                        lambda)
      norm1_residuals = get_norm1_residuals(estimated_delta,
                                            sample_covariances[[1]],
                                            sample_covariances[[2]],
                                            lambda)
      if (is.null(estimated_delta)) {
        non_zero_exact = NA
      } else {
        non_zero_exact = sum(estimated_delta != 0)
      }
      temp_dtrace_solution_path_res <- rbind(temp_dtrace_solution_path_res,
                                             cbind(lambda=lambda, 
                                                   objective_value=obj_solution_path,
                                                   non_zero_exact=non_zero_exact,
                                                   norm1_residuals=norm1_residuals,
                                                   mem_total_mb=dtrace_solution_path_mem_total_mb,
                                                   best_time=dtrace_solution_path_time))
      
      best_itr_dtrace <- NA
      temp_dtrace_res_pt <- data.frame()
      for(iteration in 1:maximum_iteration){
        l = length(dtrace_results$iteration_times)
        dtrace_reuslt_adjacency_matrix <- NULL
        if (iteration <= l) {
          dtrace_reuslt_adjacency_matrix <- dtrace_results$estimated_delta_logs[[iteration+1]]
        }
        
        temp_dtrace_performances <- rbind(temp_dtrace_performances,
                                          cbind(iteration, lambda, result_evaluator(real_differences = differential_structure,
                                                                                    estimated_differences = dtrace_reuslt_adjacency_matrix)));
        
        # Record the estimation of solution path corresponding to each tuning parameters
        obj_dtrace = get_obj_value(dtrace_reuslt_adjacency_matrix,
                                   sample_covariances[[1]],
                                   sample_covariances[[2]],
                                   lambda)
        norm1_residuals = get_norm1_residuals(dtrace_reuslt_adjacency_matrix,
                                              sample_covariances[[1]],
                                              sample_covariances[[2]],
                                              lambda)
        if (is.null(estimated_delta)) {
          dist_l1_norm = NA
          dist_frob_norm <- NA
          skeleton_change <- NA
        } else {
          dist = dtrace_reuslt_adjacency_matrix - estimated_delta
          dist_l1_norm = sum(abs(dist))
          dist_frob_norm <- norm(dist, type = "F")
          skeleton_change <- sum((dtrace_reuslt_adjacency_matrix!=0) != (estimated_delta!=0))
        }
        temp_dtrace_res_pt <- rbind(temp_dtrace_res_pt,
                                    cbind(iteration=iteration,
                                          lambda=lambda, 
                                          objective_value=obj_dtrace,
                                          objective_diff=obj_dtrace - obj_solution_path,
                                          l1_norm_approximation_error=dist_l1_norm,
                                          frob_norm_approximation_error=dist_frob_norm,
                                          skeleton_change=skeleton_change,
                                          non_zero_exact=non_zero_exact,
                                          skeleton_change_ratio=skeleton_change/non_zero_exact,
                                          norm1_residuals=norm1_residuals))
        
        
        if (is.na(best_itr_dtrace) & !is.na(skeleton_change) & skeleton_change==0)
          best_itr_dtrace <- iteration
      }
      temp_dtrace_res_pt <- cbind(temp_dtrace_res_pt,
                                  mem_total_mb=dtrace_mem$Peak_RAM_Used_MiB/maximum_iteration)
      temp_dtrace_res <- rbind(temp_dtrace_res,
                               cbind(temp_dtrace_res_pt,
                                     best_iteration=best_itr_dtrace,
                                     best_time=dtrace_results$iteration_times[best_itr_dtrace]))
      temp_dtrace_time <- c(temp_dtrace_time, dtrace_results$iteration_times);
    }
    
    dtrace_performances[[k]] <- cbind(temp_dtrace_performances, time = temp_dtrace_time);

    dtrace_solution_path_res[[k]] <- temp_dtrace_solution_path_res
    dtrace_res[[k]] <- temp_dtrace_res;
  }
  
  # Averaging Evaluation Results
  performances_colnames <- colnames(dtrace_solution_path_performances[[1]])
  
  dtrace_performances_data.frame <- bind_rows(dtrace_performances)
  mean_dtrace_performances <- aggregate(dtrace_performances_data.frame,
                                        by = list(dtrace_performances_data.frame$iteration,
                                                  dtrace_performances_data.frame$lambda),
                                        FUN = mean, na.rm=T)
  
  # Averaging objective results
  dtrace_solution_path_res_data.frame <- bind_rows(dtrace_solution_path_res)
  mean_dtrace_solution_path_res <- aggregate(. ~ lambda, data = dtrace_solution_path_res_data.frame,
                                             FUN = function(x) {
                                               c(mean = mean(x, na.rm = TRUE), 
                                                 se = sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))
                                             })
  mean_dtrace_solution_path_res_flat <- do.call(data.frame, mean_dtrace_solution_path_res)
  
  
  dtrace_res_data.frame <- bind_rows(dtrace_res)
  mean_dtrace_res <- dtrace_res_data.frame %>%
    group_by(iteration, lambda) %>%
    summarise(
      # Apply fun1 to best_time and best_iteration
      across(c(best_time, best_iteration, skeleton_change_ratio), 
             list(
               mean = ~mean(., na.rm = FALSE),
               se = ~sd(., na.rm = FALSE) / sqrt(sum(!is.na(.)))
             ),
             .names = "{.col}.{.fn}"),
      
      # Apply fun2 to norm1_residuals, objective_diff and skeleton_change_ratio
      across(c(norm1_residuals, objective_diff,
               objective_value, l1_norm_approximation_error,
               frob_norm_approximation_error, skeleton_change,
               non_zero_exact, mem_total_mb), 
             list(
               mean = ~mean(., na.rm = TRUE),
               se = ~sd(., na.rm = TRUE) / sqrt(sum(!is.na(.)))
             ),
             .names = "{.col}.{.fn}")
    )
  mean_dtrace_res_flat <- do.call(data.frame, mean_dtrace_res)
  
  # Averaging results obtained by solution path method
  set_index <- rep(1, number_of_repetition)
  mean_dtrace_solution_path_performances <- NULL;
  remaining_states <- sum(sapply(dtrace_solution_path_performances, nrow))
  
  number_of_results = number_of_repetition
  while(remaining_states > 0 && number_of_results == number_of_repetition) {
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
    temp_performance_record <- apply(temp_performance_record,2,mean, na.rm=T);
    temp_performance_record[1] <- max_lambda_value
    mean_dtrace_solution_path_performances <- rbind(mean_dtrace_solution_path_performances,
                                                    temp_performance_record)
    if(is.na(max_index))
      break;
    set_index[max_index] = set_index[max_index]+1;
    remaining_states <- remaining_states - 1;
  }
  colnames(mean_dtrace_solution_path_performances) <- colnames(dtrace_solution_path_performances[[1]])
  mean_dtrace_solution_path_performances <- head(mean_dtrace_solution_path_performances,-1)
  
  # Data shaping. Integrate all evaluation results.
  performances <- data.frame(mean_dtrace_solution_path_performances,
                             time = mean(dtrace_solution_path_times),
                             type="Solution_Path")
  colnames(performances) <- c(performances_colnames, "time", "type")
  
  res <- data.frame()
  # res_by_non_zero_exact <- data.frame()
  for(i in 1:maximum_iteration){
    temp_performances <- mean_dtrace_performances[mean_dtrace_performances$iteration == i,
                                                  c(performances_colnames, "time")]
    performances <- rbind(performances, data.frame(temp_performances,
                                                   type=paste0("D-trace, Iteration ", i)))
    
  }
  
  # Data aggregation and standardization
  names(mean_dtrace_res_flat)[names(mean_dtrace_res_flat) == "lambda.mean"] <- "lambda"
  
  missing_cols <- setdiff(names(mean_dtrace_res_flat)[-c(1:2)], names(mean_dtrace_solution_path_res_flat))
  
  mean_dtrace_res_flat$type = paste0("D-trace, Iteration ", mean_dtrace_res_flat$iteration)
  res <- rbind(res, mean_dtrace_res_flat[-c(1)])
  
  mean_dtrace_solution_path_res_flat[missing_cols] <- 0
  
  res <- rbind(res, cbind(mean_dtrace_solution_path_res_flat,
                          type="Solution_Path"))
  
  setting_results[[setting_index]] <- list(performances = performances,
                                           results = res)
}


#--------------------------------------------------------
#---------------- Visualization the Results -------------
#--------------------------------------------------------
#
# Define a color-blind safe, professional palette
# Okabe-Ito palette works well for up to 4 categories.
algorithm_colors <- c(
  "D-trace" = "#E69F00",           # Orange
  "Solution_Path" = "#56B4E9"     # Sky Blue
  # "CrossFDTL Cversion" = "#009E73", # Bluish Green
  # "Other" = "#F0E442" #Yello
)
# Hue palette
# algorithm_colors <- c(
#   "D-trace" = "#F8766D",           # Red
#   "Solution_Path" = "#00BA38",     # Green
#   "CrossFDTL Cversion" = "#619CFF" # Blue 
# )

# Define consistent labels
algorithm_labels <- c(
  "D-trace" = "APGD (Approximate)",
  "Solution_Path" = "Solution Path (Exact)"
)

# Line types for iterations
iteration_linetypes <- c(
  "1" = "solid",
  "5" = "dashed", 
  "10" = "dotted",
  "20" = "dotdash",
  "50" = "longdash"
)


optimization_results <- NULL;
for(setting_index in 1:length(setting_results)){
  setting_index_results <- setting_results[[setting_index]]$results
  optimization_results <- rbind(optimization_results,
                                cbind(settings[setting_index,], setting_index_results))
}

# rename the columns
colnames(optimization_results)[colnames(optimization_results)=="number_of_nodes"] = "d"
colnames(optimization_results)[colnames(optimization_results)=="number_of_changes"] = "s"
colnames(optimization_results)[colnames(optimization_results)=="number_of_samples"] = "m"
colnames(optimization_results)[colnames(optimization_results)=="maximum_legal_sparsity_in_solution_path"] = "c"
optimization_results <- optimization_results[order(optimization_results$lambda),]

# Illustrating the accuracy comparison of D-trace, SP D-trace and CrossFDTL using Precision-Recall graph
solution_path_index <- optimization_results$type == "Solution_Path"
non_solution_path_index <- optimization_results$type == "D-trace, Iteration 1" |
  optimization_results$type == "D-trace, Iteration 10" |
  optimization_results$type == "D-trace, Iteration 20" |
  optimization_results$type == "D-trace, Iteration 50" |
  optimization_results$type == "D-trace, Iteration 5"
row_index <-  (non_solution_path_index | solution_path_index)

optimization_results$algorithm <- sapply(str_split(optimization_results$type, ", "),"[[",1)
optimization_results$iteration <- 1
optimization_results$iteration[optimization_results$type != "Solution_Path"] <- sapply(str_split(optimization_results$type[optimization_results$type != "Solution_Path"], "Iteration "),"[[",2)

########### PLOT 1 ############
objective_diff_plot <- ggplot() +
  geom_line(mapping = aes(x = lambda, y = objective_diff.mean, color=algorithm, linetype =iteration),
            data = optimization_results[non_solution_path_index,], size = 0.5) +
  # 4. CRUCIAL: Synchronize COLOR and FILL scales
  scale_linetype_manual(values = iteration_linetypes,
                        limits = names(iteration_linetypes)) +
  scale_colour_manual(values = algorithm_colors,
                      limits = names(algorithm_labels),
                      labels = algorithm_labels) +
  scale_y_log10(
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x))
  ) +
  # scale_y_log10(
  #   limits = c(1e-10, 1),
  #   breaks = c(1e-10, 1e-8, 1e-6, 1e-4, 1e-2, 1),
  #   labels = trans_format("log10", math_format(10^.x))
  # ) +
  labs(y="Objective Value", x="Lambdas") +
  facet_grid(
    rows = vars(d), 
    labeller = label_both) +
  theme_bw() +
  # Labels
  labs(
    # title = "Convergence Profiles of Optimization Algorithms",
    x = "Regularization parameter",
    y = "Objective gap",
    linetype = "Iteration"
  ) +
  labs(color = "Method") + 
  theme(legend.position = "top",
        legend.background = element_rect(color = "black"),
        panel.grid.minor = element_blank(),
        strip.background = element_blank())
print(objective_diff_plot)

########### PLOT 2 ############
#Ratio Error Plot
valid_indices_ratio <- !is.na(optimization_results$skeleton_change_ratio.mean)
non_solution_path_index_ratio <- (optimization_results$type == "D-trace, Iteration 5" |
                                    optimization_results$type == "D-trace, Iteration 10" |
                                    optimization_results$type == "D-trace, Iteration 20" |
                                    optimization_results$type == "D-trace, Iteration 50" |
                                    optimization_results$type == "D-trace, Iteration 1") &
  valid_indices_ratio


valid_lambdas_ratio <- unique(optimization_results$lambda[non_solution_path_index_ratio])
solution_path_index_ratio <- optimization_results$type == "Solution_Path" & 
  optimization_results$lambda %in% valid_lambdas_ratio

row_index_ratio <-  (non_solution_path_index_ratio | solution_path_index_ratio)

skeleton_error_ratio_plot <- ggplot() +# 2. Line for the mean (uses COLOR aesthetic)
  geom_line(
    mapping = aes(
      x = lambda, 
      y = skeleton_change_ratio.mean, 
      color = algorithm,
      linetype = iteration
    ),
    data = optimization_results[non_solution_path_index_ratio,], 
    size = 0.5
  ) +
  scale_y_log10(
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x))
  ) +
  # scale_y_continuous(
  #   trans = pseudo_log_trans(sigma = 0.1),
  #   breaks = c(0, 0.1, 0.2, 0.5, 1, 2, 5, 10),
  #   labels = c("0", "0.1", "0.2", "0.5", "1", "2", "5", "10")
  # ) +
  # 4. CRUCIAL: Synchronize COLOR and FILL scales
  scale_linetype_manual(values = iteration_linetypes,
                        limits = names(iteration_linetypes)) +
  scale_colour_manual(values = algorithm_colors,
                      limits = names(algorithm_labels),
                      labels = algorithm_labels) +
  # 5. Labels and Faceting
  labs(
    y = "False selection rate", 
    x = "Regularization parameter",
    color = "Method",  # Legend title for color
    fill = "Method"    # Legend title for fill (this merges them!)
    # linetype = "Number of iterations" # Removed if not used elsewhere
  ) +
  facet_grid(
    rows = vars(d), 
    labeller = label_both
  ) +
  theme_bw() +
  theme(
    legend.position = "top",
    legend.background = element_rect(color = "black"),
    panel.grid.minor = element_blank(),
    strip.background = element_blank()
  )

# 6. Print the plot (do this on a separate line, not with '+')
print(skeleton_error_ratio_plot)

########### PLOT 3 ############
residual_index <- optimization_results$type == "D-trace, Iteration 5" |
                                    optimization_results$type == "D-trace, Iteration 10" |
                                    optimization_results$type == "D-trace, Iteration 20" |
                                    optimization_results$type == "D-trace, Iteration 50" |
                                    optimization_results$type == "D-trace, Iteration 1" |
                                    optimization_results$type == "Solution_Path"

residuals_plot <- ggplot() +
  # 2. Line for the mean (uses COLOR aesthetic)
  geom_line(
    mapping = aes(
      x = lambda, 
      y = norm1_residuals.mean, 
      color = algorithm,
      linetype = iteration
    ),
    data = optimization_results[residual_index,], 
    size = 0.5
  ) +
  # Machine precision reference line
  geom_hline(yintercept = 1e-12, linetype = "dotdash", 
             color = "grey50", linewidth = 0.5) +
  annotate("text", x = max(optimization_results$lambda), 
           y = 1e-11, label = "Machine precision", 
           hjust = 1, vjust = -0.5, size = 3) +
  scale_y_log10(
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x))
  ) +
  # 4. CRUCIAL: Synchronize COLOR and FILL scales
  scale_linetype_manual(values = iteration_linetypes,
                        limits = names(iteration_linetypes)) +
  scale_colour_manual(values = algorithm_colors,
                      limits = names(algorithm_labels),
                      labels = algorithm_labels) +
  # 5. Labels and Faceting
  labs(
    # title = "Exact vs Approximate Optimization: KKT Residuals",
    # subtitle = "Exact solution achieves machine precision across all λ values",
    y = "Total KKT residuals", 
    x = "Regularization parameter",
    linetype = "Iteration",
    color = "Method",  # Legend title for color
  ) +
  facet_grid(
    rows = vars(d), 
    labeller = label_both
  ) +
  theme_bw() +
  theme(
    legend.position = "top",
    legend.background = element_rect(color = "black"),
    panel.grid.minor = element_blank(),
    strip.background = element_blank()
  )

print(residuals_plot)


library(cowplot)
library(patchwork)

# 1. Extract legend
combined_legend <- get_legend(residuals_plot)

# 2. Remove legends from subplots
objective_diff_plot       <- objective_diff_plot       + theme(legend.position = "none")
skeleton_error_ratio_plot <- skeleton_error_ratio_plot + theme(legend.position = "none")
residuals_plot            <- residuals_plot            + theme(legend.position = "none")

# 3. Combine and tag ONLY this block
plots_combined <- residuals_plot + objective_diff_plot + skeleton_error_ratio_plot +
  plot_annotation(tag_levels = "A") +
  theme(
    plot.tag = element_text(size = 14, face = "bold"),
    plot.tag.position = c(0.02, 0.98),   # adjust to taste
    plot.margin = margin(5, 5, 5, 5)
  )

# 4. Wrap ONLY the legend, then stack with tagged plots
legend_plot <- wrap_elements(combined_legend)

final_plot <- legend_plot / plots_combined +
  plot_layout(heights = c(0.12, 1))  # do NOT add another plot_annotation here

print(final_plot)

ggsave("./Simulation/Results/exact_approx.pdf", plot = final_plot, width = 10, height = 5)







performances <- NULL;
for(setting_index in 1:length(setting_results)){
  perf_temp <- setting_results[[setting_index]]$performances
  sp_perf_temp <- perf_temp[perf_temp$type == "Solution_Path",][1,]
  performances <- rbind(performances,
                        cbind(settings[setting_index,],
                              sp_perf_temp))
}

# rename the columns
colnames(performances)[colnames(performances)=="number_of_nodes"] = "d"
colnames(performances)[colnames(performances)=="number_of_changes"] = "s"
colnames(performances)[colnames(performances)=="number_of_samples"] = "m"
colnames(performances)[colnames(performances)=="maximum_legal_sparsity_in_solution_path"] = "c"

# Illustrating the accuracy comparison of D-trace, SP D-trace and CrossFDTL using Precision-Recall graph
performances$algorithm <- sapply(str_split(performances$type, ", "),"[[",1)
performances$iteration <- 1
performances$iteration[performances$type != "Solution_Path"] <- sapply(str_split(performances$type[performances$type != "Solution_Path"], "Iteration "),"[[",2)

########### PLOT 4 ############
best_iteration_index <- (optimization_results$type=="D-trace, Iteration 1" |
                           optimization_results$type == "Solution_Path") &
  !is.na(optimization_results$best_iteration.se)

best_time_plot <- ggplot() +
  # 1. Ribbon for the error area (uses FILL aesthetic)
  geom_ribbon(
    mapping = aes(
      x = lambda,
      ymin = best_time.mean - best_time.se,
      ymax = best_time.mean + best_time.se,
      fill = algorithm,
    ),
    data = optimization_results[best_iteration_index,],
    alpha = 0.3,
    color = NA
  ) +
  geom_line(aes(x = lambda, y = best_time.mean - best_time.se, color = algorithm),
            data = optimization_results[best_iteration_index,],
            linetype = "dashed", alpha = 0.5) +
  geom_line(aes(x = lambda, y = best_time.mean + best_time.se, color = algorithm),
            data = optimization_results[best_iteration_index,],
            linetype = "dashed", alpha = 0.5) +
  # 2. Line for the mean (uses COLOR aesthetic)
  geom_line(
    mapping = aes(
      x = lambda, 
      y = best_time.mean,
      color = algorithm
    ),
    data = optimization_results[best_iteration_index,], 
    size = 0.5
  ) +
  scale_colour_manual(values = algorithm_colors,
                      limits = names(algorithm_labels),
                      labels = algorithm_labels) +
  scale_fill_manual(values = algorithm_colors,
                    labels = algorithm_labels,
                    guide = "none") +
  scale_y_log10(
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x))
  ) +
  # 5. Labels and Faceting
  labs(
    y = "Time (sec.)", 
    x = "Regularization parameter",
    color = "Method",  # Legend title for color
  ) +
  facet_grid(
    cols = vars(d), 
    labeller = label_both
  ) +
  theme_bw() +
  theme(
    legend.position = "top",
    legend.background = element_rect(color = "black"),
    panel.grid.minor = element_blank(),
    strip.background = element_blank()
  )

# 6. Print the plot (do this on a separate line, not with '+')
print(best_time_plot)



########### PLOT 4 ############
mem_index <- optimization_results$type == "D-trace, Iteration 1" |
  optimization_results$type == "Solution_Path"

mem_plot <- ggplot() +
  geom_ribbon(
    mapping = aes(
      x = lambda,
      ymin = mem_total_mb.mean - mem_total_mb.se,
      ymax = mem_total_mb.mean + mem_total_mb.se,
      fill = algorithm,
    ),
    data = optimization_results[mem_index,],
    alpha = 0.3,
    color = NA
  ) +
  geom_line(aes(x = lambda, y = mem_total_mb.mean - mem_total_mb.se, color = algorithm),
            data = optimization_results[mem_index,],
            linetype = "dashed", alpha = 0.5) +
  geom_line(aes(x = lambda, y = mem_total_mb.mean + mem_total_mb.se, color = algorithm),
            data = optimization_results[mem_index,],
            linetype = "dashed", alpha = 0.5) +
  geom_line(
    mapping = aes(
      x = lambda, 
      y = mem_total_mb.mean, 
      color = algorithm
    ),
    data = optimization_results[mem_index,], 
    size = 0.5
  ) +
  scale_colour_manual(values = algorithm_colors,
                      limits = names(algorithm_labels),
                      labels = algorithm_labels) +
  scale_fill_manual(values = algorithm_colors,
                    labels = algorithm_labels,
                    guide = "none") +
  scale_y_log10(
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x))
  ) +
  # 5. Labels and Faceting
  labs(
    y = "Memory (MB)", 
    x = "Regularization parameter",
    color = "Method",  # Legend title for color
  ) +
  facet_grid(
    cols = vars(d), 
    labeller = label_both
  ) +
  theme_bw() +
  theme(
    legend.position = "top",
    legend.background = element_rect(color = "black"),
    panel.grid.minor = element_blank(),
    strip.background = element_blank()
  )

# 6. Print the plot (do this on a separate line, not with '+')
print(mem_plot)

final_plot_resource <- ggarrange(mem_plot, best_time_plot, 
          common.legend = T,
          # legend="bottom",
          labels = c("A", "B"),
          ncol = 2, nrow = 1, widths = c(1,1))

ggsave("./Simulation/Results/exact_approx_resource.pdf", plot = final_plot_resource, width = 12, height = 4)
