source('../SPDtrace/Libraray.R')
source('./editedDtrace.R')
Rcpp::sourceCpp('../SPDtrace/SolutionPath.cpp')
Rcpp::sourceCpp('./CrossFDTL.cpp')

library(dplyr)
library(ggplot2)
library(stringr)
library(ggpubr)

#--------------------------------------------------------
#------------------ Initialization ----------------------
#--------------------------------------------------------
#
set.seed(1)
number_of_nodes = c(50, 100)
number_of_changes = c(20)
number_of_tuning_parameters = c(50);
maximum_legal_sparsity_in_solution_path = c(100);
number_of_samples = c(500, 1000)
number_of_sources = c(1)
maximum_iteration = c(5);
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
  
  # Experiment Setting
  number_of_nodes = settings$number_of_nodes[setting_index]
  number_of_changes = settings$number_of_changes[setting_index]
  number_of_tuning_parameters = settings$number_of_tuning_parameters[setting_index]
  maximum_legal_sparsity_in_solution_path = settings$maximum_legal_sparsity_in_solution_path[setting_index]
  number_of_samples = settings$number_of_samples[setting_index]
  number_of_sources = settings$number_of_sources[setting_index]
  maximum_iteration = settings$maximum_iteration[setting_index]
  
  dtrace_solution_path_performances <- list();
  dtrace_solution_path_times <- vector();
  dtrace_performances <- list();
  crossFDTL_performances <- list();
  
  # Generate tuning parameters
  tuning_parameters <- 0.1+1.9*(0:number_of_tuning_parameters)/(number_of_tuning_parameters)
  
  for(k in 1:number_of_repetition) {
    # Generate Models
    reference <- generate_reference_models(
      number_of_nodes = number_of_nodes,
      number_of_samples = number_of_samples*number_of_sources,
      number_of_changes = number_of_changes
    )
    
    model_precision_A = reference$precision_matrix_A
    model_precision_B = reference$precision_matrix_B
    
    differential_structure <- sign(model_precision_B - model_precision_A)
    diag(differential_structure) <- 0
    
    sample_kendall_A <- cor(reference$samples_A, method = "kendall");
    sample_kendall_B <- cor(reference$samples_B, method = "kendall");
    
    sample_pearson_A = kendall_tau_matrix_to_pearson_correlation_matrix(sample_kendall_A)
    sample_pearson_B = kendall_tau_matrix_to_pearson_correlation_matrix(sample_kendall_B)
    
    # Making positive semi definite
    sample_pearson_A = positive_semi_definite_maker(sample_pearson_A)
    sample_pearson_B = positive_semi_definite_maker(sample_pearson_B)
    
    sample_covariances <- list(sample_pearson_A, sample_pearson_B)
    
    # Differential Network Learning
    # Inference Using Solution Path Method
    dtrace_solution_path_start_time <- Sys.time()
    dtrace_solution_path_results <- inference_Dtrace_solution_path(sample_covariances[[1]], sample_covariances[[2]],
                                                                   maximum_legal_sparsity_in_solution_path)
    dtrace_solution_path_stop_time <- Sys.time()
    
    dtrace_solution_path_time <- dtrace_solution_path_stop_time - dtrace_solution_path_start_time
    dtrace_solution_path_times[k] <- as.numeric(dtrace_solution_path_time);
    dtrace_solution_path_performances[[k]] <- solution_path_performance_evaluator(solution_path = dtrace_solution_path_results,
                                                                                  differential_structure = differential_structure)
    
    # Inference Using D-trace and CrossFDTL Methods
    temp_dtrace_performances <- NULL;
    temp_dtrace_time <- vector();
    temp_crossFDTL_performances <- NULL;
    temp_crossFDTL_time <- vector();
    
    for(iterator in 1:length(tuning_parameters)) {
      lambda <- tuning_parameters[iterator]
      
      dtrace_results <- edited_Dtrace(X = sample_covariances, lambda = lambda, maxiter = maximum_iteration)
      crossFDTL_results <- CrossFDTL(CovA = sample_covariances[[1]], CovB = sample_covariances[[2]],
                                     lambda = lambda, rho = 0, maxiter = maximum_iteration)
      
      for(iteration in 1:length(dtrace_results$iteration_times)){
        dtrace_reuslt_adjacency_matrix <- dtrace_results$estimated_delta_logs[[iteration+1]]
        temp_dtrace_performances <- rbind(temp_dtrace_performances,
                                          cbind(iteration, lambda, result_evaluator(real_differences = differential_structure,
                                                                                    estimated_differences = dtrace_reuslt_adjacency_matrix)));
      }
      
      for(iteration in 1:length(crossFDTL_results$iteration_times)){
        crossFDTL_result_adjacency_matrix <- as.matrix(crossFDTL_results$estimated_delta_logs[[iteration]])
        temp_crossFDTL_performances <- rbind(temp_crossFDTL_performances,
                                             cbind(iteration, lambda, result_evaluator(real_differences = differential_structure,
                                                                                       estimated_differences = crossFDTL_result_adjacency_matrix)));
      }
      
      temp_dtrace_time <- c(temp_dtrace_time, dtrace_results$iteration_times)
      temp_crossFDTL_time <- c(temp_crossFDTL_time, crossFDTL_results$iteration_times)
    }
    
    dtrace_performances[[k]] <- cbind(temp_dtrace_performances, time  = temp_dtrace_time);
    crossFDTL_performances[[k]] <- cbind(temp_crossFDTL_performances, time  = temp_crossFDTL_time);
  }
  
  # Averaging Evaluation Results
  performances_colnames <- colnames(dtrace_solution_path_performances[[1]])
  
  dtrace_performances_data.frame <- bind_rows(dtrace_performances)
  mean_dtrace_performances <- aggregate(dtrace_performances_data.frame, by = list(dtrace_performances_data.frame$iteration,
                                                                                  dtrace_performances_data.frame$lambda),
                                        FUN = mean, na.rm=T)
  
  crossFDTL_performances_data.frame <- bind_rows(crossFDTL_performances)
  mean_crossFDTL_performances_data.frame <- aggregate(crossFDTL_performances_data.frame, by = list(crossFDTL_performances_data.frame$iteration,
                                                                                                   crossFDTL_performances_data.frame$lambda),
                                                      FUN = mean, na.rm=T)
  
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
  
  for(i in 1:maximum_iteration){
    temp_performances <- mean_dtrace_performances[mean_dtrace_performances$iteration == i,
                                                  c(performances_colnames, "time")]
    performances <- rbind(performances, data.frame(temp_performances,
                                                   type=paste0("D-trace, Iteration ", i)))
    
    temp_performances <- mean_crossFDTL_performances_data.frame[mean_crossFDTL_performances_data.frame$iteration == i,
                                                                c(performances_colnames, "time")]
    performances <- rbind(performances, data.frame(temp_performances,
                                                   type=paste0("CrossFDTL Cversion, Iteration ", i)))
  }
  
  setting_results[[setting_index]] <- list(performances = performances)
}


#--------------------------------------------------------
#---------------- Visualization the Results -------------
#--------------------------------------------------------
#
performances <- NULL;
for(setting_index in 1:length(setting_results)){
  performances <- rbind(performances,
                        cbind(settings[setting_index,], setting_results[[setting_index]]$performances))
}

# rename the columns
colnames(performances)[colnames(performances)=="number_of_nodes"] = "d"
colnames(performances)[colnames(performances)=="number_of_changes"] = "s"
colnames(performances)[colnames(performances)=="number_of_samples"] = "m"
colnames(performances)[colnames(performances)=="maximum_legal_sparsity_in_solution_path"] = "c"
performances <- performances[order(performances$NRecall, performances$NPrecision),]

# Illustrating the accuracy comparison of D-trace, SP D-trace and CrossFDTL using Precision-Recall graph
solution_path_index <- performances$type == "Solution_Path"
non_solution_path_index <- performances$type == "D-trace, Iteration 1" |
  performances$type == "D-trace, Iteration 5" |
  performances$type == "CrossFDTL Cversion, Iteration 1" |
  performances$type == "CrossFDTL Cversion, Iteration 5" |
  performances$type == "CrossFDTL, Iteration 1" |
  performances$type == "CrossFDTL, Iteration 5" 
row_index <-  (non_solution_path_index | solution_path_index)

performances$algorithm <- sapply(str_split(performances$type, ", "),"[[",1)
performances$iteration <- 1
performances$iteration[performances$type != "Solution_Path"] <- sapply(str_split(performances$type[performances$type != "Solution_Path"], "Iteration "),"[[",2)

ROC_plot <- ggplot() +
  geom_line(mapping = aes(x = NRecall, y = NPrecision, color=algorithm, linetype =iteration),
            data = performances[row_index,], size = 0.5) +
  geom_point(mapping = aes(x = NRecall, y = NPrecision, color=algorithm),
             data = performances[non_solution_path_index,]) +
  scale_colour_discrete(limits = c("D-trace", "Solution_Path", "CrossFDTL Cversion"),
                        labels = c("ADMM D-trace", "Solution Path D-trace", "CrossFDTL")) +
  labs(y="Precision", x="Recall") +
  facet_grid(rows = vars(d), cols = vars(m),
             labeller = label_both) +
  theme_bw() +
  labs(linetype = "Number of iterations") + 
  labs(color = "Method") + 
  theme(legend.position = "top",
        legend.background = element_rect(color = "black"),
        panel.grid.minor = element_blank(),
        strip.background = element_blank())
print(ROC_plot)


# Illustrating the time comparison of D-trace, SP D-trace and CrossFDTL using Precision-Recall graph
performances_times <- performances[,c("type", "m", "d", "time")]
performances_times_NSP <- performances[performances$type != "Solution_Path",]
performances_times_NSP <- aggregate(x = performances_times_NSP$time, by = list(performances_times_NSP$type,
                                                                               performances_times_NSP$d,
                                                                               performances_times_NSP$m),
                                    FUN = sum)
colnames(performances_times_NSP) <- c("type", "d", "m", "time")
performances_times_NSP <- aggregate(x = performances_times_NSP$time, by = list(performances_times_NSP$type,
                                                                               performances_times_NSP$d),
                                    FUN = mean)
colnames(performances_times_NSP) <- c("type", "d", "time")

performances_times_NSP$algorithm <- sapply(str_split(performances_times_NSP$type, ", "),"[[",1)
performances_times_NSP$iteration <- sapply(str_split(performances_times_NSP$type, "Iteration "),"[[",2)

performances_times_solution_path <- performances[performances$type == "Solution_Path",]
performances_times_solution_path <- aggregate(x = performances_times_solution_path$time, by = list(performances_times_solution_path$type,
                                                                                                 performances_times_solution_path$d),
                                             FUN = mean)
colnames(performances_times_solution_path) <- c("type", "d", "time")

iterations <- unique(performances_times_NSP$iteration)
performances_times_solution_path <- merge(performances_times_solution_path, data.frame(algorithm = "Solution_Path",
                                                                                     iteration = iterations))
times <- rbind(performances_times_solution_path, performances_times_NSP)


times_plot <- ggplot()+
  geom_line(data = times,
            mapping = aes(x = as.numeric(iteration), y = time,
                          group=algorithm, color=algorithm),
            size = 0.5) +
  geom_point(data = times[times$algorithm != "Solution_Path",],
             mapping = aes(x = as.numeric(iteration), y = time,
                           group=algorithm, color=algorithm)) +
  labs(y="Time (sec.)", x="Number of iterations") +
  scale_colour_discrete(limits = c("D-trace", "Solution_Path",
                                   "CrossFDTL Cversion"),
                        labels = c("ADMM D-trace", "Solution Path D-trace",
                                   "CrossFDTL")) +
  facet_grid(rows = vars(d), labeller = label_both) +
  theme_bw() +
  theme(legend.position = "top",
        legend.background = element_rect(color = "black"),
        legend.title = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank())
print(times_plot)


ggarrange(ROC_plot, times_plot, 
          common.legend = T,
          # legend="bottom",
          labels = c("A", "B"),
          ncol = 2, nrow = 1, widths = c(2,1))

