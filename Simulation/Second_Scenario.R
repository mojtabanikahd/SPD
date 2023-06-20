rm(list = ls())
#Set working directory
source('../SPDtrace/Libraray.R')
Rcpp::sourceCpp('../SPDtrace/SolutionPath.cpp')
library(ggplot2)



######## performance comparision in multi-source setting (multi-source)
set.seed(1)
number_of_nodes_items = c(50)
number_of_changes_items = c(20, 10, 5)
number_of_tuning_parameters_items = c(20);
maximum_legal_sparsity_in_solution_path_items = c(200);
number_of_samples_items = c(200)
number_of_sources_items = c(4)
maxiter_items = c(10);
number_of_repetition = 100;



settings <- expand.grid(number_of_nodes_items = number_of_nodes_items,
                        number_of_changes_items = number_of_changes_items,
                        number_of_tuning_parameters_items = number_of_tuning_parameters_items,
                        maximum_legal_sparsity_in_solution_path_items = maximum_legal_sparsity_in_solution_path_items,
                        number_of_samples_items = number_of_samples_items,
                        number_of_sources_items = number_of_sources_items,
                        maxiter_items = maxiter_items)

setting_results <- list();




for(setting_index in 1:nrow(settings)){

  number_of_nodes = settings$number_of_nodes_items[setting_index]
  number_of_changes = settings$number_of_changes_items[setting_index]
  number_of_tuning_parameters = settings$number_of_tuning_parameters_items[setting_index]
  maximum_legal_sparsity_in_solution_path = settings$maximum_legal_sparsity_in_solution_path_items[setting_index]
  number_of_samples = settings$number_of_samples_items[setting_index]
  number_of_sources = settings$number_of_sources_items[setting_index]
  maxiter = settings$maxiter_items[setting_index]
  
  
  single_dtrace_solution_path_performances <- list();
  bulk_dtrace_solution_path_performances <- list();
  integrated_dtrace_solution_path_performances <- list();
  
  
  for(k in 1:number_of_repetition){
    reference <- generate_reference_models(
      number_of_nodes = number_of_nodes,
      number_of_samples = number_of_samples*number_of_sources,
      type = "ScaleFree",
      number_of_changes = number_of_changes
    )
    
    model_precision_A = reference$precision_matrix_A
    model_precision_B = reference$precision_matrix_B
    
    real_diff_support <- sign(model_precision_B - model_precision_A)
    diag(real_diff_support) <- 0
    
    # Bulk data
    bulk_sample_Kendall_A <- cor(reference$samples_A, method = "kendall");
    bulk_sample_Kendall_B <- cor(reference$samples_B, method = "kendall");
    
    bulk_sample_pearson_A = kendall_tau_matrix_to_pearson_correlation_matrix(bulk_sample_Kendall_A)
    bulk_sample_pearson_B = kendall_tau_matrix_to_pearson_correlation_matrix(bulk_sample_Kendall_B)
    
    # make positive semi-definite
    bulk_sample_pearson_A <- positive_semi_definite_maker(bulk_sample_pearson_A)
    bulk_sample_pearson_B <- positive_semi_definite_maker(bulk_sample_pearson_B)
    bulk_sample_covs <- list(bulk_sample_pearson_A, bulk_sample_pearson_B)
    
    
    # Single Source Data
    single_sample_Kendall_A <- cor(reference$samples_A[1:number_of_samples,], method = "kendall");
    single_sample_Kendall_B <- cor(reference$samples_B[1:number_of_samples,], method = "kendall");
    
    single_sample_pearson_A = kendall_tau_matrix_to_pearson_correlation_matrix(single_sample_Kendall_A)
    single_sample_pearson_B = kendall_tau_matrix_to_pearson_correlation_matrix(single_sample_Kendall_B)

    # make positive semi-definite
    single_sample_pearson_A <- positive_semi_definite_maker(single_sample_pearson_A)
    single_sample_pearson_B <- positive_semi_definite_maker(single_sample_pearson_B)
    single_sample_covs <- list(single_sample_pearson_A, single_sample_pearson_B)
    
    
    # Integration
    integrated_sample_Kendall_A <- matrix(data = 0, nrow = number_of_nodes, ncol = number_of_nodes)
    integrated_sample_Kendall_B <- matrix(data = 0, nrow = number_of_nodes, ncol = number_of_nodes)
    
    sample_pearsons <- NULL
    for(s in 1:number_of_sources){
      l_index <- (s-1)*number_of_samples + 1;
      u_index <- (s)*number_of_samples;
      
      sample_Kendall_view_A <- cor(reference$samples_A[l_index:u_index,], method = "kendall");
      sample_Kendall_view_B <- cor(reference$samples_B[l_index:u_index,], method = "kendall");
      
      integrated_sample_Kendall_A = integrated_sample_Kendall_A + sample_Kendall_view_A
      integrated_sample_Kendall_B = integrated_sample_Kendall_B + sample_Kendall_view_B
      
      sample_pearsonsViewA <- kendall_tau_matrix_to_pearson_correlation_matrix(sample_Kendall_view_A)
      sample_pearsonsViewB <- kendall_tau_matrix_to_pearson_correlation_matrix(sample_Kendall_view_B)
      # make positive semi definite
      sample_pearsonsViewA <- positive_semi_definite_maker(sample_pearsonsViewA)
      sample_pearsonsViewB <- positive_semi_definite_maker(sample_pearsonsViewB)
      
      sample_pearsons <- rbind(sample_pearsons, c(list(sample_pearsonsViewA),
                                                list(sample_pearsonsViewB)))
    }
    
    integrated_sample_Kendall_A = integrated_sample_Kendall_A/number_of_sources
    integrated_sample_Kendall_B = integrated_sample_Kendall_B/number_of_sources
    
    integrated_sample_pearson_A = kendall_tau_matrix_to_pearson_correlation_matrix(integrated_sample_Kendall_A)
    integrated_sample_pearson_B = kendall_tau_matrix_to_pearson_correlation_matrix(integrated_sample_Kendall_B)
    
    # make positive semi-definite
    integrated_sample_pearson_A <- positive_semi_definite_maker(integrated_sample_pearson_A)
    integrated_sample_pearson_B <- positive_semi_definite_maker(integrated_sample_pearson_B)
    
    integrated_sample_covs <- list(integrated_sample_pearson_A, integrated_sample_pearson_B)
    
    
    
    # differential network learning
    # solution path computation
    single_dtrace_solution_path_results <- inference_Dtrace_solution_path(single_sample_covs[[1]], single_sample_covs[[2]],
                                                                    sparsityLevel = maximum_legal_sparsity_in_solution_path)

    # Performance evaluation
    single_dtrace_solution_path_performances[[k]] <- solution_path_performance_evaluator(solution_path = single_dtrace_solution_path_results,
                                                                                   differential_structure = real_diff_support)
    
    
    # differential network learning
    # solution path computation
    bulk_dtrace_solution_path_results <- inference_Dtrace_solution_path(bulk_sample_covs[[1]], bulk_sample_covs[[2]],
                                                                  sparsityLevel = maximum_legal_sparsity_in_solution_path)

    # Performance evaluation
    bulk_dtrace_solution_path_performances[[k]] <- solution_path_performance_evaluator(solution_path = bulk_dtrace_solution_path_results,
                                                                                 differential_structure = real_diff_support)
    
    
    # differential network learning
    # solution path computation
    integrated_dtrace_solution_path_results <- inference_Dtrace_solution_path(integrated_sample_covs[[1]], integrated_sample_covs[[2]],
                                                                        sparsityLevel = maximum_legal_sparsity_in_solution_path)

    # Performance evaluation
    integrated_dtrace_solution_path_performances[[k]] <- solution_path_performance_evaluator(solution_path = integrated_dtrace_solution_path_results,
                                                                                             differential_structure = real_diff_support)
  }
  
  ## Create Data
  performances_colnames <- colnames(single_dtrace_solution_path_performances[[1]])
  
  mean_single_dtrace_solution_path_performances <- aggregate_dtrace_solution_path_resuts(number_of_repetition, single_dtrace_solution_path_performances)
  mean_bulk_dtrace_solution_path_performances <- aggregate_dtrace_solution_path_resuts(number_of_repetition, bulk_dtrace_solution_path_performances)
  mean_integrated_dtrace_solution_path_performances <- aggregate_dtrace_solution_path_resuts(number_of_repetition, integrated_dtrace_solution_path_performances)

  
  ### Create the diagrams
  # Data shaping
  performances <- data.frame(mean_single_dtrace_solution_path_performances, type="Single View")
  colnames(performances) <- c(performances_colnames, "type")
  performances <- rbind(performances, data.frame(mean_bulk_dtrace_solution_path_performances, type="Overal View"))
  performances <- rbind(performances, data.frame(mean_integrated_dtrace_solution_path_performances, type="Integrated View"))
  
  setting_results[[setting_index]] <- list(performances = performances)
}


# Visualizing the results
# Accuracy Curves
performances <- NULL;
for(setting_index in 1:length(setting_results)){
  performances <- rbind(performances,
                        cbind(settings[setting_index,], setting_results[[setting_index]]$performances))
}
# rename the columns
colnames(performances)[colnames(performances)=="number_of_nodes_items"] = "d"
colnames(performances)[colnames(performances)=="number_of_changes_items"] = "s"
colnames(performances)[colnames(performances)=="number_of_samples_items"] = "m"
colnames(performances)[colnames(performances)=="maximum_legal_sparsity_in_solution_path_items"] = "c"
performances <- performances[order(performances$NFPR, performances$NTPR),]


views <- expand.grid(number_of_changes_items,
                     maximum_legal_sparsity_in_solution_path_items)
performances <- performances[order(performances$NRecall, performances$NPrecision),]

v_index = 1
v = views[v_index,]
specific_change_index <- (performances$s == as.numeric(v[1]))
specific_control_index <- (performances$c == as.numeric(v[2]))
SPD_index <- (performances$type %in% c("Single View", "Overal View", "Integrated View"))
row_index <- specific_change_index & specific_control_index & SPD_index

pdf(file = "./Results/Multi_View/precision_vs_recall.pdf", width = 5.5, height = 3)
p <- ggplot(performances[row_index,], aes(x = NRecall, y = NPrecision, color=type))+
  geom_line(size = 0.5) +
  labs(y="Precision", x="Recall") +
  scale_colour_discrete(limits = c("Single View", "Overal View", "Integrated View"),
                        labels = c("Single Dataset", "Homogeneous Dataset", "Heterogeneous Datasets")) +
  theme_bw() +
  theme(legend.position = "top",
        legend.background = element_rect(color = "black"),
        legend.title = element_blank(),
        panel.grid.minor = element_blank(),
        # remove the facet background
        strip.background = element_blank())
print(p)
dev.off()
