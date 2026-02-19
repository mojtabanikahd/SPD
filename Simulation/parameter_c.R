source('./SPDtrace/Libraray.R')

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
maximum_legal_sparsity_in_solution_path = (1:8)*20;
number_of_samples = c(500)
number_of_repetition = 100;

settings <- expand.grid(number_of_nodes = number_of_nodes,
                        number_of_changes = number_of_changes,
                        maximum_legal_sparsity_in_solution_path = maximum_legal_sparsity_in_solution_path,
                        number_of_samples = number_of_samples)

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
  maximum_legal_sparsity_in_solution_path = settings$maximum_legal_sparsity_in_solution_path[setting_index]
  number_of_samples = settings$number_of_samples[setting_index]
  
  dtrace_solution_path_performances <- list();
  dtrace_solution_path_res <- list();
  dtrace_solution_path_times <- vector();
  dtrace_solution_path_mem_total_mb <- vector();
  dtrace_solution_path_knots <- vector();
  dtrace_solution_path_res <- list();

  
  for(k in 1:number_of_repetition) {
    # Generate Models
    reference <- generate_reference_models(
      number_of_nodes = number_of_nodes,
      number_of_samples = number_of_samples,
      number_of_changes = number_of_changes
    )
    print(paste("Repetition number: (", setting_index, " : ", k, ")"))
    
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
    
    save(sample_pearson_A,
         sample_pearson_B,
         maximum_legal_sparsity_in_solution_path,
         file = "./Simulation/Data/temp_input.RData")
    system("Rscript ./Simulation/solution_path_Rscript.R")
    # dtrace_solution_path_mem <- peakRAM({
    #   dtrace_solution_path_start_time <- Sys.time()
    #   dtrace_solution_path_results <- inference_Dtrace_solution_path(sample_covariances[[1]], sample_covariances[[2]],
    #                                                                  maximum_legal_sparsity_in_solution_path)
    #   dtrace_solution_path_stop_time <- Sys.time()
    # })
    # dtrace_solution_path_time <- dtrace_solution_path_stop_time - dtrace_solution_path_start_time
    load("./Simulation/Data/temp_output.RData")
    dtrace_solution_path_mem_total_mb[k] <- dtrace_solution_path_results$peak_MB;
    dtrace_solution_path_knots[k] <- length(dtrace_solution_path_results$solution_path);
    dtrace_solution_path_times[k] <- as.numeric(dtrace_solution_path_time);
    dtrace_solution_path_performances[[k]] <- solution_path_performance_evaluator(solution_path = dtrace_solution_path_results,
                                                                                  differential_structure = differential_structure)
  }
  
  # Averaging Evaluation Results
  performances_colnames <- colnames(dtrace_solution_path_performances[[1]])

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
                             time.se = sd(dtrace_solution_path_times) / sqrt(sum(!is.na(dtrace_solution_path_times))),
                             mem = mean(dtrace_solution_path_mem_total_mb),
                             mem.se = sd(dtrace_solution_path_mem_total_mb) / sqrt(sum(!is.na(dtrace_solution_path_mem_total_mb))),
                             knot = mean(dtrace_solution_path_knots),
                             knot.se = sd(dtrace_solution_path_knots) / sqrt(sum(!is.na(dtrace_solution_path_knots))))
  # colnames(performances) <- c(performances_colnames, "time", "mem", "knot")

  setting_results[[setting_index]] <- list(performances = performances)
}


#--------------------------------------------------------
#---------------- Visualization the Results -------------
#--------------------------------------------------------
#
library(RColorBrewer)
my_colors <- brewer.pal(3, "Set1")[1:2]  # Red and blue

d_colors <- c(
  "50" = my_colors[1],
  "100" = my_colors[2]
)

# Define consistent labels
d_labels <- c(
  "50" = "50",
  "100" = "100"
)


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
row_index <- (performances$c == max(performances$c))

max_lambdas <- performances %>%
  group_by(c,d) %>%
  summarise(min_lambda = min(lambda))

points <- data.frame()
max_c <- max(unique(max_lambdas$c))
for (i in 1:nrow(max_lambdas)) {
  c_value <- max_lambdas$c[i]
  d_value <- max_lambdas$d[i]
  l <- max_lambdas$min_lambda[i]
  per <- performances[performances$c == max_c & performances$d == d_value,]
  lam <- min(per$lambda[per$lambda >= l])
  points <- rbind(points,
                  cbind(c=c_value,
                        d=d_value,
                        lambda=lam,
                        NPrecision=per$NPrecision[per$lambda==lam],
                        NRecall=per$NRecall[per$lambda==lam]))
}

ROC_plot <- ggplot() +
  geom_line(mapping = aes(x = NRecall, y = NPrecision, color=as.factor(d)),
            data = performances[row_index,], size = 0.5) +
  geom_point(mapping = aes(x = NRecall, y = NPrecision,
                           color=as.factor(d)),
             data = points) +
  geom_text(mapping = aes(x = NRecall, y = NPrecision, label = c),
            data = points,
            hjust = 1.2,  # Adjust horizontal position
            vjust = 0.5,   # Adjust vertical position
            size = 3) +    # Adjust text size
  scale_colour_manual(values = d_colors,
                      limits = names(d_labels),
                      labels = d_labels) +
  labs(y="Precision", x="Recall",
       color = "Dimension") +
  theme_bw() +
  labs(linetype = "Number of iterations") + 
  theme(legend.position = "top",
        legend.background = element_rect(color = "black"),
        panel.grid.minor = element_blank(),
        strip.background = element_blank())
print(ROC_plot)

times_plot <- ggplot()+
  # geom_ribbon(
  #   mapping = aes(
  #     x = as.numeric(c),
  #     ymin = time - time.se,
  #     ymax = time + time.se,
  #     fill = as.factor(d),
  #   ),
  #   data = performances,
  #   alpha = 0.3,
  #   color = NA
  # ) +
  # geom_line(aes(x = as.numeric(c), y = time - time.se, color = as.factor(d)),
  #           data = performances,
  #           linetype = "dashed", alpha = 0.5) +
  # geom_line(aes(x = as.numeric(c), y = time + time.se, color = as.factor(d)),
  #           data = performances,
  #           linetype = "dashed", alpha = 0.5) +
  geom_line(data = performances,
            mapping = aes(x = as.numeric(c), y = time,
                          group=as.factor(d), color=as.factor(d)),
            size = 0.5) +
  geom_point(data = performances,
            mapping = aes(x = as.numeric(c), y = time,
                          group=as.factor(d), color=as.factor(d)),
            size = 1) +
  scale_y_log10(
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x))
  ) +
  labs(y="Time (sec.)", x="c",
       color = "Dimension") +
  scale_colour_manual(values = d_colors,
                      limits = names(d_labels),
                      labels = d_labels) +
  scale_fill_manual(values = d_colors,
                    labels = d_labels,
                    guide = "none") +
  theme_bw() +
  theme(legend.position = "top",
        legend.background = element_rect(color = "black"),
        panel.grid.minor = element_blank(),
        strip.background = element_blank())
print(times_plot)


mems_plot <- ggplot()+
  # geom_ribbon(
  #   mapping = aes(
  #     x = as.numeric(c),
  #     ymin = mem - mem.se,
  #     ymax = mem + mem.se,
  #     fill = as.factor(d),
  #   ),
  #   data = performances,
  #   alpha = 0.3,
  #   color = NA
  # ) +
  # geom_line(aes(x = as.numeric(c), y = mem - mem.se, color = as.factor(d)),
  #           data = performances,
  #           linetype = "dashed", alpha = 0.5) +
  # geom_line(aes(x = as.numeric(c), y = mem + mem.se, color = as.factor(d)),
  #           data = performances,
  #           linetype = "dashed", alpha = 0.5) +
  geom_line(data = performances,
            mapping = aes(x = as.numeric(c), y = mem,
                          group=as.factor(d), color=as.factor(d)),
            size = 0.5) +
  geom_point(data = performances,
            mapping = aes(x = as.numeric(c), y = mem,
                          group=as.factor(d), color=as.factor(d)),
            size = 1) +
  scale_y_log10(
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x))
  ) +
  labs(y="Memory (MB)", x="c",
       color = "Dimension") +
  scale_colour_manual(values = d_colors,
                      limits = names(d_labels),
                      labels = d_labels) +
  scale_fill_manual(values = d_colors,
                    labels = d_labels,
                    guide = "none") +
  theme_bw() +
  theme(legend.position = "top",
        legend.background = element_rect(color = "black"),
        panel.grid.minor = element_blank(),
        strip.background = element_blank())
print(mems_plot)

knots_plot <- ggplot()+
  # geom_ribbon(
  #   mapping = aes(
  #     x = as.numeric(c),
  #     ymin = knot - knot.se,
  #     ymax = knot + knot.se,
  #     fill = as.factor(d),
  #   ),
  #   data = performances,
  #   alpha = 0.3,
  #   color = NA
  # ) +
  # geom_line(aes(x = as.numeric(c), y = knot - knot.se, color = as.factor(d)),
  #           data = performances,
  #           linetype = "dashed", alpha = 0.5) +
  # geom_line(aes(x = as.numeric(c), y = knot + knot.se, color = as.factor(d)),
  #           data = performances,
  #           linetype = "dashed", alpha = 0.5) +
  geom_line(data = performances,
            mapping = aes(x = as.numeric(c), y = knot,
                          color=as.factor(d)),
            size = 0.5) +
  geom_point(data = performances,
            mapping = aes(x = as.numeric(c), y = knot,
                          color=as.factor(d)),
            size = 1) +
  labs(y="# Knots", x="c",
       color = "Dimension") +
  scale_colour_manual(values = d_colors,
                      limits = names(d_labels),
                      labels = d_labels) +
  scale_fill_manual(values = d_colors,
                    labels = d_labels,
                    guide = "none") +
  theme_bw() +
  theme(legend.position = "top",
        legend.background = element_rect(color = "black"),
        panel.grid.minor = element_blank(),
        strip.background = element_blank())
print(knots_plot)





final_performance_plot <- ggarrange(ROC_plot, times_plot, mems_plot, knots_plot,
                                    common.legend = T,
                                    # legend="bottom",
                                    labels = c("A", "B", "C", "D"),
                                    ncol = 4, nrow = 1, widths = c(1,1,1,1))
print(final_performance_plot)
ggsave("./Simulation/Results/parameter_c.pdf", plot = final_performance_plot, width = 12, height = 4)

