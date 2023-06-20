source('../SPDtrace/Libraray.R')
Rcpp::sourceCpp('../SPDtrace/SolutionPath.cpp')

set.seed(1)

library(rmatio)
library(igraph)
library(readxl)



#--------------------------------------------------------
#------------------ Initialization ----------------------
#--------------------------------------------------------
#
number_of_random_set <- 10
random_set_ratio_of_StARS <- 0.8
maximum_legal_sparsity_in_solution_path <- 1000
threshold_of_StARS <- 0.001

# Retrieve name of important genes in platinum sensitivity mentioned by Ting Xu.
# Data retrieved from the following link.
# Link: https://github.com/XuTing95/WDtrace/tree/master/Ovarian_cancer
ov <- read.mat(filename = "./Data/OV_WDtrace.mat")
seven_important_pathway_gene_names <- unlist(ov$OV_WDtrace$Gene)

# Retrieve name of important genes in platinum sensitivity mentioned by Xiao-Fei Zhang.
# Load data. Data retrieved from the following link.
# Link: https://github.com/Zhangxf-ccnu/TDJGL
load("./Data/OV_PI3K.Akt.mTOR.rda")
PI3K_AKT_mTOR_pathway_gene_names <- colnames(OV_PI3K.Akt.mTOR$G450$Drug.Sensitivity)

# Integrate important gene names
the_most_important_gene_names <- sort(union(seven_important_pathway_gene_names, PI3K_AKT_mTOR_pathway_gene_names))

# Retrieve samples data (Gene expression + labels)
# Data retrieved from the following link.
# Link: https://github.com/Zhangxf-ccnu/TDJGL
load("./Data/OV_GeneExp.Rdata")

sensitive_cases_index <- (OV_GeneExp$Patient_drug_information[,"platinum_status"] == "Sensitive")
resistant_cases_index <- (OV_GeneExp$Patient_drug_information[,"platinum_status"] == "Resistant")

# # data grouping
# sensitive cases
sensitive_cases <- list()
sensitive_cases$G450 <- OV_GeneExp$G450[sensitive_cases_index, the_most_important_gene_names]
sensitive_cases$U133 <- OV_GeneExp$U133[sensitive_cases_index, the_most_important_gene_names]
sensitive_cases$HuEX <- OV_GeneExp$HuEx[sensitive_cases_index, the_most_important_gene_names]
# resistant cases
resistant_cases <- list()
resistant_cases$G450 <- OV_GeneExp$G450[resistant_cases_index, the_most_important_gene_names]
resistant_cases$U133 <- OV_GeneExp$U133[resistant_cases_index, the_most_important_gene_names]
resistant_cases$HuEX <- OV_GeneExp$HuEx[resistant_cases_index, the_most_important_gene_names]

number_of_nodes <- ncol(sensitive_cases$G450)



#--------------------------------------------------------
#------ Parameter Tuning (Adopted StARS mothod) ---------
#--------------------------------------------------------
#
solution_pathes_in_validation_step <- list();

# Generate validation models using applying SPDtrace method on random sub sets
for(i in 1:number_of_random_set) {
  # Calculate pearson correlation using a random sub set of sensitive cases
  pearson_of_sensitive_cases <- weighted_rank_based_pearson_correlation_estimator(datasets = sensitive_cases,
                                                                                  random_set_ratio = random_set_ratio_of_StARS)
  # Calculate pearson correlation using a random sub set of resistant cases
  pearson_of_resistant_cases <- weighted_rank_based_pearson_correlation_estimator(datasets = resistant_cases,
                                                                                  random_set_ratio = random_set_ratio_of_StARS)
  
  result <- inference_Dtrace_solution_path(pearson_of_sensitive_cases,
                                        pearson_of_resistant_cases,
                                        sparsityLevel = maximum_legal_sparsity_in_solution_path)
  result <- result$solution_path
  solution_pathes_in_validation_step[[i]] <- result
}

instability <- instability_evaluator_of_solution_pathes(solution_pathes_in_validation_step,
                                                        number_of_nodes)

selected_lambda <- min(instability$lambda[instability$instability < threshold_of_StARS])


#--------------------------------------------------------
#------------- Solution Path Inference ------------------
#--------------------------------------------------------
#
pearson_of_sensitive_cases <- weighted_rank_based_pearson_correlation_estimator(sensitive_cases)
pearson_of_resistant_cases <- weighted_rank_based_pearson_correlation_estimator(resistant_cases)

result <- inference_Dtrace_solution_path(pearson_of_sensitive_cases,
                                      pearson_of_resistant_cases,
                                      sparsityLevel = maximum_legal_sparsity_in_solution_path)
result <- result$solution_path

lambda_index <- sum(unlist(lapply(result, function(x) x$knots_lambdas > selected_lambda)))

if(lambda_index == 0) {
  warning("The identified differential network is empty!")
}

#--------------------------------------------------------
#-------- Differential Network Visualization ------------
#--------------------------------------------------------
# 
number_of_genes <- ncol(sensitive_cases$G450)
from_gene <- as.integer(result[[lambda_index]]$active_set/number_of_genes)+1
to_gene <- as.integer(result[[lambda_index]]$active_set%%number_of_genes)+1

# remove self loop
selfLoopIndicator <- (from_gene - to_gene) == 0
from_gene <- from_gene[!selfLoopIndicator]
to_gene <- to_gene[!selfLoopIndicator]

g <- graph(edges = the_most_important_gene_names[from_gene], directed = F)

deg <- degree(g, mode = "all")

pdf("./Results/DifferentialNetwork.pdf", width = 4)
print(plot(g, vertex.size=deg*5, vertex.label.cex=0.25,
           vertex.label.color="black",
           vertex.color= c(adjustcolor("green", alpha.f = 0.1),
                           adjustcolor("blue", alpha.f = 0.1),
                           adjustcolor("red", alpha.f = 0.1))[deg],
           vertex.frame.color = NA, layout =  layout_with_kk,
           edge.width = 0.5))
dev.off()


#--------------------------------------------------------
#----------------- Fisher's exact test ------------------
#--------------------------------------------------------
# 
# platinum_resistent_genes retrieved from http://ptrc-ddr.cptac-data-view.org/#/
platinum_resistent_genes <- read_excel("./Data/Platinum_Resistent_Genes.xlsx")
PR_genes <- platinum_resistent_genes$`HUGO Gene symbol`
PR0 = fisher_test(the_most_important_gene_names, names(deg)[deg>0], PR_genes)
PR1 = fisher_test(the_most_important_gene_names, names(deg)[deg>1], PR_genes)
