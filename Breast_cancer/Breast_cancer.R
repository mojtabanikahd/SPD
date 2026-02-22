source('./SPDtrace/Libraray.R')
Rcpp::sourceCpp('./SPDtrace/SolutionPath.cpp')

set.seed(1)

# 1. Load Essential Libraries
library(curatedTCGAData)
library(TCGAutils)
library(KEGGREST)
library(SummarizedExperiment)
library(rmatio)
library(igraph)
library(readxl)
library(ggplot2)



#--------------------------------------------------------
#------------------ Initialization ----------------------
#--------------------------------------------------------
#
number_of_random_set <- 10
random_set_ratio_of_StARS <- 0.8
maximum_legal_sparsity_in_solution_path <- 100
threshold_of_StARS <- 0.001

#--------------------------------------------------------
#------------------ Data Loading ------------------------
#--------------------------------------------------------
#
# 2. Download Breast Cancer Data (RNA-seq and Microarray)
# We download both assays at once into one MAE object
brca_mae <- curatedTCGAData(
  diseaseCode = "BRCA", 
  assays = c("RNASeqGene", "mRNAArray"), 
  version = "2.1.1", 
  dry.run = FALSE
)

# 3. Get KEGG Breast Cancer Pathway Genes (hsa05224)
kegg_data <- keggGet("hsa05224")[[1]]$GENE
# Extract symbols (KEGG provides them in a list where symbols are even-indexed)
gene_symbols_raw <- kegg_data[seq(2, length(kegg_data), by = 2)]
kegg_breast_genes <- sapply(strsplit(gene_symbols_raw, ";"), `[`, 1)
selected_genes <- kegg_breast_genes

# 4. Filter by PAM50 Subtype
clinical <- as.data.frame(colData(brca_mae))
lumA_ids <- rownames(clinical)[which(clinical$PAM50.mRNA == "Luminal A")]
basal_ids <- rownames(clinical)[which(clinical$PAM50.mRNA == "Basal-like")]

# 5. Create the Filtered Lists
# Function to extract and filter matrices for a given set of patient IDs
get_pathway_list <- function(mae, patient_ids, gene_list) {
  # Subset MAE for specific patients
  subset_mae <- mae[, patient_ids, ]
  
  # Extract RNA-seq matrix and filter genes
  rna_mat <- assays(subset_mae)[["BRCA_RNASeqGene-20160128"]]
  common_rna <- intersect(rownames(rna_mat), gene_list)
  
  # Extract Microarray matrix and filter genes
  micro_mat <- assays(subset_mae)[["BRCA_mRNAArray-20160128"]]
  common_genes <- intersect(rownames(micro_mat), common_rna)
  
  rna_final <- rna_mat[common_genes, ]
  micro_final <- micro_mat[common_genes, ]
  
  return(list(rnaseq = t(rna_final), microarray = t(micro_final)))
}

# Generate the final lists
LumA_cases <- get_pathway_list(brca_mae, lumA_ids, selected_genes)
Basal_cases <- get_pathway_list(brca_mae, basal_ids, selected_genes)

# 6. Verify Results
cat("Luminal A - RNA-seq patients:", nrow(LumA_cases$rnaseq), "\n")
cat("Luminal A - Microarray patients:", nrow(LumA_cases$microarray), "\n")
cat("Basal - RNA-seq patients:", nrow(Basal_cases$rnaseq), "\n")
cat("Basal - Microarray patients:", nrow(Basal_cases$microarray), "\n")

# Integrate important gene names
the_most_important_gene_names <- colnames(Basal_cases$microarray)
number_of_nodes <- ncol(Basal_cases$microarray)

#--------------------------------------------------------
#------ Parameter Tuning (Adopted StARS mothod) ---------
#--------------------------------------------------------
#
solution_pathes_in_validation_step <- list();

# Generate validation models using applying SPDtrace method on random sub sets
for(i in 1:number_of_random_set) {
  # Calculate pearson correlation using a random sub set of LumA cases
  pearson_of_LumA_cases <- weighted_rank_based_pearson_correlation_estimator(datasets = LumA_cases,
                                                                             random_set_ratio = random_set_ratio_of_StARS)
  # Calculate pearson correlation using a random sub set of Basal cases
  pearson_of_Basal_cases <- weighted_rank_based_pearson_correlation_estimator(datasets = Basal_cases,
                                                                              random_set_ratio = random_set_ratio_of_StARS)
  
  result <- inference_Dtrace_solution_path(pearson_of_LumA_cases,
                                           pearson_of_Basal_cases,
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
pearson_of_LumA_cases <- weighted_rank_based_pearson_correlation_estimator(LumA_cases)
pearson_of_Basal_cases <- weighted_rank_based_pearson_correlation_estimator(Basal_cases)

result <- inference_Dtrace_solution_path(pearson_of_LumA_cases,
                                         pearson_of_Basal_cases,
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
number_of_genes <- ncol(LumA_cases$rnaseq)
from_gene <- as.integer(result[[lambda_index]]$active_set/number_of_genes)+1
to_gene <- as.integer(result[[lambda_index]]$active_set%%number_of_genes)+1

# remove self loop
selfLoopIndicator <- (from_gene - to_gene) == 0
from_gene <- from_gene[!selfLoopIndicator]
to_gene <- to_gene[!selfLoopIndicator]

g <- graph(edges = the_most_important_gene_names[from_gene], directed = F)

deg <- degree(g, mode = "all")


pdf("./Breast_cancer/Results/DifferentialNetworkBreastCancer.pdf",
    width = 2, height = 2)
par(mar = c(0, 0, 0, 0))  # remove margins
print(plot(g, vertex.size=deg*5, vertex.label.cex=0.25,
           vertex.label.color="black",
           vertex.color= c(adjustcolor("green", alpha.f = 0.1),
                           adjustcolor("blue", alpha.f = 0.1),
                           adjustcolor("red", alpha.f = 0.1))[deg],
           vertex.frame.color = NA,
           edge.width = 0.5))
dev.off()

# number of nodes
vcount(g)
# number of edges
ecount(g)

#--------------------------------------------------------
#----------------- Enrichment Analysis ------------------
#--------------------------------------------------------
# 
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)

network_genes <- names(deg)
background_genes <- the_most_important_gene_names

network_entrez <- bitr(
  network_genes,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

background_entrez <- bitr(
  background_genes,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

library(msigdbr)

# Using C2 (Curated) gene sets - Chemical and Genetic Perturbations (CGP)
c2_cgp_genesets <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CGP")
c2_cgp_list <- split(c2_cgp_genesets$entrez_gene, c2_cgp_genesets$gs_name)

ec2cgp <- enricher(
  gene = network_entrez$ENTREZID,
  universe = background_entrez$ENTREZID,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  TERM2GENE = data.frame(
    term = c2_cgp_genesets$gs_name,
    gene = c2_cgp_genesets$entrez_gene
  )
)

dotplot(ec2cgp, showCategory = 10) + ggtitle("C2 CGP Gene Sets")