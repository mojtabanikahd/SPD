# This file is an edited version retrieved from the DiffGraph package (DiffGraph.R file)

# -------------------------- The Expand Operator ----------------------------- #
# --- Expand_n(A, rho, n) = argmin_{Theta} -n\log\det(Theta) + (rho/2)*||Theta - A||_F^2    						       #	
# --- Let A = U*D*U' 		       #	
# Set D2_{ii} = 0.5*(D_{ii} + sqrt(D_{ii}^2 + 4n/rho));            #  			 
# Expand_n(A, rho,n) = U*D2*U'					       #	
#									       

expand = function(A, rho, n){
  edecomp = eigen(A)
  D = edecomp$values
  U = edecomp$vectors
  D2 = 0.5*(D + sqrt(D^2 + 4*n/rho ))
  Theta = U %*% diag(D2) %*% t(U)
  Theta
}



### infer differential network using FGL
edited_FGL = function(X, lambda1, lambda2, covType = "pearson", weights="equal",  penalize.diagonal = FALSE, tol = 1e-5, maxiter = 500, rho = 0.1, rho.incr = 1.05, rho.max = 1e5){
  
  # number of states
  G = length(X)
  # number of genes
  p = dim(X[[1]])[2]
  # sample sizes
  n = c()
  if (is.data.frame(X[[1]])){
    for (g in seq(G)){
      X[[g]] = as.matrix(X[[g]])
      n[g] = dim(X[[g]])[1]
    }
  }
  
  # assign gene names if none exist
  if(length(dimnames(X[[1]])[[2]])==0){
    for(g in seq(G)){
      dimnames(X[[g]])[[2]]=paste("V",1:p,sep="")
    }
  }
  
  
  # obtain the cov/cor matrices
  S = list()
  if (isSymmetric(X[[1]])){
    S = X
  }else{
    try(if (covType %in% c("pearson","kendall", "spearman") == FALSE) stop("The cov/cor type you provide is not include in this package. Please use your own function to obtain the list of cov/cor and use them as the input of FGL()"))
    for (g in seq(G)){
      S[[g]] = cov.compute(X[[g]], covType)
    }  
  }
  
  # set weights
  if(length(weights)==1){
    if(weights == "equal"){
      weights = rep(1,G)
    }
  }
  if(length(weights)==1){
    if(weights == "sample.size"){
      weights = n/sum(n)
    }
  }
  
  # solve optimization model using FGL.solve()
  result = FGL.solve(S, lambda1, lambda2, weights, penalize.diagonal, tol, maxiter, rho, rho.incr, rho.max)
  
  # create differential network from the estimated difference between two precision matrices, Delta = Theta_2^{-1} - Theta_1^{-1} 
  Delta.graph = (abs(result$Delta)>1e-5)*1
  diag(Delta.graph) = 0
  Degree = apply(Delta.graph, 1, sum)
  Delta.graph.connected = Delta.graph[Degree>0, Degree>0]
  result$Delta.graph.full =  graph_from_adjacency_matrix(Delta.graph, mode = "undirected", weighted = TRUE, diag = FALSE)
  result$Delta.graph.connected =  graph_from_adjacency_matrix(Delta.graph.connected, mode = "undirected", weighted = TRUE, diag = FALSE)
  
  # create state-specific gene networks from the estimated  precision matrices, Theta
  result$Theta.graph.full = list()
  result$Theta.graph.connected = list()
  for(g in 1:G){
    Theta.graph = (abs(result$Theta[[g]])>1e-5)*1
    diag(Theta.graph) = 0
    Degree = apply(Theta.graph, 1, sum)
    Theta.graph.connected = Theta.graph[Degree>0, Degree>0]
    result$Theta.graph.full[[g]] =  graph_from_adjacency_matrix(Theta.graph, mode = "undirected", weighted = TRUE, diag = FALSE)
    result$Theta.graph.connected[[g]] =  graph_from_adjacency_matrix(Theta.graph.connected, mode = "undirected", weighted = TRUE, diag = FALSE)
    
  }
  
  result   
}





### ADMM for FGL
FGL.solve = function(S, lambda1, lambda2, n = c(1,1), penalize.diagonal = FALSE, tol = 1e-5, maxiter = 500, rho = 0.1, rho.incr = 1.05, rho.max = 1e10){
  
  # number of genes
  p = dim(S[[1]])[1]
  G = length(S)
  
  # assign gene names if none exist:
  if(length(dimnames(S[[1]])[[2]])==0){
    for(g in seq(G)){
      dimnames(S[[g]])[[1]]=paste("V",1:p,sep="")
      dimnames(S[[g]])[[2]]=paste("V",1:p,sep="")
    }
  }
  
  # initialize 
  Theta = list()
  Z = list()
  U = list()
  Theta[[1]] = diag(p)
  Theta[[2]] = diag(p)
  Z[[1]] = diag(p)
  Z[[2]] = diag(p)
  U[[1]] = matrix(0, p, p)
  U[[2]] = matrix(0, p, p)
  
  iterationTimes <- NULL
  dataLogDelta <- list()
  iteration.start_time <- Sys.time()
  # ADMM iterations
  for (i in seq(maxiter)){
    
    # update theta
    Theta.prev = Theta
    Theta[[1]] = expand(Z[[1]] - (n[1]*S[[1]]+U[[1]])/rho, rho, n[1])
    Theta[[2]] = expand(Z[[2]] - (n[2]*S[[2]]+U[[2]])/rho, rho, n[2])
    
    # update Z
    A = list()
    A[[1]] = Theta[[1]] + U[[1]]/rho
    A[[2]] = Theta[[2]] + U[[2]]/rho
    Z = flsa2(A, rho, lambda1, lambda2, penalize.diagonal)
    
    iteration.stop_time <- Sys.time()
    iteration.time <- iteration.stop_time-iteration.start_time
    iterationTimes <- c(iterationTimes, iteration.time)
    dataLogDelta <- append(dataLogDelta, list(Delta = Z[[2]] - Z[[1]]))
    
    # update U
    U[[1]] = U[[1]] + rho*(Theta[[1]] - Z[[1]])
    U[[2]] = U[[2]] + rho*(Theta[[2]] - Z[[2]])
    
    # check the convergence condition
    diff1 = sum(abs(Theta[[1]] - Theta.prev[[1]])) + sum(abs(Theta[[2]] - Theta.prev[[2]]))
    diff2 = sum(abs(Theta[[1]] - Z[[1]])) + sum(abs(Theta[[2]] - Z[[2]]))
    norm_value = sum(abs(Theta[[1]])) + sum(abs(Theta[[2]]))
    if (max(diff1, diff2) < tol*max(norm_value, 1)){
      break
    }
    
    # update rho
    rho = min(rho*rho.incr,rho.max)
    
  }
  
  # compute precision matrix difference
  Theta = Z
  Delta = Theta[[2]] - Theta[[1]]
  diag(Delta) = 0
  
  # assign gene names
  row.names(Theta[[1]]) = row.names(S[[1]])
  colnames(Theta[[1]]) = colnames(S[[1]])
  row.names(Theta[[2]]) = row.names(S[[1]])
  colnames(Theta[[2]]) = colnames(S[[1]])
  row.names(Delta) = row.names(S[[1]])
  colnames(Delta) = colnames(S[[1]])
  
  result = list(Delta = Delta, Theta = Theta,
                estimated_delta_logs = dataLogDelta,
                iteration_times=iterationTimes)
  result
}




## -------------------- fused lasso signal approximator ------------------- ##
##
##---------------- minimize rho/2*(||Z_1 - A_1||_2^2 + ||Z_2 - A_2||_2^2)+ lambda1*(||Z_1||_1 + ||Z_2||_1) + lambda2*||Z_1 - Z_2||_1 -------- ##

flsa2 <-
  function(A, rho, lam1,lam2,penalize.diagonal)  #A is a list of 2 matrices from which we apply an L2 penalty to departures
  {
    
    S1 = abs(A[[1]]-A[[2]])<=2*lam2/rho
    X1 = (A[[1]]+A[[2]])/2
    Y1 = X1
    
    S2 = (A[[1]] > A[[2]]+2*lam2/rho)
    X2 = A[[1]] - lam2/rho
    Y2 = A[[2]] + lam2/rho
    
    S3 = (A[[2]] > A[[1]]+2*lam2/rho)
    X3 = A[[1]] + lam2/rho
    Y3 = A[[2]] - lam2/rho
    
    Z = list()
    Z[[1]] = soft(S1*X1 + S2*X2 + S3*X3, lam1/rho, penalize.diagonal)
    Z[[2]] = soft(S1*Y1 + S2*Y2 + S3*Y3, lam1/rho, penalize.diagonal)
    
    Z
  }





