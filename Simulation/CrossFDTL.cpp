// ---------------------------------------------------------------
// CrossFDTL.cpp  (corrected version)
// ---------------------------------------------------------------
#include <RcppArmadillo.h>
using namespace Rcpp;

// ---------------------------------------------------------------
// Helper: element‑wise max that returns a *value* (not a reference)
// ---------------------------------------------------------------
inline double max_val(double a, double b) { return (a < b) ? b : a; }

// ---------------------------------------------------------------
// Objective function
// ---------------------------------------------------------------
double LQ(const arma::sp_mat& Delta,
          const arma::mat&    A,
          const arma::mat&    B,
          double              lambda)
{
  double Q = 0.5 * arma::trace(Delta * Delta * A) + arma::trace(Delta * B);
  return Q + lambda * arma::accu(arma::abs(Delta));
}

// ---------------------------------------------------------------
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List CrossFDTL(NumericMatrix CovA, NumericMatrix CovB,
               double lambda, double rho, int maxiter)
{
  // -----------------------------------------------------------
  // Timing
  // -----------------------------------------------------------
  auto start_time = std::chrono::system_clock::now();
  
  const int nrow = CovA.nrow(), ncol = CovA.ncol(), nVar = CovA.ncol();
  
  // -----------------------------------------------------------
  // Convert R matrices to Armadillo
  // -----------------------------------------------------------
  arma::mat CA(nVar, nVar), CB(nVar, nVar);
  for (int i = 0; i < nrow; ++i)
    for (int j = 0; j < ncol; ++j) {
      CA(i, j) = CovA(i, j);
      CB(i, j) = CovB(i, j);
    }
    
  // -----------------------------------------------------------
  // Initialise variables
  // -----------------------------------------------------------
  arma::sp_mat Delta(nVar, nVar), Delta_pre(nVar, nVar), D(nVar, nVar);
  arma::mat    Sum = arma::zeros(nVar, nVar);
  arma::mat    dQ  = arma::zeros(nVar, nVar);
  arma::mat    A   = CA + CB + rho * arma::eye(nVar, nVar);
  arma::mat    B   = arma::zeros(nVar, nVar);
  arma::mat    C   = arma::zeros(nVar, nVar);
  
  // eigen‑decomposition of A
  arma::vec eigval;
  arma::mat eigvec;
  arma::eig_sym(eigval, eigvec, A, "std");
  eigvec = arma::reverse(eigvec, 1);
  eigval = arma::reverse(eigval);
  
  // matrix C (used later)
  for (int i = 0; i < nVar; ++i)
    for (int j = 0; j < nVar; ++j)
      C(i, j) = 2.0 / (eigval(i) + eigval(j));
  
  // -----------------------------------------------------------
  // Containers for results
  // -----------------------------------------------------------
  List Rresults;
  NumericVector iterationTimes;
  List dataLogData;
  
  // -----------------------------------------------------------
  // Main outer loop
  // -----------------------------------------------------------
  double iteration_time = 0.0;   // accumulated time
  
  for (int itr = 0; itr < maxiter; ++itr) {
    
    // ---- recompute B and Sum for the current Delta -------------
    B   = 2.0 * arma::eye(nVar, nVar) + 0.5 * (Delta * (CA - CB))
    + 0.5 * ((CA - CB) * Delta);
    Sum = eigvec * ((eigvec.t() * B * eigvec) % C) * eigvec.t();
    
    // ---- inner coordinate‑descent loop -------------------------
    int iteration = 0;
    while (true) {
      
      D.zeros();
      
      B  = 0.5 * (Sum * (CB - CA)) + 0.5 * ((CB - CA) * Sum);
      dQ = 0.5 * (Delta * A + A * Delta) + B;
      dQ = dQ.t();
      
      // ---- coordinate‑wise update ---------------------------------
      for (int i = 0; i < nVar; ++i) {
        for (int j = i + 1; j < nVar; ++j) {
          if ( (std::abs(dQ(i, j)) > lambda) || (Delta(i, j) != 0) ) {
            double a = A(i, i) + A(j, j);
            
            double b = arma::as_scalar(A.row(i) * Delta.col(j))
              + arma::as_scalar(A.row(j) * Delta.col(i))
              + B(i, j) + B(j, i);
              
              double c = Delta(i, j);
              double z = c - (b / a);
              double r = 2.0 * lambda / a;
              
              // sign(z) = (0<z)-(z<0)
              D(i, j) = D(j, i) =
                ((0 < z) - (z < 0)) * max_val(0.0, std::abs(z) - r) - c;
          }
        }
      }
      
      // ---- back‑tracking line search (Armijo) --------------------
      int    k_Armijo = 0;
      const double beta = 0.5;
      arma::sp_mat Delta_tilde(nVar, nVar);
      double LQ_Delta, LQ_Delta_tilde;
      
      while (true) {
        double alpha = std::pow(beta, k_Armijo);
        Delta_tilde = Delta + alpha * D;
        
        LQ_Delta        = LQ(Delta,        A, B, lambda);
        LQ_Delta_tilde  = LQ(Delta_tilde,  A, B, lambda);
        
        double Armijo_rule_delta =
          arma::trace(dQ * D) +
          lambda * arma::accu(arma::abs(Delta + D)) -
          lambda * arma::accu(arma::abs(Delta));
        
        if ( (LQ_Delta_tilde <= LQ_Delta) ||
             (LQ_Delta_tilde <= LQ_Delta + 0.5 * alpha * Armijo_rule_delta) )
          break;
        
        ++k_Armijo;
      }
      
      // ---- accept the step ---------------------------------------
      Delta = Delta_tilde;
      
      // ---- convergence test --------------------------------------
      const double epsilon = 1e-3;
      if ( std::abs(LQ_Delta_tilde - LQ_Delta) <=
           epsilon * std::max(std::abs(LQ_Delta), std::abs(LQ_Delta_tilde)) )
        break;
      
      ++iteration;
    }   // end inner while
    
    // ---- timing for this outer iteration -------------------------
    auto stop_time   = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed = stop_time - start_time;
    iteration_time = elapsed.count();
    
    iterationTimes.push_back(iteration_time);
    dataLogData.push_back(Delta);
  }
  
  // -----------------------------------------------------------
  // Return results
  // -----------------------------------------------------------
  Rresults["delta"]                = Delta;
  Rresults["estimated_delta_logs"] = dataLogData;
  Rresults["iteration_times"]      = iterationTimes;
  
  return Rresults;
}