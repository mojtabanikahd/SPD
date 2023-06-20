# include <RcppArmadillo.h>


#include <List>
#include <cmath>

#include <chrono>
#include <ctime> 

#include <stdlib.h>


using namespace Rcpp;


const double& max (const double& a, const double& b) {
  return (a<b)?b:a;     // or: return comp(a,b)?b:a; for version (2)
}

double LQ(arma::sp_mat Delta, arma::mat A, arma::mat B, double lambda){
  double Q = 0.5*arma::trace(Delta*Delta*A) + arma::trace(Delta*B);
  return Q + lambda*arma::accu(abs(Delta));
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List CrossFDTL(NumericMatrix CovA, NumericMatrix CovB,
               double lambda, double rho, int maxiter) {
  auto start_time = std::chrono::system_clock::now();
  
  
  // numboer of row and column
  int nrow = CovA.nrow(), ncol = CovA.ncol(), nVar = CovA.ncol();
  
  // Define int iterator
  std::list<int>::iterator it;
  
  List Rresults;
  NumericVector iterationTimes;
  List dataLogData;
  
  // CA and CB are covariance matrises.
  arma::mat CA(nVar, nVar), CB(nVar,nVar);
  
  // Initialize covariance matrices
  for(int i=0; i<nrow; i++){
    for(int j=0; j<ncol; j++){
      CA.at(i,j) = CovA(i,j);
      CB.at(i,j) = CovB(i,j);
    }
  }
  
  
  
  arma::sp_mat Delta(nVar, nVar);
  arma::sp_mat Delta_pre(nVar, nVar);
  arma::sp_mat D(nVar, nVar);
  arma::mat Sum = arma::zeros(nVar, nVar);
  arma::mat dQ(nVar, nVar);
  
  
  arma::mat A = CA + CB + rho*arma::eye(nVar, nVar);
  arma::mat B(nVar, nVar);
  
  arma::vec eigval;
  arma::mat eigvec;
  
  eig_sym(eigval, eigvec, A, "std");
  eigvec = arma::reverse(eigvec,1);
  eigval = arma::reverse(eigval);
  
  arma::mat C = arma::mat(nVar, nVar);
  // Initialize matrix A
  for(int i=0; i<nrow; i++){
    for(int j=0; j<ncol; j++){
      C.at(i,j) = 2/(eigval.at(i) + eigval.at(j));
    }
  }
  
  auto stop_time = std::chrono::system_clock::now();
  
  std::chrono::duration<double> elapsed_time = stop_time-start_time;
  double iteration_time = elapsed_time.count();  
  
  for(int itr=0; itr<maxiter; itr++){
    start_time = std::chrono::system_clock::now();
    B = 2*arma::eye(nVar,nVar) + 0.5*( Delta*(CA-CB) ) + 0.5*( (CA-CB)*Delta);    
    Sum = eigvec * ((eigvec.t()*B*eigvec)%C) * eigvec.t();
    
    
    int iteration = 0;
    while(true){
      D.zeros();
      B = 0.5*( Sum*(CB-CA) ) + 0.5*( (CB-CA)*Sum );
      dQ = 0.5*(Delta*A + A*Delta) + B;
      dQ = dQ.t();
      
      for(int i=0; i<nrow; i++){
        for(int j=i+1; j<ncol; j++){
          if( (abs(dQ.at(i,j)) > lambda) | (Delta.at(i,j) != 0)){
            double a = A.at(i,i) + A.at(j,j);
            double b = (A.row(i)*Delta.col(j)).at(0,0) + (A.row(j)*Delta.col(i)).at(0,0) + B.at(i,j) + B.at(j,i);
            double c = Delta.at(i,j);
            
            double z = c-(b/a);
            double r = 2*lambda/a;
            
            // ((0<z)-(z<0)) is equal to sign(z)
            D.at(i,j) = D.at(j,i) = ((0<z)-(z<0))*max(0, abs(z)-r) - c;            
          }
        }
      }
      
      int k_Arimijo = 0;
      double beta = 0.5;
      arma::sp_mat Delta_tilde(nVar, nVar);
      
      double LQ_Delta;
      double LQ_Delta_tilde;
      
      while(TRUE){
        double alpha = pow(beta, k_Arimijo);
        Delta_tilde = Delta + alpha*D;
        
        LQ_Delta = LQ(Delta, A, B, lambda);
        LQ_Delta_tilde = LQ(Delta_tilde, A, B, lambda);
        
        double Armijo_rule_delta = arma::trace(dQ*D) + lambda*arma::accu(abs(Delta+D)) - lambda*arma::accu(abs(Delta));
        
        if(LQ_Delta_tilde <= LQ_Delta | LQ_Delta_tilde <= LQ_Delta + 0.5*alpha*Armijo_rule_delta){
          break;
        }
        k_Arimijo = k_Arimijo+1;
      }
      
      Delta = Delta_tilde;
      
      double epsilon = 0.001;
      if(abs(LQ_Delta_tilde - LQ_Delta) <= epsilon*max(abs(LQ_Delta), abs(LQ_Delta_tilde)))
        break;
      
      iteration++;
    }
    
    
    stop_time = std::chrono::system_clock::now();
    elapsed_time = stop_time-start_time;
    iteration_time = iteration_time + elapsed_time.count();
    
    iterationTimes.push_back(iteration_time);
    dataLogData.push_back(Delta);
  }
  
  
  Rresults["delta"] = Delta;
  Rresults["estimated_delta_logs"] = dataLogData;
  Rresults["iteration_times"] = iterationTimes;
  
  return Rresults;
}
