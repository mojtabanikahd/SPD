#include <RcppArmadillo.h>
#include <List>
#include <cmath>

using namespace Rcpp;

bool AreSame(double a, double b, double EPSILON=0.00000001) {
  return fabs(a - b) < EPSILON;
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List inference_Dtrace_solution_path(NumericMatrix CovA, NumericMatrix CovB, int sparsityLevel) {
  int nrow = CovA.nrow(), ncol = CovA.ncol(), nVar = CovA.ncol();
  // Define int iterator
  std::list<int>::iterator it;
  
  // CA and CB are covariance matrises. D is vectorized of CB-CA
  double **CA = new double* [nrow];
  for(int i=0; i<nrow; i++)
    CA[i] = new double [nrow];
  
  double **CB = new double* [nrow];
  for(int i=0; i<nrow; i++)
    CB[i] = new double [nrow];
  
  double *D = new double [nrow*ncol];
  NumericVector RDV(nrow*ncol);
  
  List Rresults;
  
  // Initialize the difference of covariance matrices and activation indicator
  for(int i=0; i<nrow; i++) {
    for(int j=0; j<ncol; j++) {
      CA[i][j] = CovA(i,j);
      CB[i][j] = CovB(i,j);
      
      RDV(j*nrow + i) = D[j*nrow + i] = CB[i][j] - CA[i][j];
    }
  }
  Rresults["D"] = RDV;
  
  // main loop of algorithm
  List RsolutionPathList;
  
  // Declare essential variable for main loop
  // size of active set
  int activeSetSize = 0;
  double preLambda=2;
  std::list<int> preActiveSetIndices;
  std::list<int> preActiveSetSigns;
  
  while(activeSetSize <= sparsityLevel) {
    if(preLambda == 0)
      break;
    if(activeSetSize == 0) {
      // Identify the maximum absolute difference
      double maximumValue = 0;
      
      // Indices are related to the maximum value
      std::list<int> maxIndices;
      
      // Identify the maximum value and corresponding indices.
      for(int i=0; i<nrow*ncol; i++) {
        double absD = std::abs(D[i]);
        if(absD > maximumValue){
          maximumValue = absD;
          maxIndices.clear();
          maxIndices.push_back(i);
        }
        else if(absD == maximumValue){
          maxIndices.push_back(i);
        }
      }
      
      // Retrieve sign of activated variable
      std::list<int> signs;
      for (it=maxIndices.begin(); it!=maxIndices.end(); ++it)
        signs.push_back((D[*it] > 0) ? -1 : 1);
      
      //////////////////////////////
      // Configure Output results
      /////////////////////////////
      List Rknot;
      Rknot["knots_lambdas"] = maximumValue;
      preLambda = maximumValue;

      int indexSize = maxIndices.size();
      activeSetSize += indexSize;

      // filling active set indices
      NumericVector RactiveSetIndices(indexSize);
      int i=0;
      for (it=maxIndices.begin(); it!=maxIndices.end(); ++it){
        RactiveSetIndices(i) = *it;
        i++;
        preActiveSetIndices.push_back(*it);
      }
      Rknot["active_set"] = RactiveSetIndices;
      
      // filling active set sign
      NumericVector RactiveSetSigns(indexSize);
      i = 0;
      for (it=signs.begin(); it!=signs.end(); ++it){
        RactiveSetSigns(i) = *it;
        i++;
        preActiveSetSigns.push_back(*it);
      }
      Rknot["active_set_signs"] = RactiveSetSigns;
      
      // filling active set values
      NumericVector RactiveSetValues(indexSize);
      for (i=0; i<indexSize; ++i){
        RactiveSetValues(i) = 0;
      }
      Rknot["active_set_values"] = RactiveSetValues;
      
      RsolutionPathList.push_back(Rknot);
    }
    // There are some active variables
    else {
      // Allocate the memory for active matrix inversion Computation.
      arma::mat IKA = arma::ones(activeSetSize, activeSetSize);
      
      // Define int iterator on The Gamma matrix
      std::list<int>::iterator itc;
      itc = preActiveSetIndices.begin();
      for (int i = 0; i < activeSetSize; ++i) {
        // Compute the active matrix
        it = preActiveSetIndices.begin();
        for (int j = 0; j < activeSetSize; ++j) {
          int columnA = (*itc)/nVar, columnB = (*itc)%nVar;
          int rowA = (*it)/nVar, rowB = (*it)%nVar;
          double val = CA[rowA][columnA]*CB[rowB][columnB] + CB[rowA][columnA]*CA[rowB][columnB];
          IKA.at(i,j) = val;
          ++it;
        }
        ++itc;
      }
      
      if(!IKA.is_sympd()){
        Rcout<<"The problem hasn't any unique solution at lambda = "<<preLambda<<"\n";
        break;
      } else {
        IKA = inv_sympd(IKA);
      }
      // compute IKA * DA
      double* IKADA = new double[activeSetSize];
      for(int i=0; i<activeSetSize; ++i){
        IKADA[i]=0;
        it=preActiveSetIndices.begin();
        for(int j=0; j<activeSetSize; ++j){
          IKADA[i] += IKA.at(i,j)*D[*it];
          ++it;
        }
      }
      // compute IKA * SA
      double* IKASA = new double[activeSetSize];
      for(int i=0; i<activeSetSize; ++i){
        IKASA[i]=0;
        it=preActiveSetSigns.begin();
        for(int j=0; j<activeSetSize; ++j){
          IKASA[i] += IKA.at(i,j)*(*it);
          ++it;
        }
      }
      
      // Compute the next crossing time
      double crossingLambda = 0;
      std::list<int> crossingIndices;
      it = preActiveSetIndices.begin();
      for(int i=0; i<activeSetSize; i++){
        double tempLambda = -IKADA[i]/IKASA[i];
        if(tempLambda < preLambda && !AreSame(tempLambda, preLambda)){
          if(AreSame(tempLambda, crossingLambda)){
            crossingIndices.push_back(*it);
          }else if(tempLambda > crossingLambda){
            crossingIndices.clear();
            crossingIndices.push_back(*it);
            crossingLambda = tempLambda;
          }
        }
        it++;
      }
      
      // Compute the next hitting time
      double hittingLambda = 0;
      std::list<int> hittingIndices;
      std::list<int> hittingSigns;
      
      for(int i=0; i<nVar*nVar; i++){
        bool isActivated = FALSE;
        for(it = preActiveSetIndices.begin(); it!=preActiveSetIndices.end(); ++it){
          if(i == *it)
            isActivated = TRUE;
        }
        if(isActivated == FALSE){
          double KIKADA=0, KIKASA=0;
          int j =0;
          for(itc = preActiveSetIndices.begin(); itc!=preActiveSetIndices.end(); ++itc){
            int columnA = (*itc)/nVar, columnB = (*itc)%nVar;
            int rowA = (i)/nVar, rowB = (i)%nVar;
            double val = CA[rowA][columnA]*CB[rowB][columnB] + CB[rowA][columnA]*CA[rowB][columnB];
            KIKADA += val*IKADA[j];
            KIKASA += val*IKASA[j];
            j++;
          }
          // positive sign checking
          double tempLambda = (KIKADA-D[i])/(1 -KIKASA);
          if(tempLambda < preLambda && !AreSame(tempLambda, preLambda)) {
            if(AreSame(tempLambda, hittingLambda)) {
              hittingIndices.push_back(i);
              hittingSigns.push_back(1);
            } else if(tempLambda > hittingLambda) {
              hittingIndices.clear();
              hittingSigns.clear();
              hittingIndices.push_back(i);
              hittingSigns.push_back(1);
              hittingLambda = tempLambda;
            }
          }
          
          // negative sign checking
          tempLambda = (KIKADA-D[i])/(-1 -KIKASA);
          if(tempLambda < preLambda && !AreSame(tempLambda, preLambda)) {
            if(AreSame(tempLambda, hittingLambda)) {
              hittingIndices.push_back(i);
              hittingSigns.push_back(-1);
            } else if(tempLambda > hittingLambda) {
              hittingIndices.clear();
              hittingSigns.clear();
              hittingIndices.push_back(i);
              hittingSigns.push_back(-1);
              hittingLambda = tempLambda;
            }
          }
        }
      }
      
      //////////////////////////////
      // Configure Output results
      /////////////////////////////
      
      List Rknot;
      if(hittingLambda > crossingLambda) {
        Rknot["knots_lambdas"] = hittingLambda;
        preLambda = hittingLambda;
        int indexSize = hittingIndices.size();
        
        // filling active set values
        NumericVector RactiveSetValues(activeSetSize);
        for(int i=0; i<activeSetSize; i++) {
          RactiveSetValues(i) = -2*(IKADA[i] + (preLambda*IKASA[i]));
        }
        Rknot["active_set_values"] = RactiveSetValues;
        
        
        activeSetSize += indexSize;

        // filling active set indices
        for (it=hittingIndices.begin(); it!=hittingIndices.end(); ++it){
          preActiveSetIndices.push_back(*it);
        }
        
        // filling active set Index
        NumericVector RactiveSetIndices;
        for (it=preActiveSetIndices.begin(); it!=preActiveSetIndices.end(); ++it){
          RactiveSetIndices.push_back(*it);
        }
        Rknot["active_set"] = RactiveSetIndices;
        
        // filling active set sign
        for (it=hittingSigns.begin(); it!=hittingSigns.end(); ++it) {
          preActiveSetSigns.push_back(*it);
        }
        NumericVector RactiveSetSigns;
        for (it=preActiveSetSigns.begin(); it!=preActiveSetSigns.end(); ++it) {
          RactiveSetSigns.push_back(*it);
        }
        Rknot["active_set_signs"] = RactiveSetSigns;
      } else {
        Rknot["knots_lambdas"] = crossingLambda;
        preLambda = crossingLambda;
        
        int indexSize = crossingIndices.size();
        
        
        // filling active set values
        NumericVector RactiveSetValues(activeSetSize);
        for(int i=0; i<activeSetSize; i++){
          RactiveSetValues(i) = -2*(IKADA[i] + (preLambda*IKASA[i]));
        }
        Rknot["active_set_values"] = RactiveSetValues;
        
        activeSetSize -= indexSize;

        // filling active set indices
        for (it=crossingIndices.begin(); it!=crossingIndices.end(); ++it) {
          std::list<int>::iterator itIndices = preActiveSetIndices.begin();
          std::list<int>::iterator itSigns = preActiveSetSigns.begin();
          for(itIndices = preActiveSetIndices.begin(); itIndices!=preActiveSetIndices.end(); itIndices++) {
            if(*itIndices == *it) {
              preActiveSetIndices.erase(itIndices);
              preActiveSetSigns.erase(itSigns);
            }
            itSigns++;
          }
        }
        
        // filling active set Index
        NumericVector RactiveSetIndices;
        for (it=preActiveSetIndices.begin(); it!=preActiveSetIndices.end(); ++it){
          RactiveSetIndices.push_back(*it);
        }
        Rknot["active_set"] = RactiveSetIndices;
        
        // filling active set sign
        NumericVector RactiveSetSigns;
        for (it=preActiveSetSigns.begin(); it!=preActiveSetSigns.end(); ++it){
          RactiveSetSigns.push_back(*it);
        }
        Rknot["active_set_signs"] = RactiveSetSigns;
      }
      
      RsolutionPathList.push_back(Rknot);
      
      // Release the memory
      delete [] IKADA;
    }
  }
  
  delete[] D;
  
  for (int i = 0; i<nrow; i++){
    delete [] CA[i];
  }
  delete [] CA;
  
  for (int i = 0; i<nrow; i++){
    delete [] CB[i];
  }
  delete [] CB;
  
  Rresults["solution_path"] = RsolutionPathList;
  return Rresults;
}