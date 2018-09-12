#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix pickCoordinates(unsigned Dim, unsigned N, unsigned fe, ListOf<NumericMatrix> VT, NumericMatrix U){
  NumericMatrix VTend(Dim,N);
  for(unsigned i=0; i<N; i++){
    NumericMatrix VTi = VT[i]; 
    for(unsigned j=0; j<Dim; j++){
      if(U(j,i) < 0.5){
        if(j < fe){
          VTend(j,i) = min(VTi(j,_));
        }else{
          double x = min(VTi(j,_));
          if(x < 0){
            x = 0;
          }
          VTend(j,i) = x;
        }
      }else{
        if(j < fe){
          VTend(j,i) = max(VTi(j,_));
        }else{
          double x = max(VTi(j,_));
          if(x < 0){
            x = 0;
          }
          VTend(j,i) = x;
        }
      }
    }
  }
  return VTend;
}