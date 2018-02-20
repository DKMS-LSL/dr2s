#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix cpp_hamming(NumericMatrix sequences) {
  int seqs = sequences.nrow();
  int pos = sequences.ncol();
  NumericMatrix result(seqs,seqs);
  for(int i = 0; i < seqs; ++i){
    NumericVector x = sequences(i,_);
    for (int j = 0; j < seqs; ++j){
      NumericVector y = sequences(j,_);
      double res = 0;
      for(int s = 0; s < pos; ++s){
        if (x[s] != y[s]){
          ++res;
        }
      }
      result(i,j) = res/pos;
    }
  }
  return(result);
}
