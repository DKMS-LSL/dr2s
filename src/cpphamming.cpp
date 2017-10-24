#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

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
