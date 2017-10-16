#include <Rcpp.h>
using namespace Rcpp;

//NumericMatrix PSDM (NumericMatrix consensusMatrix, CharacterMatrix sequences) {
// [[Rcpp::export]]
NumericMatrix cpp_PSDM(NumericMatrix consMat, NumericMatrix sequences) {
  int seqs = sequences.nrow();
  int pos = sequences.ncol();
  NumericMatrix result(seqs,seqs);
  for(int i = 0; i < seqs; ++i){
    NumericVector x = sequences(i,_);
    for (int j = 0; j < seqs; ++j){
      NumericVector y = sequences(j,_);
      DoubleVector r(pos);
      int n = 0;
      for(int s = 0; s < pos; ++s){
        int xs = x[s]-1;
        int ys = y[s]-1;
        if (xs == ys){
         r[s] = 0.0;
          ++n;
        } else if (xs == 4 || ys == 4){
          r[s] = 0.0;
        } else{
          r[s] = consMat(xs,s) * consMat(ys,s);
           ++n;
         }
      }
      result(i,j) = sum(r)/n;//+r[1];
    }
  }
  return(result);
}
