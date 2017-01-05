#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <queue>

// [[Rcpp::plugins(cpp11)]]

bool is_true (bool i) { return (i == true); }

// [[Rcpp::export]]
Rcpp::IntegerVector cpp_polymorphic_positions(Rcpp::NumericMatrix x, const double threshold) {
  auto nrows = x.rows();
  Rcpp::IntegerVector ppos;
  for (int i{}; i < nrows; ++i) {
    Rcpp::LogicalVector r = x.row(i) >= threshold;
    auto count = std::count_if(r.begin(), r.end(), is_true);
    if (count >= 2) {
      ppos.push_back(i + 1);
    }
  }
  return ppos;
}

bool desc (double i, double j) { return (i > j); }

// [[Rcpp::export]]
Rcpp::List cpp_top2_cols(Rcpp::NumericMatrix x) {
  auto nrows = x.rows();
  Rcpp::List out;
  Rcpp::IntegerVector i1, i2; // column indices for greatest/second greatest value
  Rcpp::NumericVector v1, v2; // greatest/second greatest value
  // in priority queues the top element is always the greatest element it contains.
  std::priority_queue<std::pair<double, int>> q;
  for (int i{}; i < nrows; ++i) {
    Rcpp::NumericVector r = x.row(i);
    // Push row into queue
    for (int j{}; j < r.size(); ++j) {
      q.push(std::pair<double, int>{r[j], j + 1});
    }
    // get first index/value pair
    i1.push_back(q.top().second);
    v1.push_back(q.top().first);
    q.pop();
    // get second index/value pair
    i2.push_back(q.top().second);
    v2.push_back(q.top().first);
    // reset priority queue
    q = std::priority_queue<std::pair<double, int>>{};
  }
  out["i1"] = i1;
  out["v1"] = v1;
  out["i2"] = i2;
  out["v2"] = v2;
  return out;
}
