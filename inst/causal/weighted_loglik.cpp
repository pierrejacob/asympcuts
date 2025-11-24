#include <Rcpp.h>
using namespace Rcpp;


//' @note Function for evaluating weighted log likelihood for a logistic regression model
//' @param beta Parameter beta (including intercept)
//' @param y Outcome variable
//' @param x Matrix with covariates
//' @param weights Vector of weights
//' @param ncol Number of columns of x
//' @param nrow Number of rows of x
// [[Rcpp::export]]
double weighted_logistic_log_likC(NumericVector beta,
                                  NumericVector y,
                                  NumericMatrix x,
                                  NumericVector weights,
                                  int ncol,
                                  int nrow) {
  double out = 0;

  for (int i = 0; i < nrow; i++) {
    double prop_score = 0.;
    double exp_prop_score = 0.;
    for (int j = 0; j < ncol; j++) {
      prop_score += x(i, j)*beta[j];
    }

    exp_prop_score = exp(prop_score);
    if (y[i] == 1){
      out += weights[i]*log(exp_prop_score / (exp_prop_score + 1.));
    } else{
      out += weights[i]*log(1. / (exp_prop_score + 1.));
    }
  }

  return out;
}

//' @note Function for evaluating weighted log likelihood for a linear regression model
//' @param beta Parameter beta (including intercept)
//' @param y Outcome variable
//' @param x Matrix with covariates
//' @param weights Vector of weights
//' @param ncol Number of columns of x
//' @param nrow Number of rows of x
//'
// [[Rcpp::export]]
double weighted_linear_regressionC(NumericVector beta,
                                   NumericVector y,
                                   NumericMatrix x,
                                   NumericVector weights,
                                   int ncol,
                                   int nrow) {
  double out = 0;
  for (int i = 0; i < nrow; i++) {
    double temp = 0;
    for (int j = 0; j < ncol; j++) {
      temp += x(i, j)*beta[j];
    }
    out += weights[i] * 0.5 * pow(y[i]-temp, 2); // don't forget the 0.5 factor, not useful for MAP but useful for Laplace approx
  }
  //out = -out;
  return -out;
}



