#include <Rcpp.h>
using namespace Rcpp;



//' @note Function for evaluation of the log likelihood of the second module in
//' the epidemiological example (up to a constant)
//' @param theta_1 Value of theta_1 from the first module
//' @param theta_2 Value of theta_2 (two-dim vector) at which theta_2 is evaluated
//' @param n_cases Vector with numbers of cases
//' @param follow_up Vector with numbers of follow up years

// [[Rcpp::export]]
double epidemiology_loglikC(NumericVector theta_1,
                            NumericVector theta_2,
                            NumericVector n_cases,
                            NumericVector follow_up){
  double res = 0;
  double log_mu, mu;

  int n2 = n_cases.length();
  int n1_check = follow_up.length();
  int n2_check = theta_1.length();

  if(n2!= n1_check){
    Rcerr << "Error\n";
    return(-9999.0);
  }

  if(n2!= n2_check){
    Rcerr << "Error\n";
    return(-9999.0);
  }

  for (int i = 0; i < n2; i++){
    // Use linear relationship, then evaluate mean of the Poisson distribution
    log_mu = theta_2[0] + theta_1[i]*theta_2[1] + follow_up[i];
    mu = exp(log_mu);
    // Evaluate log density of the Poisson distribution (up to a constant)
    res += n_cases[i]*log_mu - mu;
  }

  return res;
}

//' @note Function for evaluation of the log likelihood of the second module in
//' the epidemiological example (up to a constant)
//' @param theta_1 Value of theta_1 from the first module
//' @param theta_2 Value of theta_2 (two-dim vector) at which theta_2 is evaluated
//' @param weights Vector of WLB weights for the second module
//' @param n_cases Vector with numbers of cases
//' @param follow_up Vector with numbers of follow up years


// [[Rcpp::export]]
double epidemiology_weighted_loglikC(NumericVector theta_1,
                                     NumericVector theta_2,
                                     NumericVector weights,
                                     NumericVector n_cases,
                                     NumericVector follow_up){
  double res = 0;
  double log_mu, mu;

  int n2 = n_cases.length();
  int n1_check = follow_up.length();
  int n2_check = theta_1.length();

  if(n2!= n1_check){
    Rcerr << "Error\n";
    return(-9999.0);
  }

  if(n2!= n2_check){
    Rcerr << "Error\n";
    return(-9999.0);
  }
  for (int i = 0; i < n2; i++){
    log_mu = theta_2[0] + theta_1[i]*theta_2[1] + follow_up[i];
    mu = exp(log_mu);
    // Evaluate log density of the Poisson distribution (up to a constant)
    // and multiply by the weight
    res += weights[i]*(n_cases[i]*log_mu - mu);
  }

  return res;
}



