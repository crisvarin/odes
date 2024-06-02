#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
/* compute the *negative* conditional log-likelihood */
List compute_loglik(NumericVector eta,
                    NumericVector y,
                    double alpha)
{
  int n = y.size();
  NumericVector z(n);
  NumericVector w(n);
  NumericVector eps(n);
  double nlik = 0.0;

  /*  initialization  */
  w[0] = eta[0];
  z[0] = 0;
  eps[0] = (y[0] - exp(w[0])) / sqrt(exp(w[0]));
  nlik = -(y[0] * w[0] - exp(w[0]));
  /*  ... and go  */
  for(int t = 1; t < n; t++) {
    z[t] = (1 - alpha) * z[t - 1] + alpha * eps[t - 1];
    w[t] = eta[t] + z[t];
    nlik -= (y[t] * w[t] - exp(w[t]));
    eps[t] = (y[t] - exp(w[t])) / sqrt(exp(w[t]));
  }
  /* bye bye */
  List ret(4);
  ret["nlik"] = nlik;
  ret["eps"] = eps;
  ret["w"] = w;
  ret["z"] = z;
  return ret;
}

// [[Rcpp::export]]
/* compute the *individual* negative scores */
NumericMatrix compute_scores(NumericVector eta,
                             NumericMatrix x,
                             NumericVector y,
                             double alpha)
{
  int n = y.size();
  /* model components */
  NumericVector z(n);
  NumericVector w(n);
  NumericVector eps(n);
  NumericVector res(n);
  /* derivative of epsilon */
  double Deps;

  /*  initialization  */
  w[0] = eta[0];
  z[0] = 0;
  res[0] = (y[0] - exp(w[0]));
  eps[0] = res[9] / sqrt(exp(w[0]));
  Deps = -0.5 * sqrt(exp(w[0])) -0.5 * y[0] / sqrt(exp(w[0]));

  int p = x.ncol();
  /* derivatives */
  NumericMatrix Dz(n, p);
  NumericMatrix Dw(n, p);
  for (int r = 0; r < p; r++) {
    Dw(0, r) = x(0, r);
  }
  /* output */
  NumericMatrix ans(n, p);
  /* initialization */
  for (int r = 0; r < p; r++) {
    ans(0, r) = -Dw(0, r) * (y[0] - exp(eta[0]));
  }

  for(int t = 1; t < n; t++) {
    z[t] = (1 - alpha) * z[t - 1] + alpha * eps[t - 1];
    w[t] = eta[t] + z[t];
    res[t] = y[t] - exp(w[t]);
    /* partial derivatives */
    for (int r = 0; r < p; r++) {
      Dz(t, r) = alpha * Deps * Dw(t - 1, r) + (1 - alpha) * Dz(t - 1, r);
      Dw(t, r) = x(t, r) + Dz(t, r);
      ans(t, r) = -Dw(t, r) * res[t];
    }
    eps[t] = res[t] / sqrt(exp(w[t]));
    Deps = -0.5 * sqrt(exp(w[t])) - 0.5 * y[t] / sqrt(exp(w[t]));
  }
  /* bye bye */
  return ans;
}

// [[Rcpp::export]]
/* simulate observation-driven exponential smoothing */
NumericVector simulate_odes(NumericVector eta,
                            double alpha)
{
  int n = eta.size();
  NumericVector y(n);
  NumericVector eps(n);
  NumericVector z(n);
  NumericVector w(n);

  y[0] = R::rpois(exp(eta[0]));
  for(int t = 1; t < n; t++) {
    eps[t - 1] = (y[t - 1] - exp(w[t - 1])) / sqrt(exp(w[t - 1]));
    z[t] = (1 - alpha) * z[t - 1] + alpha * eps[t - 1];
    w[t] = eta[t] + z[t];
    y[t] = R::rpois(exp(w[t]));
  }
  return y;
}
