#include <Rcpp.h>
#include "find_lambda.h"
using namespace Rcpp;

// Define sum_lambda function
struct SumLambda
{
  NumericMatrix Z;
  int n, m, d;
  Function r;

  SumLambda(NumericMatrix Z_, int n_, int m_, int d_, Function r_) : Z(Z_), n(n_), m(m_), d(d_), r(r_) {}

  double operator()(double lambda) const
  {
    Environment base = Environment::namespace_env("base");
    Function do_call = base["do.call"];

    double sum = 0.0;
    for (int k = 0; k < (n + m); k++)
    {
      List Zk(d);
      for (int j = 0; j < d; j++)
      {
        Zk[j] = Z(k, j);
      }
      double r_val = as<double>(do_call(r, Zk));
      sum += 1.0 / (n + m * lambda * r_val);
    }
    return sum - 1.0;
  }
};

// Bisection method
double bisection(SumLambda &func, double lower, double upper, double tol, int max_iter)
{
  double mid, f_lower, f_mid;
  f_lower = func(lower);

  for (int i = 0; i < max_iter; i++)
  {
    mid = (lower + upper) / 2.0;
    f_mid = func(mid);

    if (std::abs(f_mid) < tol || (upper - lower) / 2.0 < tol)
    {
      return mid;
    }

    if (f_lower * f_mid < 0)
    {
      upper = mid;
    }
    else
    {
      lower = mid; 
      f_lower = f_mid;
    }
  }
  return mid;
}

// [[Rcpp::export]]
double compute_lambda_star_C(Rcpp::NumericMatrix X, Rcpp::NumericMatrix Y,
                             Rcpp::Function r)
{
  int n = X.nrow();
  int m = Y.nrow();
  int d = X.ncol();

  // Combine X and Y into one matrix Z
  Rcpp::NumericMatrix Z(n + m, d);

  // Copy data from X and Y into Z
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < d; j++)
    {
      Z(i, j) = X(i, j);
    }
  }
  for (int i = 0; i < m; i++)
  {
    for (int j = 0; j < d; j++)
    {
      Z(n + i, j) = Y(i, j);
    }
  }

  // Find the root
  SumLambda sum_lambda(Z, n, m, d, r);
  double lower_bound = 0.0, upper_bound = 100.0, tol = std::pow(std::numeric_limits<double>::epsilon(), 4);

  return bisection(sum_lambda, lower_bound, upper_bound, tol, 100);
}
