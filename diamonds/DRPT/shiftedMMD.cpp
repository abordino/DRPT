#include <Rcpp.h>
#include "shiftedMMD.h"
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
      // Create a list of columns for r_func
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

// Bisection root finding
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
      return mid; // Root found within tolerance
    }

    if (f_lower * f_mid < 0)
    {
      upper = mid; // Root is in lower half
    }
    else
    {
      lower = mid; // Root is in upper half
      f_lower = f_mid;
    }
  }
  return mid; // Return best approximation after max iterations
}

// Define shifted MMD
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

  // Root finding setup
  SumLambda sum_lambda(Z, n, m, d, r);
  double lower_bound = 0.0, upper_bound = 100.0, tol = std::pow(std::numeric_limits<double>::epsilon(), 4);

  // Perform root finding using manual bisection method
  return bisection(sum_lambda, lower_bound, upper_bound, tol, 100);
}

// [[Rcpp::export]]
double compute_mmd_C(NumericMatrix X, NumericMatrix Y,
                     Function r_func, Function k_func)
{
  int n = X.nrow();
  int m = Y.nrow();
  double tau = n / m;
  int d = X.ncol();
  double lambda_star = compute_lambda_star_C(X, Y, r_func);

  Environment base = Environment::namespace_env("base");
  Function do_call = base["do.call"];

  // First term
  double first_term = 0.0;
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < n; j++)
    {
      if (i != j)
      {
        // Create a list of columns for r_func
        List Xi(d);
        List Xj(d);
        for (int k = 0; k < d; k++)
        {
          Xi[k] = X(i, k);
          Xj[k] = X(j, k);
        }

        double r_i = as<double>(do_call(r_func, Xi));
        double r_j = as<double>(do_call(r_func, Xj));
        double k_ij = as<double>(k_func(X(i, _), X(j, _)));

        double den = (tau + lambda_star * r_i) * (tau + lambda_star * r_j);
        double num = k_ij * r_i * r_j * lambda_star * lambda_star;

        first_term += num / den;
      }
    }
  }
  first_term /= (n * (n - 1));

  // Second term
  double second_term = 0.0;
  for (int i = 0; i < m; i++)
  {
    for (int j = 0; j < m; j++)
    {
      if (i != j)
      {
        // Create a list of columns for r_func
        List Yi(d);
        List Yj(d);
        for (int k = 0; k < d; k++)
        {
          Yi[k] = Y(i, k);
          Yj[k] = Y(j, k);
        }

        double r_i = as<double>(do_call(r_func, Yi));
        double r_j = as<double>(do_call(r_func, Yj));
        double k_ij = as<double>(k_func(Y(i, _), Y(j, _)));

        double den = (tau + lambda_star * r_i) * (tau + lambda_star * r_j);
        second_term += k_ij / den;
      }
    }
  }
  second_term /= (m * (m - 1));

  // Mixed term
  double mixed_term = 0.0;
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < m; j++)
    {
      // Create a list of columns for r_func
      List Xi(d);
      List Yj(d);
      for (int k = 0; k < d; k++)
      {
        Xi[k] = X(i, k);
        Yj[k] = Y(j, k);
      }

      double r_i = as<double>(do_call(r_func, Xi));
      double r_j = as<double>(do_call(r_func, Yj));
      double k_ij = as<double>(k_func(X(i, _), Y(j, _)));

      double den = (tau + lambda_star * r_i) * (tau + lambda_star * r_j);
      double num = k_ij * r_i * lambda_star;

      mixed_term += num / den;
    }
  }
  mixed_term /= (n * m);

  // Compute and return shifted MMD
  return first_term + second_term - 2 * mixed_term;
}
