#include <Rcpp.h>
#include "star_sampler.h"

using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List star_sampler_C(const Rcpp::NumericMatrix X,
                          const Rcpp::NumericMatrix Y,
                          int S, int H, Rcpp::Function r_func)
{

  int n = X.nrow();
  int m = Y.nrow();
  int h = std::min(n, m);
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

  // Save original copy of Z as Z^{(0)}
  Rcpp::NumericMatrix Z_original = clone(Z);

  // Function to swap elements based on r() ratio - Algorithm 1
  auto swap_based_on_odds_ratio = [&](std::vector<int> &ind_n, std::vector<int> &ind_m)
  {
    int len = ind_n.size();

    // Create a list of columns for r_func
    List r_n_args(d);
    List r_m_args(d);
    for (int j = 0; j < d; j++)
    {
      r_n_args[j] = Z(Rcpp::Range(0, n - 1), Rcpp::Range(j, j));     
      r_m_args[j] = Z(Rcpp::Range(n, n + m - 1), Rcpp::Range(j, j));
    }

    Environment base = Environment::namespace_env("base");
    Function do_call = base["do.call"];

    // Call r_func with extracted columns
    NumericVector r_n = as<NumericVector>(do_call(r_func, r_n_args));
    NumericVector r_m = as<NumericVector>(do_call(r_func, r_m_args));

    for (int i = 0; i < len; i++)
    {
      if (ind_n[i] < n && ind_m[i] >= n && ind_m[i] < (n + m))
      {
        double odds_ratio = r_n[ind_n[i]] / r_m[ind_m[i] - n];
        if (R::runif(0, 1) < odds_ratio / (1 + odds_ratio))
        {
          for (int j = 0; j < d; j++)
          {
            std::swap(Z(ind_n[i], j), Z(ind_m[i], j));
          }
        }
      }
    }
  };

  // **Main Sampling Loop** - Generate Z^ast
  for (int j = 0; j < S; j++)
  {
    std::vector<int> ind_n = Rcpp::as<std::vector<int>>(Rcpp::sample(n, h, false));
    std::vector<int> ind_m = Rcpp::as<std::vector<int>>(Rcpp::sample(m, h, false));

    for (int i = 0; i < h; i++)
    {
      ind_m[i] += n; // Shift `m` indices to match Y range
    }

    swap_based_on_odds_ratio(ind_n, ind_m);
  }

  // Store intermediate `Z_star`
  Rcpp::NumericMatrix Z_star = clone(Z);

  // **Preallocate List** to store the Z^{(h)}'s
  Rcpp::List exch_Z(H + 1);
  exch_Z[0] = Z_original; // Store original data in the first position

  for (int b = 1; b <= H; b++)
  {
    for (int j = 0; j < S; j++)
    {
      std::vector<int> ind_n = Rcpp::as<std::vector<int>>(Rcpp::sample(n, h, false));
      std::vector<int> ind_m = Rcpp::as<std::vector<int>>(Rcpp::sample(m, h, false));

      for (int i = 0; i < h; i++)
      {
        ind_m[i] += n; // Shift `m` indices to match Y range
      }

      swap_based_on_odds_ratio(ind_n, ind_m);
    }
    exch_Z[b] = clone(Z); // Store results into Z^{(h)}
    Z = clone(Z_star);    // Reset Z to Z^star for next iteration
  }

  return exch_Z;
}