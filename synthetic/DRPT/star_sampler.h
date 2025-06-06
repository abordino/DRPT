#ifndef STAR_SAMPLER_H
#define STAR_SAMPLER_H

#include <Rcpp.h>

Rcpp::List star_sampler_C(const Rcpp::NumericMatrix X,
                          const Rcpp::NumericMatrix Y,
                          int S, int H, Rcpp::Function r_func);

#endif // STAR_SAMPLER_H