#ifndef SHIFTED_MMD_H
#define SHIFTED_MMD_H

#include <Rcpp.h>

double compute_mmd_C(Rcpp::NumericMatrix X, Rcpp::NumericMatrix Y,
                     Rcpp::Function r_func, Rcpp::Function k_func);

#endif