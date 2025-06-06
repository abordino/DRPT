#ifndef FIND_LAMBDA_H
#define FIND_LAMBDA_H

#include <Rcpp.h>

double compute_lambda_star_C(Rcpp::NumericMatrix X, Rcpp::NumericMatrix Y, Rcpp::Function r);

#endif