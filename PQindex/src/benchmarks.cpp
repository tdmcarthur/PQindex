#include <RcppArmadillo.h>
#include <Rcpp.h>
// [[Rcpp::depends(RcppProgress)]]
// #include <progress.hpp>
using namespace Rcpp;
using namespace arma;
using namespace std;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
arma::sp_mat sp_mat_mult(arma::sp_mat x) {
 // return(log(conv_to<mat>::from(x * x.t()))) ;
 //  return(conv_to<mat>::from(x * x.t())) ;
     return(x * x.t()) ;
}

// [[Rcpp::export]]
arma::mat reg_mat_mult(arma::mat x) {
//  return(log(x * x.t()));
  return(x * x.t());
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R

*/
