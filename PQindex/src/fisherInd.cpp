#include <RcppArmadillo.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//



// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export]]
arma::mat fisherInd (arma::mat Q, arma::mat P, int base_period) {
  // Fisher index; direct approach; yes transitivity normalization
  //int M = Q.n_rows ;
  //int K = Q.n_cols ;
  //arma::mat check_nonneg_zeros(M, K) ;
  //check_nonneg_zeros.zeros() ;

  //arma::mat Q_nonneg = Q >= check_nonneg_zeros.zeros()
  //bool all_nonneg_check_Q = arma::as_scalar( arma::all(arma::all(Q_nonneg, 0), 1) ) ;
  // bool all_nonneg_check_Q = arma::as_scalar( arma::all(arma::all(Q >= 0, 0), 1) ) ;
  // bool all_nonneg_check_P = arma::as_scalar( arma::all(arma::all(P >= 0, 0), 1) )  ;
  //if ( ! (all_nonneg_check_Q && all_nonneg_check_P) ) {
  //  stop("All elements of the quantity and price matrices must be non-negative.") ;
  //}

  // Eventually I want to check that all elemets are posi, but I failed to figure out how to do it
  // Probably better to wrap those checks within R code.


  int M = Q.n_rows ;
  // double M_double = M ;
  arma::mat x(M, M) ;
  arma::vec I_row(M);
  arma::vec I_col(M);
  arma::vec ret(M);
  I_row.ones() ;
  I_col.ones() ;
  ret.ones() ;

  for (int col=0; col<M; col++) {
    double Q_ind_L = dot(P.row(base_period), Q.row(col)) / dot(P.row(base_period), Q.row(base_period)) ;
    double Q_ind_P = dot(P.row(col), Q.row(col)) / dot(P.row(col), Q.row(base_period)) ;
    I_row(col) = sqrt(Q_ind_L * Q_ind_P) ;
  }

  for (int col=0; col<M; col++) {

    for (int row=0; row<M; row++) {
      double Q_ind_L = dot(P.row(row), Q.row(col)) / dot(P.row(row), Q.row(row)) ;
      double Q_ind_P = dot(P.row(col), Q.row(col)) / dot(P.row(col), Q.row(row)) ;
      I_col(row) = sqrt(Q_ind_L * Q_ind_P) ;
    }

    arma::vec interm_vec_prod = I_row % I_col ;

      // Rcout << "The value is " << sum(interm_vec_prod) << std::endl ;
      if ( ! arma::is_finite(interm_vec_prod) ) {
        // Note that arma::is_finite checks the whole vector for any non-finite values
        stop("NaNs produced in quantity index. Check the quantity and price matrix inputs. Quantity indices must be positive, so the product of quantities and prices must be positive in all cases.") ;
        // Throws error in the case of zeros in the computed quantity index, which
        // will only occur in the case of all zeros in quantity
        // Thanks to explanation here: http://gallery.rcpp.org/articles/intro-to-exceptions/
      }
      ret(col) = exp(mean(log(interm_vec_prod))) ;
  }

  return(ret) ;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
timesTwo(42)
*/
