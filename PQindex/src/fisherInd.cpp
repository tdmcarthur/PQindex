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

  base_period = base_period - 1 ;
  // Since C++ is zero-indexed
  int M = Q.n_rows ;
  // double M_double = M ;
  // arma::mat x(M, M) ;
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




// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat fisherIndfast (arma::mat Q,         arma::mat P,
                         arma::mat Q_consol,  arma::mat P_consol,
                         arma::mat Q_freq,    arma::mat P_freq) {
  // so P_freq has to be a row matrix
  // so Q_freq has to be a column matrix
  int base_period = 0 ;
  // Since C++ is zero-indexed
  int M = Q.n_rows ;
  double M_dbl = Q.n_rows ;
  // double M_double = M ;
  // arma::mat x(M, M) ;
  arma::vec I_row(M);
  arma::vec I_col(M);
  arma::vec ret(M);
  I_row.ones() ;
  I_col.ones() ;
  ret.ones() ;
  inplace_trans(P_consol) ;
  inplace_trans(Q_consol) ;

  for (int col=0; col<M; col++) {
    double Q_ind_L = dot(P.row(base_period), Q.row(col)) / dot(P.row(base_period), Q.row(base_period)) ;
    double Q_ind_P = dot(P.row(col), Q.row(col)) / dot(P.row(col), Q.row(base_period)) ;
    I_row(col) = (Q_ind_L * Q_ind_P) ;
    // NOTE: important change fom the original fn: This is not sqrt since
    // that is done below with (1 / (2 * M_dbl))
  }
  // Above is the computation of the non-transitive indices. I wrote this earlier so it is inefficient
  double top_row = sum(log(I_row)) ;

  arma::vec revenue = log( sum(Q % P, 1) ) ;
  double sum_revenue = sum(revenue) ;
  // stop("test") ;
  for (int col=0; col<M; col++) {
    // stop("test2") ;
    // arma::mat test = P.row(col) * Q_consol ;
    // test.print("test:") ;

    //for (int row=0; row<M; row++) {
      // double Q_ind_L = dot(P.row(row), Q.row(col)) / dot(P.row(row), Q.row(row)) ;
      // double Q_ind_P = dot(P.row(col), Q.row(col)) / dot(P.row(col), Q.row(row)) ;
      // these two parts are a ok: double Q_ind_L = dot(P.row(row), Q.row(col)) ;
      //    double Q_ind_P = 1 / dot(P.row(col), Q.row(row)) ;
      // I_col(row) = log(Q_ind_L * Q_ind_P) ;
      // This one is ok too: double Q_ind_L = 1 / dot(P.row(row), Q.row(row)) ;
      //double Q_ind_L = dot(P.row(col), Q.row(col)) ;
      //I_col(row) = log(Q_ind_L ) ;
      //double Q_ind_L = dot(P.row(row), Q.row(col)) / dot(P.row(row), Q.row(row)) ;
      //double Q_ind_P = dot(P.row(col), Q.row(col)) / dot(P.row(col), Q.row(row)) ;
      //I_col(row) = (Q_ind_L * Q_ind_P) ;
    //}

    ret(col) = (1 / (2 * M_dbl)) * (
        top_row +
        M_dbl * arma::as_scalar(revenue(col))  -
        sum_revenue +
        arma::as_scalar( sum(log( Q.row(col) * P_consol) % P_freq )) -
        arma::as_scalar( sum(log( P.row(col) * Q_consol) % Q_freq ))
    ) ;
    // "For matrix M, return the sum of elements in each column (dim=0), or each row (dim=1) "
  }

  return(exp(ret)) ;
}














// [[Rcpp::export]]
arma::mat fisherIndfullmat (arma::mat Q, arma::mat P) {
  // Fisher index; direct approach; yes transitivity normalization

  int M = Q.n_rows ;
  // double M_double = M ;
  // arma::mat x(M, M) ;
  arma::vec I_row(M);
  arma::vec I_col(M);
  arma::mat ret(M,M);
  I_row.ones() ;
  I_col.ones() ;
  ret.ones() ;

  for (int col=0; col<M; col++) {

    for (int row=0; row<M; row++) {
      double Q_ind_L_1 = dot(P.row(col), Q.row(row)) / dot(P.row(col), Q.row(col)) ;
      double Q_ind_P_1 = dot(P.row(row), Q.row(row)) / dot(P.row(row), Q.row(col)) ;
      I_row(col) = sqrt(Q_ind_L_1 * Q_ind_P_1) ;

      double Q_ind_L_2 = dot(P.row(row), Q.row(col)) / dot(P.row(row), Q.row(row)) ;
      double Q_ind_P_2 = dot(P.row(col), Q.row(col)) / dot(P.row(col), Q.row(row)) ;
      I_col(row) = sqrt(Q_ind_L_2 * Q_ind_P_2) ;


    arma::vec interm_vec_prod = I_row % I_col ;

      // Rcout << "The value is " << sum(interm_vec_prod) << std::endl ;
      if ( ! arma::is_finite(interm_vec_prod) ) {
        // Note that arma::is_finite checks the whole vector for any non-finite values
        stop("NaNs produced in quantity index. Check the quantity and price matrix inputs. Quantity indices must be positive, so the product of quantities and prices must be positive in all cases.") ;
        // Throws error in the case of zeros in the computed quantity index, which
        // will only occur in the case of all zeros in quantity
        // Thanks to explanation here: http://gallery.rcpp.org/articles/intro-to-exceptions/
      }
      ret(row, col) = exp(mean(log(interm_vec_prod))) ;
    }
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
