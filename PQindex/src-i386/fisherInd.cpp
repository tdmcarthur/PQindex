#include <RcppArmadillo.h>
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



// arma::mat matmult_cpp(SEXP Xr, const arma::mat Y) {
//    if (Rf_isS4(Xr)) {
//        if(Rf_inherits(Xr, "dgCMatrix")) {
//            return matmult_sp(as<arma::sp_mat>(Xr), Y) ;
//        } ;



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


   // inplace_trans(P_freq) ;

    ret(col) = (1 / (2 * M_dbl)) * (
        top_row +
        M_dbl * arma::as_scalar(revenue(col))  -
        sum_revenue  +
        arma::as_scalar( log( Q.row(col) * P_consol) * P_freq ) -
        arma::as_scalar( log( P.row(col) * Q_consol) * Q_freq )
    ) ;
    // "For matrix M, return the sum of elements in each column (dim=0), or each row (dim=1) "
  }

  return(exp(ret)) ;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat testfn (arma::mat Q_consol, arma::uvec Q_ind) {
 // vec v = vec(Q_ind.n_elem) ;
//  v.ones() ;
  mat Q_consol_redef = Q_consol.rows(Q_ind) ;
  return(Q_consol_redef) ;
}

// arma::mat fisherIndfastest (arma::sp_mat  Q_consol,  arma::sp_mat P_consol,
 //                           arma::sp_mat  Q_freq,    arma::sp_mat P_freq,
   //                         arma::uvec Q_ind,     arma::uvec P_ind) {

// arma::mat fisherIndfastest (arma::mat  Q_consol,  arma::mat P_consol,
   //                         arma::mat  Q_freq,    arma::mat P_freq,
     //                       arma::uvec Q_ind,     arma::uvec P_ind) {

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat fisherIndfastest (arma::sp_mat  Q_consol,  arma::sp_mat  P_consol,
                            arma::mat     Q_freq,    arma::mat     P_freq,
                            arma::uvec    Q_ind,     arma::uvec    P_ind) {


  uword base_period = 0 ;
  double M_dbl = Q_ind.n_elem ;
 // arma::mat x(5, 5) ;
//  mat Q_consol_dense = conv_to<mat>::from(Q_consol);
//  mat P_consol_dense = conv_to<mat>::from(P_consol);
  arma::wall_clock timer;
  timer.tic();
  arma::vec revenue = log( sum( conv_to<mat>::from(Q_consol).rows(Q_ind) %
                                conv_to<mat>::from(P_consol).rows(P_ind), 1) ) ;
  cout << "time taken = " << timer.toc()  << endl;
  double sum_revenue = sum(revenue) ;

//  inplace_trans(P_consol) ;
  arma::vec Q_x_P_consol = log( conv_to<mat>::from(Q_consol * P_consol.t() )) * P_freq;
  //arma::vec Q_x_P_consol = log( Q_consol * P_consol.t()) * P_freq ;
  // arma::vec Q_x_P_consol = Q_x_P_consol_inter(Q_ind) ;

//  inplace_trans(Q_consol) ;
 // inplace_trans(P_consol) ;
  // arma::as_scalar( sum(log( Q.row(col) * P_consol) * P_freq ))
  arma::vec P_x_Q_consol = log( conv_to<mat>::from( P_consol * Q_consol.t() )) * Q_freq.t() ;
  // arma::vec P_x_Q_consol = P_x_Q_consol_inter(P_ind) ;

   double top_row  = (-1) * M_dbl * arma::as_scalar(revenue(base_period))  +
        sum_revenue -
        Q_x_P_consol(base_period) +
        P_x_Q_consol(base_period) ;

  arma::vec ret = (1 / (2 * M_dbl)) * (
        top_row +
        M_dbl * revenue  -
        sum_revenue  +
        Q_x_P_consol(Q_ind) -
        P_x_Q_consol(P_ind)
    ) ;
    // "For matrix M, return the sum of elements in each column (dim=0), or each row (dim=1) "



  return(exp(ret)) ;
}



// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::sp_mat matmult_cpp(const arma::sp_mat X, const arma::sp_mat Y) {
// arma::sp_mat ret =  ;
 return(X * Y) ;
}



// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat fisherIndfaster (arma::mat  Q_consol,  arma::mat P_consol,
                           arma::mat  Q_freq,    arma::mat P_freq,
                           arma::uvec Q_ind,     arma::uvec P_ind) {
  // so P_freq has to be a row matrix
  // so Q_freq has to be a column matrix
  uword base_period = 0 ;
  // Since C++ is zero-indexed
  // int M = Q_ind.n_elem ;
  double M_dbl = Q_ind.n_elem ;
  // double M_double = M ;
  // arma::mat x(M, M) ;
  //arma::vec I_row(M);
  //arma::vec I_col(M);
 // arma::vec ret(M);
  // I_row.ones() ;
  // I_col.ones() ;
  // ret.ones() ;


 // for (int col=0; col<M; col++) {
 //   double Q_ind_L = dot(P.row(base_period), Q.row(col)) / dot(P.row(base_period), Q.row(base_period)) ;
 //   double Q_ind_P = dot(P.row(col), Q.row(col)) / dot(P.row(col), Q.row(base_period)) ;
 //   I_row(col) = (Q_ind_L * Q_ind_P) ;
    // NOTE: important change fom the original fn: This is not sqrt since
    // that is done below with (1 / (2 * M_dbl))
  // }
  // Above is the computation of the non-transitive indices. I wrote this earlier so it is inefficient
 // double top_row = sum(log(I_row)) ;


  // mat Q_consol_redef = Q_consol.each_row(Q_ind) + v;
  arma::vec revenue = log( sum(Q_consol.rows(Q_ind) % P_consol.rows(P_ind), 1) ) ;
 // arma::vec revenue = log( sum(Q_consol.rows(Q_ind) % v, 1) ) ;
  double sum_revenue = sum(revenue) ;

  inplace_trans(P_consol) ;
  // arma::as_scalar( sum(log( Q.row(col) * P_consol) * P_freq ))
  arma::vec Q_x_P_consol = log( Q_consol * P_consol) * P_freq ;
  //arma::vec Q_x_P_consol = Q_x_P_consol_inter(Q_ind) ;
 // inplace_trans(P_consol) ;
//  inplace_trans(Q_consol) ;
//   vec test = Q_freq.t();
 // stop("test") ;

  inplace_trans(Q_consol) ;
  inplace_trans(P_consol) ;
  // arma::as_scalar( sum(log( Q.row(col) * P_consol) * P_freq ))
  arma::vec P_x_Q_consol = log( P_consol * Q_consol) * Q_freq.t() ;
//  arma::vec P_x_Q_consol = P_x_Q_consol_inter(P_ind) ;

  // arma::mat P_x_Q_consol = Q_freq * log(  Q_consol * P_consol.cols(P_ind))   ;
  // Q_freq *
   // arma::as_scalar( sum(log( P.row(col) * Q_consol) % Q_freq ))
  // stop("test") ;
  //inplace_trans(P_x_Q_consol) ;


   double top_row  = (-1) * M_dbl * arma::as_scalar(revenue(base_period))  +
        sum_revenue -
        Q_x_P_consol(base_period) +
        P_x_Q_consol(base_period) ;

    // "For matrix M, return the sum of elements in each column (dim=0), or each row (dim=1) "
  // }


  // stop("test") ;
//  for (uword col=0; col<M; col++) {
//    ret(col) = (1 / (2 * M_dbl)) * (
//        top_row +
//        M_dbl * arma::as_scalar(revenue(col))  -
//        sum_revenue +
//        Q_x_P_consol(Q_ind(col)) -
//        P_x_Q_consol(P_ind(col))
//    ) ;
    // "For matrix M, return the sum of elements in each column (dim=0), or each row (dim=1) "
//  }


  arma::vec ret = (1 / (2 * M_dbl)) * (
        top_row +
        M_dbl * revenue  -
        sum_revenue  +
        Q_x_P_consol(Q_ind) -
        P_x_Q_consol(P_ind)
    ) ;
    // "For matrix M, return the sum of elements in each column (dim=0), or each row (dim=1) "



  return(exp(ret)) ;
}











// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat fisherIndfasterold (arma::mat Q_consol,  arma::mat P_consol,
                           arma::mat Q_freq,    arma::mat P_freq,
                           arma::uvec Q_ind,     arma::uvec P_ind) {
  // so P_freq has to be a row matrix
  // so Q_freq has to be a column matrix
  uword base_period = 0 ;
  // Since C++ is zero-indexed
  int M = Q_ind.n_elem ;
  double M_dbl = Q_ind.n_elem ;
  // double M_double = M ;
  // arma::mat x(M, M) ;
  arma::vec I_row(M);
  arma::vec I_col(M);
  arma::vec ret(M);
  I_row.ones() ;
  I_col.ones() ;
  ret.ones() ;


 // for (int col=0; col<M; col++) {
 //   double Q_ind_L = dot(P.row(base_period), Q.row(col)) / dot(P.row(base_period), Q.row(base_period)) ;
 //   double Q_ind_P = dot(P.row(col), Q.row(col)) / dot(P.row(col), Q.row(base_period)) ;
 //   I_row(col) = (Q_ind_L * Q_ind_P) ;
    // NOTE: important change fom the original fn: This is not sqrt since
    // that is done below with (1 / (2 * M_dbl))
  // }
  // Above is the computation of the non-transitive indices. I wrote this earlier so it is inefficient
 // double top_row = sum(log(I_row)) ;


  // mat Q_consol_redef = Q_consol.each_row(Q_ind) + v;
  arma::vec revenue = log( sum(Q_consol.rows(Q_ind) % P_consol.rows(P_ind), 1) ) ;
 // arma::vec revenue = log( sum(Q_consol.rows(Q_ind) % v, 1) ) ;
  double sum_revenue = sum(revenue) ;

  inplace_trans(P_consol) ;

  arma::vec Q_x_P_consol = log( Q_consol.rows(Q_ind) * P_consol) * P_freq ;
  inplace_trans(P_consol) ;
  inplace_trans(Q_consol) ;
  arma::vec P_x_Q_consol = log( P_consol.rows(P_ind) * Q_consol) * Q_freq ;
   // arma::as_scalar( sum(log( P.row(col) * Q_consol) % Q_freq ))

   double top_row  = (-1) * M_dbl * arma::as_scalar(revenue(base_period))  +
        sum_revenue -
        Q_x_P_consol(Q_ind(base_period)) +
        P_x_Q_consol(P_ind(base_period)) ;

    // "For matrix M, return the sum of elements in each column (dim=0), or each row (dim=1) "
  // }


  // stop("test") ;
  for (uword col=0; col<M; col++) {
    ret(col) = (1 / (2 * M_dbl)) * (
        top_row +
        M_dbl * arma::as_scalar(revenue(col))  -
        sum_revenue +
        Q_x_P_consol(Q_ind(col)) -
        P_x_Q_consol(P_ind(col))
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
if ( F) {
library(data.table)
library(Matrix)

consol.matrix <- function(x) {
  if (!is.data.table(x)) x <- as.data.table(x)
  x.ret <- x[, .(.N), by = names(x)]
  N.ret <- matrix(x.ret$N, ncol = 1)
  x.ret[, N := NULL]
  list(mat = as.matrix(x.ret), freq = N.ret)
}



set.seed(100)
# n.col <- 100; n.row = 40000
# With these params, fastest index fn get 77 secs. Faster index fn gets 320 secs (4 times faster):
# n.row.fact <- 20000 ; real.rows.factor = 20 ; n.col <- 400;
# With the below, I have fastest 0.13; faster 0.185; naive 18.4 secs :
# n.row.fact <- 1000 ; real.rows.factor = 5 ; n.col <- 100;
# With below, I get fastest 0.013; faster 0.014; naive 112.533:
n.row.fact <- 100 ; real.rows.factor = 100 ; n.col <- 100;
# n.row.fact <- 10 ; real.rows.factor = 2 ; n.col <- 4;
n.row = real.rows.factor; n.row = n.row * n.row.fact
n.real.rows = n.row / real.rows.factor
P.mat <- matrix(runif(n.real.rows*n.col), ncol = n.col, nrow = n.row, byrow = TRUE )
P.mat <- rbind(P.mat[-1, ], P.mat[1, ])
#P.mat <- rbind(P.mat[1, ], P.mat[1, ], P.mat[2, ], P.mat[2, ])
#P.mat <- matrix(runif(n.col*n.row), nrow = n.row )
# Q.mat <- matrix(runif(n.col*n.row), ncol = n.col)
Q.mat <- matrix(runif(n.real.rows*n.col), ncol = n.col, nrow = n.row, byrow = TRUE )
Q.mat[, 10:ncol(Q.mat)] <- 0
# Making the matrix sparse


Q.mat.consol <- consol.matrix(Q.mat)
P.mat.consol <- consol.matrix(P.mat)

if (F) {
  print( system.time( fisherInd.ret <- fisherInd(Q.mat, P.mat, 1) ) )
}

if (F) {
print(system.time(
fisherIndfast.ret <-
  fisherIndfast(Q = Q.mat, P = P.mat,
                Q_consol = Q.mat.consol$mat,
                P_consol = P.mat.consol$mat,
                Q_freq = Q.mat.consol$freq,
                P_freq = P.mat.consol$freq ) # t(P.mat.consol$freq ))
))
}


if (T) {
print(system.time(
fisherIndfaster.ret <- fisherIndfaster(Q_consol = Q.mat.consol$mat,
                P_consol = P.mat.consol$mat,
                Q_freq = t(Q.mat.consol$freq),
                #Q_freq = Q.mat.consol$freq,
                P_freq = P.mat.consol$freq,
                Q_ind = rep((1:n.real.rows) - 1, real.rows.factor),
                P_ind = rep((1:n.real.rows) - 1, real.rows.factor))
                # P_ind = c(rep((1:n.real.rows) - 1, real.rows.factor)[-1], 0))
))
}



if (T) {
print(system.time(
fisherIndfastest.ret <- fisherIndfastest(
              # Q_consol = Q.mat.consol$mat,
              # P_consol = P.mat.consol$mat,
                Q_consol = Matrix(Q.mat.consol$mat, sparse = TRUE),
                P_consol = Matrix(P.mat.consol$mat, sparse = TRUE),
                Q_freq = t(Q.mat.consol$freq),
                #Q_freq = Q.mat.consol$freq,
                P_freq = P.mat.consol$freq,
                Q_ind = rep((1:n.real.rows) - 1, real.rows.factor),
                P_ind = rep((1:n.real.rows) - 1, real.rows.factor))
                # P_ind = c(rep((1:n.real.rows) - 1, real.rows.factor)[-1], 0))
))
}


try( print(summary(fisherInd.ret - fisherIndfast.ret)) )

try( print(summary(fisherIndfaster.ret - fisherIndfast.ret)) )

try( print(summary(fisherIndfastest.ret - fisherIndfaster.ret)) )

try( print(summary(fisherIndfastest.ret - fisherInd.ret)) )


}
*/
