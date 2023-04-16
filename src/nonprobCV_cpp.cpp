#include <RcppArmadillo.h>
#include <Rcpp.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;


// [[Rcpp::export]]
arma::vec fit_nonprobsvy_rcpp(arma::mat X,
                             arma::vec R,
                             arma::vec weights,
                             std::string method_selection,
                             std::string h,
                             double lambda,
                             int maxit,
                             double eps,
                             bool warn = false
) {
  int p = X.n_cols;
  arma::vec init_theta(p, arma::fill::zeros);
  arma::vec par0 = init_theta;
  arma::vec LAMBDA(p, arma::fill::zeros);
  arma::vec par;

  Environment nonprobsvy_env = Environment::namespace_env("nonprobsvy");

  Rcpp::Function u_theta_cpp = nonprobsvy_env["u_theta"];
  Rcpp::Function u_theta_der_cpp = nonprobsvy_env["u_theta_der"];
  Rcpp::Function q_lambda_cpp = nonprobsvy_env["q_lambda"];

  //Rcpp::Function u_theta_cpp("u_theta");
  //Rcpp::Function u_theta_der_cpp("u_theta_der");
  //Rcpp::Function q_lambda_cpp("q_lambda");
  //Rcpp::Function abs("abs");
  //Rcpp::Function sum("sum");

  int it = 0;
  for (int jj = 1; jj <= maxit; jj++) {
    it++;
    if (warn && it == maxit) {
      Rcpp::warning("Convergence not obtained in %d iterations of fitting algorithm for variables selection", maxit);
      break;
    }

    Function u_theta_output = u_theta_cpp(R,
                                          X,
                                          weights,
                                          method_selection,
                                          h);
    SEXP u_theta_outputt = u_theta_output(par0);
    arma::vec u_theta0v(Rcpp::as<arma::vec>(u_theta_outputt));

    Function u_theta_der_output = u_theta_der_cpp(R,
                                                  X,
                                                  weights,
                                                  method_selection,
                                                  h);
    SEXP u_theta_der_outputt = u_theta_der_output(par0);
    arma::mat u_theta0_derv(Rcpp::as<arma::mat>(u_theta_der_outputt));

    SEXP q_lambda_output = q_lambda_cpp(par0, lambda);
    LAMBDA = arma::abs(Rcpp::as<arma::vec>(q_lambda_output)) / (eps + arma::abs(par0)); // use abs instead of f for consistency

    // Fix problem with theta, theta_der
    //Rcout << "The value of x is: " << arma::reshape(u_theta0v, p, p) << "\n";
    par = par0 + arma::inv(arma::reshape(u_theta0_derv, p, p) + arma::diagmat(LAMBDA)) * (u_theta0v - arma::diagmat(LAMBDA) * par0);

    if (arma::sum(arma::abs(par - par0)) < eps) break;
    if (arma::sum(arma::abs(par - par0)) > 1000) break;

    par0 = par;
  }

   par(arma::find(arma::abs(par) < 0.001)).zeros(); // use abs instead of f for consistency
   arma::vec theta_est = par;


  return theta_est;
}

