#include <RcppArmadillo.h>
#include <Rcpp.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

double loss_theta(const vec& par,
                  const vec& R,
                  const mat& X,
                  const vec& weights,
                  const std::string& method_selection,
                  const std::string& h) { // TODO add weights

  Environment nonprobsvy_env = Environment::namespace_env("nonprobsvy");

  Rcpp::Function get_method = nonprobsvy_env["get_method"];
  List method = get_method(method_selection);

  Function inv_link = method["make_link_inv"];

  vec eta_pi = X * par;
  vec ps = as<vec>(inv_link(eta_pi));
  vec R_rand = 1 - R;

  // Preallocate temporary matrix
  mat temp;

  uvec loc_nons = find(R == 1);
  uvec loc_rand = find(R == 0);

  double N_nons = sum(1/ps);
  double N_rand = sum(weights);

  // Calculate the loss using the appropriate method
  double loss;
  if (h == "1") {
    temp = X.each_col() % (R % weights / ps / N_nons - R_rand % weights / N_rand);
    loss = accu(square(sum(temp, 0)));
  } else if (h == "2") {
    temp = X.each_col() % (R % weights / N_nons - R_rand % weights % ps / N_rand);
    loss = accu(square(sum(temp, 0)));
  } else {
    loss = 0;
  }

  return loss;
}


arma::vec u_theta(const arma::vec& par,
                  const arma::vec& R,
                  const arma::mat& X,
                  const arma::vec& weights,
                  const std::string& method_selection,
                  const std::string& h,
                  Nullable<int> N = R_NilValue,
                  Nullable<NumericVector> pop_totals = R_NilValue,
                  Nullable<double> pop_size = R_NilValue) { // TODO add weights

  Environment nonprobsvy_env = Environment::namespace_env("nonprobsvy");

  Rcpp::Function get_method = nonprobsvy_env["get_method"];
  List method = get_method(method_selection);

  Function inv_link = method["make_link_inv"];

  vec eta_pi = X * par;
  vec ps = as<vec>(inv_link(eta_pi));
  vec R_rand = 1 - R;
  double N_nons = sum(1/ps);

  vec eq;
  mat temp;
  int h_ = stoi(h);
  if (pop_totals.isNull()) {
    switch(h_) {
    case 1:
      temp = X.each_col() % (R/ps % weights - R_rand % weights);
      eq = sum(temp, 0).t() / N_nons;
      break;
    case 2:
      temp = X.each_col() % (R % weights - R_rand % weights % ps);
      eq = sum(temp, 0).t() / N_nons;
      break;
    }
  } else {
    vec total_pop = join_cols(vec{N_nons}, as<vec>(pop_totals));
    temp = X.each_col() % (R/ps % weights);
    eq = (sum(temp, 0).t() - total_pop) / N_nons;
  }

  return eq;
}


arma::mat u_theta_der(const arma::vec& par,
                      const arma::vec& R,
                      const arma::mat& X,
                      const arma::vec& weights,
                      const std::string& method_selection,
                      const std::string& h,
                      Nullable<int> N = R_NilValue,
                      Nullable<arma::vec> pop_totals = R_NilValue) { // TODO add weights

  Environment nonprobsvy_env = Environment::namespace_env("nonprobsvy");

  Rcpp::Function get_method = nonprobsvy_env["get_method"];
  List method = get_method(method_selection);

  Function inv_link = method["make_link_inv"];
  Function inv_link_der = method["make_link_inv_der"];

  //int p = X0.n_cols;
  arma::vec eta_pi = X * par;
  arma::vec ps = as<arma::vec>(inv_link(eta_pi));
  arma::vec R_rand = 1 - R;
  arma::vec psd;

  if (method_selection == "probit") {
    inv_link_der = method["make_link_inv_der"];
    psd = as<arma::vec>(inv_link_der(eta_pi));
  }

  int n = X.n_rows;
  int p = X.n_cols;
  arma::mat mxDer(p, p, arma::fill::zeros);
  double N_nons = sum(1/ps);

  arma::rowvec X_row;
  arma::mat temp;

  if (h == "1" || !pop_totals.isNull()) {
    if (method_selection == "logit") {
      for(int i = 0; i < n; i++) {
        X_row = X.row(i);
        temp = R(i) * weights(i) * (1-ps(i))/ps(i) * X_row.t();
        mxDer += temp * X_row;
      }
    } else if (method_selection == "cloglog") {
      for(int i = 0; i < n; i++) {
        X_row = X.row(i);
        temp = R(i) * weights(i) * (1-ps(i))/pow(ps(i), 2) * exp(eta_pi(i)) * X_row.t();
        mxDer += temp * X_row;
      }
    } else if (method_selection == "probit") {
      for(int i = 0; i < n; i++) {
        X_row = X.row(i);
        temp = R(i) * weights(i) * psd(i)/pow(ps(i), 2) * X_row.t();
        mxDer += temp * X_row;
      }
    }
  } else if (h == "2") {
    if (method_selection == "logit") {
      for(int i = 0; i < n; i++) {
        X_row = X.row(i);
        temp = R_rand(i) * weights(i) * ps(i)/(exp(eta_pi(i)) + 1) * X_row.t();
        mxDer += temp * X_row;
      }
    } else if (method_selection == "cloglog") {
      for(int i = 0; i < n; i++) {
        X_row = X.row(i);
        temp = R_rand(i) * weights(i) * (1-ps(i)) * exp(eta_pi(i)) * X_row.t();
        mxDer += temp * X_row;
      }
    } else if (method_selection == "probit") {
      for(int i = 0; i < n; i++) {
        X_row = X.row(i);
        temp = R_rand(i) * weights(i) * psd(i) * X_row.t();
        mxDer += temp * X_row;
      }
    }
    else {
      Rcpp::stop("Unknown method selection");
    }
  } else {
    Rcpp::stop("Unknown h selection");
  }

  return mxDer / N_nons;

}

arma::vec q_lambda_cpp(const arma::vec& par, double lambda, double a = 3.7) {

  arma::vec penaltyd(par.size(), arma::fill::zeros);

  double abs_par_i;
  for (arma::vec::iterator it = penaltyd.begin() + 1; it != penaltyd.end(); ++it) {
    abs_par_i = std::abs(*it);
    if (abs_par_i < lambda) {
      *it = lambda;
    } else {
      double tmp = ((a * lambda) - abs_par_i) / (a - 1);
      *it = tmp * (tmp > 0);
    }
  }

  // No penalty on the intercept
  penaltyd(0) = 0;

  return penaltyd;
}


// pass by reference where possible to avoid unnecessary copying
arma::vec fit_nonprobsvy_rcpp(const arma::mat& X,
                              const arma::vec& R,
                              const arma::vec& weights,
                              const std::string& method_selection,
                              const std::string& h,
                              double lambda,
                              int maxit,
                              double eps,
                              bool warn = false
) { // TODO add weights
  int p = X.n_cols;
  // preallocate large vectors
  arma::vec par0(p, arma::fill::zeros);
  arma::vec LAMBDA(p, arma::fill::zeros);
  arma::vec par(p, arma::fill::zeros);

  int it = 0;
  for (int jj = 1; jj <= maxit; jj++) {
    it++;
    if (warn && it == maxit) {
      Rcpp::warning("Convergence not obtained in %d iterations of fitting algorithm for variables selection", maxit);
      break;
    }
    // avoid unnecessary computations by saving results
    arma::vec u_theta0v = u_theta(par0, R, X, weights, method_selection, h);
    arma::mat u_theta0_derv = u_theta_der(par0, R, X, weights, method_selection, h);
    arma::vec q_lambda_output = q_lambda_cpp(par0, lambda);

    LAMBDA = arma::abs(q_lambda_output) / (eps + arma::abs(par0));
    // use efficient Armadillo functions
    par = par0 + solve(arma::reshape(u_theta0_derv, p, p) + arma::diagmat(LAMBDA), u_theta0v - arma::diagmat(LAMBDA) * par0);

    if (arma::sum(arma::abs(par - par0)) < eps) break;
    if (arma::sum(arma::abs(par - par0)) > 1000) break;

    par0 = par;
  }

  par(arma::find(arma::abs(par) < 0.001)).zeros();
  return par;
}


// [[Rcpp::export]]
Rcpp::List cv_nonprobsvy_rcpp(const arma::mat& X,
                              const arma::vec& R,
                              const arma::vec& weights_X,
                              const std::string& method_selection,
                              const std::string& h,
                              int maxit,
                              double eps,
                              double lambda_min,
                              int nlambda,
                              int nfolds,
                              double lambda = -1) { // TODO add weights

  Environment nonprobsvy_env = Environment::namespace_env("nonprobsvy");
  Rcpp::Function setup_lambda_cpp = nonprobsvy_env["setup_lambda"];

  Rcpp::Function logit = nonprobsvy_env["logit"];
  Rcpp::Function cloglog = nonprobsvy_env["cloglog"];
  Rcpp::Function probit = nonprobsvy_env["probit"];

  arma::vec weights;
  arma::vec loss_theta_av(nlambda, arma::fill::zeros);
  arma::vec R_ = R;
  arma::mat X_ = X;

  if(lambda == -1) {
    arma::uvec loc_nons = find(R == 1);
    arma::uvec loc_rand = find(R == 0);
    arma::mat X_nons = join_rows(X.rows(loc_nons), weights_X(loc_nons), R(loc_nons));
    arma::mat X_rand = join_rows(X.rows(loc_rand), weights_X(loc_rand), R(loc_rand));

    SEXP lambdas = setup_lambda_cpp(X, R, weights_X, method_selection, lambda_min, nlambda);
    arma::vec lambdas1(Rcpp::as<arma::vec>(lambdas));

    arma::uvec shuffle_nons = arma::shuffle(arma::linspace<arma::uvec>(0, X_nons.n_rows-1, X_nons.n_rows));
    arma::uvec shuffle_rand = arma::shuffle(arma::linspace<arma::uvec>(0, X_rand.n_rows-1, X_rand.n_rows));

    arma::uvec folds_nons = arma::randi<arma::uvec>(X_nons.n_rows, arma::distr_param(0, nfolds-1));
    arma::uvec folds_rand = arma::randi<arma::uvec>(X_rand.n_rows, arma::distr_param(0, nfolds-1));

    arma::uvec sample_nons = arma::shuffle(arma::linspace<arma::uvec>(0, nfolds-1, nfolds));
    arma::uvec sample_rand = arma::shuffle(arma::linspace<arma::uvec>(0, nfolds-1, nfolds));

    for(int i = 0; i < nlambda; i++) {
      lambda = lambdas1(i);
      arma::vec loss_theta_vec(nfolds, arma::fill::zeros);

      for(int j = 0; j < nfolds; j++) {
        arma::uvec idx_nons = find(folds_nons != sample_nons(j));
        arma::mat X_nons_train = X_nons.rows(idx_nons);
        arma::mat X_nons_test = X_nons.rows(find(folds_nons == sample_nons(j)));

        arma::uvec idx_rand = find(folds_rand != sample_rand(j));
        arma::mat X_rand_train = X_rand.rows(idx_rand);
        arma::mat X_rand_test = X_rand.rows(find(folds_rand == sample_rand(j)));

        arma::mat X_train = arma::join_cols(X_rand_train, X_nons_train);
        arma::mat X_test = arma::join_cols(X_rand_test, X_nons_test);
        int ncols = X_test.n_cols;

        arma::uvec idxx = arma::regspace<arma::uvec>(0, ncols - 3);

        arma::vec theta_est = fit_nonprobsvy_rcpp(X_train.cols(idxx),
                                                  X_train.col(ncols - 1),
                                                  X_train.col(ncols - 2),
                                                  method_selection,
                                                  h,
                                                  lambda,
                                                  maxit,
                                                  eps);

        if (arma::any(theta_est == 0)) {
          idxx.shed_rows(arma::find(theta_est == 0));
        }
        arma::mat X_testloss = X_test.cols(idxx);
        arma::vec R_testloss = X_test.col(ncols - 1);
        arma::vec weights_testloss = X_test.col(ncols - 2);
        arma::vec par = theta_est(idxx);

        double loss = loss_theta(par, R_testloss, X_testloss, weights_testloss, method_selection, h);
        loss_theta_vec(j) = loss;
      }
      loss_theta_av(i) = mean(loss_theta_vec);
    }
    lambda = lambdas1(loss_theta_av.index_min());
  }

  arma::vec theta = fit_nonprobsvy_rcpp(X_,
                                        R_,
                                        weights_X,
                                        method_selection,
                                        h,
                                        lambda,
                                        maxit,
                                        eps);

  arma::uvec theta_selected = find(theta != 0);

  return List::create(_["theta_est"] = theta,
                      _["theta_selected"] = theta_selected,
                      _["min"] = loss_theta_av.min(),
                      _["lambda"] = lambda);
}
