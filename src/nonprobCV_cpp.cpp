#include <RcppArmadillo.h>
#include <Rcpp.h>
//#include <Eigen/Dense>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
using namespace std;
//using namespace Eigen;

inline double loss_theta(const vec& par,
                         const vec& R,
                         const mat& X,
                         const vec& weights,
                         const std::string& method_selection,
                         const int& gee_h_fun,
                         const uvec& idx,
                         Nullable<arma::vec> pop_totals) { // TODO add weights

  Environment nonprobsvy_env = Environment::namespace_env("nonprobsvy");
  Rcpp::Function method_ps = nonprobsvy_env["method_ps"];
  List method = method_ps(method_selection);

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
  if (pop_totals.isNull()) {
    if (gee_h_fun == 1) {
      temp = X.each_col() % (R % weights / ps / N_nons - R_rand % weights / N_rand);
      loss = accu(square(sum(temp, 0)));
    } else if (gee_h_fun == 2) {
      temp = X.each_col() % (R % weights / N_nons - R_rand % weights % ps / N_rand);
      loss = accu(square(sum(temp, 0)));
    } else {
      loss = 0;
    }
  } else {
   // vec total_pop = join_cols(vec{N_nons}, as<vec>(pop_totals));
    temp = X.each_col() % (R % weights / ps / N_nons);
    //vec colSums_result = sum(temp, 0);

    // Calculate (colSums(weights * X/pi) - pop_totals)^2
    vec diff_squared = square(sum(temp, 0).t() - as<vec>(pop_totals)(idx) / N_rand);
    //loss = accu(square(sum(temp, 0)));

    // Calculate the sum of squared differences
    loss = accu(diff_squared);
  }
  return loss;
}


inline arma::vec u_theta(const arma::vec& par,
                         const arma::vec& R,
                         const arma::mat& X,
                         const arma::vec& weights,
                         const std::string& method_selection,
                         const int& h,
                         Nullable<arma::vec> pop_totals,
                         Nullable<double> pop_size = R_NilValue,
                         Nullable<int> N = R_NilValue) { // TODO add weights

  Environment nonprobsvy_env = Environment::namespace_env("nonprobsvy");
  Rcpp::Function method_ps = nonprobsvy_env["method_ps"];
  List method = method_ps(method_selection);

  Function inv_link = method["make_link_inv"];

  vec eta_pi = X * par;
  vec ps = as<vec>(inv_link(eta_pi));
  vec R_rand = 1 - R;
  double N_nons = sum(1/ps);

  vec eq;
  mat temp;
  if (pop_totals.isNull()) {
    switch(h) {
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
    //vec total_pop = join_cols(vec{N_nons}, as<vec>(pop_totals));
    temp = X.each_col() % (R/ps % weights);
    eq = (sum(temp, 0).t() - as<vec>(pop_totals)) / N_nons;
  }
  return eq;
}

arma::mat u_theta_der(const arma::vec& par,
                      const arma::vec& R,
                      const arma::mat& X,
                      const arma::vec& weights,
                      const std::string& method_selection,
                      const int& gee_h_fun,
                      Nullable<arma::vec> pop_totals,
                      Nullable<int> N = R_NilValue) { // TODO add weights

  Environment nonprobsvy_env = Environment::namespace_env("nonprobsvy");
  Rcpp::Function method_ps = nonprobsvy_env["method_ps"];
  List method = method_ps(method_selection);

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

  if (gee_h_fun == 1 || !pop_totals.isNull()) {
    if (method_selection == "logit") {
      for(int i = 0; i < n; i++) {
        X_row = X.row(i);
        temp = R(i) * weights(i) * (1-ps(i))/ps(i) * X_row.t();
        mxDer += temp * X_row;
      }
      //mxDer = X.t() * X;
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
  } else if (gee_h_fun == 2) {
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

arma::vec q_lambda_cpp(const arma::vec& par,
                       double lambda,
                       const std::string& penalty,
                       double a) { // TODO add a to control

  arma::vec penaltyd(par.size(), arma::fill::zeros);

  if (penalty == "SCAD") {
    double abs_par_i;
    for (arma::vec::iterator it = penaltyd.begin() + 1; it != penaltyd.end(); ++it) {
      abs_par_i = std::abs(*it);
      if (abs_par_i <= lambda) {
        *it = lambda;
      } else {
        double tmp = ((a * lambda) - abs_par_i) / (a - 1);
        *it = tmp * (tmp > 0);
      }
    }
  } else if (penalty == "lasso") {
    penaltyd = lambda * arma::sign(par);
  } else if (penalty == "MCP") {
    for (std::size_t i = 0; i < par.size(); i++) { // int i = 0; i < par.size(); i++
      if (std::abs(par[i]) <= a*lambda) {
        if (par[i] < 0) {
          penaltyd[i] = - (lambda - std::abs(par[i]) / a);
        }
        else if (par[i] > 0) {
          penaltyd[i] = (lambda - std::abs(par[i]) / a);
        } else {
          penaltyd[i] = 0;
        }
      } else {
        penaltyd[i] = 0;
      }
    }
  }
  return penaltyd;
}

// HybridNonLinearSolver (?)
// pass by reference where possible to avoid unnecessary copying
arma::vec fit_nonprobsvy_rcpp(const arma::mat& X,
                              const arma::vec& R,
                              const arma::vec& weights,
                              const std::string& method_selection,
                              const int& gee_h_fun,
                              double lambda,
                              int maxit,
                              double eps,
                              const std::string& penalty,
                              double a,
                              Nullable<arma::vec> pop_totals,
                              bool warn = false
) {
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

    arma::vec u_theta0v = u_theta(par0, R, X, weights, method_selection, gee_h_fun, pop_totals);
    arma::mat u_theta0_derv = u_theta_der(par0, R, X, weights, method_selection, gee_h_fun, pop_totals);
    arma::vec q_lambda_output = q_lambda_cpp(par0, lambda, penalty, a);

    LAMBDA = arma::abs(q_lambda_output) / (eps + arma::abs(par0)); // TODO  q_lambda_output instead of arma::abs(q_lambda_output)
    //LAMBDA = arma::abs(q_lambda_output);
    // LAMBDA = q_lambda_output;
    par = par0 + inv(arma::reshape(u_theta0_derv, p, p) + arma::diagmat(LAMBDA)) * (u_theta0v - arma::diagmat(LAMBDA) * par0);

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
                              const int& gee_h_fun,
                              int maxit,
                              double eps,
                              double lambda_min,
                              int nlambda,
                              int nfolds,
                              const std::string& penalty,
                              double a,
                              Nullable<arma::vec> pop_totals,
                              bool verbose,
                              double lambda = -1) { // TODO add weights

  Environment nonprobsvy_env = Environment::namespace_env("nonprobsvy");
  Rcpp::Function setup_lambda_cpp = nonprobsvy_env["setup_lambda"];

  //Rcpp::Function logit = nonprobsvy_env["logit_model_nonprobsvy"];
  //Rcpp::Function cloglog = nonprobsvy_env["cloglog_model_nonprobsvy"];
  //Rcpp::Function probit = nonprobsvy_env["probit_model_nonprobsvy"];

  arma::vec weights;
  arma::vec loss_theta_av(nlambda);
  const arma::vec& R_ = R;
  const arma::mat& X_ = X;

  if(lambda == -1) {
    arma::uvec loc_nons = find(R == 1);
    arma::uvec loc_rand = find(R == 0);
    const arma::mat& X_nons = join_rows(X.rows(loc_nons), weights_X(loc_nons), R(loc_nons));
    const arma::mat& X_rand = join_rows(X.rows(loc_rand), weights_X(loc_rand), R(loc_rand));

    SEXP lambdas = setup_lambda_cpp(X, R, weights_X, method_selection, lambda_min, nlambda, pop_totals);
    arma::vec lambdas1(Rcpp::as<arma::vec>(lambdas));

    //arma::uvec shuffle_nons = arma::shuffle(arma::linspace<arma::uvec>(0, X_nons.n_rows-1, X_nons.n_rows));
    //arma::uvec shuffle_rand = arma::shuffle(arma::linspace<arma::uvec>(0, X_rand.n_rows-1, X_rand.n_rows));

    arma::uvec folds_nons = arma::randi<arma::uvec>(X_nons.n_rows, arma::distr_param(0, nfolds-1));
    arma::uvec folds_rand = arma::randi<arma::uvec>(X_rand.n_rows, arma::distr_param(0, nfolds-1));

    arma::uvec sample_nons = arma::shuffle(arma::linspace<arma::uvec>(0, nfolds-1, nfolds));
    arma::uvec sample_rand = arma::shuffle(arma::linspace<arma::uvec>(0, nfolds-1, nfolds));

    arma::field<arma::vec> loss_theta_fld(nfolds, nlambda);
    //#pragma omp parallel for
    for(int j = 0; j < nfolds; j++) {
      if (verbose) {
        wcout << "Starting CV fold #" << j+1 << endl;
      }
      arma::uvec idx_nons = find(folds_nons != sample_nons(j));
      const arma::mat& X_nons_train = X_nons.rows(idx_nons);
      const arma::mat& X_nons_test = X_nons.rows(find(folds_nons == sample_nons(j)));

      arma::uvec idx_rand = find(folds_rand != sample_rand(j));
      const arma::mat& X_rand_train = X_rand.rows(idx_rand);
      const arma::mat& X_rand_test = X_rand.rows(find(folds_rand == sample_rand(j)));

      // Randomize the columns (features) in the training data
      // arma::uvec col_indices = arma::shuffle(arma::regspace<arma::uvec>(0, X_nons_train.n_cols - 1));
      // arma::mat X_nons_train_randomized = X_nons_train.cols(col_indices);
      // arma::mat X_rand_train_randomized = X_rand_train.cols(col_indices);
      // arma::mat X_nons_test_randomized = X_nons_test.cols(col_indices);
      // arma::mat X_rand_test_randomized = X_rand_test.cols(col_indices);

      const arma::mat& X_train = arma::join_cols(X_rand_train, X_nons_train);
      const arma::mat& X_test = arma::join_cols(X_rand_test, X_nons_test);
      int ncols = X_test.n_cols;
      arma::uvec idxx = arma::regspace<arma::uvec>(0, ncols - 3);

      //#pragma omp parallel for
      for(int i = 0; i < nlambda; i++) {
        // lambda = lambdas1(i);
        //arma::vec loss_theta_vec(nfolds, arma::fill::zeros);
        arma::vec theta_est = fit_nonprobsvy_rcpp(X_train.cols(idxx),
                                                  X_train.col(ncols - 1),
                                                  X_train.col(ncols - 2),
                                                  method_selection,
                                                  gee_h_fun,
                                                  lambdas1(i),
                                                  maxit,
                                                  eps,
                                                  penalty,
                                                  a,
                                                  pop_totals);
        // cout << theta_est << "\n";

        const arma::mat& X_testloss = X_test.cols(arma::find(theta_est != 0));
        const arma::vec& R_testloss = X_test.col(ncols - 1);
        const arma::vec& weights_testloss = X_test.col(ncols - 2);
        const arma::vec& par = theta_est(arma::find(theta_est != 0));

        double loss = loss_theta(par,
                                 R_testloss,
                                 X_testloss,
                                 weights_testloss,
                                 method_selection,
                                 gee_h_fun,
                                 arma::find(theta_est != 0),
                                 pop_totals);
        loss_theta_fld(j, i) = loss;
      }
      //loss_theta_av(i) = mean(loss_theta_vec);
    }

    arma::vec loss_theta_vec(nfolds);
    // Vector to store means, one for each field
    for (int i = 0; i < nlambda; i++) {
      // arma::vec loss_theta_vec(nfolds);
      for (int j = 0; j < nfolds; j++) {
        loss_theta_vec(j) = loss_theta_fld(j, i)(0);
      }
      loss_theta_av(i) = mean(loss_theta_vec);
    }
    lambda = lambdas1(loss_theta_av.index_min());
  }
  arma::vec theta = fit_nonprobsvy_rcpp(X_,
                                        R_,
                                        weights_X,
                                        method_selection,
                                        gee_h_fun,
                                        lambda,
                                        maxit,
                                        eps,
                                        penalty,
                                        a,
                                        pop_totals);

  arma::uvec theta_selected = find(theta != 0);

  return List::create(_["theta_est"] = theta,
                      _["theta_selected"] = theta_selected,
                      _["min"] = loss_theta_av.min(),
                      _["lambda"] = lambda,
                      _["cv_error"] = loss_theta_av
                      );
}
