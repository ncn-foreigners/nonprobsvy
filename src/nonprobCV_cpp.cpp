#include <RcppArmadillo.h>
#include <Rcpp.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

double loss_theta(arma::vec par,
                  arma::vec R,
                  arma::mat X,
                  arma::vec weights,
                  std::string method_selection,
                  std::string h) {

  Environment nonprobsvy_env = Environment::namespace_env("nonprobsvy");

  Rcpp::Function get_method = nonprobsvy_env["get_method"];
  List method = get_method(method_selection);

  Function inv_link = method["make_link_inv"];

  arma::mat X0 = X;
  arma::vec theta = par;
  arma::vec R_rand = 1 - R;

  // Call the R function to get the inverse link function
  arma::vec eta_pi = X0 * theta;
  arma::vec ps = as<arma::vec>(inv_link(eta_pi));

  // Subset indices for random and non-random samples
  arma::uvec loc_nons = find(R == 1);
  arma::uvec loc_rand = find(R == 0);

  double N_nons = sum(1/ps);
  double N_rand = sum(weights);

  // Calculate the loss using the appropriate method
  double loss;
  if (h == "1") {
    arma::mat temp = X0.each_col() % (R / ps / N_nons - R_rand % weights / N_rand);
    loss = arma::accu(arma::sum(temp, 1) % arma::sum(temp, 1));
  } else if (h == "2") {
    arma::mat temp = X0.each_col() % (R / N_nons - R_rand % weights % ps / N_rand);
    loss = arma::accu(arma::sum(temp, 1) % arma::sum(temp, 1));
  } else {
    loss = 0;
  }

  return loss;
}

arma::vec u_theta(arma::vec par,
                  arma::vec R,
                  arma::mat X,
                  arma::vec weights,
                  std::string method_selection,
                  std::string h,
                  Nullable<int> N = R_NilValue,
                  Nullable<NumericVector> pop_totals = R_NilValue,
                  Nullable<double> pop_size = R_NilValue) {

  Environment nonprobsvy_env = Environment::namespace_env("nonprobsvy");

  Rcpp::Function get_method = nonprobsvy_env["get_method"];
  List method = get_method(method_selection);

  Function inv_link = method["make_link_inv"];

  arma::vec theta = par;
  arma::mat X0 = X;
  arma::vec eta_pi = X0 * theta;
  arma::vec ps = as<arma::vec>(inv_link(eta_pi));
  arma::vec R_rand = 1 - R;
  double N_nons = sum(1/ps);

  arma::vec eq;
  arma::mat temp;
  int h_ = stoi(h);
  if (pop_totals.isNull()) {
    switch(h_) {
    case 1:
      temp = X0.each_col() % (R/ps - R_rand % weights);
      eq = sum(temp, 0).t() / N_nons;
      break;
    case 2:
      temp = X0.each_col() % (R - R_rand % weights % ps);
      eq = sum(temp, 0).t() / N_nons;
      break;
    }
  } else { // to consider
    arma::vec total_pop = arma::join_cols(arma::vec{N_nons}, as<arma::vec>(pop_totals));
    eq = (sum(X0.each_col() % (R/ps), 0).t() - total_pop) / N_nons;
  }

  return eq;
}

arma::mat u_theta_der(arma::vec par,
                      arma::vec R,
                      arma::mat X,
                      arma::vec weights,
                      std::string method_selection,
                      std::string h,
                      Nullable<int> N = R_NilValue,
                      Nullable<arma::vec> pop_totals = R_NilValue) {

  Environment nonprobsvy_env = Environment::namespace_env("nonprobsvy");

  Rcpp::Function get_method = nonprobsvy_env["get_method"];
  List method = get_method(method_selection);

  Function inv_link = method["make_link_inv"];

  arma::mat X0 = X;
  arma::vec theta = par;
  //int p = X0.n_cols;
  arma::vec eta_pi = X0 * theta;
  SEXP pss = inv_link(eta_pi);
  arma::vec ps(Rcpp::as<arma::vec>(pss));
  arma::vec R_rand = 1 - R;
  arma::vec psd;
  if (method_selection == "probit") {
    Function inv_link_der = method["make_link_inv_der"];
    SEXP psdd = inv_link_der(eta_pi);
    psd = Rcpp::as<arma::vec>(psdd);
  }

  int n = X0.n_rows;
  int p = X0.n_cols;
  arma::mat mxDer(p, p, arma::fill::zeros);
  double N_nons = sum(1/ps);


  if (h == "1" || !pop_totals.isNull()) {
    if (method_selection == "logit") {
      for(int i = 0; i < n; i++) {
        arma::mat temp = R(i) * (1-ps(i))/ps(i) * X0.row(i).t();
        mxDer += temp * X0.row(i);
      }
    } else if (method_selection == "cloglog") {
        for(int i = 0; i < n; i++) {
          arma::mat temp = R(i) * (1-ps(i))/pow(ps(i), 2) * exp(eta_pi(i)) * X0.row(i).t();
          mxDer += temp * X0.row(i);
      }
    } else if (method_selection == "probit") {
        for(int i = 0; i < n; i++) {
          arma::mat temp = R(i) * psd(i)/pow(ps(i), 2) * X0.row(i).t();
          mxDer += temp * X0.row(i);
      }
    }
    //arma::mat mxDer = (X0.each_col() % (R % (1-ps)/ps)).t() * X0; OK
    //return mxDer.t() / arma::sum(1/ps);
  } else if (h == "2") {
    if (method_selection == "logit") {
      for(int i = 0; i < n; i++) {
        arma::mat temp = R_rand(i) * weights(i) * ps(i)/(exp(eta_pi(i)) + 1) * X0.row(i).t();
        mxDer += temp * X0.row(i);
      }
    } else if (method_selection == "cloglog") {
      for(int i = 0; i < n; i++) {
        arma::mat temp = R_rand(i) * weights(i) * (1-ps(i)) * exp(eta_pi(i)) * X0.row(i).t();
        mxDer += temp * X0.row(i);
      }
    } else if (method_selection == "probit") {
      for(int i = 0; i < n; i++) {
        arma::mat temp = R_rand(i) * weights(i) * psd(i) * X0.row(i).t();
        mxDer += temp * X0.row(i);
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

arma::vec q_lambda_cpp(arma::vec par, double lambda, double a = 3.7) {

  arma::vec abs_par = arma::abs(par);
  arma::vec penaltyd(par.size());

  penaltyd.elem(arma::find(abs_par < lambda)).fill(lambda);

  arma::vec tmp = ((a * lambda) - abs_par.elem(arma::find(abs_par >= lambda))) / (a - 1);
  penaltyd.elem(arma::find(abs_par >= lambda)) = tmp % (tmp > 0);

  // No penalty on the intercept
  penaltyd(0) = 0;

  return penaltyd;
}

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
  //double N = sum(weights);
  //int n = X.n_rows;
  arma::vec init_theta(p, arma::fill::zeros);
  arma::vec par0 = init_theta;
  arma::vec LAMBDA(p, arma::fill::zeros);
  arma::vec par;

  //Environment nonprobsvy_env = Environment::namespace_env("nonprobsvy");
  //Rcpp::Function q_lambda_cpp = nonprobsvy_env["q_lambda"]; // to replace on Cpp function

  int it = 0;
  for (int jj = 1; jj <= maxit; jj++) {
    it++;
    if (warn && it == maxit) {
      Rcpp::warning("Convergence not obtained in %d iterations of fitting algorithm for variables selection", maxit);
      break;
    }

    arma::vec u_theta0v = u_theta(par0,
                                  R,
                                  X,
                                  weights,
                                  method_selection,
                                  h);

    arma::mat u_theta0_derv = u_theta_der(par0,
                                          R,
                                          X,
                                          weights,
                                          method_selection,
                                          h);

    //SEXP q_lambda_output = q_lambda_cpp(par0, lambda);
    //LAMBDA = arma::abs(Rcpp::as<arma::vec>(q_lambda_output)) / (eps + arma::abs(par0));
    arma::vec q_lambda_output = q_lambda_cpp(par0, lambda);
    LAMBDA = arma::abs(q_lambda_output) / (eps + arma::abs(par0));
    //Rcout << "The value of lambda is: " << q_lambda_output << "\n";

    par = par0 + arma::inv(arma::reshape(u_theta0_derv, p, p) + arma::diagmat(LAMBDA)) * (u_theta0v - arma::diagmat(LAMBDA) * par0);

    if (arma::sum(arma::abs(par - par0)) < eps) break;
    if (arma::sum(arma::abs(par - par0)) > 1000) break;

    par0 = par;
  }

   par(arma::find(arma::abs(par) < 0.001)).zeros();
   arma::vec theta_est = par;


  return theta_est;
}



// [[Rcpp::export]]
Rcpp::List cv_nonprobsvy_rcpp(arma::mat X,
                              arma::vec R,
                              arma::vec weights_X,
                              std::string method_selection,
                              std::string h,
                              int maxit,
                              double eps,
                              double lambda_min,
                              int nlambda,
                              int nfolds,
                              double lambda = -1) {

  Environment nonprobsvy_env = Environment::namespace_env("nonprobsvy");
  Rcpp::Function setup_lambda_cpp = nonprobsvy_env["setup_lambda"];

  arma::vec weights;
  arma::vec loss_theta_av(nlambda);
  arma::vec R_ = R;
  arma::mat X_ = X;

  if(lambda == -1) {
    arma::uvec loc_nons = find(R == 1);
    arma::uvec loc_rand = find(R == 0);
    arma::mat X_nons = join_rows(X.rows(loc_nons), weights_X(loc_nons), R(loc_nons));
    arma::mat X_rand = join_rows(X.rows(loc_rand), weights_X(loc_rand), R(loc_rand));

    SEXP lambdas = setup_lambda_cpp(X,
                                    R,
                                    weights_X,
                                    method_selection,
                                    lambda_min,
                                    nlambda);
    arma::vec lambdas1(Rcpp::as<arma::vec>(lambdas));

    // shuffle rows of X_nons and X_rand
    arma::uvec shuffle_nons = arma::shuffle(arma::linspace<arma::uvec>(0, X_nons.n_rows-1, X_nons.n_rows));
    arma::uvec shuffle_rand = arma::shuffle(arma::linspace<arma::uvec>(0, X_rand.n_rows-1, X_rand.n_rows));
    X_nons = X_nons.rows(shuffle_nons);
    X_rand = X_rand.rows(shuffle_rand);

    // generate folds for cross-validation
    arma::uvec folds_nons = arma::randi<arma::uvec>(X_nons.n_rows, arma::distr_param(0, nfolds-1));
    arma::uvec folds_rand = arma::randi<arma::uvec>(X_rand.n_rows, arma::distr_param(0, nfolds-1));

    // randomly sample nfolds folds without replacement
    arma::uvec sample_nons = arma::shuffle(arma::linspace<arma::uvec>(0, nfolds-1, nfolds));
    arma::uvec sample_rand = arma::shuffle(arma::linspace<arma::uvec>(0, nfolds-1, nfolds));


    //arma::vec loss_theta_av(nlambda);
    for(int i = 0; i < nlambda; i++) {
      lambda = lambdas1(i);
      arma::vec loss_theta_vec(nfolds);
    #pragma omp parallel for schedule(dynamic)
      for(int j = 0; j < nfolds; j++) {
        // train data for X_nons
        arma::uvec idx_nons = find(folds_nons != sample_nons(j));
        arma::mat X_nons_train = X_nons.rows(idx_nons);
        // test data for X_nons
        arma::mat X_nons_test = X_nons.rows(find(folds_nons == sample_nons(j)));

        // train data for X_rand
        arma::uvec idx_rand = find(folds_rand != sample_rand(j));
        arma::mat X_rand_train = X_rand.rows(idx_rand);
        // test data for X_rand
        arma::mat X_rand_test = X_rand.rows(find(folds_rand == sample_rand(j)));


        arma::mat X_train = arma::join_cols(X_rand_train, X_nons_train);
        arma::mat X_test = arma::join_cols(X_rand_test, X_nons_test);
        int ncols = X_test.n_cols;

        arma::uvec idxx = arma::regspace<arma::uvec>(0, ncols - 3);

        arma::vec theta_est = fit_nonprobsvy_rcpp(X = X_train.cols(idxx),
                                              R = X_train.col(ncols - 1),
                                              weights = X_train.col(ncols - 2),
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
    int min_idx = 0;
    double min_val = loss_theta_av[0];

    for (int i = 1; i < nlambda; i++) {
      if (loss_theta_av[i] <= min_val) {
        min_idx = i;
        min_val = loss_theta_av[i];
      }
    }
    lambda = lambdas1(min_idx);

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
  //Rcout << "The value of w is: " << theta_selected << "\n";

  return List::create(_["theta_est"] = theta,
                      _["theta_selected"] = theta_selected,
                      _["min"] = loss_theta_av.min(),
                      _["lambda"] = lambda);
}

