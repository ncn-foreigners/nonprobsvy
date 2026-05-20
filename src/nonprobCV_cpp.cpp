#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <limits>
//#include <Eigen/Dense>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
using namespace std;
//using namespace Eigen;

enum class LinkMethod {
  logit,
  probit,
  cloglog
};

enum class PenaltyType {
  SCAD,
  lasso,
  MCP
};

struct EstimatingSystem {
  arma::vec score;
  arma::mat derivative;
};

inline LinkMethod parse_link_method(const std::string& method_selection) {
  if (method_selection == "logit") return LinkMethod::logit;
  if (method_selection == "probit") return LinkMethod::probit;
  if (method_selection == "cloglog") return LinkMethod::cloglog;
  Rcpp::stop("Unknown method selection");
}

inline PenaltyType parse_penalty_type(const std::string& penalty) {
  if (penalty == "SCAD") return PenaltyType::SCAD;
  if (penalty == "lasso") return PenaltyType::lasso;
  if (penalty == "MCP") return PenaltyType::MCP;
  Rcpp::stop("Unknown penalty");
}

inline double clamp_probability_scalar(double ps) {
  const double eps = std::numeric_limits<double>::epsilon();
  const double upper = 1.0 - eps;

  if (std::isnan(ps) ||
      ps == -std::numeric_limits<double>::infinity() ||
      ps < eps) {
    return eps;
  }
  if (ps == std::numeric_limits<double>::infinity() || ps > upper) {
    return upper;
  }
  return ps;
}

inline double stable_cloglog_rate(double eta) {
  const double eps = std::numeric_limits<double>::epsilon();
  const double lower = std::log(eps);
  const double upper = std::log(-std::log(eps));
  return std::exp(std::min(std::max(eta, lower), upper));
}

inline double stable_logit(double eta) {
  if (eta >= 0) {
    const double exp_neg_eta = std::exp(-eta);
    return 1.0 / (1.0 + exp_neg_eta);
  }

  const double exp_eta = std::exp(eta);
  return exp_eta / (1.0 + exp_eta);
}

inline void fill_link_values(const arma::vec& eta,
                             const LinkMethod method_selection,
                             arma::vec& ps,
                             arma::vec& ps_der,
                             arma::vec& cloglog_rate) {
  const arma::uword n = eta.n_elem;
  ps.set_size(n);
  ps_der.set_size(n);
  cloglog_rate.reset();

  switch (method_selection) {
  case LinkMethod::logit:
    for (arma::uword i = 0; i < n; ++i) {
      ps(i) = clamp_probability_scalar(stable_logit(eta(i)));
      ps_der(i) = ps(i) * (1.0 - ps(i));
    }
    break;
  case LinkMethod::probit:
    for (arma::uword i = 0; i < n; ++i) {
      ps(i) = clamp_probability_scalar(R::pnorm5(eta(i), 0.0, 1.0, 1, 0));
      ps_der(i) = R::dnorm4(eta(i), 0.0, 1.0, 0);
    }
    break;
  case LinkMethod::cloglog:
    cloglog_rate.set_size(n);
    for (arma::uword i = 0; i < n; ++i) {
      cloglog_rate(i) = stable_cloglog_rate(eta(i));
      ps(i) = clamp_probability_scalar(-std::expm1(-cloglog_rate(i)));
      ps_der(i) = cloglog_rate(i) * std::exp(-cloglog_rate(i));
      if (!std::isfinite(ps_der(i))) ps_der(i) = 0.0;
    }
    break;
  }
}

inline arma::vec inv_link_cpp(const arma::vec& eta,
                              const LinkMethod method_selection) {
  arma::vec ps;
  arma::vec ps_der;
  arma::vec cloglog_rate;
  fill_link_values(eta, method_selection, ps, ps_der, cloglog_rate);
  return ps;
}

inline double loss_theta(const arma::vec& par,
                         const arma::mat& X,
                         const arma::vec& R,
                         const arma::vec& weights,
                         const LinkMethod method_selection,
                         const int& gee_h_fun,
                         const arma::uvec& selected,
                         const arma::vec& pop_totals,
                         const bool has_pop_totals) {
  if (selected.n_elem == 0) {
    return 0.0;
  }

  arma::vec eta_pi(X.n_rows, arma::fill::zeros);
  for (arma::uword k = 0; k < selected.n_elem; ++k) {
    const arma::uword col = selected(k);
    eta_pi += X.col(col) * par(col);
  }

  const arma::vec ps = inv_link_cpp(eta_pi, method_selection);
  const arma::vec R_rand = 1.0 - R;
  const arma::vec inv_ps = 1.0 / ps;

  const double N_nons = arma::sum(inv_ps);
  const double N_rand = arma::sum(weights);
  arma::vec row_factor;

  if (!has_pop_totals) {
    if (gee_h_fun == 1) {
      row_factor = R % weights % inv_ps / N_nons - R_rand % weights / N_rand;
    } else if (gee_h_fun == 2) {
      row_factor = R % weights / N_nons - R_rand % weights % ps / N_rand;
    } else {
      return 0.0;
    }
  } else {
    row_factor = R % weights % inv_ps / N_nons;
  }

  double loss = 0.0;
  for (arma::uword k = 0; k < selected.n_elem; ++k) {
    const arma::uword col = selected(k);
    double col_balance = arma::dot(X.col(col), row_factor);
    if (has_pop_totals) {
      col_balance -= pop_totals(col) / N_rand;
    }
    loss += col_balance * col_balance;
  }

  return loss;
}

inline EstimatingSystem estimating_system(const arma::vec& par,
                                          const arma::mat& X,
                                          const arma::vec& R,
                                          const arma::vec& R_rand,
                                          const arma::vec& weights,
                                          const LinkMethod method_selection,
                                          const int& gee_h_fun,
                                          const arma::vec& pop_totals,
                                          const bool has_pop_totals) {
  const arma::vec eta_pi = X * par;
  arma::vec ps;
  arma::vec ps_der;
  arma::vec cloglog_rate;
  fill_link_values(eta_pi, method_selection, ps, ps_der, cloglog_rate);

  const arma::vec inv_ps = 1.0 / ps;
  const arma::vec inv_ps_sq = arma::square(inv_ps);
  const double N_nons = arma::sum(inv_ps);

  arma::vec score_weight;
  if (!has_pop_totals) {
    switch (gee_h_fun) {
    case 1:
      score_weight = R % weights % inv_ps - R_rand % weights;
      break;
    case 2:
      score_weight = R % weights - R_rand % weights % ps;
      break;
    default:
      Rcpp::stop("Unknown h selection");
    }
  } else {
    score_weight = R % weights % inv_ps;
  }

  EstimatingSystem out;
  out.score = X.t() * score_weight;
  if (has_pop_totals) {
    out.score -= pop_totals;
  }
  out.score /= N_nons;

  arma::vec coeff;
  if (gee_h_fun == 1 || has_pop_totals) {
    switch (method_selection) {
    case LinkMethod::logit:
      coeff = R % weights % (1.0 - ps) % inv_ps;
      break;
    case LinkMethod::cloglog:
      coeff = R % weights % (1.0 - ps) % inv_ps_sq % cloglog_rate;
      break;
    case LinkMethod::probit:
      coeff = R % weights % ps_der % inv_ps_sq;
      break;
    }
  } else if (gee_h_fun == 2) {
    switch (method_selection) {
    case LinkMethod::logit:
      coeff = R_rand % weights % ps % (1.0 - ps);
      break;
    case LinkMethod::cloglog:
      coeff = R_rand % weights % (1.0 - ps) % cloglog_rate;
      break;
    case LinkMethod::probit:
      coeff = R_rand % weights % ps_der;
      break;
    }
  } else {
    Rcpp::stop("Unknown h selection");
  }

  out.derivative = X.t() * (X.each_col() % coeff);
  out.derivative /= N_nons;

  return out;
}

arma::vec q_lambda_cpp(const arma::vec& par,
                       double lambda,
                       const PenaltyType penalty,
                       double a) {
  arma::vec penaltyd(par.size(), arma::fill::zeros);

  switch (penalty) {
  case PenaltyType::SCAD: {
    double abs_par_i;
    for (std::size_t i = 0; i < par.size(); i++) {
      abs_par_i = std::abs(par[i]);
      if (abs_par_i <= lambda) {
        penaltyd[i] = lambda;
      } else {
        double tmp = ((a * lambda) - abs_par_i) / (a - 1);
        penaltyd[i] = std::max(tmp, 0.0);
      }
    }
    break;
  }
  case PenaltyType::lasso:
    penaltyd = lambda * arma::sign(par);
    break;
  case PenaltyType::MCP:
    for (std::size_t i = 0; i < par.size(); i++) {
      if (std::abs(par[i]) <= a * lambda) {
        if (par[i] < 0) {
          penaltyd[i] = -(lambda - std::abs(par[i]) / a);
        } else if (par[i] > 0) {
          penaltyd[i] = lambda - std::abs(par[i]) / a;
        } else {
          penaltyd[i] = 0;
        }
      } else {
        penaltyd[i] = 0;
      }
    }
    break;
  }

  return penaltyd;
}

// [[Rcpp::export]]
arma::vec q_lambda_test_cpp(const arma::vec& par,
                            double lambda,
                            const std::string& penalty,
                            double a) {
  return q_lambda_cpp(par, lambda, parse_penalty_type(penalty), a);
}

arma::uvec r_sample_zero_based(int population_size, int sample_size, bool replace) {
  Rcpp::IntegerVector sampled = Rcpp::sample(population_size, sample_size, replace);
  arma::uvec out = Rcpp::as<arma::uvec>(sampled);
  out -= 1;
  return out;
}

// HybridNonLinearSolver (?)
// pass by reference where possible to avoid unnecessary copying
arma::vec fit_nonprobsvy_rcpp(const arma::mat& X,
                              const arma::vec& R,
                              const arma::vec& weights,
                              const LinkMethod method_selection,
                              const int& gee_h_fun,
                              double lambda,
                              int maxit,
                              double eps,
                              const PenaltyType penalty,
                              double a,
                              const arma::vec& pop_totals,
                              const bool has_pop_totals,
                              bool warn = false
) {
  const int p = X.n_cols;
  const arma::vec R_rand = 1.0 - R;

  arma::vec par0(p, arma::fill::zeros);
  arma::vec LAMBDA(p, arma::fill::zeros);
  arma::vec par(p, arma::fill::zeros);

  bool converged = false;
  for (int jj = 1; jj <= maxit; jj++) {

    EstimatingSystem system = estimating_system(
      par0,
      X,
      R,
      R_rand,
      weights,
      method_selection,
      gee_h_fun,
      pop_totals,
      has_pop_totals
    );
    arma::vec q_lambda_output = q_lambda_cpp(par0, lambda, penalty, a);

    LAMBDA = arma::abs(q_lambda_output) / (eps + arma::abs(par0));

    arma::mat system_matrix = system.derivative;
    system_matrix.diag() += LAMBDA;
    arma::vec rhs = system.score - LAMBDA % par0;
    arma::vec step(p, arma::fill::zeros);
    bool solved = arma::solve(step, system_matrix, rhs);
    if (!solved) {
      step = arma::pinv(system_matrix) * rhs;
    }
    par = par0 + step;

    if (arma::sum(arma::abs(par - par0)) < eps) { converged = true; break; }
    if (arma::sum(arma::abs(par - par0)) > 1000) break;

    par0 = par;
  }
  if (warn && !converged) {
    Rcpp::warning("Convergence not obtained in %d iterations of fitting algorithm for variable selection", maxit);
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
                              double lambda = -1) {

  Environment nonprobsvy_env = Environment::namespace_env("nonprobsvy");
  Rcpp::Function setup_lambda_cpp = nonprobsvy_env["setup_lambda"];

  const LinkMethod link_method = parse_link_method(method_selection);
  const PenaltyType penalty_type = parse_penalty_type(penalty);
  const bool has_pop_totals = !pop_totals.isNull();
  const arma::vec pop_totals_vec = has_pop_totals ? Rcpp::as<arma::vec>(pop_totals) : arma::vec();

  arma::vec loss_theta_av(nlambda, arma::fill::zeros);
  const arma::vec& R_ = R;
  const arma::mat& X_ = X;

  if (lambda == -1) {
    arma::uvec loc_nons = arma::find(R == 1);
    arma::uvec loc_rand = arma::find(R == 0);

    SEXP lambdas = setup_lambda_cpp(X, R, weights_X, method_selection, lambda_min, nlambda, pop_totals);
    arma::vec lambdas1(Rcpp::as<arma::vec>(lambdas));

    arma::uvec folds_nons = r_sample_zero_based(nfolds, static_cast<int>(loc_nons.n_elem), true);
    arma::uvec folds_rand = r_sample_zero_based(nfolds, static_cast<int>(loc_rand.n_elem), true);

    arma::uvec sample_nons = r_sample_zero_based(nfolds, nfolds, false);
    arma::uvec sample_rand = r_sample_zero_based(nfolds, nfolds, false);

    arma::mat loss_theta_fld(nfolds, nlambda, arma::fill::zeros);
    //#pragma omp parallel for
    for (int j = 0; j < nfolds; j++) {
      if (verbose) {
        wcout << "Starting CV fold #" << j + 1 << endl;
      }

      arma::uvec idx_nons_train = arma::find(folds_nons != sample_nons(j));
      arma::uvec idx_nons_test = arma::find(folds_nons == sample_nons(j));
      arma::uvec idx_rand_train = arma::find(folds_rand != sample_rand(j));
      arma::uvec idx_rand_test = arma::find(folds_rand == sample_rand(j));

      arma::uvec train_rows = arma::join_cols(
        loc_rand.elem(idx_rand_train),
        loc_nons.elem(idx_nons_train)
      );
      arma::uvec test_rows = arma::join_cols(
        loc_rand.elem(idx_rand_test),
        loc_nons.elem(idx_nons_test)
      );

      arma::mat X_train = X.rows(train_rows);
      arma::vec R_train = R.elem(train_rows);
      arma::vec weights_train = weights_X.elem(train_rows);

      arma::mat X_test = X.rows(test_rows);
      arma::vec R_test = R.elem(test_rows);
      arma::vec weights_test = weights_X.elem(test_rows);

      //#pragma omp parallel for
      for (int i = 0; i < nlambda; i++) {
        arma::vec theta_est = fit_nonprobsvy_rcpp(
          X_train,
          R_train,
          weights_train,
          link_method,
          gee_h_fun,
          lambdas1(i),
          maxit,
          eps,
          penalty_type,
          a,
          pop_totals_vec,
          has_pop_totals
        );

        const arma::uvec selected = arma::find(theta_est != 0);
        double loss = loss_theta(
          theta_est,
          X_test,
          R_test,
          weights_test,
          link_method,
          gee_h_fun,
          selected,
          pop_totals_vec,
          has_pop_totals
        );
        loss_theta_fld(j, i) = loss;
      }
    }

    for (int i = 0; i < nlambda; i++) {
      loss_theta_av(i) = mean(loss_theta_fld.col(i));
    }
    lambda = lambdas1(loss_theta_av.index_min());
  }

  arma::vec theta = fit_nonprobsvy_rcpp(
    X_,
    R_,
    weights_X,
    link_method,
    gee_h_fun,
    lambda,
    maxit,
    eps,
    penalty_type,
    a,
    pop_totals_vec,
    has_pop_totals,
    verbose
  );

  arma::uvec theta_selected = arma::find(theta != 0);

  return List::create(_["theta_est"] = theta,
                      _["theta_selected"] = theta_selected,
                      _["min"] = loss_theta_av.min(),
                      _["lambda"] = lambda,
                      _["cv_error"] = loss_theta_av
  );
}
