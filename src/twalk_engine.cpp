// Compiled t-walk engine.
//
// This file is a line-by-line translation of the reference implementation in
// R/twalk_engine.R (Christen & Fox, 2010). Random numbers are drawn from R's
// generator in exactly the same order as the reference implementation, so
// both engines produce identical chains for the same seed; the test suite
// checks this property.
//
// The engine is resumable: its state lives in an external pointer. When the
// user's log posterior raises an R error, the error unwinds out of
// `twalk_engine_run()` (Rcpp's unwind protection keeps the C++ state
// consistent), the R driver records it, calls `twalk_engine_reject_pending()`
// and resumes sampling. This avoids wrapping every evaluation in `tryCatch()`.

#include <Rcpp.h>
#include <cfloat>
#include <cmath>
#include <vector>

using namespace Rcpp;

namespace {

enum Kernel { TRAVERSE = 0, WALK = 1, BLOW = 2, HOP = 3 };

struct TwalkEngine {
  int n;
  int n_iter;
  int thin;
  int n_store;
  int iter;
  int n_evals;

  std::vector<double> x, xp;
  double lp_x, lp_xp;

  double p_phi, aw, at;
  double cum_probs[4];
  bool store_companion;

  NumericMatrix out_x, out_xp;
  NumericVector out_lp_x, out_lp_xp;
  RObject names;

  int proposed[4], accepted[4], moved[4];

  // State of an iteration interrupted while evaluating the log posterior.
  bool pending;
};

TwalkEngine* get_engine(SEXP ptr) {
  XPtr<TwalkEngine> engine(ptr);
  if (engine.get() == NULL) {
    stop("Invalid t-walk engine state.");
  }
  return engine.get();
}

// Same operations, in the same order, as `simulate_beta()` in R.
double simulate_beta(double at) {
  if (unif_rand() < (at - 1.0) / (2.0 * at)) {
    return std::exp(1.0 / (at + 1.0) * std::log(unif_rand()));
  }
  return std::exp(1.0 / (1.0 - at) * std::log(unif_rand()));
}

// -log g(h | a, b) for the Blow (centre b) and Hop (centre a) proposals,
// matching `log_density_blow()` and `log_density_hop()` in R. Only called
// with n_phi > 0: proposals that select no coordinate are identity moves.
double neg_log_density(const std::vector<double>& h,
                       const std::vector<double>& a,
                       const std::vector<double>& b,
                       const std::vector<double>& phi,
                       int n_phi, bool hop) {
  double sigma = 0.0;
  for (std::size_t j = 0; j < h.size(); ++j) {
    double s = phi[j] * std::fabs(b[j] - a[j]);
    if (s > sigma) sigma = s;
  }
  if (hop) sigma = sigma / 3.0;
  if (sigma < DBL_EPSILON) sigma = DBL_EPSILON;

  const std::vector<double>& centre = hop ? a : b;
  long double ss = 0.0;
  for (std::size_t j = 0; j < h.size(); ++j) {
    if (phi[j] != 0.0) {
      double d = h[j] - centre[j];
      ss += d * d;
    }
  }

  return (n_phi / 2.0) * std::log(2.0 * M_PI) + n_phi * std::log(sigma) +
    0.5 * static_cast<double>(ss) / (sigma * sigma);
}

void store_iteration(TwalkEngine* e) {
  e->iter += 1;
  if (e->iter % e->thin != 0) {
    return;
  }

  int row = e->iter / e->thin - 1;
  for (int j = 0; j < e->n; ++j) {
    e->out_x(row, j) = e->x[j];
  }
  e->out_lp_x[row] = e->lp_x;

  if (e->store_companion) {
    for (int j = 0; j < e->n; ++j) {
      e->out_xp(row, j) = e->xp[j];
    }
    e->out_lp_xp[row] = e->lp_xp;
  }
}

double checked_log_posterior(TwalkEngine* e, SEXP call,
                             const std::vector<double>& h) {
  NumericVector params(h.begin(), h.end());
  if (!e->names.isNULL()) {
    params.attr("names") = e->names;
  }
  SETCADR(call, params);

  e->pending = true;
  Shield<SEXP> value(Rcpp_fast_eval(call, R_GlobalEnv));
  e->pending = false;
  e->n_evals += 1;

  int type = TYPEOF(value);
  if ((type != REALSXP && type != INTSXP) || Rf_length(value) != 1 ||
      Rf_isFactor(value)) {
    stop("`log_posterior` must return a single numeric value.");
  }

  double lp;
  if (type == INTSXP) {
    int v = INTEGER(value)[0];
    lp = (v == NA_INTEGER) ? NA_REAL : static_cast<double>(v);
  } else {
    lp = REAL(value)[0];
  }

  if (ISNAN(lp) || lp == R_PosInf) {
    stop("`log_posterior` must return a finite value or -Inf outside the support.");
  }

  return lp;
}

// One t-walk transition. Mirrors `.engine_r_step()` in R/twalk_engine.R.
void step(TwalkEngine* e, SEXP call) {
  const int n = e->n;

  double u_kernel = unif_rand();
  int kernel = HOP;
  for (int k = 0; k < 3; ++k) {
    if (u_kernel < e->cum_probs[k]) {
      kernel = k;
      break;
    }
  }
  bool move_xp = unif_rand() < 0.5;

  // `a` is the point being moved and `b` the other point: the paper's
  // proposals (x, h(x', x)) and (h(x, x'), x') respectively.
  std::vector<double>& a = move_xp ? e->xp : e->x;
  const std::vector<double>& b = move_xp ? e->x : e->xp;
  double& lp_a = move_xp ? e->lp_xp : e->lp_x;

  std::vector<double> h(n), phi(n);
  int n_phi = 0;
  double beta = 0.0;

  if (kernel == TRAVERSE) {
    beta = simulate_beta(e->at);
    for (int j = 0; j < n; ++j) {
      phi[j] = (unif_rand() < e->p_phi) ? 1.0 : 0.0;
      n_phi += (phi[j] != 0.0);
      h[j] = (phi[j] != 0.0) ? b[j] + beta * (b[j] - a[j]) : a[j];
    }
  } else if (kernel == WALK) {
    std::vector<double> u(n);
    for (int j = 0; j < n; ++j) u[j] = unif_rand();
    for (int j = 0; j < n; ++j) {
      phi[j] = (unif_rand() < e->p_phi) ? 1.0 : 0.0;
      n_phi += (phi[j] != 0.0);
    }
    const double aw = e->aw;
    for (int j = 0; j < n; ++j) {
      double z = (aw / (1.0 + aw)) * (aw * (u[j] * u[j]) + 2.0 * u[j] - 1.0);
      z = z * phi[j];
      h[j] = a[j] + (a[j] - b[j]) * z;
    }
  } else {
    for (int j = 0; j < n; ++j) {
      phi[j] = (unif_rand() < e->p_phi) ? 1.0 : 0.0;
      n_phi += (phi[j] != 0.0);
    }
    double sigma = 0.0;
    for (int j = 0; j < n; ++j) {
      double s = phi[j] * std::fabs(b[j] - a[j]);
      if (s > sigma) sigma = s;
    }
    if (kernel == HOP) sigma = sigma / 3.0;
    if (sigma < DBL_EPSILON) sigma = DBL_EPSILON;

    std::vector<double> noise(n);
    for (int j = 0; j < n; ++j) noise[j] = norm_rand() * phi[j];

    for (int j = 0; j < n; ++j) {
      if (kernel == BLOW) {
        h[j] = b[j] * phi[j] + sigma * noise[j] + a[j] * (1.0 - phi[j]);
      } else {
        h[j] = (phi[j] != 0.0) ? a[j] + sigma * noise[j] : a[j];
      }
    }
  }

  e->proposed[kernel] += 1;

  double alpha = 0.0;
  double lp_h = lp_a;

  if (n_phi == 0) {
    // Identity proposal: accepted with probability one, no evaluation needed.
    alpha = 1.0;
  } else {
    bool degenerate = false;
    if (kernel != TRAVERSE) {
      degenerate = true;
      for (int j = 0; j < n; ++j) {
        if (h[j] != b[j]) {
          degenerate = false;
          break;
        }
      }
    }

    if (!degenerate) {
      lp_h = checked_log_posterior(e, call, h);

      if (R_FINITE(lp_h)) {
        double log_alpha = lp_h - lp_a;
        if (kernel == TRAVERSE) {
          log_alpha = log_alpha + (n_phi - 2) * std::log(beta);
        } else if (kernel == BLOW || kernel == HOP) {
          bool hop = (kernel == HOP);
          double w1 = neg_log_density(h, a, b, phi, n_phi, hop);
          double w2 = neg_log_density(a, h, b, phi, n_phi, hop);
          log_alpha = log_alpha + (w1 - w2);
        }
        alpha = (log_alpha >= 0.0) ? 1.0 : std::exp(log_alpha);
      }
    }
  }

  if (unif_rand() < alpha) {
    bool changed = false;
    for (int j = 0; j < n; ++j) {
      if (h[j] != a[j]) {
        changed = true;
        break;
      }
    }
    a = h;
    lp_a = lp_h;
    e->accepted[kernel] += 1;
    if (changed) e->moved[kernel] += 1;
  }

  store_iteration(e);
}

}  // namespace

// [[Rcpp::export]]
SEXP twalk_engine_new(NumericVector x0, NumericVector xp0,
                      double lp_x0, double lp_xp0,
                      int n_iter, int thin, double p_phi,
                      double aw, double at, NumericVector cum_probs,
                      bool store_companion, SEXP names) {
  TwalkEngine* e = new TwalkEngine();
  e->n = x0.size();
  e->n_iter = n_iter;
  e->thin = thin;
  e->n_store = n_iter / thin;
  e->iter = 0;
  e->n_evals = 0;
  e->x.assign(x0.begin(), x0.end());
  e->xp.assign(xp0.begin(), xp0.end());
  e->lp_x = lp_x0;
  e->lp_xp = lp_xp0;
  e->p_phi = p_phi;
  e->aw = aw;
  e->at = at;
  for (int k = 0; k < 4; ++k) {
    e->cum_probs[k] = cum_probs[k];
    e->proposed[k] = 0;
    e->accepted[k] = 0;
    e->moved[k] = 0;
  }
  e->store_companion = store_companion;
  e->out_x = NumericMatrix(e->n_store, e->n);
  e->out_lp_x = NumericVector(e->n_store);
  if (store_companion) {
    e->out_xp = NumericMatrix(e->n_store, e->n);
    e->out_lp_xp = NumericVector(e->n_store);
  }
  e->names = names;
  e->pending = false;

  XPtr<TwalkEngine> ptr(e, true);
  return ptr;
}

// [[Rcpp::export]]
int twalk_engine_run(SEXP ptr, Function log_posterior, int n_steps) {
  TwalkEngine* e = get_engine(ptr);
  if (e->pending) {
    stop("The t-walk engine has a pending evaluation; reject it first.");
  }

  Language call(log_posterior, R_NilValue);
  int last = e->iter + n_steps;
  if (last > e->n_iter) last = e->n_iter;

  RNGScope rng_scope;
  while (e->iter < last) {
    if (e->iter % 1024 == 0) {
      checkUserInterrupt();
    }
    step(e, call);
  }

  return e->iter;
}

// Complete an iteration whose log-posterior evaluation raised an error: the
// proposal is treated as outside the support (alpha = 0). The acceptance
// uniform is still drawn, as it is for any proposal with alpha = 0.
// [[Rcpp::export]]
int twalk_engine_reject_pending(SEXP ptr) {
  TwalkEngine* e = get_engine(ptr);
  if (!e->pending) {
    stop("The t-walk engine has no pending evaluation.");
  }
  e->pending = false;
  {
    RNGScope rng_scope;
    unif_rand();
  }
  store_iteration(e);
  return e->iter;
}

// [[Rcpp::export]]
bool twalk_engine_pending(SEXP ptr) {
  return get_engine(ptr)->pending;
}

// [[Rcpp::export]]
int twalk_engine_iter(SEXP ptr) {
  return get_engine(ptr)->iter;
}

// [[Rcpp::export]]
List twalk_engine_result(SEXP ptr) {
  TwalkEngine* e = get_engine(ptr);
  int n_done = e->iter / e->thin;

  NumericMatrix samples = e->out_x;
  NumericVector lp = e->out_lp_x;
  if (n_done == 0) {
    samples = NumericMatrix(0, e->n);
    lp = NumericVector(0);
  } else if (n_done < e->n_store) {
    samples = e->out_x(Range(0, n_done - 1), _);
    lp = e->out_lp_x[Range(0, n_done - 1)];
  }

  RObject companion = R_NilValue;
  RObject companion_lp = R_NilValue;
  if (e->store_companion) {
    if (n_done < e->n_store) {
      if (n_done == 0) {
        companion = NumericMatrix(0, e->n);
        companion_lp = NumericVector(0);
      } else {
        NumericMatrix cs = e->out_xp(Range(0, n_done - 1), _);
        NumericVector clp = e->out_lp_xp[Range(0, n_done - 1)];
        companion = cs;
        companion_lp = clp;
      }
    } else {
      companion = e->out_xp;
      companion_lp = e->out_lp_xp;
    }
  }

  IntegerVector proposed(4), accepted(4), moved(4);
  for (int k = 0; k < 4; ++k) {
    proposed[k] = e->proposed[k];
    accepted[k] = e->accepted[k];
    moved[k] = e->moved[k];
  }

  return List::create(
    _["samples"] = samples,
    _["companion_samples"] = companion,
    _["log_posterior"] = lp,
    _["companion_log_posterior"] = companion_lp,
    _["proposed"] = proposed,
    _["accepted"] = accepted,
    _["moved"] = moved,
    _["x"] = NumericVector(e->x.begin(), e->x.end()),
    _["xp"] = NumericVector(e->xp.begin(), e->xp.end()),
    _["lp_x"] = e->lp_x,
    _["lp_xp"] = e->lp_xp,
    _["iter"] = e->iter,
    _["n_evals"] = e->n_evals
  );
}
