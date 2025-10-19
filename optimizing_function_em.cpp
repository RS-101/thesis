// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <algorithm>
#include <vector>
#include <cmath>
using namespace Rcpp;

struct ModelData {
  // --- counts (ints) ---
  int J, C, K_tilde, U, N_tilde, M, W, N1_obs_of_T_star, U_pos_obs_of_T_star,
  N, N_star, M_mark, I, K, I_mark;
  
  // --- scalars (double) ---
  double s_max, R_max, e_star_max;
  
  // --- vectors ---
  NumericVector s_j, L_c, t_c, e_k, L_u, t_u, t_m_in_N_tilde,
  L_m, R_m_vec, t_m, E_star, T_star,
  L_bar, R_bar, s_j_c, s_j_full, Q, Q_i_mark, A_union;
  
  // counts per unique time (int vectors)
  IntegerVector d_n;
  arma::rowvec c_k;
  
  // --- 2-col interval matrices ---
  NumericMatrix A_m, A_u, A_c, full_A_m, Q_i;
  
  // --- alpha and beta matrices --- 
  NumericMatrix alpha_ij;  // I_mark x (J+C)
  NumericMatrix beta_im;   // I x M'
  
  
  // ---------- (18) alpha: I' x (J+C) ----------
  void cal_alpha() {
    int I = Q_i.nrow();
    int JpC = s_j_full.size();
    
    alpha_ij = Rcpp::NumericMatrix(I_mark, JpC);
    for (int i = 0; i < I_mark; ++i) {
      double left = (i < I) ? Q_i(i,0) : Q_i_mark[i - I];
      for (int j = 0; j < JpC; ++j)
        alpha_ij(i,j) = (s_j_full[j] <= left) ? 1.0 : 0.0;
    }
  }
  
  // ---------- (19) beta: I x M' ----------
  void cal_beta() {
    int I = Q_i.nrow();
    int M_mark = full_A_m.nrow();
    beta_im = Rcpp::NumericMatrix(I, M_mark);
    
    for (int i = 0; i < I; ++i) {
      double Li = Q_i(i,0), Ri = Q_i(i,1);
      for (int m = 0; m < M_mark; ++m) {
        double Lm = full_A_m(m,0), Rm = full_A_m(m,1);
        beta_im(i,m) = (Lm <= Li && Ri <= Rm) ? 1.0 : 0.0;
      }
    }
  }
  
  // ctor from named R list and calculates alpha and beta
  ModelData(const List& x) {
    // ints
    J  = as<int>(x["J"]);
    C  = as<int>(x["C"]);
    K_tilde = as<int>(x["K_tilde"]);
    U  = as<int>(x["U"]);
    N_tilde = as<int>(x["N_tilde"]);
    M  = as<int>(x["M"]);
    W  = as<int>(x["W"]);
    N1_obs_of_T_star = as<int>(x["N1_obs_of_T_star"]);
    U_pos_obs_of_T_star = as<int>(x["U_pos_obs_of_T_star"]);
    N  = as<int>(x["N"]);
    N_star = as<int>(x["N_star"]);
    M_mark = as<int>(x["M_mark"]);
    I  = as<int>(x["I"]);
    K  = as<int>(x["K"]);
    I_mark = as<int>(x["I_mark"]);
    
    // scalars
    s_max = as<double>(x["s_max"]);
    R_max = as<double>(x["R_max"]);
    e_star_max = as<double>(x["e_star_max"]);
    
    // vectors
    s_j = as<NumericVector>(x["s_j"]);
    L_c = as<NumericVector>(x["L_c"]);
    t_c = as<NumericVector>(x["t_c"]);
    e_k = as<NumericVector>(x["e_k"]);
    L_u = as<NumericVector>(x["L_u"]);
    t_u = as<NumericVector>(x["t_u"]);
    t_m_in_N_tilde = as<NumericVector>(x["t_m_in_N_tilde"]);
    L_m = as<NumericVector>(x["L_m"]);
    R_m_vec = as<NumericVector>(x["R_m"]);
    t_m = as<NumericVector>(x["t_m"]);
    E_star = as<NumericVector>(x["E_star"]);
    T_star = as<NumericVector>(x["T_star"]);
    L_bar = as<NumericVector>(x["L_bar"]);
    R_bar = as<NumericVector>(x["R_bar"]);
    s_j_c = as<NumericVector>(x["s_j_c"]);
    s_j_full = as<NumericVector>(x["s_j_full"]);
    Q = as<NumericVector>(x["Q"]);
    Q_i_mark = as<NumericVector>(x["Q_i_mark"]);
    A_union = as<NumericVector>(x["A_union"]);
    
    // integer vectors
    c_k = Rcpp::as<arma::rowvec>(x["c_k"]);   // coerces integer
    d_n = as<IntegerVector>(x["d_n"]);
    
    // matrices (2 columns) — the important part
    A_m      = as<NumericMatrix>(x["A_m"]);
    A_u      = as<NumericMatrix>(x["A_u"]);
    A_c      = as<NumericMatrix>(x["A_c"]);
    full_A_m = as<NumericMatrix>(x["full_A_m"]);
    Q_i      = as<NumericMatrix>(x["Q_i"]);
  
    cal_alpha();
    cal_beta();

    // checks
    if ((int)s_j.size() != J) stop("length(s_j) != J");
    if ((int)c_k.size() != K) stop("length(c_k) != K");
    if ((int)E_star.size() != K) stop("length(E_star) != K");
    if (Q_i.ncol() != 2) stop("Q_i must have 2 columns");
    if (A_m.ncol() != 2 || A_u.ncol() != 2 || A_c.ncol() != 2 || full_A_m.ncol() != 2)
      stop("All A_* matrices must have 2 columns");
  }
};

// [[Rcpp::export]]
SEXP make_model_data(List x) {
  return XPtr<ModelData>(new ModelData(x), true);
}

// ---------- workspace + helpers ----------

struct Workspace {
  // parameters updated by EM
  Rcpp::NumericVector lambda_n;  // length N
  Rcpp::NumericVector z_i;       // length I_mark
  
  Rcpp::NumericVector t_sorted;   // sorted t_n
  Rcpp::NumericVector log_prefix; // prefix sum of log(1 - lambda_n) over sorted order
  Rcpp::NumericVector zero_prefix;// prefix count of zeros in (1 - lambda_n)
  
  // E-step pieces
  Rcpp::NumericMatrix mu_mi;     // M x I
  Rcpp::NumericMatrix mu_bar_ji; // J x I_mark
  Rcpp::NumericMatrix eta_ui;    // U x I_mark
  Rcpp::NumericMatrix gamma_ci;  // C x I_mark
  
  // aggregation for lambda
  Rcpp::NumericMatrix rho_mn;    // M x N
  Rcpp::NumericMatrix pi_un;     // U x N
  Rcpp::NumericMatrix sigma_cn;  // C x N
};

// Precompute sorted order and prefix arrays.
// - Sort indices by md.T_star
// - t_sorted[i] = sorted T_star
// - log_prefix[i+1] = sum_{k<=i, lambda_k!=1} log(1 - lambda_k)
// - zero_prefix[i+1] = count_{k<=i} [lambda_k == 1]
void setup_prod(const ModelData& md, Workspace& ws) {
  const R_xlen_t N = md.T_star.size();
  if (ws.lambda_n.size() != N) {
    stop("lambda_n length (%d) must match T_star length (%d).", 
         (int)ws.lambda_n.size(), (int)N);
  }
  
  // build sorted index of T_star
  std::vector<int> idx(N);
  for (R_xlen_t i = 0; i < N; ++i) idx[i] = static_cast<int>(i);
  std::sort(idx.begin(), idx.end(),
            [&](int a, int b){ return md.T_star[a] < md.T_star[b]; });
  
  ws.t_sorted  = NumericVector(N);
  ws.log_prefix = NumericVector(N + 1);
  ws.zero_prefix = NumericVector(N + 1);
  
  ws.log_prefix[0]  = 0.0;
  ws.zero_prefix[0] = 0.0;
  
  // Treat lambda exactly equal to 1 as a zero factor (product becomes 0 if present in interval)
  const double TOL = 0.0; // set to, e.g., 1e-15 if you want tolerance
  for (R_xlen_t i = 0; i < N; ++i) {
    int k = idx[i];
    const double t   = md.T_star[k];
    const double lam = ws.lambda_n[k];
    
    ws.t_sorted[i] = t;
    
    const bool is_one = std::abs(lam - 1.0) <= TOL;
    ws.zero_prefix[i + 1] = ws.zero_prefix[i] + (is_one ? 1.0 : 0.0);
    
    // Only add to log-prefix if not zero; if zero, keep running sum unchanged.
    if (!is_one) {
      const double one_minus = 1.0 - lam;
      if (one_minus <= 0.0) {
        // Defensive: log undefined. You can choose to stop() or handle differently.
        stop("Encountered 1 - lambda_n <= 0 at sorted index %d (lambda=%.17g).", (int)i, lam);
      }
      ws.log_prefix[i + 1] = ws.log_prefix[i] + std::log(one_minus);
    } else {
      ws.log_prefix[i + 1] = ws.log_prefix[i];
    }
  }
}

// Evaluate product_{t_n* in (L, R)} (1 - lambda_n).
// We implement (L, R) with L open, R open by using:
//   i = upper_bound(t_sorted, L)   -> first index with t > L
//   j = lower_bound(t_sorted, R)   -> first index with t >= R
// If any lambda==1 in (i, j-1), product is 0.
// Else product = exp(log_prefix[j] - log_prefix[i]).
double evaluate(const Workspace& ws, double L, double R) {
  const R_xlen_t N = ws.t_sorted.size();
  if (ws.log_prefix.size() != N + 1 || ws.zero_prefix.size() != N + 1)
    stop("Workspace not initialized. Call setup_prod first.");
  
  if (R <= L || N == 0) return 1.0;
  
  // Binary searches on sorted t
  const double* begin = &ws.t_sorted[0];
  const double* end   = begin + N;
  
  const R_xlen_t i = std::upper_bound(begin, end, L) - begin; // t > L
  const R_xlen_t j = std::lower_bound(begin, end, R) - begin; // t >= R
  
  if (j <= i) return 1.0; // empty interval
  
  const double zeros_in_range = ws.zero_prefix[j] - ws.zero_prefix[i];
  if (zeros_in_range > 0.0) return 0.0;
  
  const double log_prod = ws.log_prefix[j] - ws.log_prefix[i];
  return std::exp(log_prod);
}

// Move to the next/previous representable double (1 ULP)
inline double next_double(double x) {
  return std::nextafter(x, std::numeric_limits<double>::infinity());
}
inline double prev_double(double x) {
  return std::nextafter(x, -std::numeric_limits<double>::infinity());
}



static inline void na0_inplace(Rcpp::NumericMatrix m) {
  for (R_xlen_t i = 0; i < m.nrow(); ++i)
    for (R_xlen_t j = 0; j < m.ncol(); ++j)
      if (!R_finite(m(i,j))) m(i,j) = 0.0;
}


// membership of a point x in (L,R) with open/closed flags
static inline bool in_interval(
    double x,
    double L, double R,
    bool L_open, bool R_open) {
  bool left  = L_open ? (x >  L) : (x >= L);
  bool right = R_open ? (x <  R) : (x <= R);
  return left && right;
}

// [a,b] (closed) subset-of (L,R) (with flags on the outer interval)
static inline bool closed_subset_of(
    double a, double b, 
    double L, double R, 
    bool L_open, bool R_open) {
  bool left  = L_open ? (a >  L) : (a >= L);
  bool right = R_open ? (b <  R) : (b <= R);
  return left && right;
}

// vector<bool> : is each Q_i row a subset of (L,R) with flags
static inline Rcpp::LogicalVector Qi_subset_indicator(
    const Rcpp::NumericMatrix& Q_i,
    double L, double R,
    bool L_open, bool R_open) {
  R_xlen_t I = Q_i.nrow();
  Rcpp::LogicalVector res(I);
  for (R_xlen_t i = 0; i < I; ++i) {
    double a = Q_i(i,0);
    double b = Q_i(i,1);
    res[i] = closed_subset_of(a, b, L, R, L_open, R_open);
  }
  return res;
}

// Everything at once
void calc_all(const ModelData& md, Workspace& ws) { 
  // Define z = z_1, ... , z_I of lenght I 
  // Maybe include in above  
  // Define z = z_I+1, ... , z_I' of length I' - I
  // Define λ = λ_1, ..., λ_N of length N
  
  // For λ we need results form the calculation of z.
  // Thus we first loop over I and then N
  Rcpp::NumericVector z(md.I_mark, 0.0);
  Rcpp::NumericVector lambda(md.N, 0.0);
  
  
  // We need matrices to storage values
  // using armadillo allows for fast col sums
  arma::mat mu_mi(md.M, md.I, arma::fill::zeros);
  arma::mat mu_bar_ji(md.J, md.I_mark, arma::fill::zeros);
  arma::mat eta_ui(md.U, md.I_mark, arma::fill::zeros);
  arma::mat gamma_ci(md.C, md.I_mark, arma::fill::zeros);
  
  
  const int L = std::min(md.M, 
                         std::min(md.U,
                                  std::min(md.C, md.J)));
  
  
  // These loops over I and I_mark we calculate the matrices. 
  // I cannot find a better way than to store these as matrices
  
  // loop of I 
  for(int i = 0; i < md.I; ++i) {
    // main fused range (branch-free body)
    for (int l = 0; l < L; ++l) {
      
      // Not normalized only numerator
      mu_mi(l,i) = md.beta_im(i,l) * ws.z_i[i] * evaluate(ws, md.Q_i(i,1), prev_double(md.R_m_vec[l]));
      
      // Not normalized only numerator
      mu_bar_ji(l,i) = md.alpha_ij(i,l) * ws.z_i[i];
      
      
      eta_ui(l,i) = ws.lambda_n[md.M + l] * 
        evaluate(ws, md.Q_i(l,i), md.t_u[l]) * 
        md.beta_im(i, md.M + l) * 
        ws.z_i[i];
      
      
      gamma_ci(l,i) = md.alpha_ij(i, md.J + l) *
        ws.z_i[i] +
        evaluate(ws, md.Q_i(l,i), md.t_c[l]) * 
        md.beta_im(i, md.W + l) * 
        ws.z_i[i];
    }
    
    // tails
    for (int l = L; l < md.M; ++l) {
      // Not normalized only numerator
      mu_mi(l,i) = md.beta_im(i,l) * ws.z_i[i] * evaluate(ws, md.Q_i(l,i), md.R_m_vec[l]);
    } 
    for (int l = L; l < md.J; ++l) {
      // Not normalized only numerator
      mu_bar_ji(l,i) = md.alpha_ij(i,l) * ws.z_i[i];
    }
    for (int l = L; l < md.U; ++l) {
      eta_ui(l,i) = ws.lambda_n[md.M + l] * 
        evaluate(ws, md.Q_i(l,i), md.t_u[l]) * 
        md.beta_im(i, md.M + l) * 
        ws.z_i[i];
    }
    for (int l = L; l < md.C; ++l) {
      gamma_ci(l,i) = md.alpha_ij(i, md.J + l) *
        ws.z_i[i] +
        evaluate(ws, md.Q_i(l,i), md.t_c[l]) * 
        md.beta_im(i, md.W + l) * 
        ws.z_i[i];
    }
  }
  // Extra for I mark
  for(int i = md.I; i < md.I_mark; ++i) {
    // main fused range (branch-free body)
    for (int l = 0; l < L; ++l) {
      // Not normalized only numerator
      mu_bar_ji(l,i) = md.alpha_ij(i,l) * ws.z_i[i];
      
      eta_ui(l,i) = md.t_u[l] == md.E_star[i-md.I] ? ws.z_i[i] : 0;
  
      gamma_ci(l,i) = md.alpha_ij(i, md.J + l) * ws.z_i[i];
    } 
    
    // tails
    for (int l = L; l < md.J; ++l) {
      // Not normalized only numerator
      mu_bar_ji(l,i) = md.alpha_ij(i,l) * ws.z_i[i];
    }
    for (int l = L; l < md.U; ++l) {
      eta_ui(l,i) = md.t_u[l] == md.E_star[i-md.I] ? ws.z_i[i] : 0;
    }
    for (int l = L; l < md.C; ++l) {
      gamma_ci(l,i) = md.alpha_ij(i, md.J + l) * ws.z_i[i];
    }
  }
  
  // Normalizing by rows (by summing over I) ONLY WORK FOR POSITIVE VALUES
  mu_mi = arma::normalise(mu_mi, 1, 0); // L1-normalize each row
    
  mu_bar_ji = arma::normalise(mu_bar_ji, 1, 0); // L1-normalize each row
  
  eta_ui = arma::normalise(eta_ui, 1, 0); // L1-normalize each row
  
  gamma_ci = arma::normalise(gamma_ci, 1, 0); // L1-normalize each row
  
  arma::rowvec base = arma::sum(mu_bar_ji + eta_ui + gamma_ci, 0);
  arma::rowvec new_z = base;
  
  new_z.head(md.I) += arma::sum(mu_mi, 0);
  new_z.subvec(md.I, md.I_mark - 1) += md.c_k;
  new_z /= md.N_star;

  
  
  // omega_un = I(t_M+u == t_n_star) sum_i=1 ^ I eta_ui
  // t_n_star is in T_star
  double rowsums_omega_n = 0.0;
  double rowsums_sigma_n = 0.0;
  double rowsums_pi_n = 0.0;
  double rowsums_rho_n = 0.0;
  
  // Outer loop over N
  for(int n = 0; n < md.N; ++n) {
    
    rowsums_omega_n = 0.0;
    
    rowsums_sigma_n = 0.0;
    rowsums_pi_n = 0.0;
    rowsums_rho_n = 0.0;
    
    
      lambda[n] = n < md.N1_obs_of_T_star ? md.d_n[n] : 0;
      
      lambda[n] = lambda[n] + 0; 

  }
}

