// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <algorithm>
#include <vector>
#include <cmath>
using namespace Rcpp;

struct ModelData {
  // --- counts (ints) ---
  int J, C, K_tilde, U, N_tilde, M, W, N1_obs_of_T_star, U_pos_obs_of_T_star,
  N, N_star, M_mark, I, K, I_mark, L;
  
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
  
  
  IntegerVector lambda_M_plus_u;
  
  // lambda_M_p_u <- lambda_n[T_star == t_u[u]]
  void calc_lambda_M_plus_u() {
    lambda_M_plus_u = IntegerVector(U, NA_INTEGER);
    for(int u = 0; u < U; ++u){
      // lambda_M_p_u <- lambda_n[T_star == t_u[u]]
      const double val = t_u[u];
      auto it = std::find_if(T_star.begin(), T_star.end(), [&](double x){ return std::abs(x - val) < 1e-9; });
      int index = (it != T_star.end()) ? std::distance(T_star.begin(), it) : -1;
      lambda_M_plus_u[u] = index;
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
    
    calc_lambda_M_plus_u();
    
    L = std::min(M, std::min(U, std::min(C, J)));
  
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
  arma::rowvec lambda_n;  // length N
  arma::rowvec z_i;       // length I_mark
  
  Rcpp::NumericVector t_sorted;   // sorted t_n
  Rcpp::NumericVector log_prefix; // prefix sum of log(1 - lambda_n) over sorted order
  Rcpp::NumericVector zero_prefix;// prefix count of zeros in (1 - lambda_n)
  
  // E-step pieces
  arma::mat mu_mi;
  arma::mat mu_bar_ji;
  arma::mat eta_ui;
  arma::mat gamma_ci;
  
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

void print_summary(const ModelData& md, const Workspace& ws) {
  // --- Debug: sizes of core scalars/lengths ---------------------------------
  Rcpp::Rcout << "[calc_all] sizes:"
              << " I=" << md.I
              << " I_mark=" << md.I_mark
              << " N=" << md.N
              << " M=" << md.M
              << " J=" << md.J
              << " U=" << md.U
              << " C=" << md.C
              << " W=" << md.W
              << " T_star=" << md.T_star.size()
              << " N_star=" << md.N_star
              << " N1_obs_of_T_star=" << md.N1_obs_of_T_star
              << std::endl;
  
  Rcpp::Rcout << "[calc_all] vector lengths:"
              << " z=" << ws.z_i.size()
              << " lambda=" << ws.lambda_n.size()
              << std::endl;
  
  // Input vectors used inside loops
  Rcpp::Rcout << "[calc_all] input vector lengths:"
              << " ws.z_i=" << ws.z_i.size()
              << " ws.lambda_n=" << ws.lambda_n.size()
              << " R_m_vec=" << md.R_m_vec.size()
              << " t_u=" << md.t_u.size()
              << " t_c=" << md.t_c.size()
              << " E_star=" << md.E_star.size()
              << " c_k=" << md.c_k.size()
              << " d_n=" << md.d_n.size()
              << std::endl;
  
  Rcpp::Rcout << "[calc_all] matrix shapes (rows x cols):"
              << " ws.mu_mi=" << ws.mu_mi.n_rows << "x" << ws.mu_mi.n_cols
              << " ws.mu_bar_ji=" << ws.mu_bar_ji.n_rows << "x" << ws.mu_bar_ji.n_cols
              << " ws.eta_ui=" << ws.eta_ui.n_rows << "x" << ws.eta_ui.n_cols
              << " ws.gamma_ci=" << ws.gamma_ci.n_rows << "x" << ws.gamma_ci.n_cols
              << " beta_im=" << md.beta_im.nrow()  << "x" << md.beta_im.ncol()
              << " alpha_ij=" << md.alpha_ij.nrow() << "x" << md.alpha_ij.ncol()
              << " Q_i=" << md.Q_i.nrow() << "x" << md.Q_i.ncol()
              << std::endl;
  
  Rcpp::Rcout << "[calc_all] L (min of M,U,C,J)=" << md.L << std::endl;
  
  
  // Loop ranges summary (avoid printing at every iteration)
  Rcpp::Rcout << "[calc_all] loop ranges:"
              << " i_main: 0.." << (md.I > 0 ? md.I - 1 : -1)
              << " i_extra(I_mark): " << md.I << ".." << (md.I_mark > 0 ? md.I_mark - 1 : -1)
              << " l_main: 0.." << (md.L > 0 ? md.L - 1 : -1)
              << " tails: M:" << md.L << ".." << (md.M > 0 ? md.M - 1 : -1)
              << " J:" << md.L << ".." << (md.J > 0 ? md.J - 1 : -1)
              << " U:" << md.L << ".." << (md.U > 0 ? md.U - 1 : -1)
              << " C:" << md.L << ".." << (md.C > 0 ? md.C - 1 : -1)
              << std::endl;
  
  Rcpp::Rcout << "[calc_all] Final z_i" << ws.z_i << std::endl;
}


void run_em_once(const ModelData& md, Workspace& ws) { 
  // --- Resets temps ----------------------------------------------------------
  ws.mu_mi.zeros(md.M, md.I);
  ws.mu_bar_ji.zeros(md.J, md.I_mark);
  ws.eta_ui.zeros(md.U, md.I_mark);
  ws.gamma_ci.zeros(md.C, md.I_mark);

  // --- Loops over I ----------------------------------------------------------
  setup_prod(md,ws);
  
  // --- Loops over I ----------------------------------------------------------
  for(int i = 0; i < md.I; ++i) {
    for (int l = 0; l < md.L; ++l) {
      // Not normalized only numerator
      ws.mu_mi(l,i) = md.beta_im(i,l) * 
        ws.z_i[i] * 
        evaluate(ws, md.Q_i(i,1), next_double(md.R_m_vec[l]));
      
      ws.mu_bar_ji(l,i) = md.alpha_ij(i,l) * ws.z_i[i];

      ws.eta_ui(l,i) = ws.lambda_n[md.lambda_M_plus_u[l]] *
        evaluate(ws, md.Q_i(i,1), md.t_u[l]) *
        md.beta_im(i, md.M + l) *
        ws.z_i[i];

      ws.gamma_ci(l,i) = md.alpha_ij(i, md.J + l) *
        ws.z_i[i] +
        evaluate(ws, md.Q_i(i,1), next_double(md.t_c[l])) *
        md.beta_im(i, md.W + l) *
        ws.z_i[i];
    }

    // tails
    for (int l = md.L; l < md.M; ++l) {
      ws.mu_mi(l,i) = md.beta_im(i,l) * 
        ws.z_i[i] * 
        evaluate(ws, md.Q_i(i,1), next_double(md.R_m_vec[l]));
    }
    for (int l = md.L; l < md.J; ++l) {
      ws.mu_bar_ji(l,i) = md.alpha_ij(i,l) * ws.z_i[i];
    }
    for (int l = md.L; l < md.U; ++l) {
      ws.eta_ui(l,i) = ws.lambda_n[md.lambda_M_plus_u[l]] *
        evaluate(ws, md.Q_i(i,1), md.t_u[l]) *
        md.beta_im(i, md.M + l) *
        ws.z_i[i];
    }
    for (int l = md.L; l < md.C; ++l) {
      ws.gamma_ci(l,i) = md.alpha_ij(i, md.J + l) *
        ws.z_i[i] +
        evaluate(ws, md.Q_i(i,1), next_double(md.t_c[l])) *
        md.beta_im(i, md.W + l) *
        ws.z_i[i];
    }
  }
  
  // --- Extra for I_mark ------------------------------------------------------
  for(int i = md.I; i < md.I_mark; ++i) {
    for (int l = 0; l < md.L; ++l) {
      ws.mu_bar_ji(l,i) = md.alpha_ij(i,l) * ws.z_i[i];
      ws.eta_ui(l,i) = md.t_u[l] == md.E_star[i-md.I] ? ws.z_i[i] : 0;
      ws.gamma_ci(l,i) = md.alpha_ij(i, md.J + l) * ws.z_i[i];
    } 
    for (int l = md.L; l < md.J; ++l) {
      ws.mu_bar_ji(l,i) = md.alpha_ij(i,l) * ws.z_i[i];
    }
    for (int l = md.L; l < md.U; ++l) {
      ws.eta_ui(l,i) = md.t_u[l] == md.E_star[i-md.I] ? ws.z_i[i] : 0;
    }
    for (int l = md.L; l < md.C; ++l) {
      ws.gamma_ci(l,i) = md.alpha_ij(i, md.J + l) * ws.z_i[i];
    }
  }
  
  // --- Normalizes rows -------------------------------------------------------
  ws.mu_mi = arma::normalise(ws.mu_mi, 1, 1);
  ws.mu_bar_ji = arma::normalise(ws.mu_bar_ji, 1, 1);
  ws.eta_ui = arma::normalise(ws.eta_ui, 1, 1);
  ws.gamma_ci = arma::normalise(ws.gamma_ci, 1, 1);
  
  // --- Loops over N ----------------------------------------------------------
  for(int n = 0; n < md.N; ++n) {
    ws.lambda_n[n] = n < md.N1_obs_of_T_star ? md.d_n[n] : 0;
    ws.lambda_n[n] = ws.lambda_n[n] + 0;
  }
  
  
  
  
  
  // --- Calculates new z ------------------------------------------------------
  arma::rowvec base = arma::sum(ws.mu_bar_ji, 0);
  base += arma::sum(ws.eta_ui, 0);
  base += arma::sum(ws.gamma_ci, 0);
  arma::rowvec new_z = base;
  
  new_z.head(md.I) += arma::sum(ws.mu_mi, 0);
  new_z.subvec(md.I, md.I_mark - 1) += md.c_k;
  new_z /= md.N_star;
  
  ws.z_i = new_z;
}


// [[Rcpp::export]]
Rcpp::List em_fit(SEXP md_ptr,
                  Rcpp::Nullable<Rcpp::NumericVector> z_init = R_NilValue,
                  Rcpp::Nullable<Rcpp::NumericVector> lambda_init = R_NilValue,
                  int max_iter = 1) {
  
  Rcpp::XPtr<ModelData> p(md_ptr);
  const ModelData& md = *p;
  
  const arma::uword N       = static_cast<arma::uword>(md.T_star.size()); // or md.N
  const arma::uword I_mark  = static_cast<arma::uword>(md.I_mark);     // adapt if different
  
  Workspace ws;
  

  ws.z_i.set_size(I_mark);
  ws.z_i.fill(1.0 / static_cast<double>(I_mark));
  
  Rcpp::Rcout << "[em_fit] Initial ws.z_i" << ws.z_i << std::endl;

  ws.lambda_n.set_size(N);
  ws.lambda_n.fill(0.5);

  run_em_once(md, ws);
  
  print_summary(md, ws);
  
  return Rcpp::List::create(
    Rcpp::_["z_i"]      = ws.z_i,
    Rcpp::_["lambda_n"] = ws.lambda_n,
    Rcpp::_["alpha_ij"] = md.alpha_ij,
    Rcpp::_["beta_im"]  = md.beta_im,
    Rcpp::_["mu_mi"]    = ws.mu_mi,
    Rcpp::_["ws.mu_bar_ji"]= ws.mu_bar_ji,
    Rcpp::_["ws.eta_ui"]   = ws.eta_ui,
    Rcpp::_["ws.gamma_ci"] = ws.gamma_ci,
    Rcpp::_["rho_mn"]   = ws.rho_mn,
    Rcpp::_["pi_un"]    = ws.pi_un,
    Rcpp::_["sigma_cn"] = ws.sigma_cn
  );
}
