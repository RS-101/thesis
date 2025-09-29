#include <Rcpp.h>
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
  IntegerVector c_k, d_n;
  
  // --- 2-col interval matrices ---
  NumericMatrix A_m, A_u, A_c, full_A_m, Q_i;
  
  // ctor from named R list
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
    c_k = as<IntegerVector>(x["c_k"]);
    d_n = as<IntegerVector>(x["d_n"]);
    
    // matrices (2 columns) — the important part
    A_m      = as<NumericMatrix>(x["A_m"]);
    A_u      = as<NumericMatrix>(x["A_u"]);
    A_c      = as<NumericMatrix>(x["A_c"]);
    full_A_m = as<NumericMatrix>(x["full_A_m"]);
    Q_i      = as<NumericMatrix>(x["Q_i"]);
    
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

static inline double mean_num(const Rcpp::NumericVector& v) {
  if (v.size() == 0) return NA_REAL;
  return Rcpp::mean(Rcpp::na_omit(v));
}
static inline double mean_int(const Rcpp::IntegerVector& v) {
  if (v.size() == 0) return NA_REAL;
  Rcpp::NumericVector nv = Rcpp::as<Rcpp::NumericVector>(v);
  return Rcpp::mean(Rcpp::na_omit(nv));
}

// [[Rcpp::export]]
Rcpp::List model_data_summary(SEXP xp) {
  Rcpp::XPtr<ModelData> p(xp);
  return Rcpp::List::create(
    // counts
    Rcpp::_["J"]=p->J, Rcpp::_["C"]=p->C, Rcpp::_["K_tilde"]=p->K_tilde, Rcpp::_["U"]=p->U,
    Rcpp::_["N_tilde"]=p->N_tilde, Rcpp::_["M"]=p->M, Rcpp::_["W"]=p->W, Rcpp::_["N"]=p->N,
    Rcpp::_["N_star"]=p->N_star, Rcpp::_["M_mark"]=p->M_mark, Rcpp::_["I"]=p->I, Rcpp::_["K"]=p->K,
            Rcpp::_["I_mark"]=p->I_mark,
            
            // sizes
            Rcpp::_["len(s_j)"]=p->s_j.size(), Rcpp::_["len(E_star)"]=p->E_star.size(),
            Rcpp::_["nrow(Q_i)"]=p->Q_i.nrow(),
            
            // means (NumericVector)
            Rcpp::_["mean_s_j"]        = mean_num(p->s_j),
            Rcpp::_["mean_L_c"]        = mean_num(p->L_c),
            Rcpp::_["mean_t_c"]        = mean_num(p->t_c),
            Rcpp::_["mean_e_k"]        = mean_num(p->e_k),
            Rcpp::_["mean_L_u"]        = mean_num(p->L_u),
            Rcpp::_["mean_t_u"]        = mean_num(p->t_u),
            Rcpp::_["mean_t_m_in_N_tilde"] = mean_num(p->t_m_in_N_tilde),
            Rcpp::_["mean_L_m"]        = mean_num(p->L_m),
            Rcpp::_["mean_R_m"]        = mean_num(p->R_m_vec),
            Rcpp::_["mean_t_m"]        = mean_num(p->t_m),
            Rcpp::_["mean_E_star"]     = mean_num(p->E_star),
            Rcpp::_["mean_T_star"]     = mean_num(p->T_star),
            Rcpp::_["mean_L_bar"]      = mean_num(p->L_bar),
            Rcpp::_["mean_R_bar"]      = mean_num(p->R_bar),   // will be Inf if R_bar includes Inf
            Rcpp::_["mean_s_j_c"]      = mean_num(p->s_j_c),
            Rcpp::_["mean_s_j_full"]   = mean_num(p->s_j_full),
            Rcpp::_["mean_Q"]          = mean_num(p->Q),
            Rcpp::_["mean_Q_i_mark"]   = mean_num(p->Q_i_mark),
            Rcpp::_["mean_A_union"]    = mean_num(p->A_union),
            
            // means (IntegerVector)
            Rcpp::_["mean_c_k"]        = mean_int(p->c_k),
            Rcpp::_["mean_d_n"]        = mean_int(p->d_n)
  );
}

// ---------- workspace + helpers ----------

struct Workspace {
  // parameters updated by EM
  Rcpp::NumericVector lambda_n;  // length N
  Rcpp::NumericVector z_i;       // length I_mark
  
  // E-step pieces
  Rcpp::NumericMatrix alpha_ij;  // I_mark x (J+C)
  Rcpp::NumericMatrix beta_im;   // I x M'
  Rcpp::NumericMatrix mu_mi;     // M x I
  Rcpp::NumericMatrix mu_bar_ji; // J x I_mark
  Rcpp::NumericMatrix eta_ui;    // U x I_mark
  Rcpp::NumericMatrix gamma_ci;  // C x I_mark
  
  // aggregation for lambda
  Rcpp::NumericMatrix rho_mn;    // M x N
  Rcpp::NumericMatrix pi_un;     // U x N
  Rcpp::NumericMatrix sigma_cn;  // C x N
};

// --- small utilities ---
static inline double safe_mean(const Rcpp::NumericVector& x) {
  if (x.size() == 0) return NA_REAL;
  return Rcpp::mean(Rcpp::na_omit(x));
}

static inline void na0_inplace(Rcpp::NumericMatrix m) {
  for (R_xlen_t i = 0; i < m.nrow(); ++i)
    for (R_xlen_t j = 0; j < m.ncol(); ++j)
      if (!R_finite(m(i,j))) m(i,j) = 0.0;
}

static inline Rcpp::NumericVector colSums_num(const Rcpp::NumericMatrix& m) {
  Rcpp::NumericVector res(m.ncol());
  for (R_xlen_t j = 0; j < m.ncol(); ++j) {
    double s = 0.0;
    for (R_xlen_t i = 0; i < m.nrow(); ++i) {
      double v = m(i,j);
      if (R_finite(v)) s += v;
    }
    res[j] = s;
  }
  return res;
}

static inline Rcpp::NumericVector rowSums_num(const Rcpp::NumericMatrix& m) {
  Rcpp::NumericVector res(m.nrow());
  for (R_xlen_t i = 0; i < m.nrow(); ++i) {
    double s = 0.0;
    for (R_xlen_t j = 0; j < m.ncol(); ++j) {
      double v = m(i,j);
      if (R_finite(v)) s += v;
    }
    res[i] = s;
  }
  return res;
}

// membership of a point x in (L,R) with open/closed flags
static inline bool in_interval(double x, double L, double R, bool L_open, bool R_open) {
  bool left  = L_open ? (x >  L) : (x >= L);
  bool right = R_open ? (x <  R) : (x <= R);
  return left && right;
}

// [a,b] (closed) subset-of (L,R) (with flags on the outer interval)
static inline bool closed_subset_of(double a, double b, double L, double R, bool L_open, bool R_open) {
  bool left  = L_open ? (a >  L) : (a >= L);
  bool right = R_open ? (b <  R) : (b <= R);
  return left && right;
}

// product over t* in interval (L,R) with flags
static inline double product_over_t_stars_interval(double L, double R, bool L_open, bool R_open,
                                                   const Rcpp::NumericVector& T_star,
                                                   const Rcpp::NumericVector& lambda_n) {
  double prod = 1.0;
  for (R_xlen_t k = 0; k < T_star.size(); ++k) {
    double t = T_star[k];
    if (in_interval(t, L, R, L_open, R_open)) {
      prod *= (1.0 - lambda_n[k]);
    }
  }
  return prod;
}

// vector<bool> : is each Q_i row a subset of (L,R) with flags
static inline Rcpp::LogicalVector Qi_subset_indicator(const Rcpp::NumericMatrix& Q_i,
                                                      double L, double R, bool L_open, bool R_open) {
  R_xlen_t I = Q_i.nrow();
  Rcpp::LogicalVector res(I);
  for (R_xlen_t i = 0; i < I; ++i) {
    double a = Q_i(i,0);
    double b = Q_i(i,1);
    res[i] = closed_subset_of(a, b, L, R, L_open, R_open);
  }
  return res;
}

// ---------- (18) alpha: I' x (J+C) ----------
void cal_alpha(const ModelData& md, Workspace& ws) {
  const Rcpp::NumericMatrix& Q_i = md.Q_i;
  const Rcpp::NumericVector& Q_i_mark = md.Q_i_mark;
  const Rcpp::NumericVector& s_j_full = md.s_j_full;
  
  int I = Q_i.nrow();
  int I_mark = md.I_mark; // I + K
  int JpC = s_j_full.size();
  
  ws.alpha_ij = Rcpp::NumericMatrix(I_mark, JpC);
  for (int i = 0; i < I_mark; ++i) {
    double left = (i < I) ? Q_i(i,0) : Q_i_mark[i - I];
    for (int j = 0; j < JpC; ++j)
      ws.alpha_ij(i,j) = (s_j_full[j] <= left) ? 1.0 : 0.0;
  }
}

// ---------- (19) beta: I x M' ----------
void cal_beta(const ModelData& md, Workspace& ws) {
  const Rcpp::NumericMatrix& Q_i = md.Q_i;
  const Rcpp::NumericMatrix& full_A_m = md.full_A_m;
  
  int I = Q_i.nrow();
  int M_mark = full_A_m.nrow();
  ws.beta_im = Rcpp::NumericMatrix(I, M_mark);
  
  for (int i = 0; i < I; ++i) {
    double Li = Q_i(i,0), Ri = Q_i(i,1);
    for (int m = 0; m < M_mark; ++m) {
      double Lm = full_A_m(m,0), Rm = full_A_m(m,1);
      ws.beta_im(i,m) = (Lm <= Li && Ri <= Rm) ? 1.0 : 0.0;
    }
  }
}

// ---------- μ_{mi} : M x I ----------
void cal_mu_MI(const ModelData& md, Workspace& ws) {
  const int I = md.Q_i.nrow();
  const int M = md.A_m.nrow();
  ws.mu_mi = Rcpp::NumericMatrix(M, I);
  
  // r_i are Q_i[,2], R_m are A_m[,2]
  Rcpp::NumericVector r_i(I), Rm(M);
  for (int i = 0; i < I; ++i) r_i[i] = md.Q_i(i,1);
  for (int m = 0; m < M; ++m) Rm[m] = md.A_m(m,1);
  
  // beta_im restricted to first M columns
  if (ws.beta_im.nrow() != I || ws.beta_im.ncol() < M) Rcpp::stop("beta_im not ready");
  
  for (int m = 0; m < M; ++m) {
    // precompute products for each i over (r_i, Rm[m]] with left open, right closed
    Rcpp::NumericVector prod_res(I);
    for (int i = 0; i < I; ++i)
      prod_res[i] = product_over_t_stars_interval(r_i[i], Rm[m], true, false, md.T_star, ws.lambda_n);
    
    // denominator
    double denom = 0.0;
    for (int i = 0; i < I; ++i)
      denom += ws.beta_im(i,m) * ws.z_i[i] * prod_res[i];
    
    for (int i = 0; i < I; ++i) {
      double num = ws.beta_im(i,m) * ws.z_i[i] * prod_res[i];
      ws.mu_mi(m,i) = (num != 0.0 && denom != 0.0) ? (num / denom) : 0.0;
    }
  }
  na0_inplace(ws.mu_mi);
}

// ---------- \bar{μ}_{ji} : J x I' ----------
void cal_mu_bar_JI_mark(const ModelData& md, Workspace& ws) {
  int I_mark = md.I_mark;
  int J = md.J;
  
  if (ws.alpha_ij.nrow() != I_mark || ws.alpha_ij.ncol() < J)
    Rcpp::stop("alpha_ij not ready");
  
  ws.mu_bar_ji = Rcpp::NumericMatrix(J, I_mark);
  // prod = alpha[ ,1:J] * z_i (row-wise multiply by z_i), then normalize by colSums
  Rcpp::NumericMatrix prod(I_mark, J);
  for (int i = 0; i < I_mark; ++i)
    for (int j = 0; j < J; ++j)
      prod(i,j) = ws.alpha_ij(i,j) * ws.z_i[i];
  
  // col sums over i (i.e., sums for each j)
  Rcpp::NumericVector cs(J);
  for (int j = 0; j < J; ++j) {
    double s = 0.0;
    for (int i = 0; i < I_mark; ++i) s += prod(i,j);
    cs[j] = s;
  }
  
  // transpose + sweep: result size J x I_mark
  for (int j = 0; j < J; ++j) {
    double d = cs[j];
    for (int i = 0; i < I_mark; ++i)
      ws.mu_bar_ji(j,i) = (d != 0.0) ? (prod(i,j) / d) : 0.0;
  }
  na0_inplace(ws.mu_bar_ji);
}

// ---------- η_{ui} : U x I' ----------
void cal_eta_UI_mark(const ModelData& md, Workspace& ws) {
  int U = md.A_u.nrow();
  int I = md.Q_i.nrow();
  int I_mark = md.I_mark;
  int M = md.M;
  
  ws.eta_ui = Rcpp::NumericMatrix(U, I_mark);
  Rcpp::NumericVector r_i(I);
  for (int i = 0; i < I; ++i) r_i[i] = md.Q_i(i,1);
  
  for (int u = 0; u < U; ++u) {
    double tu = md.t_u[u];
    
    // λ at t_u[u]
    double lambda_Mu = 0.0;
    // find index in T_star exactly equal to t_u[u]
    for (int n = 0; n < md.T_star.size(); ++n)
      if (md.T_star[n] == tu) { lambda_Mu = ws.lambda_n[n]; break; }
      
      // prod over (r_i, t_u[u]) with both open
      Rcpp::NumericVector prod_res(I);
      for (int i = 0; i < I; ++i)
        prod_res[i] = product_over_t_stars_interval(r_i[i], tu, true, true, md.T_star, ws.lambda_n);
      
      // denom
      double denom_head = 0.0;
      for (int i = 0; i < I; ++i)
        denom_head += prod_res[i] * ws.beta_im(i, M + u) * ws.z_i[i];
      denom_head *= lambda_Mu;
      
      // tail term: sum over i>I where E_star == t_u[u]
      double denom_tail = 0.0;
      for (int k = 0; k < md.E_star.size(); ++k)
        if (md.E_star[k] == tu) denom_tail += ws.z_i[I + k];
        
        double denom = denom_head + denom_tail;
        if (denom == 0.0) {
          // all zeros
          for (int i = 0; i < I_mark; ++i) ws.eta_ui(u,i) = 0.0;
          continue;
        }
        
        for (int i = 0; i < I_mark; ++i) {
          if (i < I) {
            ws.eta_ui(u,i) = lambda_Mu * prod_res[i] * ws.beta_im(i, M + u) * ws.z_i[i] / denom;
          } else {
            int k = i - I; // tail index into E_star
            ws.eta_ui(u,i) = (md.E_star[k] == tu ? ws.z_i[i] / denom : 0.0);
          }
        }
  }
  na0_inplace(ws.eta_ui);
}

// ---------- γ_{ci} : C x I' ----------
void cal_gamma_CI_mark(const ModelData& md, Workspace& ws) {
  int I = md.Q_i.nrow();
  int I_mark = md.I_mark;
  int C = md.A_c.nrow();
  int W = md.W;
  int J = md.J;
  
  ws.gamma_ci = Rcpp::NumericMatrix(C, I_mark);
  if (C == 0) return;
  
  Rcpp::NumericVector r_i(I);
  for (int i = 0; i < I; ++i) r_i[i] = md.Q_i(i,1);
  
  for (int c = 0; c < C; ++c) {
    double tc = md.t_c[c];
    
    // prod over (r_i, t_c[c]] with left open, right closed
    Rcpp::NumericVector prod_res(I);
    for (int i = 0; i < I; ++i)
      prod_res[i] = product_over_t_stars_interval(r_i[i], tc, true, false, md.T_star, ws.lambda_n);
    
    // denom: head (I part) + alpha tail column (J+c)
    double denom = 0.0;
    for (int i = 0; i < I; ++i)
      denom += prod_res[i] * ws.beta_im(i, W + c) * ws.z_i[i];
    
    // α_{i,(J+c)} z_i for all i up to I_mark
    for (int i = 0; i < I_mark; ++i)
      denom += ws.alpha_ij(i, J + c) * ws.z_i[i];
    
    if (denom == 0.0) {
      for (int i = 0; i < I_mark; ++i) ws.gamma_ci(c,i) = 0.0;
      continue;
    }
    
    for (int i = 0; i < I_mark; ++i) {
      if (i < I) {
        double term1 = ws.alpha_ij(i, J + c) * ws.z_i[i] / denom;
        double term2 = (prod_res[i] * ws.beta_im(i, W + c) * ws.z_i[i]) / denom;
        ws.gamma_ci(c,i) = term1 + term2;
      } else {
        ws.gamma_ci(c,i) = ws.alpha_ij(i, J + c) * ws.z_i[i] / denom;
      }
    }
  }
  na0_inplace(ws.gamma_ci);
}

// ---------- ρ_{mn} : M x N ----------
void cal_rho_MN(const ModelData& md, Workspace& ws) {
  int M = md.A_m.nrow();
  int N = md.T_star.size();
  int I = md.Q_i.nrow();
  ws.rho_mn = Rcpp::NumericMatrix(M, N);
  
  for (int n = 0; n < N; ++n) {
    double tn = md.T_star[n];
    for (int m = 0; m < M; ++m) {
      if (md.t_m[m] >= tn) {
        // interval [A_m[m,1], tn) : left closed, right open
        Rcpp::LogicalVector ind = Qi_subset_indicator(md.Q_i, md.A_m(m,0), tn, false, true);
        double s = 0.0;
        for (int i = 0; i < I; ++i) if (ind[i]) s += ws.mu_mi(m,i);
        ws.rho_mn(m,n) = s;
      } else {
        ws.rho_mn(m,n) = 0.0;
      }
    }
  }
  na0_inplace(ws.rho_mn);
}

// ---------- π_{un} : U x N ----------
void cal_pi_UN(const ModelData& md, Workspace& ws) {
  int U = md.A_u.nrow();
  int N = md.T_star.size();
  int I = md.Q_i.nrow();
  ws.pi_un = Rcpp::NumericMatrix(U, N);
  
  for (int n = 0; n < N; ++n) {
    double tn = md.T_star[n];
    for (int u = 0; u < U; ++u) {
      if (md.t_u[u] >= tn) {
        Rcpp::LogicalVector ind = Qi_subset_indicator(md.Q_i, md.A_u(u,0), tn, false, true);
        double s = 0.0;
        for (int i = 0; i < I; ++i) if (ind[i]) s += ws.eta_ui(u,i);
        ws.pi_un(u,n) = s;
      } else {
        ws.pi_un(u,n) = 0.0;
      }
    }
  }
  na0_inplace(ws.pi_un);
}

// ---------- σ_{cn} : C x N ----------
void cal_sigma_CN(const ModelData& md, Workspace& ws) {
  int C = md.A_c.nrow();
  int N = md.T_star.size();
  int I = md.Q_i.nrow();
  ws.sigma_cn = Rcpp::NumericMatrix(C, N);
  if (C == 0) return;
  
  for (int n = 0; n < N; ++n) {
    double tn = md.T_star[n];
    for (int c = 0; c < C; ++c) {
      if (md.t_c[c] >= tn) {
        Rcpp::LogicalVector ind = Qi_subset_indicator(md.Q_i, md.A_c(c,0), tn, false, true);
        double s = 0.0;
        for (int i = 0; i < I; ++i) if (ind[i]) s += ws.gamma_ci(c,i);
        ws.sigma_cn(c,n) = s;
      } else {
        ws.sigma_cn(c,n) = 0.0;
      }
    }
  }
  na0_inplace(ws.sigma_cn);
}

// ---------- (23): z (head, i ≤ I) ----------
void e_23(const ModelData& md, Workspace& ws) {
  // z_head over i=1..I
  int I = md.Q_i.nrow();
  if (ws.mu_mi.ncol() != I) Rcpp::stop("mu_mi not ready");
  if (ws.mu_bar_ji.ncol() != md.I_mark) Rcpp::stop("mu_bar_ji not ready");
  if (ws.eta_ui.ncol() != md.I_mark) Rcpp::stop("eta_ui not ready");
  if (ws.gamma_ci.ncol() != md.I_mark && md.A_c.nrow() > 0) Rcpp::stop("gamma_ci not ready");
  
  Rcpp::NumericVector z_head(I, 0.0);
  
  // colSums over I for each matrix, then take first I entries
  Rcpp::NumericVector s1 = colSums_num(ws.mu_mi);
  Rcpp::NumericVector s2 = colSums_num(ws.mu_bar_ji);
  Rcpp::NumericVector s3 = colSums_num(ws.eta_ui);
  Rcpp::NumericVector s4 = (md.A_c.nrow() > 0) ? colSums_num(ws.gamma_ci) : Rcpp::NumericVector(md.I_mark, 0.0);
  
  for (int i = 0; i < I; ++i) {
    double num = s1[i] + s2[i] + s3[i] + s4[i];
    ws.z_i[i] = num / md.N_star;
  }
}

// ---------- (24): z (tail, i > I) ----------
void e_24(const ModelData& md, Workspace& ws) {
  int I = md.Q_i.nrow();
  int I_mark = md.I_mark;
  int tail = I_mark - I;
  if (tail <= 0) return;
  
  // slice columns I..I_mark-1
  Rcpp::NumericMatrix mu_bar_tail(md.J, tail);
  Rcpp::NumericMatrix eta_tail(ws.eta_ui.nrow(), tail);
  Rcpp::NumericMatrix gamma_tail(ws.gamma_ci.nrow(), tail);
  
  for (int j = 0; j < md.J; ++j)
    for (int k = 0; k < tail; ++k)
      mu_bar_tail(j,k) = ws.mu_bar_ji(j, I + k);
  
  for (int u = 0; u < ws.eta_ui.nrow(); ++u)
    for (int k = 0; k < tail; ++k)
      eta_tail(u,k) = ws.eta_ui(u, I + k);
  
  for (int c = 0; c < ws.gamma_ci.nrow(); ++c)
    for (int k = 0; k < tail; ++k)
      gamma_tail(c,k) = ws.gamma_ci(c, I + k);
  
  Rcpp::NumericVector cks = Rcpp::clone(Rcpp::as<Rcpp::NumericVector>(md.c_k));// sums
  if (cks.size() != tail) Rcpp::stop("c_k length must equal I_mark - I");
  Rcpp::NumericVector s2 = colSums_num(mu_bar_tail);
  Rcpp::NumericVector s3 = colSums_num(eta_tail);
  Rcpp::NumericVector s4 = colSums_num(gamma_tail);
  
  for (int k = 0; k < tail; ++k) {
    double num = cks[k] + s2[k] + s3[k] + s4[k];
    ws.z_i[I + k] = num / md.N_star;
  }
}

// ---------- (25): λ_n ----------
void e_25(const ModelData& md, Workspace& ws) {
  int N = md.T_star.size();
  if (ws.rho_mn.ncol() != N || ws.pi_un.ncol() != N || ws.sigma_cn.ncol() != N)
    Rcpp::stop("rho/pi/sigma not ready");
  
  // denom
  Rcpp::NumericVector denom(N, 0.0);
  Rcpp::NumericVector add;
  
  add = colSums_num(ws.rho_mn);
  for (int n = 0; n < N; ++n) denom[n] += add[n];
  add = colSums_num(ws.pi_un);
  for (int n = 0; n < N; ++n) denom[n] += add[n];
  add = colSums_num(ws.sigma_cn);
  for (int n = 0; n < N; ++n) denom[n] += add[n];
  
  // Σ_u 1(t_{M+u}=t*_n) Σ_i η_{ui}
  Rcpp::NumericVector sum_eta_u = rowSums_num(ws.eta_ui); // length U
  Rcpp::NumericVector eta_term(N, 0.0);
  for (int u = 0; u < md.t_u.size(); ++u) {
    double tu = md.t_u[u];
    // find index n where T_star[n] == tu
    for (int n = 0; n < N; ++n)
      if (md.T_star[n] == tu) eta_term[n] += sum_eta_u[u];
  }
  
  // I(n <= N1) d_n
  Rcpp::NumericVector num(N, 0.0);
  for (int n = 0; n < N; ++n) {
    double ind = (n < md.N1_obs_of_T_star) ? 1.0 : 0.0; // n is 0-based here
    num[n] = ind * md.d_n[n] + eta_term[n];
    ws.lambda_n[n] = (denom[n] != 0.0) ? (num[n] / denom[n]) : 0.0;
  }
}

// ---------- EM orchestrator (void) ----------
void em_estimate(const ModelData& md, Workspace& ws, int max_iter, double tol, bool verbose) {
  const int I      = md.Q_i.nrow();
  const int I_mark = md.I_mark;
  const int N      = md.T_star.size();
  
  // initialize if empty
  if (ws.z_i.size() != I_mark) {
    ws.z_i = Rcpp::NumericVector(I_mark, 1.0 / I_mark);
  }
  if (ws.lambda_n.size() != N) {
    ws.lambda_n = Rcpp::runif(N); // simple init, user can overwrite before call
  }
  
  for (int iter = 0; iter < max_iter; ++iter) {
    Rcpp::NumericVector z_prev = Rcpp::clone(ws.z_i);
    Rcpp::NumericVector l_prev = Rcpp::clone(ws.lambda_n);
    
    // E-step pieces
    cal_alpha(md, ws);
    cal_beta(md, ws);
    
    cal_mu_MI(md, ws);
    cal_mu_bar_JI_mark(md, ws);
    cal_eta_UI_mark(md, ws);
    cal_gamma_CI_mark(md, ws);
    
    cal_rho_MN(md, ws);
    cal_pi_UN(md, ws);
    cal_sigma_CN(md, ws);
    
    // M-step
    e_23(md, ws);
    e_24(md, ws);
    e_25(md, ws);
    
    // convergence
    double dz = 0.0, dl = 0.0;
    for (int i = 0; i < I_mark; ++i) dz = std::max(dz, std::abs(ws.z_i[i] - z_prev[i]));
    for (int n = 0; n < N; ++n)      dl = std::max(dl, std::abs(ws.lambda_n[n] - l_prev[n]));
    
    if (verbose) Rcpp::Rcout << "iter " << (iter+1) << ": max|Δz|=" << dz << " max|Δλ|=" << dl << "\n";
    if (R_finite(dz) && R_finite(dl) && std::max(dz, dl) < tol) break;
  }
}

// ---------- R interface: run EM and return workspace ----------
/*
 em_fit(
 md_ptr,             // externalptr from make_model_data()
 z_init = NULL,      // optional numeric
 lambda_init = NULL, // optional numeric
 max_iter = 200,
 tol = 1e-8,
 verbose = FALSE
 )
 Returns a list with z_i, lambda_n and all intermediate matrices.
 */
// [[Rcpp::export]]
Rcpp::List em_fit(SEXP md_ptr,
                  Rcpp::Nullable<Rcpp::NumericVector> z_init = R_NilValue,
                  Rcpp::Nullable<Rcpp::NumericVector> lambda_init = R_NilValue,
                  int max_iter = 200,
                  double tol = 1e-8,
                  bool verbose = false) {
  Rcpp::XPtr<ModelData> p(md_ptr);
  const ModelData& md = *p;
  
  Workspace ws;
  
  // initialize if provided
  if (z_init.isNotNull())  ws.z_i = Rcpp::NumericVector(z_init.get());
  if (lambda_init.isNotNull()) ws.lambda_n = Rcpp::NumericVector(lambda_init.get());
  
  em_estimate(md, ws, max_iter, tol, verbose);
  
  return Rcpp::List::create(
    Rcpp::_["z_i"]      = ws.z_i,
    Rcpp::_["lambda_n"] = ws.lambda_n,
    Rcpp::_["alpha_ij"] = ws.alpha_ij,
    Rcpp::_["beta_im"]  = ws.beta_im,
    Rcpp::_["mu_mi"]    = ws.mu_mi,
    Rcpp::_["mu_bar_ji"]= ws.mu_bar_ji,
    Rcpp::_["eta_ui"]   = ws.eta_ui,
    Rcpp::_["gamma_ci"] = ws.gamma_ci,
    Rcpp::_["rho_mn"]   = ws.rho_mn,
    Rcpp::_["pi_un"]    = ws.pi_un,
    Rcpp::_["sigma_cn"] = ws.sigma_cn
  );
}


// [[Rcpp::export]]
Rcpp::List model_data_to_list(SEXP xp) {
  Rcpp::XPtr<ModelData> p(xp);
  const ModelData& md = *p;
  
  return Rcpp::List::create(
    // ints
    Rcpp::_["J"] = md.J,
    Rcpp::_["C"] = md.C,
    Rcpp::_["K_tilde"] = md.K_tilde,
    Rcpp::_["U"] = md.U,
    Rcpp::_["N_tilde"] = md.N_tilde,
    Rcpp::_["M"] = md.M,
    Rcpp::_["W"] = md.W,
    Rcpp::_["N1_obs_of_T_star"] = md.N1_obs_of_T_star,
    Rcpp::_["U_pos_obs_of_T_star"] = md.U_pos_obs_of_T_star,
    Rcpp::_["N"] = md.N,
    Rcpp::_["N_star"] = md.N_star,
    Rcpp::_["M_mark"] = md.M_mark,
    Rcpp::_["I"] = md.I,
    Rcpp::_["K"] = md.K,
    Rcpp::_["I_mark"] = md.I_mark,
    
    // scalars
    Rcpp::_["s_max"] = md.s_max,
    Rcpp::_["R_max"] = md.R_max,
    Rcpp::_["e_star_max"] = md.e_star_max,
    
    // vectors
    Rcpp::_["s_j"]      = Rcpp::clone(md.s_j),
    Rcpp::_["L_c"]      = Rcpp::clone(md.L_c),
    Rcpp::_["t_c"]      = Rcpp::clone(md.t_c),
    Rcpp::_["e_k"]      = Rcpp::clone(md.e_k),
    Rcpp::_["L_u"]      = Rcpp::clone(md.L_u),
    Rcpp::_["t_u"]      = Rcpp::clone(md.t_u),
    Rcpp::_["t_m_in_N_tilde"] = Rcpp::clone(md.t_m_in_N_tilde),
    Rcpp::_["L_m"]      = Rcpp::clone(md.L_m),
    Rcpp::_["R_m"]      = Rcpp::clone(md.R_m_vec), // expose as R_m
    Rcpp::_["t_m"]      = Rcpp::clone(md.t_m),
    Rcpp::_["E_star"]   = Rcpp::clone(md.E_star),
    Rcpp::_["T_star"]   = Rcpp::clone(md.T_star),
    Rcpp::_["L_bar"]    = Rcpp::clone(md.L_bar),
    Rcpp::_["R_bar"]    = Rcpp::clone(md.R_bar),
    Rcpp::_["s_j_c"]    = Rcpp::clone(md.s_j_c),
    Rcpp::_["s_j_full"] = Rcpp::clone(md.s_j_full),
    Rcpp::_["Q"]        = Rcpp::clone(md.Q),
    Rcpp::_["Q_i_mark"] = Rcpp::clone(md.Q_i_mark),
    Rcpp::_["A_union"]  = Rcpp::clone(md.A_union),
    
    // integer vectors
    Rcpp::_["c_k"] = Rcpp::clone(md.c_k),
    Rcpp::_["d_n"] = Rcpp::clone(md.d_n),
    
    // matrices (2 columns)
    Rcpp::_["A_m"]      = Rcpp::clone(md.A_m),
    Rcpp::_["A_u"]      = Rcpp::clone(md.A_u),
    Rcpp::_["A_c"]      = Rcpp::clone(md.A_c),
    Rcpp::_["full_A_m"] = Rcpp::clone(md.full_A_m),
    Rcpp::_["Q_i"]      = Rcpp::clone(md.Q_i)
  );
}





