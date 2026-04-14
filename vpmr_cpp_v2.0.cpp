#include <Rcpp.h>
#include <algorithm>
#include <cstdint>
using namespace Rcpp;

// Helper: classify a pair (same logic as before)
static inline int classify_code(double o_i, double o_j, double p_i, double p_j) {
  if (o_i == o_j)      return (p_i == p_j ? 4 : 2);   // Both tied vs Obs tied
  if (p_i == p_j)      return 3;                      // Pred tied only
  if ((o_i > o_j && p_i > p_j) || (o_j > o_i && p_j > p_i)) return 0; // Concordant
  return 1; // Discordant
}

// [[Rcpp::export]]
std::string classify_pair(int i, int j,
                          const NumericVector& obs,
                          const NumericVector& pred) {
  int ii = i - 1;
  int jj = j - 1;
  double o_i = obs[ii], o_j = obs[jj];
  double p_i = pred[ii], p_j = pred[jj];

  int code = classify_code(o_i, o_j, p_i, p_j);
  switch (code) {
    case 0: return "Concordant";
    case 1: return "Discordant";
    case 2: return "ObsTied";
    case 3: return "PredTied";
    default: return "BothTied"; // 4
  }
}

// [[Rcpp::export]]
DataFrame generate_pairs(CharacterVector ids,
                         NumericVector obs,
                         NumericVector pred) {
  int n = ids.size();
  std::int64_t n64 = static_cast<std::int64_t>(n);
  std::int64_t M64 = n64 * (n64 - 1) / 2; // number of pairs, 64-bit

  // Reserve using size_t (may throw if too large for memory)
  std::vector<std::string> id1; id1.reserve(static_cast<size_t>(M64));
  std::vector<std::string> id2; id2.reserve(static_cast<size_t>(M64));
  std::vector<std::string> pt;  pt.reserve(static_cast<size_t>(M64));

  for (int i = 0; i < n - 1; ++i) {
    double o_i = obs[i], p_i = pred[i];
    for (int j = i + 1; j < n; ++j) {
      double o_j = obs[j], p_j = pred[j];
      id1.push_back(as<std::string>(ids[i]));
      id2.push_back(as<std::string>(ids[j]));
      int code = classify_code(o_i, o_j, p_i, p_j);
      switch (code) {
        case 0: pt.push_back("Concordant"); break;
        case 1: pt.push_back("Discordant"); break;
        case 2: pt.push_back("ObsTied");    break;
        case 3: pt.push_back("PredTied");   break;
        default: pt.push_back("BothTied");   break;
      }
    }
  }

  return DataFrame::create(
    _["ID1"]      = id1,
    _["ID2"]      = id2,
    _["PairType"] = pt
  );
}

// [[Rcpp::export]]
double rec_CIndex(CharacterVector ids,
                  NumericVector obs,
                  NumericVector pred,
                  std::string type) {
  int n = ids.size();
  if (n < 2) Rcpp::stop("Need at least 2 observations.");
  std::int64_t n64 = static_cast<std::int64_t>(n);
  std::int64_t M = n64 * (n64 - 1) / 2; // total number of pairs, 64-bit
  if (M <= 0) Rcpp::stop("No comparable pairs or overflow.");

  // Global counts (64-bit)
  std::int64_t C = 0, D = 0, Ot = 0, Pt_ct = 0, Bt = 0;

  // Single pass over all pairs
  for (int i = 0; i < n - 1; ++i) {
    double o_i = obs[i], p_i = pred[i];
    for (int j = i + 1; j < n; ++j) {
      double o_j = obs[j], p_j = pred[j];
      int code = classify_code(o_i, o_j, p_i, p_j);
      switch (code) {
        case 0: ++C;    break;
        case 1: ++D;    break;
        case 2: ++Ot;   break;
        case 3: ++Pt_ct;break;
        case 4: ++Bt;   break;
      }
    }
  }

  // Optional consistency check (cheap):
  if (C + D + Ot + Pt_ct + Bt != M) {
    Rcpp::stop("Internal error: pair counts do not sum to M.");
  }

  // Validate index type
  std::vector<std::string> valid = {
    "C type1", "C type2", "C type3",
    "kendall", "goodman", "C type4", "somer"
  };
  if (std::find(valid.begin(), valid.end(), type) == valid.end()) {
    Rcpp::stop("Invalid type specified. Choose one of: C type1, C type2, C type3, kendall, goodman, C type4, somer.");
  }

  // Compute the requested C-index
  double result = NA_REAL;
  if (type == "C type1") {
    double denom = static_cast<double>(M - Ot - Bt);
    if (denom == 0.0) Rcpp::stop("Division by zero: M - Ot - Bt equals zero.");
    result = (static_cast<double>(C) + 0.5 * static_cast<double>(Pt_ct)) / denom;
  } else if (type == "C type2") {
    double denom = static_cast<double>(M - Ot);
    if (denom == 0.0) Rcpp::stop("Division by zero: M - Ot equals zero.");
    result = (static_cast<double>(C) + 0.5 * static_cast<double>(Pt_ct + Bt)) / denom;
  } else if (type == "C type3") {
    double denom = static_cast<double>(M);
    if (denom == 0.0) Rcpp::stop("Division by zero: M equals zero.");
    result = (static_cast<double>(C) + 0.5 * static_cast<double>(Pt_ct + Bt)) / denom;
  } else if (type == "kendall") {
    double denom = static_cast<double>(M);
    if (denom == 0.0) Rcpp::stop("Division by zero: M equals zero.");
    result = (static_cast<double>(C) - static_cast<double>(D)) / denom;
  } else if (type == "goodman") {
    double denom = static_cast<double>(C + D);
    if (denom == 0.0) Rcpp::stop("Division by zero: C + D equals zero.");
    result = (static_cast<double>(C) - static_cast<double>(D)) / denom;
  } else if (type == "C type4") {
    double denom = static_cast<double>(C + D + Pt_ct);
    if (denom == 0.0) Rcpp::stop("Division by zero: C + D + Pt equals zero.");
    result = (static_cast<double>(C) + 0.5 * static_cast<double>(Pt_ct)) / denom;
  } else { // "somer"
    double denom = static_cast<double>(C + D + Pt_ct);
    if (denom == 0.0) Rcpp::stop("Division by zero: C + D + Pt equals zero.");
    double c4 = (static_cast<double>(C) + 0.5 * static_cast<double>(Pt_ct)) / denom;
    result = 2.0 * c4 - 1.0;
  }
  return result;
}

// [[Rcpp::export]]
NumericVector jk_rec_CIndex(CharacterVector ids,
                            NumericVector obs,
                            NumericVector pred,
                            std::string type) {
  int n = ids.size();
  if (n < 2) Rcpp::stop("Need at least 2 observations.");
  std::int64_t n64 = static_cast<std::int64_t>(n);
  std::int64_t M = n64 * (n64 - 1) / 2; // total pairs
  if (M <= 0) Rcpp::stop("No comparable pairs or overflow.");

  // Global totals (64-bit)
  std::int64_t C_tot = 0, D_tot = 0, Ot_tot = 0, Pt_tot = 0, Bt_tot = 0;
  // Per-ID counts (64-bit)
  std::vector<std::int64_t> C_id(n, 0), D_id(n, 0), Ot_id(n, 0), Pt_id(n, 0), Bt_id(n, 0);

  // Single pass over pairs accumulating global and per-ID counts
  for (int i = 0; i < n - 1; ++i) {
    double o_i = obs[i], p_i = pred[i];
    for (int j = i + 1; j < n; ++j) {
      double o_j = obs[j], p_j = pred[j];
      int code = classify_code(o_i, o_j, p_i, p_j);
      switch (code) {
        case 0: ++C_tot;  ++C_id[i];  ++C_id[j];  break;
        case 1: ++D_tot;  ++D_id[i];  ++D_id[j];  break;
        case 2: ++Ot_tot; ++Ot_id[i]; ++Ot_id[j]; break;
        case 3: ++Pt_tot; ++Pt_id[i]; ++Pt_id[j]; break;
        case 4: ++Bt_tot; ++Bt_id[i]; ++Bt_id[j]; break;
      }
    }
  }

  // Consistency check
  if (C_tot + D_tot + Ot_tot + Pt_tot + Bt_tot != M) {
    Rcpp::stop("Internal error: pair counts do not sum to M.");
  }

  // Validate type
  std::vector<std::string> valid = {
    "C type1", "C type2", "C type3",
    "kendall", "goodman", "C type4", "somer"
  };
  if (std::find(valid.begin(), valid.end(), type) == valid.end()) {
    Rcpp::stop("Invalid type specified. Choose one of: C type1, C type2, C type3, kendall, goodman, C type4, somer.");
  }

  NumericVector jk_vals(n);
  jk_vals.attr("names") = ids;
  std::int64_t pairs_per_id = n64 - 1; // each row participates in n-1 pairs

  // Compute each jackknife index in O(1) using 64-bit arithmetic
  for (int k = 0; k < n; ++k) {
    std::int64_t Ck  = C_tot  - C_id[k];
    std::int64_t Dk  = D_tot  - D_id[k];
    std::int64_t Otk = Ot_tot - Ot_id[k];
    std::int64_t Ptk = Pt_tot - Pt_id[k];
    std::int64_t Btk = Bt_tot - Bt_id[k];
    std::int64_t Nk  = M      - pairs_per_id;

    if (Nk <= 0) Rcpp::stop("Jackknife: non-positive number of pairs after removal.");

    double val;
    if (type == "C type1") {
      double denom = static_cast<double>(Nk - Otk - Btk);
      if (denom == 0.0) Rcpp::stop("Division by zero (jk): Nk - Ot - Bt equals zero.");
      val = (static_cast<double>(Ck) + 0.5 * static_cast<double>(Ptk)) / denom;
    } else if (type == "C type2") {
      double denom = static_cast<double>(Nk - Otk);
      if (denom == 0.0) Rcpp::stop("Division by zero (jk): Nk - Ot equals zero.");
      val = (static_cast<double>(Ck) + 0.5 * static_cast<double>(Ptk + Btk)) / denom;
    } else if (type == "C type3") {
      double denom = static_cast<double>(Nk);
      if (denom == 0.0) Rcpp::stop("Division by zero (jk): Nk equals zero.");
      val = (static_cast<double>(Ck) + 0.5 * static_cast<double>(Ptk + Btk)) / denom;
    } else if (type == "kendall") {
      double denom = static_cast<double>(Nk);
      if (denom == 0.0) Rcpp::stop("Division by zero (jk): Nk equals zero.");
      val = (static_cast<double>(Ck) - static_cast<double>(Dk)) / denom;
    } else if (type == "goodman") {
      double denom = static_cast<double>(Ck + Dk);
      if (denom == 0.0) Rcpp::stop("Division by zero (jk): C + D equals zero.");
      val = (static_cast<double>(Ck) - static_cast<double>(Dk)) / denom;
    } else if (type == "C type4") {
      double denom = static_cast<double>(Ck + Dk + Ptk);
      if (denom == 0.0) Rcpp::stop("Division by zero (jk): C + D + Pt equals zero.");
      val = (static_cast<double>(Ck) + 0.5 * static_cast<double>(Ptk)) / denom;
    } else { // "somer"
      double denom = static_cast<double>(Ck + Dk + Ptk);
      if (denom == 0.0) Rcpp::stop("Division by zero (jk): C + D + Pt equals zero.");
      double c4k = (static_cast<double>(Ck) + 0.5 * static_cast<double>(Ptk)) / denom;
      val = 2.0 * c4k - 1.0;
    }

    jk_vals[k] = val;
  }

  return jk_vals;
}

// [[Rcpp::export]]
Rcpp::NumericVector gjk_rec_partition(CharacterVector ids,
                                         NumericVector obs,
                                         NumericVector pred,
                                         std::string type,
                                         int d,
                                         bool d_is_percent = false) {
  // Respect R's RNG + set.seed()
  Rcpp::RNGScope scope;

  const int n = ids.size();
  if (n < 2) Rcpp::stop("Need at least 2 observations.");
  if (d_is_percent) {
    // Interpret 'd' as a percentage in [1, 99], e.g. 1 => ~1% of n
    if (d <= 0 || d >= 100) {
      Rcpp::stop("When d_is_percent=TRUE, d must be in 1..99 (percent of n).");
    }
    double frac = static_cast<double>(d) / 100.0;
    d = std::max(1, std::min(n - 1, static_cast<int>(std::round(frac * n))));
  } else {
    if (d <= 0 || d >= n) {
      Rcpp::stop("d must be between 1 and n-1.");
    }
  }

  // Number of non-overlapping groups (last group may be smaller)
  const int G = (n + d - 1) / d; // ceil(n / d)
  if (G < 2) {
    Rcpp::stop("Resulting number of groups G must be at least 2.");
  }

  // Random permutation of indices 0..n-1 (reproducible via set.seed in R)
  // Rcpp::sample returns 1..n; convert to 0-based
  Rcpp::IntegerVector perm1 = Rcpp::sample(n, n, false);
  std::vector<int> perm(n);
  for (int i = 0; i < n; ++i) perm[i] = perm1[i] - 1;

  // Precompute full-sample statistic (often useful downstream)
  // (You might not need this here, but it's cheap relative to O(n^2) replicates
  //  and sometimes used in pseudo-value construction in R.)
  double theta_full = rec_CIndex(ids, obs, pred, type);

  // Storage for replicate estimates (one per dropped group)
  Rcpp::NumericVector reps(G);

  // Scratch containers for subset
  std::vector<std::string> ids_keep;
  std::vector<double>      obs_keep, pred_keep;
  ids_keep.reserve(n);
  obs_keep.reserve(n);
  pred_keep.reserve(n);

  for (int g = 0; g < G; ++g) {
    // Determine which permuted indices belong to the g-th dropped block
    const int start = g * d;
    const int end   = std::min(n, (g + 1) * d); // exclusive
    // Mark to drop
    std::vector<char> drop(n, 0);
    for (int t = start; t < end; ++t) drop[ perm[t] ] = 1;

    // Build kept subset
    ids_keep.clear();
    obs_keep.clear();
    pred_keep.clear();
    ids_keep.reserve(n - (end - start));
    obs_keep.reserve(n - (end - start));
    pred_keep.reserve(n - (end - start));

    for (int i = 0; i < n; ++i) {
      if (!drop[i]) {
        ids_keep.push_back( Rcpp::as<std::string>(ids[i]) );
        obs_keep.push_back( obs[i] );
        pred_keep.push_back( pred[i] );
      }
    }

    // Safety: ensure at least 2 remain (rec_CIndex requires >= 2)
    if (ids_keep.size() < 2) {
      Rcpp::stop("Grouped jackknife replicate left fewer than 2 observations.");
    }

    // Convert kept vectors back to Rcpp types
    Rcpp::CharacterVector ids_k(ids_keep.begin(), ids_keep.end());
    Rcpp::NumericVector   obs_k(obs_keep.begin(),  obs_keep.end());
    Rcpp::NumericVector   pred_k(pred_keep.begin(), pred_keep.end());

    // Compute the statistic on the kept subset
    reps[g] = rec_CIndex(ids_k, obs_k, pred_k, type);
  }

  // Attach a few helpful attributes for downstream CI code (optional)
  reps.attr("theta_full")   = theta_full;
  reps.attr("n")            = n;
  reps.attr("d")            = d;
  reps.attr("G")            = G;
  reps.attr("d_is_percent") = d_is_percent;
  reps.attr("method")       = "gjk_rec_partition";

  return reps;
}

// [[Rcpp::export]]
Rcpp::NumericVector gjk_rec_overlap(CharacterVector ids,
                                           NumericVector obs,
                                           NumericVector pred,
                                           std::string type,
                                           int d,
                                           int n_reps,
                                           bool d_is_percent = false) {
  // Respect R's RNG + set.seed()
  Rcpp::RNGScope scope;

  const int n = ids.size();
  if (n < 2) Rcpp::stop("Need at least 2 observations.");

  // Interpret d (count or percent)
  if (d_is_percent) {
    if (d <= 0 || d >= 100) {
      Rcpp::stop("When d_is_percent=TRUE, d must be in 1..99 (percent of n).");
    }
    double frac = static_cast<double>(d) / 100.0;
    d = std::max(1, std::min(n - 1, static_cast<int>(std::round(frac * n))));
  } else {
    if (d <= 0 || d >= n) {
      Rcpp::stop("d must be between 1 and n-1.");
    }
  }

  if (n_reps < 2) {
    Rcpp::stop("n_reps must be at least 2.");
  }

  // Validate index type (keep consistent with rec_CIndex)
  std::vector<std::string> valid = {
    "C type1", "C type2", "C type3",
    "kendall", "goodman", "C type4", "somer"
  };
  if (std::find(valid.begin(), valid.end(), type) == valid.end()) {
    Rcpp::stop("Invalid type specified. Choose one of: C type1, C type2, C type3, kendall, goodman, C type4, somer.");
  }

  // Full-sample statistic (useful downstream)
  double theta_full = rec_CIndex(ids, obs, pred, type);

  // Output vector
  Rcpp::NumericVector reps(n_reps);

  // Scratch containers reused across replicates
  std::vector<char> drop(n, 0);
  std::vector<std::string> ids_keep; ids_keep.reserve(n);
  std::vector<double>      obs_keep; obs_keep.reserve(n);
  std::vector<double>      pred_keep; pred_keep.reserve(n);

  for (int r = 0; r < n_reps; ++r) {
    // Fresh random set of d distinct indices (0-based)
    Rcpp::IntegerVector drop_idx = Rcpp::sample(n, d, false);
    std::fill(drop.begin(), drop.end(), 0);
    for (int t = 0; t < d; ++t) {
      int j = drop_idx[t] - 1; // Rcpp::sample returns 1..n
      drop[j] = 1;
    }

    // Build kept subset
    ids_keep.clear();
    obs_keep.clear();
    pred_keep.clear();
    ids_keep.reserve(n - d);
    obs_keep.reserve(n - d);
    pred_keep.reserve(n - d);

    for (int i = 0; i < n; ++i) {
      if (!drop[i]) {
        ids_keep.push_back(Rcpp::as<std::string>(ids[i]));
        obs_keep.push_back(obs[i]);
        pred_keep.push_back(pred[i]);
      }
    }

    if (ids_keep.size() < 2) {
      Rcpp::stop("Grouped jackknife replicate left fewer than 2 observations.");
    }

    // Convert back to Rcpp and evaluate the statistic
    Rcpp::CharacterVector ids_k(ids_keep.begin(), ids_keep.end());
    Rcpp::NumericVector   obs_k(obs_keep.begin(),  obs_keep.end());
    Rcpp::NumericVector   pred_k(pred_keep.begin(), pred_keep.end());

    reps[r] = rec_CIndex(ids_k, obs_k, pred_k, type);
  }

  // Helpful attributes
  reps.attr("theta_full")   = theta_full;
  reps.attr("n")            = n;
  reps.attr("d")            = d;
  reps.attr("n_reps")       = n_reps;
  reps.attr("d_is_percent") = d_is_percent;
  reps.attr("method")       = "gjk_rec_overlap";

  return reps;
}

// Optional helper: accept pre-drawn drop sets from R for reproducibility / parallel pipelines.
// drop_sets must be an IntegerMatrix of dimension n_reps x d with 1-based indices.
// [[Rcpp::export]]
Rcpp::NumericVector gjk_rec_overlap_sets(CharacterVector ids,
                                                          NumericVector obs,
                                                          NumericVector pred,
                                                          std::string type,
                                                          Rcpp::IntegerMatrix drop_sets) {
  const int n = ids.size();
  if (n < 2) Rcpp::stop("Need at least 2 observations.");
  const int n_reps = drop_sets.nrow();
  const int d      = drop_sets.ncol();
  if (n_reps < 2) Rcpp::stop("n_reps must be at least 2.");
  if (d <= 0 || d >= n) Rcpp::stop("d must be between 1 and n-1.");

  std::vector<std::string> valid = {
    "C type1", "C type2", "C type3",
    "kendall", "goodman", "C type4", "somer"
  };
  if (std::find(valid.begin(), valid.end(), type) == valid.end()) {
    Rcpp::stop("Invalid type specified. Choose one of: C type1, C type2, C type3, kendall, goodman, C type4, somer.");
  }

  // Full-sample statistic
  double theta_full = rec_CIndex(ids, obs, pred, type);

  // Precompute counts
  const std::int64_t n64 = static_cast<std::int64_t>(n);
  const std::int64_t M   = n64 * (n64 - 1) / 2;

  std::int64_t C_tot = 0, D_tot = 0, Ot_tot = 0, Pt_tot = 0, Bt_tot = 0;
  std::vector<std::int64_t> C_id(n, 0), D_id(n, 0), Ot_id(n, 0), Pt_id(n, 0), Bt_id(n, 0);

  for (int i = 0; i < n - 1; ++i) {
    double o_i = obs[i], p_i = pred[i];
    for (int j = i + 1; j < n; ++j) {
      double o_j = obs[j], p_j = pred[j];
      int code = classify_code(o_i, o_j, p_i, p_j);
      switch (code) {
        case 0: ++C_tot;  ++C_id[i];  ++C_id[j];  break;
        case 1: ++D_tot;  ++D_id[i];  ++D_id[j];  break;
        case 2: ++Ot_tot; ++Ot_id[i]; ++Ot_id[j]; break;
        case 3: ++Pt_tot; ++Pt_id[i]; ++Pt_id[j]; break;
        case 4: ++Bt_tot; ++Bt_id[i]; ++Bt_id[j]; break;
      }
    }
  }
  if (C_tot + D_tot + Ot_tot + Pt_tot + Bt_tot != M) {
    Rcpp::stop("Internal error: pair counts do not sum to M.");
  }

  Rcpp::NumericVector reps(n_reps);

  // Scratch
  std::vector<int> S; S.reserve(d);

  for (int r = 0; r < n_reps; ++r) {
    // Convert 1-based indices to 0-based, sanity-check range
    S.clear();
    S.reserve(d);
    for (int t = 0; t < d; ++t) {
      int idx1 = drop_sets(r, t);
      if (idx1 < 1 || idx1 > n) Rcpp::stop("drop_sets contains out-of-range index.");
      S.push_back(idx1 - 1);
    }

    // Remove anything touching S via per-ID totals
    std::int64_t Ck  = C_tot, Dk  = D_tot, Otk = Ot_tot, Ptk = Pt_tot, Btk = Bt_tot;
    for (int u = 0; u < d; ++u) {
      int k = S[u];
      Ck  -= C_id[k];
      Dk  -= D_id[k];
      Otk -= Ot_id[k];
      Ptk -= Pt_id[k];
      Btk -= Bt_id[k];
    }

    // Add back internal(S) pairs (fix double subtraction)
    for (int a = 0; a < d - 1; ++a) {
      int i = S[a];
      double o_i = obs[i], p_i = pred[i];
      for (int b = a + 1; b < d; ++b) {
        int j = S[b];
        double o_j = obs[j], p_j = pred[j];
        int code = classify_code(o_i, o_j, p_i, p_j);
        switch (code) {
          case 0: ++Ck;  break;
          case 1: ++Dk;  break;
          case 2: ++Otk; break;
          case 3: ++Ptk; break;
          case 4: ++Btk; break;
        }
      }
    }

    // Nk for the kept set
    const std::int64_t nk = static_cast<std::int64_t>(n - d);
    const std::int64_t Nk = nk * (nk - 1) / 2;

    // Metric on kept counts
    double val;
    if (type == "C type1") {
      double denom = static_cast<double>(Nk - Otk - Btk);
      if (denom == 0.0) Rcpp::stop("Division by zero (jk): Nk - Ot - Bt equals zero.");
      val = (static_cast<double>(Ck) + 0.5 * static_cast<double>(Ptk)) / denom;
    } else if (type == "C type2") {
      double denom = static_cast<double>(Nk - Otk);
      if (denom == 0.0) Rcpp::stop("Division by zero (jk): Nk - Ot equals zero.");
      val = (static_cast<double>(Ck) + 0.5 * static_cast<double>(Ptk + Btk)) / denom;
    } else if (type == "C type3") {
      double denom = static_cast<double>(Nk);
      if (denom == 0.0) Rcpp::stop("Division by zero (jk): Nk equals zero.");
      val = (static_cast<double>(Ck) + 0.5 * static_cast<double>(Ptk + Btk)) / denom;
    } else if (type == "kendall") {
      double denom = static_cast<double>(Nk);
      if (denom == 0.0) Rcpp::stop("Division by zero (jk): Nk equals zero.");
      val = (static_cast<double>(Ck) - static_cast<double>(Dk)) / denom;
    } else if (type == "goodman") {
      double denom = static_cast<double>(Ck + Dk);
      if (denom == 0.0) Rcpp::stop("Division by zero (jk): C + D equals zero.");
      val = (static_cast<double>(Ck) - static_cast<double>(Dk)) / denom;
    } else if (type == "C type4") {
      double denom = static_cast<double>(Ck + Dk + Ptk);
      if (denom == 0.0) Rcpp::stop("Division by zero (jk): C + D + Pt equals zero.");
      val = (static_cast<double>(Ck) + 0.5 * static_cast<double>(Ptk)) / denom;
    } else { // "somer"
      double denom = static_cast<double>(Ck + Dk + Ptk);
      if (denom == 0.0) Rcpp::stop("Division by zero (jk): C + D + Pt equals zero.");
      double c4k = (static_cast<double>(Ck) + 0.5 * static_cast<double>(Ptk)) / denom;
      val = 2.0 * c4k - 1.0;
    }

    reps[r] = val;
  }

  // Attributes
  reps.attr("theta_full")   = theta_full;
  reps.attr("n")            = n;
  reps.attr("d")            = d;
  reps.attr("n_reps")       = n_reps;
  reps.attr("d_is_percent") = false; // here d is a fixed column count in drop_sets
  reps.attr("method")       = "gjk_rec_overlap_sets";

  return reps;
}
