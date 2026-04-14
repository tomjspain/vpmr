# vpmr

**vpmr** provides C++ implementations (via Rcpp) of concordance-based discrimination indices for recurrent event prediction models, along with standard and grouped jackknife resampling for uncertainty quantification.

## Requirements

- R with the `Rcpp` package installed
- A C++ compiler (provided automatically on most systems with Rtools / Xcode / build-essential)

## Getting started

```r
library(Rcpp)
sourceCpp("vpmr_cpp_v2.1.cpp")

# Load example dataset
con <- url("https://raw.githubusercontent.com/tomjspain/vpmr/main/example%20df.rds", "rb")
df <- readRDS(con)
close(con)

# Point estimate
df_kendall <- rec_CIndex(df$ID, df$ObsCount, df$PredCount, type = "kendall")

# Leave-one-out jackknife values
df_kendall_jk <- jk_rec_CIndex(df$ID, df$ObsCount, df$PredCount, type = "kendall")
```

---

## Index types

All functions that compute a discrimination index accept a `type` argument. The options are:

| `type` | Formula | Notes |
|---|---|---|
| `"C type1"` | (C + 0.5·Pt) / (M − Ot − Bt) | Excludes obs-tied and both-tied pairs from denominator |
| `"C type2"` | (C + 0.5·(Pt + Bt)) / (M − Ot) | Excludes obs-tied pairs only |
| `"C type3"` | (C + 0.5·(Pt + Bt)) / M | All pairs in denominator |
| `"C type4"` | (C + 0.5·Pt) / (C + D + Pt) | Denominator is non-obs-tied comparable pairs |
| `"kendall"` | (C − D) / M | Kendall's τ-b analogue |
| `"goodman"` | (C − D) / (C + D) | Goodman–Kruskal γ |
| `"somer"` | 2·C type4 − 1 | Somers' D rescaling of C type4 |

Where: C = concordant pairs, D = discordant pairs, Ot = obs-tied only, Pt = pred-tied only, Bt = both tied, M = total pairs.

---

## Functions

### `rec_CIndex(ids, obs, pred, type)`

Computes a single discrimination index over all pairs in the dataset.

**Arguments**

- `ids` - character vector of observation identifiers
- `obs` - numeric vector of observed event counts
- `pred` - numeric vector of predicted event counts
- `type` - index type (see table above)

**Returns** a single numeric value.

```r
rec_CIndex(df$ID, df$ObsCount, df$PredCount, type = "kendall")
```

---

### `jk_rec_CIndex(ids, obs, pred, type)`

Leave-one-out jackknife. Computes the discrimination index on the dataset with each observation removed in turn. Uses an O(n) update formula - pair counts for the full sample are computed once, then each leave-one-out index is derived in O(1) by subtracting the removed observation's pair contributions.

**Arguments** - same as `rec_CIndex`.

**Returns** a named numeric vector of length n (one value per observation). The names correspond to `ids`.

```r
jk_vals <- jk_rec_CIndex(df$ID, df$ObsCount, df$PredCount, type = "kendall")
```

---

### `gjk_rec_CIndex_partition(ids, obs, pred, type, d, d_is_percent = FALSE)`

Grouped jackknife using a **random non-overlapping partition**. Randomly permutes the data and divides it into non-overlapping blocks of size `d` (ceiling division, so the last block may be smaller). The index is computed on the data with each block omitted in turn, yielding G = ⌈n/d⌉ replicates.

**Arguments**

- `ids`, `obs`, `pred`, `type` - as above
- `d` - group size (integer). If `d_is_percent = TRUE`, interpreted as a percentage of n (e.g. `d = 10` means ~10% of n)
- `d_is_percent` - logical, default `FALSE`

**Returns** a numeric vector of length G with attributes: `theta_full`, `n`, `d`, `G`, `d_is_percent`, `method`.

Reproducible via `set.seed()` in R.

```r
set.seed(42)
reps <- gjk_rec_CIndex_partition(df$ID, df$ObsCount, df$PredCount,
                                  type = "kendall", d = 10)
```

---

### `gjk_rec_CIndex_overlap(ids, obs, pred, type, d, n_reps, d_is_percent = FALSE)`

Grouped jackknife using **random overlapping drop sets**. For each of `n_reps` replicates, a fresh random subset of `d` observations is drawn (without replacement within each replicate, but replicates may overlap) and omitted, and the index is computed on the remainder.

**Arguments**

- `ids`, `obs`, `pred`, `type` - as above
- `d` - drop size per replicate. If `d_is_percent = TRUE`, interpreted as a percentage of n
- `n_reps` - number of replicates (must be ≥ 2)
- `d_is_percent` - logical, default `FALSE`

**Returns** a numeric vector of length `n_reps` with attributes: `theta_full`, `n`, `d`, `n_reps`, `d_is_percent`, `method`.

Reproducible via `set.seed()` in R.

```r
set.seed(42)
reps <- gjk_rec_CIndex_overlap(df$ID, df$ObsCount, df$PredCount,
                                type = "kendall", d = 10, n_reps = 200)
```

---

### `gjk_rec_CIndex_overlap_sets(ids, obs, pred, type, drop_sets)`

A fast variant of `gjk_rec_CIndex_overlap` that accepts **pre-drawn drop sets** from R, intended for reproducibility across sessions or use within parallel pipelines where the RNG is managed externally. Uses an O(1)-per-replicate update via precomputed per-observation pair counts, correcting for double-subtraction of internal pairs within each drop set.

**Arguments**

- `ids`, `obs`, `pred`, `type` - as above
- `drop_sets` - an integer matrix of dimensions `n_reps × d`, where each row contains the 1-based indices of observations to drop for that replicate. All indices must be distinct within a row and in the range [1, n].

**Returns** a numeric vector of length `n_reps` with attributes: `theta_full`, `n`, `d`, `n_reps`, `d_is_percent` (always `FALSE`), `method`.

```r
set.seed(42)
drop_sets <- do.call(rbind, replicate(200, sample(nrow(df), 10), simplify = FALSE))
reps <- gjk_rec_CIndex_overlap_sets(df$ID, df$ObsCount, df$PredCount,
                                     type = "kendall", drop_sets = drop_sets)
```

---

### `classify_pair(i, j, obs, pred)` *(diagnostic)*

Classifies a single pair of observations as Concordant, Discordant, ObsTied, PredTied, or BothTied. Uses 1-based indexing.

```r
classify_pair(1, 2, df$ObsCount, df$PredCount)
```

---

### `generate_pairs(ids, obs, pred)` *(diagnostic)*

Enumerates all n(n−1)/2 pairs and returns a data frame with columns `ID1`, `ID2`, and `PairType`. Useful for small datasets to verify pair classifications. **Not recommended for large n** due to memory usage.

```r
pairs_df <- generate_pairs(df$ID, df$ObsCount, df$PredCount)
```

---

## Notes

- All pair counts use 64-bit integers internally to avoid overflow for large n.
- Internal consistency checks verify that pair counts sum to M = n(n−1)/2 after each full pass.
- The grouped jackknife functions attach `theta_full` (the full-sample statistic) as an attribute on the returned vector, which can be retrieved with `attr(reps, "theta_full")`.
