# Test Comparison Summary: Two Implementations of Estimators

## Overview
This document summarizes the comparison between two implementations of the NPMLE estimators:
- **Implementation 1**: `calc_F_and_hazards()` in `helper_functions.R`
- **Implementation 2**: Individual functions (`calc_F_12`, `calc_F_13`, `calc_A_12`, `calc_A_13`, `calc_A_23`) in `estimation_of_A_to_test_against.R`

## Test Results

### Perfect Matches (Difference = 0)
The following estimators match perfectly between both implementations:
- ✅ **F13**: Cumulative distribution function for transition 1→3
- ✅ **Lambda13**: Cumulative hazard for transition 1→3
- ✅ **Lambda23**: Cumulative hazard for transition 2→3

### Discrepancies Found
The following estimators show differences:
- ⚠️ **F12**: Max difference ranges from 0.01 to 0.13 depending on sample size
- ⚠️ **Lambda12**: Max difference ranges from 0.22 to 0.78 depending on sample size

## Root Cause Analysis

### The Problem
The difference in Lambda12 and F12 stems from a **conceptual bug** in Implementation 1 (`calc_F_and_hazards()`), not a numerical precision issue.

### Technical Details

**Lambda12 Formula:**
```
Lambda12(s) = sum_{i: r_i <= s} [z_i / (1 - F(l_i-))]
```

where F(l_i-) is the total cumulative distribution just before l_i.

**What Implementation 1 Does (INCORRECT):**
```r
F12_at_l_i <- step_cdf(Q_i[,1]-1e-6, times = Q_i[1:I, 2], masses = z_i[1:I])
denom12 <- 1 - F12_at_l_i  # Using F12 only
term12  <- ifelse(denom12 > 0, z_i[1:I] / denom12, 0)
A12 <- step_cdf(grid_points, times = Q_i[,2], term12)
```

**What Implementation 2 Does (CORRECT):**
```r
F_hat <- calc_F(z_i, Q_i, Q_i_mark, I, I_mark)  # F_total = F12 + F13
F_hat_at_l <- F_hat((Q_i[,1]-1e-10))
denom12 <- 1 - F_hat_at_l  # Using F_total
term12 <- ifelse(denom12 > 0, z_i[1:I] / denom12, 0)
A12 <- step_cdf(grid_points, times = Q_i[,2], term12)
```

### Why This Matters

In the competing risks framework:
- The **risk set** at time l_i consists of individuals who have not yet experienced ANY event
- Once someone experiences event 2 OR event 3, they leave the risk set
- Therefore, the denominator should be `1 - F_total(l_i-)` where `F_total = F12 + F13`
- Using only `1 - F12(l_i-)` incorrectly assumes that people who experienced event 3 are still at risk

## Numerical Impact

### Across Different Sample Sizes

| Scenario | n | Max |F12| diff | Max |Lambda12| diff |
|----------|---|-------------|------------------|
| Small | 50 | 0.134 | 0.223 |
| Medium | 100 | 0.049 | 0.784 |
| Large | 250 | 0.011 | 0.741 |

The Lambda12 differences are particularly large (up to 0.78), which could significantly impact inference and conclusions.

## Conclusion

**Implementation 2 (Individual functions) is THEORETICALLY CORRECT**

**Implementation 1 (calc_F_and_hazards) has a BUG** that needs to be fixed.

## Recommendation

Fix `calc_F_and_hazards()` in `helper_functions.R`:

### Required Changes:
1. **Lambda12 calculation**: Change denominator from `1 - F12(l_i-)` to `1 - F_total(l_i-)`
2. **Lambda13 calculation**: Similarly use `1 - F_total(e_k-)` instead of `1 - F13(e_k-)`

### Code Fix:
```r
# Current (WRONG):
F12_at_l_i <- step_cdf(Q_i[,1]-1e-6, times = Q_i[1:I, 2], masses = z_i[1:I])
denom12 <- 1 - F12_at_l_i

# Should be (CORRECT):
F_total_at_l_i <- step_cdf(Q_i[,1]-1e-6, times = Q_i[1:I, 2], masses = z_i[1:I]) +
                  step_cdf(Q_i[,1]-1e-6, times = Q_i_mark, masses = z_i[(I+1):I_mark])
denom12 <- 1 - F_total_at_l_i
```

## Files Generated

Test scripts:
- `frydman/test_estimator_implementations.R` - Basic comparison
- `frydman/detailed_implementation_comparison.R` - Detailed analysis with plots
- `frydman/investigate_lambda12_difference.R` - Root cause analysis

Plots:
- `frydman/comparison_F12.png`
- `frydman/comparison_F13.png`
- `frydman/comparison_Lambda12.png`

Clean reference implementation:
- `frydman/estimation_of_A_functions_only.R`
