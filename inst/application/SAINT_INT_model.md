# SAINTexpress-int model notes

This note describes the intensity-mode model implemented in `R/saint_int.R`.
It is meant as a code-level explanation of what is fitted to what, not as a
formal re-derivation of the SAINT publication.

## Data objects

The intensity engine uses three SAINT input tables:

- `inter`: one row per observed prey measurement in one IP replicate, with
  `ipId`, `baitId`, `preyId`, and `quant`.
- `prey`: prey annotation, at minimum `preyId` and `preyGeneId`.
- `bait`: IP replicate annotation, with `ipId`, `baitId`, and `CorT`, where
  `T` is a test bait replicate and `C` is a control replicate.

The model is fitted after these transformations:

1. Intensities are transformed as `log(quant)`.
2. Missing prey measurements are represented as `NA`.
3. Test and control log intensities are centered and scaled together using all
   finite test and control values.
4. Self-interactions are masked: if `preyGeneId == baitId`, the corresponding
   test measurements are set to `NA`.

After this, the key observed matrices are:

- `ctrl_mat[prey, control_ip]`: normalized log intensity in control IPs.
- `test_mat_DATA[prey, test_ip]`: normalized log intensity in test IPs.
- `test_mat[prey, bait]`: a list matrix where each cell contains the replicate
  vector for one prey-bait pair.

## Biological question encoded by the model

For each prey-bait replicate measurement, SAINTexpress-int asks whether the
measurement is better explained as:

- a background or false interaction for that prey; or
- a true interaction for that prey and bait.

The latent indicator is:

```text
Z[prey, bait, replicate] = TRUE   if the replicate is treated as a true signal
Z[prey, bait, replicate] = FALSE  if the replicate is treated as background
```

The model therefore fits two intensity distributions per prey:

- a prey-specific background distribution learned mostly from controls; and
- a prey-specific true-signal distribution learned from high test intensities.

It also fits an MRF-style prior over the `Z` indicators. Without topology, this
prior is only a global preference for true versus false labels. With topology,
linked prey can increase the prior probability of a true interaction for a
prey-bait pair.

## Control/background model

For each prey, the background mean is called `eta`.

Controls are used to estimate:

```text
eta[prey]      = expected normalized log intensity under background
sd_false[prey] = background standard deviation
```

The control model handles missing controls by estimating a global missing-value
threshold:

```text
t = mean of observed control values for prey rows that have at least one NA
```

If this is not finite, it falls back to the mean of all finite control values.
This `t` is used later as the imputed log intensity for missing replicate
values during likelihood and score calculations.

For each prey:

- if no control values are detected, `eta = t` and `sd_false = sd_NA`;
- if at least three control values are detected, `eta` is the observed control
  mean and `sd_false` is the MLE standard deviation of the detected controls;
- if only one or two control values are detected, `eta` is the observed control
  mean and `sd_false` is predicted from a fitted relationship between control
  mean and log standard deviation.

The relationship used for sparse controls is:

```text
log(sd) ~= intercept + slope * mean_control
```

fitted across complete control rows. If there are too few complete control rows,
or the control means have no variance, the implementation uses a median
standard deviation fallback. All background standard deviations are lower
bounded by `sd_NA`.

## True-signal model

The true-signal mean is modeled as a prey-specific shift above background:

```text
mu_true[prey] = eta[prey] + d[prey]
```

The shift `d` is estimated from test intensities for the same prey:

```text
d[prey] = max(mean(test values above eta) - eta, log(4))
```

If no test values are above `eta`, the mean term falls back to zero before the
`max()` is applied, so the shift is still at least `log(4)`.

Initial latent states are assigned by a simple threshold:

```text
Z = TRUE  when observed test intensity > eta + log(5)
Z = FALSE for missing values and lower observed values
```

The true-signal standard deviation is then estimated from test values currently
assigned `Z = TRUE`:

- if at least two true-assigned test values exist for the prey, `sd_true` is
  their sample standard deviation;
- otherwise `sd_true` falls back to `sd_false`;
- `sd_true` is also lower bounded by `sd_NA`.

## Observation likelihood

For a single observed or imputed normalized log intensity `y`, the model uses
two normal densities:

```text
background density: Normal(eta[prey], sd_false[prey])
true density:       Normal(eta[prey] + d[prey], sd_true[prey])
```

The implementation truncates the value used in each density so extreme values
do not dominate:

```text
true component uses       min(y, eta + d)
background component uses max(y, eta)
```

The normal log density is additionally capped at five standard deviations from
the mean inside `.saint_int_quant_log_pdf()`. This is a numerical guard copied
from the SAINTexpress behavior.

Missing test values are imputed as `t` for likelihood calculations. During
final replicate scoring, missing test values get posterior probability zero.

## MRF prior over latent states

The model includes a prior term over `Z`. In the no-topology case, the true
state weight is:

```text
mrf_true  = exp(beta1)
mrf_false = exp(beta0)
beta0     = 0
```

The parameter actually fitted is `beta1`. Larger `beta1` means the model is
more willing to assign replicate states to true interaction.

When topology is provided through `p2p_mapping`, the true state weight becomes:

```text
mrf_true = exp(beta1 + gamma * gsum[prey, bait])
```

where:

```text
gsum[prey, bait] = sum over linked prey of mean(Z[linked_prey, bait, ])
```

So `gamma` controls how much evidence from linked prey increases the prior
weight for a true state. `gamma` is constrained to be non-negative. If the
topology model does not improve over `gamma = 0`, the implementation resets
`gamma` to zero.

In the current `runSaint(engine = "r")` path, no topology input is exposed, so
the practical default is the no-topology model.

## What is optimized

The implementation alternates between two steps:

1. Update the latent state matrix `Z`.
2. Refit MRF parameters (`beta1`, and optionally `gamma`).

This is an ICM-style procedure, not a full posterior sampler.

The no-topology MRF objective for `beta1` is:

```text
sum over prey-bait cells:
  log(prod over replicates P(Z_rep | beta1)) / number_of_replicates
```

The topology objective is the same structure, but uses:

```text
P(Z_rep | beta1, gamma, linked prey Z states)
```

The optimizer is either:

- base R `stats::optim(..., method = "L-BFGS-B")`; or
- optional `nloptr` COBYLA when `optimizer = "nloptr"` and `nloptr` is
  installed.

Parameter bounds are:

```text
beta1 in [-15, 15]
gamma in [0, 10]
```

The ICM loop runs for at most 15 iterations and stops early when the likelihood
improvement is very small.

## Latent-state update

For each prey-bait-replicate state, the algorithm compares the local
log-likelihood before and after flipping `Z`.

If flipping `Z` improves the local log-likelihood, the flipped value is kept.
Otherwise the old value is restored.

This gives a hard assignment of each replicate to true or background. It is the
main reason this implementation is closer to coordinate optimization than to a
fully Bayesian posterior integration.

## Final scores

After fitting, the code computes replicate-level posterior-like probabilities:

```text
true_log  = log(mrf_true)  + log true density
false_log = log(mrf_false) + log background density

score = exp(true_log) / (exp(true_log) + exp(false_log))
```

The calculation is performed on the log scale for numerical stability.

For each prey-bait pair:

- `AvgP` is the mean of replicate scores, or the mean of the top `R` scores if
  there are more than `R` replicates.
- `MaxP` is the maximum replicate score.
- `OddsScore` is the smallest replicate log odds among the sorted replicate
  odds.
- missing replicate intensities get replicate score zero.

With no topology, `TopoAvgP`, `TopoMaxP`, and topological odds are identical to
the ordinary scores.

`SaintScore` is:

```text
max(AvgP, TopoAvgP)
```

which is equal to `AvgP` in the no-topology path.

## BFDR

BFDR is computed from the distribution of `AvgP` values over all prey-bait
matrix cells, not only the rows emitted in `list.txt`.

For a prey-bait pair with score `avg_p`:

```text
higher = all AvgP values greater than avg_p
BFDR   = 0                         if there are no higher values
BFDR   = 1 - mean(higher AvgP)      otherwise
```

This makes BFDR a monotone summary of how high a score sits relative to the
full fitted prey-bait score matrix.

## Output scale

The model is fitted on globally normalized log intensities. The final list
output is back-transformed with `exp()`:

- `Intensity` reports per-replicate test intensities, using `"."` for missing
  values.
- `IntensitySum` and `AvgIntensity` use imputed `exp(t)` for missing test
  values.
- `ctrlIntensity` reports control intensities, again using `"."` for missing
  values.
- `FoldChange` is:

```text
AvgIntensity / exp(eta[prey])
```

So the denominator is the prey-specific fitted background intensity after
back-transformation.

## Short summary

SAINTexpress-int fits a two-component prey-specific intensity model:

```text
background: Normal(eta, sd_false)
true:       Normal(eta + d, sd_true)
```

The control IPs determine the background distribution. High test intensities
determine the true-signal shift and variance. A latent binary matrix `Z`
assigns each replicate to background or true signal. An MRF prior, usually with
only `beta1` in the current R interface, regularizes those assignments. Final
SAINT scores are posterior-like probabilities derived from the fitted density
ratio and MRF prior.
