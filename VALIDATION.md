# Packaging checks

Checks performed on 17 November 2026 using the supplied `scripts.zip` and `compact.zip`.

## Supplied data and model

- All four CSVs are byte-for-byte identical to the supplied compact data. Their SHA-256 checksums are recorded in `data/compact/SHA256SUMS`.
- The supplied `attribution_20260821.stan` is byte-for-byte unchanged.
- All 26 packaged R scripts pass R's syntax parser.
- The Stan model passes standalone stanc3 syntax and semantic checking (compiler identifier `cb801dc`). This is not an HMC fit or a C++ executable build.

## Checks with the actual CSVs

The packaged `scripts/check_compact_data.R` was run for both diagnosis scopes. It verified the complete ISO calendar, eight age groups, four pathogen series and validation grids, then successfully prepared both Stan data lists.

| Scope | Study weeks | RTI diagnosis-category count sum | Validation rows |
|---|---:|---:|---:|
| All diagnoses | 365 | 515,069 | 1,460 |
| Main diagnosis | 365 | 289,355 | 1,460 |

Both scopes include 2015-W53. Virology is identical between scopes. Neither original hospital extracts nor the virology workbook are needed for these checks.

## Changed validation and plotting code

Using the actual confirmed-count CSVs and a synthetic posterior-summary interface, stage 13 completed for both scopes using the supplied 2,000-replicate, 13-week block-bootstrap settings. The compact reader used by script 16 returned the same confirmed counts. This checks data flow and computation, not the fitted scientific results.

The changed beta prior/posterior plot rendered using synthetic beta draws. Its prior quantiles were checked against `qgamma(..., shape = 2, rate = 2)`.

Synthetic posterior summaries, test figures and other generated test outputs are not included in the release package.

## Environment and limits

Checks used R 4.5.0 with readr 2.1.5, dplyr 1.1.4, tidyr 1.3.1, ggplot2 3.5.1, posterior 1.6.1, cmdstanr 0.9.0, loo 2.8.0, cowplot 1.1.3 and withr 3.0.2. These record the checking environment; they are not a lockfile of the original analysis environment.

MAP, HMC and LOO were not rerun. Saved posterior chains were not supplied, so the full eight-figure manuscript export was not executed against actual fitted posteriors. The original hospital/virology extracts were not supplied, so the original preprocessing could not be independently regenerated. Optional scripts 17–18 still require the original virological workbook, as documented in the README.
