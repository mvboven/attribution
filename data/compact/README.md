# Compact analysis inputs

These four CSVs contain the aggregate inputs supplied for the public attribution analysis. Source exports were prepared by `scripts/1_read_data.R` and `scripts/13_prepare_validation_data.R`, separately for each diagnosis scope. The CSV bytes are preserved in this package.

| File | Scope | Rows | Columns |
|---|---|---:|---:|
| `attribution_all_diagnoses.csv` | All diagnoses | 365 | 15 |
| `attribution_main_diagnosis.csv` | Main diagnosis | 365 | 15 |
| `validation_confirmed_all_diagnoses.csv` | All diagnoses | 1,460 | 3 |
| `validation_confirmed_main_diagnosis.csv` | Main diagnosis | 1,460 | 3 |

## Admission and virology files

One row represents one ISO week in 2013–2019. ISO weeks start on Monday; the first study week starts on 31 December 2012. The 365 study weeks include 2015-W53.

| Column | Definition |
|---|---|
| `year`, `week` | ISO year and week number |
| `0-2m`, `3-5m`, `6-11m` | Weekly RTI diagnosis-category counts in infants of these ages in months |
| `1-2y`, `3-4y` | Weekly counts in children of these ages in years |
| `65-74y`, `75-84y`, `85+y` | Weekly counts in older adults of these ages in years |
| `VIRINFA`, `VIRINFB` | Weekly influenza A and B laboratory detections |
| `VIRRS` | Weekly RSV laboratory detections |
| `VIRHMPV` | Weekly HMPV laboratory detections |
| `VIRPIV` | Weekly combined PIV laboratory detections: types 1–4 plus untyped/non-typeable categories |

The model uses influenza A+B (`VIRINF = VIRINFA + VIRINFB`), RSV, HMPV and combined PIV. Virology is identical in the two diagnosis scopes.

Hospital cells suppressed as `1-4` were replaced by integer draws from 1–4 by the supplied preprocessing, using seed 20240610. The compact files retain those realised counts, so rereading them does not draw replacements again. Original age groups were combined into the eight model groups. Absent hospital week/age cells were explicitly filled with zero under the supplied settings. Weekly hospital counts sum the non-COVID background and pathogen-specific single/coinfection categories; this sum is not a deduplicated count of individuals with coinfections.

The resulting sums across all study weeks and age groups are 515,069 RTI diagnosis-category counts for all diagnoses and 289,355 for main diagnosis.

## Validation files

One row represents a week/pathogen combination. The grid contains 365 weeks × four pathogens.

| Column | Definition |
|---|---|
| `week` | ISO Monday, stored as `YYYY-MM-DD` |
| `component` | `Influenza`, `RSV`, `HMPV` or `Parainfluenza` |
| `confirmed` | Sum of linked pathogen-specific single and coinfection categories across the eight study age groups |

In the linked-data preprocessing, suppressed `1-4` cells were replaced by their midpoint 2.5. Decimal counts in these aggregate files are therefore intentional. Absent week/pathogen combinations were filled with zero in the original preparation.

The compact validation counts already combine the selected age groups. Use the supplied eight groups when reproducing this validation. Age-specific linked counts and separate PIV subtype series cannot be reconstructed from these CSVs.

Fitted intensities and bootstrap results are not stored in these inputs. Stage 13 recomputes them from the current fitted model. Its validation estimand is the sum of age-specific posterior median pathogen intensities, compared with confirmed counts; it is distinct from allocating observed admissions to pathogens.

## Integrity

`SHA256SUMS` records each CSV's SHA-256 checksum. Run `scripts/check_compact_data.R` from the project root to check schemas, calendars, observed totals and validation grids and to prepare the Stan inputs for both scopes without fitting.
