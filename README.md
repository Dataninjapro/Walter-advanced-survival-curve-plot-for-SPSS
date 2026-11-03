# Walter: Advanced Survival Curve Plot for SPSS

Kaplan-Meier survival figures generated directly from SPSS Syntax using the SPSS R integration. Version 7 uses `survival` and standard R graphics; no `ggplot2`, `survminer`, extension bundle, or package installation is required when `survival` is already available.

## Quick start

1. Download [SPSS_KM_Clinic_v7_LineStyle.sps](SPSS_KM_Clinic_v7_LineStyle.sps).
2. Open your dataset in IBM SPSS Statistics and make it the active dataset.
3. Open the `.sps` file in the SPSS Syntax Editor.
4. Review the **USER SETTINGS** block at the top, especially variable names, event coding, group order and analysis options.
5. Run the **entire file**, including `BEGIN PROGRAM R.` and `END PROGRAM.`.
6. View the figure and its generated caption in the SPSS Output Viewer.

For detailed setup, settings, interpretation and troubleshooting, see the [Hebrew user guide](docs/USAGE.he.md).

## Features

- Kaplan-Meier curves by group, confidence ribbons and censoring marks.
- A number-at-risk table aligned to the time axis.
- Overall log-rank test and conditional significance wording.
- Optional unadjusted Cox HR with confidence interval for exactly two groups, with the comparison direction stated explicitly.
- A polished caption assembled from enabled analysis results.
- Distinct group colors, with either all-solid lines or different line types.
- Explicit high-resolution PNG submission to the SPSS Viewer.
- Optional PNG and vector PDF export.
- Detailed R tables hidden by default; errors and model diagnostic messages remain visible.

## Requirements and validation

- IBM SPSS Statistics with working R integration (`BEGIN PROGRAM R`). SPSS 31/32 are the intended target versions.
- The `survival` package must load in the R runtime used by SPSS. The script does not install packages.
- Windows was used during development. The user's preflight reported R 4.6.1, `survival` 3.8-6 and `spssstatistics` 32.0.1.0.
- Plot generation and analysis options were exercised with synthetic data; the user confirmed successful plotting from their SPSS dataset. This is not a claim of validation on every SPSS version or operating system.

## Expected example data

One row per person, ordinary right-censored follow-up starting at time zero:

| Variable | Meaning |
| --- | --- |
| `id` | Unique, nonmissing participant identifier |
| `clinic` | Group code, e.g. 0 and 1 |
| `status` | 1 = event; 0 = censored |
| `survtime` | Nonnegative numeric follow-up time in days |

`prison` and `methdose` are not used. All variable names and event/censor codes can be changed in the settings block. No participant data or real-data figures are included in this repository.

## Current v7 defaults

```r
SHOW_POLISHED_CAPTION <- TRUE
VERBOSE_OUTPUT <- FALSE
SHOW_LOGRANK <- TRUE
SHOW_SIGNIFICANCE_TEXT <- TRUE
SHOW_COX_HR <- TRUE
ALL_LINES_SOLID <- TRUE
WIDTH <- 3
HEIGHT <- 3
PNG_DPI <- 300
FONT_PT <- 6
CENSOR_SIZE <- 0.8
CENSOR_LWD <- 0.8
CURVE_LWD <- 0.8
```

The current dimensions produce **900 x 900 pixels**. The inline DPI comment retained in the original v7 file refers to an earlier 4 x 4 inch profile; actual dimensions are always `WIDTH * PNG_DPI` by `HEIGHT * PNG_DPI`. A larger tested viewing profile is 4 x 4 inches, font 8, and line widths 1.0.

To export, create a destination directory and set, for example:

```r
OUTPUT_DIR <- "C:/Results"
FILE_PREFIX <- "KM_clinic_run01"
```

Use forward slashes on Windows. Existing output files are not overwritten: choose a new prefix. An empty `OUTPUT_DIR` displays the image in the Viewer without permanent export. The Viewer can scale displayed images independently of the source resolution.

## Analysis scope

- **Unweighted analysis:** SPSS `WEIGHT` is not applied.
- **Pooled analysis:** SPSS `SPLIT FILE` groups are not analyzed separately.
- Active SPSS case filters are honored during data retrieval. The script does not change these SPSS settings.
- Missing time, status or group values are excluded. Invalid times/status codes or missing/duplicate IDs among complete analysis cases stop the run. Set `VERBOSE_OUTPUT <- TRUE` to inspect counts.
- Supports 2-6 observed groups. With more than two groups, set `SHOW_COX_HR <- FALSE`; log-rank is an overall test, not a pairwise comparison.
- HR is unadjusted and compares the second group to the first. `GROUP_ORDER` controls this order. No covariates are included.
- Cox proportional-hazards diagnostics are calculated when a usable HR is available. A global diagnostic P < 0.05 produces a visible message, even with verbose output off. Model warnings or nonfinite HR confidence limits cause HR omission.
- Confidence intervals use log-log intervals for KM and the configured `CONF_LEVEL` for KM/Cox. Significance wording uses the unrounded log-rank P and `SIGNIFICANCE_ALPHA`. P values below 0.001 are displayed as `P < 0.001`.
- Number at risk is counted immediately before each listed time (`time >= t`), with no delayed entry. No competing-risks, recurrent-event or adjusted-curve analysis is implemented.
- `TIME_MAX` and `TIME_BREAKS` control display only; they do not truncate the data used in the statistical models.

## Repository contents

- [Version 7 Syntax](SPSS_KM_Clinic_v7_LineStyle.sps): the current local v7 file, preserved as supplied.
- [Hebrew user guide](docs/USAGE.he.md): operation, settings and troubleshooting.

Future changes can be tracked through GitHub commits, issues and pull requests.
