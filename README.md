# halfrule

[![check](https://github.com/Marwan-Alqumbaey/halfrule/actions/workflows/check.yml/badge.svg)](https://github.com/Marwan-Alqumbaey/halfrule/actions/workflows/check.yml)

Which patients should a randomised trial enrol?

Wide eligibility criteria recruit faster, but they bring in patients who gain less from the treatment. `halfrule` finds the eligibility set that is best for a chosen aim. It is available for R and Python with the same functions and results.

## The rules

With equally variable patients:

| Question | Answer |
|:---|:---|
| Should a group be added? | Yes, if its effect is at least half the average effect of the patients already enrolled. Smaller effects can still help if the group is large; `subgroup_threshold()` gives the exact threshold. |
| Should an enrolled group be dropped? | Yes, if its effect is below half the average effect of all enrolled patients, the group included. |
| Where is the best cut-off? | Enrol patients in order of expected effect. At the best cut-off, the last patient enrolled has half the average effect of those enrolled. |
| What if part of the cost is paid per enrolled patient? | The cut-off rises to `1 - w/2` of the average, where `w` is the share of the cost that scales with screening or time. |
| What if the aim is health gain for future patients? | The cut-off falls to about one quarter of the average at 80% power. |

The general forms allow for differences in outcome variance, for model-based analyses (odds ratios, rate ratios, hazard ratios) and for unadjusted analyses. Only the relative effects of patient types are needed.

## Install

R:

```r
# install.packages("remotes")
remotes::install_github("Marwan-Alqumbaey/halfrule", subdir = "r")
```

Python:

```bash
pip install "git+https://github.com/Marwan-Alqumbaey/halfrule#subdirectory=python"
```

## Example

The benefit of a treatment falls linearly with delay and ends at 6.3 hours. Patients arrive evenly over 8 hours. How late should the trial enrol, and what does a 6-hour window cost?

R:

```r
library(halfrule)

time_window(function(t) 1 - t / 6.3, upper = 8)$close
#> 4.2

delay <- seq(0.005, 7.995, by = 0.01)
relative_duration(effect = 1 - delay / 6.3, included = delay <= 6)
#> 1.13   (the trial takes 13% longer for the same power)

# Marker-negative patients are 75% of those screened. Are they worth enrolling
# if their effect is 40% of the effect in marker-positive patients?
0.4 > subgroup_threshold(0.75 / 0.25)
#> TRUE
```

Python:

```python
import numpy as np
from halfrule import relative_duration, subgroup_threshold, time_window

time_window(lambda t: 1 - t / 6.3, upper=8)["close"]
# 4.2

delay = np.arange(0.005, 8, 0.01)
relative_duration(1 - delay / 6.3, included=delay <= 6)
# 1.13

0.4 > subgroup_threshold(0.75 / 0.25)
# True
```

## Functions

| Function | Purpose |
|:---|:---|
| `eligibility()` | Best eligibility set for power, cost or expected health gain |
| `time_window()` | Best treatment time window |
| `relative_duration()` | Duration of a trial with a given set, relative to the best set |
| `design_power()` | Power of a trial with a given set |
| `subgroup_threshold()` | Exact threshold for adding or dropping a whole group |
| `rule_threshold()` | Threshold at the best set, in closed form |
| `inclusion_bound()` | Upper bound on the cost of enrolling everyone |
| `power_elasticity()` | Elasticity of power |

## Repository layout

- `r/`: R package
- `python/`: Python package
- `analysis/`: scripts, results and figures for the accompanying paper

## Reproducing the paper

From the `analysis` folder, with the R package installed:

```bash
Rscript 00_exhaustive_checks.R
Rscript 01_published_examples.R
Rscript 02_simulation.R
IST3_PATH=/path/to/datashare_aug2015.sas7bdat Rscript 03_ist3.R
Rscript 04_figures.R
```

Every change to the repository is checked automatically: the R package is checked, the Python tests are run, and scripts `00` to `02` are rerun and compared with the published results. Script `03` needs the `haven` package and the IST-3 data, which are not included here; download them from <https://doi.org/10.7488/ds/1350>. Script `01` uses the `BioPET` package for one comparison and skips it if the package is not installed. The results are written to `analysis/output` and the figures to `analysis/figures`.

## Citation

Alqumbaey M, Nasher YM. How wide should trial eligibility be? A half rule for enrolling patients whose benefit varies. See `CITATION.cff`.

## Licence

MIT
