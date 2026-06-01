# Crime Severity Network — Analysis Code and Data

> **The Moral Architecture of Crime Severity Perception: A Network Psychometric Approach**
>
> Herman E. Elgueta & Beatriz Pérez Sánchez
>
> *Manuscript under review*

---

## Contents

| Path | Description |
|---|---|
| `scripts/network_analysis.R` | Main analysis: estimates EBICglasso networks, community detection, centrality, and all manuscript figures |
| `data/polychoric_gravedad.csv` | 15 × 15 polychoric correlation matrix of *gravedad* (seriousness) ratings — input for the primary network (Figure 2) |
| `data/spearman_composite.csv` | 15 × 15 Spearman correlation matrix of the composite severity index — input for the sensitivity network (Figure 3) |

---

## Reproducing the analyses

### Requirements

```r
install.packages(c("tidyverse", "bootnet", "qgraph", "igraph",
                   "psych", "scales", "gridExtra"))
```

### Run

Set working directory to `scripts/` and run:

```r
source("network_analysis.R")
```

Figures 2–4 are produced directly from the correlation matrices and require no individual-level data. Figure 1 (descriptive boxplots) and bootstrap supplementary figures (S1–S2) require the individual-level dataset, which is available from the corresponding author under a data use agreement.

---

## Data

**Sample:** 274 community residents, Temuco, Chile (2017). Mean age 35.5 years (SD = 15.7); 52.6% female; three socioeconomic sampling zones.

**15 crimes assessed:** Abortion, Child Sexual Abuse, Murder, Fraud, Euthanasia, Tax Evasion, Child Abuse, Piracy, Bribery, Robbery, Terrorism, Drug Trafficking, Vandalism, Rape, Partner Violence.

Individual-level data are not distributed here (ethics restriction). Contact herman.elgueta@umag.cl to request access under a data use agreement.

---

## Citation

> Elgueta, H. E., & Pérez Sánchez, B. (manuscript under review). The moral architecture of crime severity perception: A network psychometric approach.

---

## License

[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)
