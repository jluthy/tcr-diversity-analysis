# 🧬 TCR Diversity Analysis: Apheresis to Drug Product

A reproducible analysis pipeline demonstrating **T-cell receptor (TCR) repertoire diversity changes** during cell therapy manufacturing, with **multi-omic integration** of ELISpot functional immune response data. This project examines how the starting material (apheresis PBMCs) compares to the final drug product in terms of clonal diversity, expansion patterns, and clinical response associations.

> **Note:** All data in this repository is **synthetic** and generated for portfolio demonstration purposes. CDR3 sequences are randomly generated and do not correspond to real patient data.

---

## 📊 Analysis Reports

Pre-rendered HTML reports are available for immediate viewing:

| Notebook | Description | View |
|:---|:---|:---|
| **01 — Data Overview** | Dataset characterization, sample summary, V gene usage, quality checks | [View Report](https://jluthy.github.io/tcr-diversity-analysis/reports/01-jal-data-overview.html) |
| **02 — Diversity Analysis** | Shannon entropy, clonality, richness, D50 — Apheresis vs Product comparison | [View Report](https://jluthy.github.io/tcr-diversity-analysis/reports/02-jal-diversity-analysis.html) |
| **03 — Clonal Expansion** | Top expanded clonotypes, cumulative dominance curves, shared clone tracking | [View Report](https://jluthy.github.io/tcr-diversity-analysis/reports/03-jal-clonal-expansion.html) |
| **04 — edgeR Differential Abundance** | Negative binomial GLM, volcano plots, significant expander/contractor identification | [View Report](https://jluthy.github.io/tcr-diversity-analysis/reports/04-jal-edger-differential-abundance.html) |
| **05 — Clonal Tracking & AUC** | Longitudinal clonotype trajectories, time-integrated exposure, persistence analysis | [View Report](https://jluthy.github.io/tcr-diversity-analysis/reports/05-jal-clonal-tracking-auc.html) |
| **06 — Clinical Response Integration** | Integrated feature matrix, correlation heatmap, composite scoring | [View Report](https://jluthy.github.io/tcr-diversity-analysis/reports/06-jal-clinical-response-integration.html) |
| **07 — Multi-Omic Integration (TCR × ELISpot)** | Functional immune response correlation, multi-omic composite score | [View Report](https://jluthy.github.io/tcr-diversity-analysis/reports/07-jal-multiomics-elispot-integration.html) |

---

## 🔬 Key Findings

1. **Manufacturing reduces diversity by 45–60%** — Shannon entropy drops significantly from Apheresis to Product (paired Wilcoxon p = 0.031).
2. **Clonality increases 2–3×** in the Product, reflecting oligoclonal expansion during manufacturing.
3. **D50 collapses to single digits** in Product samples, indicating extreme dominance by a few clonotypes.
4. **edgeR identifies hundreds of differentially abundant clonotypes** per patient, with significant expanders representing candidate therapeutic clones.
5. **CR patients show the most focused Products and strongest post-infusion expansion**, with higher AUC and persistence ratios than PD patients.
6. **ELISpot functional responses correlate with TCR expansion** — patients with higher clonal AUC show stronger IFN-γ secretion against tumor antigens.
7. **The multi-omic composite score (TCR + ELISpot) provides the best separation** of clinical response groups, supporting a multi-modal immune monitoring approach.

---

## 🚀 Getting Started

### Prerequisites

- **R** ≥ 4.5.0
- **renv** for dependency management

### Setup

```bash
git clone https://github.com/jluthy/tcr-diversity-analysis.git
cd tcr-diversity-analysis
```

Open R in the project root and restore dependencies:

```r
renv::restore()
```

### Run the Pipeline

```bash
make all        # Full pipeline: generate data → render reports
make data       # Generate synthetic data only
make report     # Render R Markdown notebooks only
make clean      # Remove all generated outputs
```

---

## 📁 Repository Structure

```
├── Makefile                          # Reproducible pipeline commands
├── README.md
├── LICENSE
├── renv.lock                         # R environment lockfile
│
├── data/
│   ├── raw/                          # Synthetic source data
│   │   ├── synthetic_tcr_repertoires.csv
│   │   └── patient_metadata.csv
│   └── processed/                    # Computed metrics
│       ├── diversity_metrics.csv
│       ├── top_clonotypes.csv
│       ├── vgene_usage.csv
│       ├── longitudinal_clonal_tracking.csv
│       ├── clonal_auc_metrics.csv
│       ├── elispot_data.csv
│       └── elispot_summary.csv
│
├── src/R/
│   ├── generate_synthetic_data.R     # Data generation script
│   ├── diversity_metrics.R           # Shannon, clonality, D50, Gini-Simpson
│   └── visualization_helpers.R       # Custom ggplot2 theme & plot functions
│
├── notebooks/
│   ├── 01-jal-data-overview.Rmd
│   ├── 02-jal-diversity-analysis.Rmd
│   ├── 03-jal-clonal-expansion.Rmd
│   ├── 04-jal-edger-differential-abundance.Rmd
│   ├── 05-jal-clonal-tracking-auc.Rmd
│   ├── 06-jal-clinical-response-integration.Rmd
│   └── 07-jal-multiomics-elispot-integration.Rmd
│
├── reports/                          # Pre-rendered HTML reports
│   ├── 01-jal-data-overview.html
│   ├── 02-jal-diversity-analysis.html
│   ├── 03-jal-clonal-expansion.html
│   ├── 04-jal-edger-differential-abundance.html
│   ├── 05-jal-clonal-tracking-auc.html
│   ├── 06-jal-clinical-response-integration.html
│   └── 07-jal-multiomics-elispot-integration.html
│
└── www/
    └── report-style.css              # Custom report styling
```

---

## 🧪 Data Description

### TCR Repertoire Data

| Field | Description |
|:---|:---|
| `patient_id` | 6 patients (PT-001 through PT-006) |
| `clinical_response` | CR (n=2), PR (n=2), PD (n=2) |
| `sample_type` | Apheresis (starting material) or Product (drug product) |
| `cdr3_aa` | Randomly generated CDR3β amino acid sequences |
| `v_gene` / `j_gene` | TRBV and TRBJ gene segments (IMGT nomenclature) |
| `clone_count` | Read counts following power-law distribution |
| `clone_fraction` | Proportion within sample |

### Longitudinal Tracking Data

| Field | Description |
|:---|:---|
| `timepoint` | Apheresis, Product, Wk2, Wk4, Wk8, Wk12 |
| `is_product_clone` | Whether clonotype was present in manufactured Product |
| `clone_rank_in_product` | Rank by frequency in the Product sample |

### ELISpot Data

| Field | Description |
|:---|:---|
| `antigen` | WT1, PRAME, Survivin, NY-ESO-1, MAGE-A3 |
| `cytokine` | IFN-γ, IL-2, Granzyme B |
| `spots_per_1e6` | Raw spot count per 10⁶ cells |
| `spots_corrected` | Background-subtracted spot count |

---

## 🛠️ Methods & Tools

- **R** with tidyverse, edgeR, immunarch, vegan
- **Diversity metrics:** Shannon entropy, clonality index (1 − H/log N), Gini-Simpson, D50, richness
- **Differential abundance:** edgeR negative binomial GLM with fixed BCV for unreplicated designs
- **Clonal kinetics:** Trapezoidal AUC, persistence ratio, peak timing analysis
- **Functional assays:** ELISpot (IFN-γ, IL-2, Granzyme B) against tumor-associated antigens
- **Multi-omic integration:** Z-score normalized composite scoring across TCR + ELISpot modalities
- **Statistical testing:** Paired Wilcoxon, Kruskal-Wallis, Spearman correlation
- **Visualization:** ggplot2 with custom dark theme, Canvas-based interactive charts in HTML reports

---

## 📄 License

This project is released under the MIT License.

---

## 👤 Author

**Joshua Luthy** — Computational Immunologist & Bioinformatics Scientist

- 🌐 [jluthy.github.io](https://jluthy.github.io)
- 📧 jluthy123@gmail.com
- 🐙 [github.com/jluthy](https://github.com/jluthy)
