# 🧬 TCR Diversity Analysis: Apheresis to Drug Product

A reproducible analysis pipeline demonstrating **T-cell receptor (TCR) repertoire diversity changes** during cell therapy manufacturing. This project examines how the starting material (apheresis PBMCs) compares to the final drug product in terms of clonal diversity, expansion patterns, and clinical response associations.

> **Note:** All data in this repository is **synthetic** and generated for portfolio demonstration purposes. CDR3 sequences are randomly generated and do not correspond to real patient data.

---

## 📊 Analysis Reports

Pre-rendered HTML reports are available for immediate viewing:

| Notebook | Description | View |
|:---|:---|:---|
| **01 — Data Overview** | Dataset characterization, sample summary, V gene usage, quality checks | [View Report](reports/01-jal-data-overview.html) |
| **02 — Diversity Analysis** | Shannon entropy, clonality, richness, D50 — Apheresis vs Product comparison | [View Report](reports/02-jal-diversity-analysis.html) |
| **03 — Clonal Expansion** | Top expanded clonotypes, cumulative dominance curves, shared clone tracking | [View Report](reports/03-jal-clonal-expansion.html) |

---

## 🔬 Key Findings

1. **Manufacturing reduces diversity by 45–60%** — Shannon entropy drops significantly from Apheresis to Product (paired Wilcoxon p = 0.031).
2. **Clonality increases 2–3×** in the Product, reflecting oligoclonal expansion during manufacturing.
3. **D50 collapses to single digits** in Product samples, indicating extreme dominance by a few clonotypes.
4. **Complete responders (CR) show the most focused Products** — lower diversity and higher clonality compared to PD patients, suggesting efficient clonal selection may be associated with clinical benefit.

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
│
├── data/
│   ├── raw/                          # Synthetic TCR clonotype data
│   │   ├── synthetic_tcr_repertoires.csv
│   │   └── patient_metadata.csv
│   └── processed/                    # Computed metrics
│       ├── diversity_metrics.csv
│       ├── top_clonotypes.csv
│       └── vgene_usage.csv
│
├── src/R/
│   ├── generate_synthetic_data.R     # Data generation script
│   ├── diversity_metrics.R           # Shannon, clonality, D50, Gini-Simpson
│   └── visualization_helpers.R       # Custom ggplot2 theme & plot functions
│
├── notebooks/
│   ├── 01-jal-data-overview.Rmd
│   ├── 02-jal-diversity-analysis.Rmd
│   └── 03-jal-clonal-expansion.Rmd
│
├── reports/                          # Pre-rendered HTML reports
│   ├── 01-jal-data-overview.html
│   ├── 02-jal-diversity-analysis.html
│   └── 03-jal-clonal-expansion.html
│
└── www/
    └── report-style.css              # Custom report styling
```

---

## 🧪 Synthetic Data Description

| Field | Description |
|:---|:---|
| `patient_id` | 6 patients (PT-001 through PT-006) |
| `clinical_response` | CR (n=2), PR (n=2), PD (n=2) |
| `sample_type` | Apheresis (starting material) or Product (drug product) |
| `cdr3_aa` | Randomly generated CDR3β amino acid sequences |
| `v_gene` / `j_gene` | TRBV and TRBJ gene segments (IMGT nomenclature) |
| `clone_count` | Read counts following power-law distribution |
| `clone_fraction` | Proportion within sample |

The data is designed to recapitulate known biology:
- **Apheresis:** High diversity (8,000–15,000 unique clonotypes, Shannon H > 7)
- **Product:** Low diversity (800–3,000 clonotypes, Shannon H 3–5), with CR patients showing the most focused repertoires

---

## 🛠️ Methods & Tools

- **R** with tidyverse, immunarch, vegan, edgeR
- **Diversity metrics:** Shannon entropy, clonality index (1 − H/log N), Gini-Simpson, D50, richness
- **Statistical testing:** Paired Wilcoxon signed-rank test
- **Visualization:** ggplot2 with custom theme, Canvas-based interactive charts in HTML reports

---

## 📄 License

This project is released under the MIT License.

---

## 👤 Author

**Joshua Luthy** — Computational Immunologist & Bioinformatics Scientist

- 🌐 [jluthy.github.io](https://jluthy.github.io)
- 📧 jluthy123@gmail.com
- 🐙 [github.com/jluthy](https://github.com/jluthy)
