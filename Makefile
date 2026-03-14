# ============================================================================
# Makefile — TCR Repertoire Diversity Analysis
# Author: Joshua Luthy
# ============================================================================
#
# Usage:
#   make data       Generate synthetic TCR dataset
#   make report     Knit all R Markdown notebooks to HTML
#   make all        Run full pipeline (data → reports)
#   make clean      Remove generated outputs
#
# Prerequisites:
#   R >= 4.5.0, renv::restore() completed
# ============================================================================

.PHONY: all data report clean

RSCRIPT := Rscript

# ---------- Data generation ----------
DATA_RAW    := data/raw/synthetic_tcr_repertoires.csv data/raw/patient_metadata.csv
DATA_PROC   := data/processed/diversity_metrics.csv data/processed/top_clonotypes.csv

data: $(DATA_RAW)

$(DATA_RAW): src/R/generate_synthetic_data.R
	@echo ">> Generating synthetic TCR data..."
	$(RSCRIPT) src/R/generate_synthetic_data.R

# ---------- Report generation ----------
NOTEBOOKS := $(wildcard notebooks/*.Rmd)
REPORTS   := $(patsubst notebooks/%.Rmd,reports/%.html,$(NOTEBOOKS))

report: $(REPORTS)

reports/%.html: notebooks/%.Rmd $(DATA_RAW) $(DATA_PROC) www/report-style.css
	@echo ">> Knitting $<..."
	$(RSCRIPT) -e "rmarkdown::render('$<', output_dir = 'reports/')"

# ---------- Full pipeline ----------
all: data report
	@echo ">> Pipeline complete."

# ---------- Clean ----------
clean:
	@echo ">> Cleaning generated files..."
	rm -f data/raw/*.csv
	rm -f data/processed/*.csv
	rm -f reports/*.html
	rm -f reports/figures/*
