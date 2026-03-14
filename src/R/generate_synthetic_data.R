# ============================================================================
# generate_synthetic_data.R
# Generates synthetic TCR repertoire data for portfolio demonstration
# Author: Joshua Luthy
# ============================================================================

library(tidyverse)

set.seed(42)

# --- Configuration -----------------------------------------------------------
n_patients     <- 6
sample_types   <- c("Apheresis", "Product")
clinical_resp  <- c("CR", "CR", "PR", "PR", "PD", "PD")  
patient_ids    <- paste0("PT-", sprintf("%03d", 1:n_patients))

# Realistic TRBV / TRBJ gene segments
trbv_genes <- c(
  "TRBV5-1", "TRBV6-5", "TRBV7-2", "TRBV7-9", "TRBV9", "TRBV10-3",
  "TRBV11-2", "TRBV12-3", "TRBV12-4", "TRBV13", "TRBV14", "TRBV15",
  "TRBV18", "TRBV19", "TRBV20-1", "TRBV24-1", "TRBV25-1", "TRBV27",
  "TRBV28", "TRBV29-1", "TRBV30"
)

trbj_genes <- c(
  "TRBJ1-1", "TRBJ1-2", "TRBJ1-3", "TRBJ1-4", "TRBJ1-5", "TRBJ1-6",
  "TRBJ2-1", "TRBJ2-2", "TRBJ2-3", "TRBJ2-4", "TRBJ2-5", "TRBJ2-7"
)

# Amino acid alphabet for CDR3 sequences
aa <- c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")

# --- Helper functions --------------------------------------------------------

generate_cdr3 <- function(n = 1, len_range = c(12, 18)) {
  sapply(1:n, function(x) {
    len <- sample(len_range[1]:len_range[2], 1)
    # CDR3 typically starts with C and ends with F
    paste0("C", paste0(sample(aa, len - 2, replace = TRUE), collapse = ""), "F")
  })
}

generate_repertoire <- function(n_clonotypes, total_reads, skew = 1.5) {
  # Power-law-ish distribution for clone sizes
  raw_freq <- (1:n_clonotypes)^(-skew) + runif(n_clonotypes, 0, 0.001)
  freq     <- raw_freq / sum(raw_freq)
  counts   <- as.integer(round(freq * total_reads))
  counts[counts < 1] <- 1
  return(counts)
}

# --- Generate data -----------------------------------------------------------

all_data <- list()

for (i in seq_along(patient_ids)) {
  pid  <- patient_ids[i]
  resp <- clinical_resp[i]
  
  # Apheresis: high diversity (many clonotypes, even distribution)
  n_aph_clones   <- sample(8000:15000, 1)
  aph_reads      <- sample(80000:150000, 1)
  aph_counts     <- generate_repertoire(n_aph_clones, aph_reads, skew = 1.2)
  
  # Product: low diversity (fewer clonotypes, skewed toward dominant clones)
  # CR patients have more focused products; PD patients are more diffuse
  product_skew <- case_when(
    resp == "CR" ~ 2.5,
    resp == "PR" ~ 2.0,
    resp == "PD" ~ 1.5
  )
  n_prod_clones  <- sample(800:3000, 1)
  prod_reads     <- sample(50000:100000, 1)
  prod_counts    <- generate_repertoire(n_prod_clones, prod_reads, skew = product_skew)
  
  # Generate shared clonotypes (some apheresis clones appear in product)
  n_shared       <- min(n_prod_clones, sample(200:600, 1))
  shared_cdr3    <- generate_cdr3(n_shared)
  shared_v       <- sample(trbv_genes, n_shared, replace = TRUE)
  shared_j       <- sample(trbj_genes, n_shared, replace = TRUE)
  
  # Apheresis-only clonotypes
  n_aph_only     <- n_aph_clones - n_shared
  aph_only_cdr3  <- generate_cdr3(n_aph_only)
  aph_only_v     <- sample(trbv_genes, n_aph_only, replace = TRUE)
  aph_only_j     <- sample(trbj_genes, n_aph_only, replace = TRUE)
  
  # Product-only clonotypes (new expansions)
  n_prod_only    <- n_prod_clones - n_shared
  prod_only_cdr3 <- generate_cdr3(max(n_prod_only, 0))
  prod_only_v    <- sample(trbv_genes, max(n_prod_only, 0), replace = TRUE)
  prod_only_j    <- sample(trbj_genes, max(n_prod_only, 0), replace = TRUE)
  
  # Build apheresis dataframe
  aph_df <- tibble(
    patient_id       = pid,
    clinical_response = resp,
    sample_type      = "Apheresis",
    clonotype_id     = paste0(pid, "_CLN_", sprintf("%05d", 1:n_aph_clones)),
    cdr3_aa          = c(shared_cdr3, aph_only_cdr3),
    v_gene           = c(shared_v, aph_only_v),
    j_gene           = c(shared_j, aph_only_j),
    clone_count      = aph_counts[1:n_aph_clones],
    clone_fraction   = clone_count / sum(clone_count)
  )
  
  # Build product dataframe
  prod_all_cdr3 <- c(shared_cdr3, prod_only_cdr3)
  prod_all_v    <- c(shared_v, prod_only_v)
  prod_all_j    <- c(shared_j, prod_only_j)
  
  prod_df <- tibble(
    patient_id       = pid,
    clinical_response = resp,
    sample_type      = "Product",
    clonotype_id     = paste0(pid, "_CLN_", sprintf("%05d", 1:n_prod_clones)),
    cdr3_aa          = prod_all_cdr3[1:n_prod_clones],
    v_gene           = prod_all_v[1:n_prod_clones],
    j_gene           = prod_all_j[1:n_prod_clones],
    clone_count      = prod_counts[1:n_prod_clones],
    clone_fraction   = clone_count / sum(clone_count)
  )
  
  all_data[[paste0(pid, "_aph")]]  <- aph_df
  all_data[[paste0(pid, "_prod")]] <- prod_df
}

# Combine all
tcr_data <- bind_rows(all_data)

# --- Write outputs -----------------------------------------------------------
write_csv(tcr_data, "data/raw/synthetic_tcr_repertoires.csv")

# Patient metadata
patient_meta <- tibble(
  patient_id        = patient_ids,
  clinical_response = clinical_resp,
  age               = sample(35:68, n_patients),
  sex               = sample(c("M", "F"), n_patients, replace = TRUE),
  disease           = "AML",
  prior_lines       = sample(1:4, n_patients, replace = TRUE)
)

write_csv(patient_meta, "data/raw/patient_metadata.csv")

cat("Synthetic data generated successfully.\n")
cat(sprintf("  Total rows: %d\n", nrow(tcr_data)))
cat(sprintf("  Patients:   %d\n", n_patients))
cat(sprintf("  Samples:    %d\n", length(unique(paste(tcr_data$patient_id, tcr_data$sample_type)))))
