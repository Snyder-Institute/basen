# BASEN

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Input-Kraken2](https://img.shields.io/badge/input-Kraken2-4B8BBE)](https://ccb.jhu.edu/software/kraken2/)
[![Application-Long--read%20metagenomics](https://img.shields.io/badge/application-long--read%20metagenomics-2E8B57)](#overview)

> A Bracken alternative for long-read metagenomic data

**B**ase-level **A**bundance estimation with **S**pecies-assigned **E**vidence using **N**anopore

![Logo](.github/imgs/BASEN_logo_black.png#gh-dark-mode-only)
![Logo](.github/imgs/BASEN_logo_white.png#gh-light-mode-only)

---

## Overview

BASEN is an R package for **genome-size-aware, species-level abundance estimation** from **Kraken2** long-read classification outputs. It was developed for long-read metagenomic workflows, particularly those based on custom reference databases, where conventional count-based summaries can be difficult to interpret because of differences in read length, sequencing depth, taxonomic assignment characteristics, and reference genome size.

Rather than relying on read counts alone, BASEN aggregates **species-assigned k-mer evidence** from Kraken2 read-level output, normalizes that signal by **effective reference genome length**, and scales the normalized values within each sample to produce **species-level relative abundance estimates**.

BASEN also provides a **read-level compositional quality-control workflow** based on k-mer assignment profiles. These diagnostics help identify ambiguous or low-information reads before downstream abundance interpretation.

## Features

BASEN provides two main components:

1. **Species-level abundance estimation**
   - parses Kraken2 `*.kraken` and `*.report` files
   - aggregates species-assigned k-mer evidence
   - normalizes by effective reference genome length
   - reports coverage proxies and relative abundance

2. **Read-level compositional quality control**
   - profiles per-read k-mer composition
   - quantifies assigned and unassigned k-mer proportions
   - computes Shannon diversity of taxonomic k-mer assignments
   - supports filtering of low-information reads before abundance estimation

## Installation

### Requirements

- R >= 4.3.0
- Dependencies:
  - `data.table`
  - `dplyr`
  - `ggplot2` (for QC plotting)

### Install from GitHub

```r
# install.packages("devtools")
devtools::install_github("Snyder-Institute/basen")
```

## Workflow summary

### Inputs

BASEN expects:

- **Kraken2 read-level output** (`*.kraken`)
- **Kraken2 report files** (`*.report`)
- **Genome statistics table** containing:
  - column 1: TaxID
  - column 2: genome size in base pairs

Genome statistics can be derived from reference FASTA files using the helper script:

```bash
sh get_genome_length_from_fasta.sh *.fa > genome_stats.txt
```

### Core workflow

  1. **Classify reads with Kraken2**  
     Run Kraken2 against a custom reference database and generate `*.kraken` and `*.report` files.

  2. **Profile read-level k-mer composition (optional but recommended)**  
     Compute per-read QC metrics and identify low-information reads.

  3. **Filter reads based on compositional QC (optional but recommended)**  
     Retain reads with acceptable unassigned-k-mer proportion and Shannon diversity.

  4. **Estimate species-level abundance**  
     - identifies species-level taxa from the Kraken2 report
     - extracts species-assigned k-mer evidence from classified reads
     - sums that evidence by species
     - normalizes by effective reference genome length
     - scales within each sample to generate relative abundance

  5. **Optional downstream normalization**  
     BASEN can also compute size factors from coverage-proxy matrices and generate normalized matrices for downstream analyses such as ordination, PCoA, clustering, or heatmaps.

## Primary outputs

### `kraken_relative_abundance()`

The main BASEN output is a long-format table with columns including:

- `sample_id`
- `name`
- `taxid`
- `genome_size_bp`
- `bases_assigned`
- `reads_assigned`
- `coverage_proxy`
- `relative_abundance`

### `profile_kmer_composition()`

The QC workflow returns read-level metrics including:

- `total_kmers`
- `assigned_kmers`
- `unassigned_kmers`
- `assigned_perc`
- `unassigned_perc`
- `kmers_diversity`

### Relative abundance and normalized matrices

The `relative_abundance` column returned by `kraken_relative_abundance()` represents the species-level relative abundance within each sample and is the primary BASEN abundance output.

BASEN also provides optional utilities for calculating size factors from `coverage_proxy` values with `calc_size_factors()` and for generating normalized matrices from those size factors. These normalized matrices are useful for downstream cross-sample analyses, such as ordination, clustering, heatmaps, or other comparative analyses across samples.

---

## Method summary

For each species $begin:math:text$s$end:math:text$, BASEN computes:

### Species-assigned k-mer evidence

```math
B_s = \sum_r k_{r,s}
```

where $begin:math:text$k\_\{r\,s\}$end:math:text$ is the number of k-mers assigned to species $begin:math:text$s$end:math:text$ in read $begin:math:text$r$end:math:text$.

### Effective genome length

```math
G_s = \text{genome\_size\_bp} - k + 1
```

where $begin:math:text$k$end:math:text$ is the Kraken2 k-mer length (35 by default).

### Coverage proxy

```math
C_s = \frac{B_s}{G_s}
```

### Relative abundance

```math
RA_s = \frac{C_s}{\sum_i C_i}
```

This approach is conceptually related to length-normalization strategies used in transcriptomics, but is adapted here for **species-level abundance estimation from Kraken2 long-read classifications**.

---

## Example

```r
library(basen)

results <- kraken_relative_abundance(
  report_dir = report_dir,
  kraken_dir = kraken_dir,
  genome_stats_file = genome_stats_file
)
```

### Current scope and planned extension

BASEN currently implements:

- **Species-level abundance estimation**
- **Read-level compositional QC**
- **Optional filtering of low-information reads**
- **Optional downstream size-factor normalization of coverage-proxy matrices**

Genus-level summarization is planned, but robust roll-up from lower-rank assignments is not yet implemented. Meaningful aggregation across taxonomic ranks requires explicit taxonomy-aware mapping, and this functionality is being reserved for future development.

---

## Documentation

The package vignette provides a full end-to-end BASEN workflow, including:

  - Preparation of Kraken2 outputs
  - Genome statistics generation
  - Read-level QC
  - Filtering
  - Abundance estimation
  - Optional downstream normalization

## Citation

If you use BASEN in your work, please cite the BASEN manuscript once available.

## License

MIT License

Copyright (c) 2026 The Bioinformatics Hub at the Snyder Institute
