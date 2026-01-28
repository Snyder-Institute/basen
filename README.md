# BASEN
> Bracken alternatives for the long-read sequencing data

**B**ase-level **A**bundance estimation with **S**pecies-assigned **E**vidence using **N**anopore

![Logo](.github/imgs/BASEN_logo_black.png#gh-dark-mode-only)
![Logo](.github/imgs/BASEN_logo_white.png#gh-light-mode-only)

---

## Overview
Comparing taxonomic abundance across metagenomic samples is challenging, particularly for long-read sequencing data, where sequencing depth, classification resolution, and database composition vary substantially. Raw read counts and sample-specific proportional summaries are not directly comparable across samples or taxonomic ranks, often introducing technical bias and obscuring true biological variation.

BASEN addresses these limitations by converting taxonomic classification outputs into base-level relative abundance, normalizing assigned bases by reference genome length. This rank-agnostic metric enables robust, genome-aware comparison of taxonomic composition across samples and taxonomic levels, improving interpretability and reproducibility in Nanopore-based metagenomic analyses.

## Features
BASEN estimates species-level relative abundance from long-read metagenomic data by:

1. Restricting classification to species-level assignments
2. Counting only species-assigned k-mers (base-level evidence)
3. Normalizing by reference genome length
4. Computing relative abundance across species


## Installation
### Requirement
 * R > 4.3.0
 * Dependencies (alphabetical order):
   * data.table
   * dplyr

### Install via GitHub
```
# install.packages("devtools")
devtools::install_github("Snyder-Institute/basen")
```

## Workflow summary
### Input
  * Long-read metagenomic sequencing data (Nanopore FASTQ)
  * Custom Kraken2 database with TaxID-annotated reference genomes
  * Reference genome FASTA files for genome size estimation

### Step 1. Taxonomic classification
  * Run Kraken2 on each sample using the custom database
  * Generate:
    * `.kraken` files (per-read taxonomic assignments and k-mer evidence)
    * `.report` files (taxonomy summaries across ranks)

### Step 2. Genome size normalization
  * Extract genome sizes (base pairs) from reference FASTA files
  * Create a genome statistics table linking TaxID to genome length

### Step 3. Base-level signal aggregation
  * Identify species-level TaxIDs from Kraken2 reports
  * For each read:
    * Extract species-specific k-mer counts
    * Sum assigned bases per species
    * Normalize assigned bases by reference genome size to compute a coverage proxy

### Step 4. Relative abundance estimation
  * Convert coverage proxies into relative abundance values by normalization within each sample
  * Produce rank-agnostic, genome-size–normalized abundance estimates

### Output
  * Per-sample species-level relative abundance tables
  * Optional exported files for parallel processing
  * Combined long-format table or wide abundance matrix for downstream analysis

### Key advantages
  * Accounts for variable read length in long-read sequencing
  * Reduces bias from genome size differences
  * Enables consistent comparison across samples and taxonomic ranks

## Summary

BASEN estimates **species-level relative abundance** by quantifying species-assigned k-mers from long-read sequencing data, normalizing by reference genome length, and scaling across detected species to obtain genome-size–corrected relative abundances.

---

## Workflows in details
### Step 1. Input data

```text
Kraken2 output
        |-- *.kraken   (read-level k-mer assignments)
        |-- *.report   (taxonomy summary)
```

Inputs:

- Long-read sequencing data with variable read length
- Kraken2 classification output against a custom reference database
- Reference genome statistics (genome length per species)

### Step 2. Select species-level taxa

From the Kraken2 report file, retain only taxa annotated at the species rank.

This step defines the set of species included in downstream abundance estimation.

### Step 3. Parse k-mer assignments per read

Each classified read contains a k-mer assignment string of the form:

```text
taxid1:count1 taxid2:count2 ... 0:count_unassigned
```

Example:

```text
1256908:79 0:41 1256908:2 0:11 1256908:69 0:1
```

Only k-mers assigned to species taxids are counted. Unassigned k-mers (`taxid = 0`) are excluded.

### Step 4. Estimate species-assigned bases

For each read *r*:

- Let $`k_r`$ be the k-mers where Kraken assigns the read to species _s_

Define:

```math
B_s = \sum_r k_{r,s}
```

Where:

- $`B_s`$ is the total number of species-assigned bases (k-mers) for species *s*

Text-based schematic:

```text
Read 1: [#][#][#][ ][ ]  -> species s
Read 2: [#][#][ ][#][#]  -> species s
Read 3: [#][ ][ ][#][ ]  -> species s
```

- $`B_s`$ = total number of [#] boxes aligned to species _s_

## Step 5. Incorporate reference genome length

For each species *s*:

- Let $G_s$  represent the effective reference genome length, defined as the total length of the reference genome minus the k-mer size used by Kraken2 (35 bp by default)

Text-based schematic:

```text
Reference genome (species s):
|--------------------------------------|
<------------  G_s bases  ------------->
```

The reference genome length can be derived directly from FASTA files using the auxiliary script `get_genome_length_from_fasta.sh`, located in the inst/scripts directory.

### Step 6. Normalize by genome length (coverage proxy)

To correct for genome size differences, compute a genome-size–normalized coverage proxy.

```math
C_s = B_s / G_s
```

Where:

- $`C_s`$ approximates the average genome coverage of species _s_

This step is conceptually analogous to length normalization in TPM-style RNA-seq metrics.

### Step 7. Compute relative abundance

Normalize the coverage proxy across all detected species.

```math
RA_s = C_s / \sum_{i=1}^n C_i
```

Where:

- $`RA_s`$ is the relative abundance of species _s_

## Final output

For each species *s*, BASEN reports:

- $B_s$  species-assigned bases (k-mers)
- $G_s$  reference genome length
- $C_s$  genome-size–normalized coverage proxy
- $RA_s$ relative abundance

## Usage
  * For a complete description of the workflow and implementation details, refer to the vignette.

```
library(basen)
resultsL <- kraken_relative_abundance(report_dir = report_dir, genome_stats_file = genome_stats_file)
```
