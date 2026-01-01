#!/usr/bin/env bash

## ===============================
## Usage:
##   sh get_genome_length_from_fasta.sh *.fa > genome_stats.txt
##
## Output (tab-delimited):
##   TaxID  TotalBases  Accession  Description
## ===============================

for fasta in "$@"; do

    ## -------------------------------
    ## Extract header (first FASTA line)
    ## -------------------------------
    header=$(grep -m 1 "^>" "$fasta")

    ## Remove leading ">"
    header=${header#>}

    ## Genome name: before first "|"
    accession=$(echo "$header" | cut -d'|' -f1)

    ## TaxID: extract number after kraken:taxid|
    taxid=$(echo "$header" | sed -n 's/.*kraken:taxid|\([0-9]\+\).*/\1/p')

    ## Description: everything after TaxID
    description=$(echo "$header" | sed -n 's/.*kraken:taxid|[0-9]\+[[:space:]]\+//p')

    ## -------------------------------
    ## Count total bases
    ## -------------------------------
    total_bases=$(awk '
        BEGIN { sum = 0 }
        !/^>/ {
            gsub(/[^ACGTNacgtn]/, "", $0)
            sum += length($0)
        }
        END { print sum }
    ' "$fasta")

    ## -------------------------------
    ## Output row
    ## -------------------------------
    echo -e "${taxid}\t${total_bases}\t${accession}\t${description}"

done
