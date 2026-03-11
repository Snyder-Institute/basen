if (base::getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    "assigned_kmers", "assigned_perc", "bases_assigned", "coverage_proxy",
    "genome_size_bp", "kmer_string", "kmers_diversity", "name",
    "rank_code", "read_id", "reads_assigned", "relative_abundance",
    "status", "taxid", "total_kmers", "unassigned_kmers", "unassigned_perc",
    "..required_cols"
  ))
}