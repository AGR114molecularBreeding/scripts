### A collection of scripts that we use regularly.  

  * [tidy_names](tidy_names/README.md) : Read a FASTA file containing isoforms, where the headers (previously taken from a GFF file) include the IDs for the protein, mRNA, and gene. The script generates a new FASTA file with cleaned headers, where only the gene ID is shown. In cases where there are isoforms, only the longest sequence is included.
  * [family_expansion](family_expansion/analysis.R) : Read a FASTA file containing CDS, tidy the names for pairwise alignments (e.g., BLAST), and identify family members at different identity and coverage thresholds.
* [read_coverage](read_coverage/README.md): Compute the read coverage at each position in the genome based on aligned reads (from a .BAM file). Plot the read coverage and gene annotations within a specified genomic region of interest. Obtains the read counts (i.e., the number of reads assigned to each gene) for all genes in the provided GFF file.
