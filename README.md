CRISPR Target Finder
Video Demo: https://www.youtube.com/watch?v=m_TxH70sfSU
Description:
One of the most innovative technologies in contemporary biotechnology is CRISPR-Cas9, which allows researchers to modify genomes with previously unheard-of accuracy. Finding appropriate target sites, specifically 20-nucleotide guide RNAs (gRNAs) followed by the NGG PAM sequence, is a crucial step in creating CRISPR experiments.

By identifying every possible CRISPR target site by scanning a DNA sequence or FASTA file, this Python application automates that procedure. The gRNA sequence, PAM, start/end positions, and GC content—a crucial factor in CRISPR efficiency—are reported for every site.

Additionally, the program creates a visual depiction of the locations of PAM sites along the DNA sequence, which facilitates the analysis of distribution patterns.

This project integrates the fundamental ideas of bioinformatics, data analysis, and scientific visualization into a single, integrated tool by combining Biopython, regex, pandas, and matplotlib.

Features:
Two Input Choices: Directly enter a raw DNA sequence. Your sequence of interest should be included in a FASTA file.

Identification of CRISPR Sites: looks for legitimate 20 bp gRNAs followed by NGG PAM sites using regex.

Calculation of GC Content: calculates the GC% for every gRNA, which is crucial for the stability and effectiveness of Cas9.

Exporting Data: Crispr_targets.csv is a structured CSV file that is saved for further analysis.

Visualization: plots the locations of the PAM sites throughout the sequence and saves the image (crispr_targets_plot.png).

Handling Errors: gracefully handles invalid inputs, missing files, and empty results.
