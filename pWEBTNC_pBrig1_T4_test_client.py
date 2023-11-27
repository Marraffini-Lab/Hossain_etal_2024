import NGS_phage_DNA_infected_cells

phage = NGS_phage_DNA_infected_cells.phage_genome("T4_WT.fa")
print(phage)

# Generate fastq files containing only good quality reads for the two input fastq files
NGS_phage_DNA_infected_cells.fastq_quality_filter(["S1_pWEBTNC_T4_L001_R1_001.fastq", "S2_pBrig1_T4_L001_R1_001.fastq"],
                                                  ["good_reads_S1_fwd.fastq", "good_reads_S2_fwd.fastq"], 10)

print("Getting there...sorry it's slow")

# Produce a plot of the normalized read counts per base against the phage genome position
# For pWEB-TNC (control) vs pBrig1 (targeting)
# For each input fastq file, you also get two csv files
# One csv file with the list of aligned reads and one with the read counts dictionary
NGS_phage_DNA_infected_cells.count_plotter(["good_reads_S1_fwd.fastq", "good_reads_S2_fwd.fastq"], "T4_WT.fa",
                                           "Brig1_targeting_T4_8mins.png",
                                           "Brig1 targeting - T4 8mins post infection",
                                           ["red", "blue"], ["pWEB-TNC", "pBrig1"])

print("Done!")
print("Thank you!")
