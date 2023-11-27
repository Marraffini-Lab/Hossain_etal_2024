from Bio.Seq import Seq
from Bio import SeqIO
import pandas as pd
import matplotlib.pyplot as plt


# Function takes as input a fasta file of a phage genome
# and returns a string of only the DNA sequence of the phage
def phage_genome(filename):
    # Open fasta file as read-only
    f = open(filename, "r")

    # Read the file into a string
    # and eliminate new line characters
    # Creates a list of DNA sequences for each line of the fasta file
    s = f.read().split("\n")

    # Join the list entries into a single string
    # Converts the entire phage genomic DNA sequence into a single string
    genome_string = "".join(s[1:len(s)])

    # Close the file
    f.close()

    return genome_string


# Function takes as input a list of fastq files containing next-gen sequencing reads
# and returns corresponding fastq files containing only high-quality reads
# The name of the output fastq files are provided by the user as a list in the function inputs
# A high-quality read is one in which all bases have a quality score above a certain threshold
# The quality-score threshold is set by the program user
# In general, a phred quality score greater than 10 is considered good quality
def fastq_quality_filter(input_files, output_files, quality_score_threshold):
    # Loop through each fastq file
    for i in range(0, len(input_files)):
        # Parse the original fastq file using BioPython SeqIO
        input_seq_iterator = SeqIO.parse(input_files[i], "fastq")

        # Use a generator expression to avoid creating the entire list in memory
        # Only commit a read to the generator expression if its minimum base quality >= threshold
        good_reads = (record for record in input_seq_iterator
                      if min(record.letter_annotations["phred_quality"]) >= quality_score_threshold)

        # Use the SeqIO writer to write a fastq file containing only the high quality reads
        SeqIO.write(good_reads, output_files[i], "fastq")


# Function returns a list of all the sequence reads in a fastq file
def total_reads(fastq_file):
    # Open fastq file as read-only
    f = open(fastq_file, "r")

    # Create a list of sequence reads from all the entries in the fastq file
    # First, read the entire fastq file into a string
    # Split the string into a list
    # where each list entry is a line in the fastq file
    # The second line of the fastq file contains the first sequence read
    # and every fourth line from that point on contains a subsequent sequence read
    fastq_seq_list = f.read().split("\n")[1::4]

    # Close the file
    f.close()

    return fastq_seq_list


# Function takes as input a fastq file for deep-seq reads
# and a fasta file for the reference phage genome
# and returns a list containing fastq reads that align to the reference genome
# Each entry in the output list is a tuple
# Each tuple contains an aligned read, its start position in the phage genome and end position
# Function also writes a csv file with the list of aligned reads
# and their corresponding start and end positions in the phage genome
# Name of the output csv file is specified by the user as the third input in the function
def aligner(fastq_file, genome_fasta_file, output_csv_file):
    # Convert the reference phage genome fasta file into a DNA sequence string
    genome = phage_genome(genome_fasta_file)

    # Get the list of sequence reads from the fastq file
    fastq_reads = total_reads(fastq_file)

    # Initiate an empty list for aligned fastq reads
    matches = []

    # Iterate through the list of sequence reads
    for entry in fastq_reads:

        # Search through the reference phage genome
        # and if the first 15 bp of the read is found in the genome
        # then the read aligns to the phage genome
        # Record the position in the phage genome of the first base of the read
        start_base = genome.find(entry[0:15])

        # Record the position in the phage genome of the last base of the read
        end_base = start_base + len(entry)

        # If the recorded end position is greater than the length of the phage genome
        # record the end position as the final base position of the phage genome
        # (This is a minor caveat of this program since it does not
        # take into account circularization of the phage genome during replication)
        # (So the first and final 75 bp positions in the final per base read counts
        # of the phage genome will not be accurate)
        if end_base >= len(genome):
            end_base = len(genome) - 1

        # If the read aligns, add the read to the list of aligned reads
        # Record the sequence, genome start position and end position of the read
        if start_base != -1:
            matches.append([entry, start_base, end_base])

        # If the read fails to align
        # check the reverse complement of the read to see whether that aligns
        else:

            # Similar logic here as above
            # First generate the reverse complement
            # and check whether the read now aligns
            entry_seq = Seq(entry)
            entry_rev_complement = str(entry_seq.reverse_complement())
            start_base_rev = genome.find(entry_rev_complement[0:15])
            end_base_rev = start_base_rev + len(entry_rev_complement)
            if end_base_rev >= len(genome):
                end_base_rev = len(genome) - 1

            # If the reverse complement aligns
            # add the reverse complement read to the list of aligned reads
            if start_base_rev != -1:
                matches.append([entry_rev_complement, start_base_rev, end_base_rev])

    # Convert this list of tuples to a pandas dataframe
    df_phage_reads = pd.DataFrame(matches, columns=["Phage Read", "Start Position (bp)", "End Position (bp)"])

    # Write the dataframe to a csv file
    df_phage_reads.to_csv(output_csv_file)

    # Return the list of phage reads
    return matches


# Function takes as input a fastq file for deep-seq reads
# and a fasta file for the reference phage genome
# and returns a dictionary of the normalized number of reads per base
# for every bp position in the reference phage genome
# Function also writes a csv file with the number of reads per base
# for every bp position in the reference phage genome
# Name of the output csv file is specified by the user as the third input in the function without .csv
# Since this function calls on the aligner function you get two csv files
# One csv file with the list of aligned reads and one with the read counts dictionary
def read_counter(fastq_file, genome_fasta_file, output_csv_file):
    # Call the aligner function to get the list of aligned phage reads
    phage_reads = aligner(fastq_file, genome_fasta_file, output_csv_file + "_reads.csv")

    # Call the total reads function to get the total number of reads
    number_reads = float(len(total_reads(fastq_file)))

    # Convert the reference phage genome fasta file into a DNA sequence string
    genome = phage_genome(genome_fasta_file)

    # Initiate an empty dictionary for the number of reads per genome position
    read_counts_dictionary = {}

    # Populate the dictionary with null values (i.e. reads = 0)
    # for every bp in the phage genome
    # This ensures that if a certain base pair is not covered by any read,
    # that bp position exists in the dictionary and is assigned read count = 0
    for i in range(0, len(genome)):
        read_counts_dictionary[i] = 0

    # Iterate through each tuple in the list of aligned phage reads
    for read in phage_reads:

        # Extract the start position and end position for each read
        start_bp = read[1]
        end_bp = read[2]

        # Iterate through values from the start position to the end position
        for i in range(start_bp, end_bp + 1):
            # For each genome position covered, increment the read count by 1
            read_counts_dictionary[i] += (1 / number_reads) * 1000000

    # Convert dictionary to a pandas dataframe
    df_read_counts = pd.DataFrame(read_counts_dictionary.items(), columns=["Genome Position (bp)", "Read Count"])

    # Write the dataframe to a csv file
    df_read_counts.to_csv(output_csv_file + "_read_counts.csv")

    # Return the dictionary
    return read_counts_dictionary


# Function takes as input a list of fastq files for deep-seq reads
# and a fasta file for the reference phage genome
# and returns a png file with a plot of the number of reads per base
# for every bp position in the reference phage genome for each input fastq file
# Name of the output png file is specified by the user as the third input in the function
# The user also has to specify the plot title, and lists of marker colors and legend labels
# Since this function calls on the read_counter function you also get two csv files per input fastq file
# One csv file with the list of aligned reads and one with the read counts dictionary
def count_plotter(fastq_files, genome_fasta_file, plot_filename, plot_title, marker_color, legend_label):
    # Loop through each fastq file
    for i in range(0, len(fastq_files)):

        # Get the read counts per base dictionary for the fastq file
        read_counts_dictionary = read_counter(fastq_files[i], genome_fasta_file, legend_label[i])

        # Initiate empty lists for the phage genome position and the read count per base
        genome_position_list = []
        read_count_list = []

        # Loop through the read counts dictionary and add to both lists
        for base in read_counts_dictionary:
            genome_position_list.append(base)
            read_count_list.append(read_counts_dictionary[base])

        # Plot a scatter plot of the read count for each bp in the phage genome
        # The marker color is specified by the user in the function inputs
        plt.scatter(genome_position_list, read_count_list,
                    marker="o", c=marker_color[i], alpha=1, edgecolor="face", s=3)

    # Set the range for the x and y axes
    plt.axis([-10, 169000, 0, 60])

    # Add x and y labels
    plt.xlabel("Genome Position (bp)")
    plt.ylabel("Normalized Read Count")

    # Add figure legend
    # The list of legend labels is specified by the user in the function inputs
    plt.legend(legend_label)

    # Add the title of the plot
    # The title is specified by the user in the function inputs
    plt.title(plot_title)

    # Save the plot
    plt.savefig(plot_filename)
