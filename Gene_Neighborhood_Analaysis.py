from Bio import SeqIO
from Bio.Blast import NCBIXML
import re
from Bio import Entrez

# Enter Entrez email
Entrez.email = ''
Entrez.tool = 'Biopython'


# NEIGHBORHOOD ANALYSIS V2, CFB 11/25/23

# collects protein accession numbers from a blast results XML that are within a specified expectation value
# uses these protein accession numbers to obtain nucleotide entries associated with the identical protein group
# searches CDSs from associated nucleotide entries for gene neighbors upstream and downstream the original protein accession
# writes the gene neighborhood for each blast result and other blast hit information in an output .TSV file

neighbor_number = int(input('Define number of neighbors: '))
accession_string = ''

numbers_list = []
for i in range(0, (neighbor_number + 1)):
    numbers_list.append(str(i))
    numbers_list.append(str(-i))
numbers_list = list(set(numbers_list))

for i in range(0, len(numbers_list)):
    numbers_list[i] = int(numbers_list[i])

res1 = []
res2 = []
for i in numbers_list:
    if(i > 0):
        res1.append(i)
    else:
        res2.append(i)
res1.sort()
res2.sort()

numbers_list = (res2 + res1)



E_VALUE_THRESH = 1
blast_results_list = []
results_accession_numbers = []

# opens an XML file from a blast output and collects protein accession numbers for hits
for record in NCBIXML.parse(open('/Users/christianbaca/Desktop/Python/brig1.xml')):
    blast_results_list.append(record)


for align in blast_results_list[0].alignments:
    for hsp in align.hsps:
        if hsp.expect < E_VALUE_THRESH:
            if 'pdb|' in align.title:
                continue
            x = re.search("\|.*?\|", str(align.title))
            accession_string = x.group()[1:-1]
            results_accession_numbers.append(accession_string)


# writes neighborhood results file with first row containing blast hit neighbor position
# followed by all the neighbors upstream and downstream
with open ('neighborhood_results.tsv', 'w') as results_file:
    for i in numbers_list:
        results_file.write(str(i) + '\t')
    results_file.write('\n')
    for accession_number in results_accession_numbers:
        new_handle = Entrez.efetch(db="ipg", id=(accession_number)[:-2], retmode="xml")

        read_file = Entrez.read(new_handle)

        larger_dict = read_file['IPGReport']

        smaller_dict = ((read_file['IPGReport'])['ProteinList'])
        smaller_list = list(smaller_dict)
        dict = smaller_list[0]
        organism_nucleotide_info_string = dict['CDSList']
        match = re.search(r"'accver': '[\w.]+'", str(organism_nucleotide_info_string))
        rough_nucleotide_accession = match.group()
        split_str = rough_nucleotide_accession.split(sep=':')
        cleaner_str = split_str[1]
        final_nucleotide_accession = cleaner_str.replace("'", "")

        count = 0

        master_index_list = []
        protein_name_list = []
        accessions_list = []
        protein_name = ''
        accession_name = ''

        new_nucleotide_handle = Entrez.efetch(db="nucleotide", id=final_nucleotide_accession, rettype='gb',retmode="text", seq_start = 1)
        results = (SeqIO.read(new_nucleotide_handle, 'genbank'))
        for feature in results.features:
            print(feature.type)
            if feature.type == "CDS" or feature.type == 'tRNA':
                count += 1
                master_index_list.append(count)
                if 'protein_id' in feature.qualifiers:
                    accessions_list.append((feature.qualifiers)['protein_id'])
                else:
                    accessions_list.append('no protein id')
                if 'product' in feature.qualifiers:
                    protein_name_list.append((feature.qualifiers)['product'])
                else:
                    protein_name_list.append('no description')

        for i in range(0, len(accessions_list)):
            accessions_list[i] = str((accessions_list[i])[0])

        for i in range(0, len(protein_name_list)):
            protein_name_list[i] = str((protein_name_list[i])[0])

        if accession_number in accessions_list:
            blast_hit_index = accessions_list.index(accession_number)

        final_neighbor_accession_list = []
        final_neighbor_description_list = []

        for i in range((blast_hit_index - (neighbor_number)), blast_hit_index):
            if i not in range(len(accessions_list)):
                final_neighbor_accession_list.append('no neighbor')
            else:
                final_neighbor_accession_list.append(accessions_list[i])
            if i not in range(len(protein_name_list)):
                final_neighbor_description_list.append('no neighbor')
            else:
                final_neighbor_description_list.append(protein_name_list[i])

        for i in range(blast_hit_index, (blast_hit_index + (neighbor_number + 1))):
            if i not in range(len(accessions_list)):
                final_neighbor_accession_list.append('no neighbor')
            else:
                final_neighbor_accession_list.append(accessions_list[i])
            if i not in range(len(protein_name_list)):
                final_neighbor_description_list.append('no neighbor')
            else:
                final_neighbor_description_list.append(protein_name_list[i])

        for neighbor in final_neighbor_description_list:
            results_file.write(neighbor + '\t')
        results_file.write('\n')

# writes a file containing information on the blast hits with the first row specifying the information
# for each column followed by the information for each of the blast hits

with open('blast_hits_info.tsv', 'w') as info_file:
    info_file.write(f'Organism\tBlast Hit Description\tBlast Hit Accession\tBlast Hit Percent Identity\tBlast Hit Amino Acid Length\n')
    info_list = []
    for align in blast_results_list[0].alignments:
        for hsp in align.hsps:
            if hsp.expect < E_VALUE_THRESH:
                if 'pdb|' in align.title:
                    continue

                percent_id_str = str(100 * (hsp.identities / hsp.align_length))

                full_id = (align.hit_id)

                match = re.search("\|.*?\|", full_id)
                match = match.group()[1:-1]


                full_definition = (align.hit_def)
                new_match = re.search("\[.*?\]", full_definition)
                new_match = new_match.group()[1:-1]


                amino_acid_length = (align.length)
                # prints the blast hit protein description
                def_match = re.search("^[\s\w\:\,\&\;\.\]\/\-\(\)]+\[", full_definition)
                protein_definition = def_match.group()[:-1]
                info_list = [new_match, protein_definition, match, percent_id_str, amino_acid_length]

                for i in info_list:
                    info_file.write(str(i) + '\t')
                info_file.write('\n')

# writes a file containing both the blast hit information and all the neighbors in one file

blast_hits_info_list = []
neighborhood_results_list = []
with open('Brig1_complete_neighborhood_results.tsv', 'w') as output_file:
    with open('blast_hits_info.tsv', 'r') as file:
        with open('neighborhood_results.tsv', 'r') as second_file:
            for line in file:
                blast_hits_info_list.append(line.replace('\n', ''))

            for new_line in second_file:
                neighborhood_results_list.append(new_line.replace('\n', ''))

            neighborhood_results_list[0] = '\t' + neighborhood_results_list[0]

    output_file.write('Blast Hit #\t')


    for i in range(0, len(neighborhood_results_list)):
        if i == 0:
            output_file.write(blast_hits_info_list[i] + neighborhood_results_list[i] + '\n')
        else:
            output_file.write(str(i) + '\t' + blast_hits_info_list[i] + neighborhood_results_list[i] + '\n')






















