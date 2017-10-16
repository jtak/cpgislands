import argparse
import sys
import numpy as np
from Bio import SeqIO
from Bio import Entrez

"""
Creates a parser for command line arguments
"""
def get_parser():
    parser = argparse.ArgumentParser(description='Find CpG-islands in a nucleotide sequence. '
                                                 'Either -f or --genbank is required.')
    parser.add_argument('-f', '--file', help='Input file (FASTA)')
    parser.add_argument('--genbank', help='Sequence ID to be fetched from GenBank')
    parser.add_argument('--email', help='Email-address for Entrez. ')
    #parser.add_argument('--plot', help='Plot the locations of CpG sites', action="store_true")
    parser.add_argument('-v', '--verbose', action="store_true", help="increase output verbosity")
    parser.add_argument('-l', '--length', type=int, default=200, help="Set the minimum CpG island lenght. Default: 200")
    parser.add_argument('--cgthreshold', type=float, default=0.5,
                        help="The minimum CG-percent required in CpG-islands. Default: 0.5")
    parser.add_argument('--oethreshold', type=float, default=0.6,
                        help="The minimum observed to expected CpG-site count ratio in CpG-islands. Default: 0.6")
    return parser


def parse_args():
    return get_parser().parse_args()


args = parse_args()
verbose = args.verbose

"""
A helper function for printing additional output if -v is set
"""
def print_verbose(s):
    if(verbose):
        print(s)

"""
Read sequences from a fasta-file
"""
def read_seqs(fasta):
    seqs = []
    for seq_record in SeqIO.parse(fasta, "fasta"):
        seqs.append(seq_record)
    return seqs

"""
Fetch sequence from GenBank. 
"""
def fetch_genbank(id, email):
    if email is not None:
        Entrez.email = email
    handle = Entrez.efetch(db="nucleotide", \
                               rettype="gb", retmode="text", id=id)
    seq_record = SeqIO.read(handle, "gb")
    handle.close()
    return seq_record

"""
Creates sum lists for C- and G-bases, and CpG-sites.
"""
def create_sum_lists(seq):
    cpg_count = 0
    c_count = 0
    g_count = 0

    n_elements = len(seq)

    # A sum list of cpg-sites.
    # The i:th element is the number of cpg sites encountered up to that index.
    cpg_list = np.zeros((n_elements,), np.int)
    # Sum lists for C and G
    c_list = np.zeros((n_elements,), np.int)
    g_list = np.zeros((n_elements,), np.int)

    # Count values for sum lists
    for i in range(1, len(seq)):
        current_nuc = seq[i-1]
        next_nuc = seq[i]
        if current_nuc == 'C':
            c_count += 1
            if next_nuc == 'G':
                cpg_count += 1
        if current_nuc == 'G':
            g_count += 1
        c_list[i] = c_count
        g_list[i] = g_count
        cpg_list[i] = cpg_count


    return c_list, g_list, cpg_list


"""
Check if the given range is a CpG-island.
If verbose input is set, prints the statistics that the decision is based on.
"""
def is_cpg_island(c_list, g_list, cpg_list, start, stop, cg_perc_threshold, obs_to_exp_ratio_threshold):
    island_size = stop - start
    if island_size == 0 or stop >= len(cpg_list):
        return False

    c_count = c_list[stop] - c_list[start]
    g_count = g_list[stop] - g_list[start]
    cpg_count = cpg_list[stop] - cpg_list[start]

    c_percent = c_count / island_size
    g_percent = g_count / island_size
    cg_percent = c_percent + g_percent

    exp_cpg_count = (c_count * g_count) / island_size
    obs_to_exp_ratio = 0
    if cpg_count > 0:
     obs_to_exp_ratio= cpg_count / exp_cpg_count

    result = cg_percent >= cg_perc_threshold and obs_to_exp_ratio >= obs_to_exp_ratio_threshold

    print_verbose("Island: " + str(start) + " - " + str(stop))
    print_verbose("C: " + str(c_count))
    print_verbose("G: " + str(g_count))
    print_verbose("CG%: " + str(cg_percent))
    print_verbose("Obs: " + str(cpg_count))
    print_verbose("Exp: " + str(exp_cpg_count))
    print_verbose("Obs to Exp ratio: " + str(obs_to_exp_ratio))
    print_verbose("Result: " + str(result))

    return result

"""
Splits the sequence to segments of island_length and checks if the segments are CpG-islands.
Returns CpG-islands as a list of tuples: (start-index, end-index).

Returns a list of cpg-islands
"""
def find_cpg_islands(seq, island_length, cg_perc_threshold, obs_to_exp_ratio_threshold):
    c_sum_list, g_sum_list, cpg_sum_list = create_sum_lists(seq)

    islands = []
    start = 0
    stop = island_length
    while start < len(seq) - 1:
        if is_cpg_island(c_sum_list, g_sum_list, cpg_sum_list, start, stop, cg_perc_threshold, obs_to_exp_ratio_threshold):
            islands.append((start, stop))
        start = stop
        stop += island_length
        stop = min(len(seq) -1, stop)

    return islands

"""
Prints help message and exits the program
"""
def print_help_and_exit():
    parser = get_parser()
    parser.print_help()
    sys.exit(0)

file = args.file
genbank_id = args.genbank

# Only allow a file OR a genbank id
if file is not None and genbank_id is not None:
    print_help_and_exit()

# If neither a file nor a GenBank id are given, print help and exit
if file is None and genbank_id is None:
    print_help_and_exit()

# A list of islands found per sequence in case the fasta file is a multifasta.
# value at index i is the number of CpG-islands in sequence number i
num_islands = []
sequences = []
island_length = args.length

cg_perc_threshold=args.cgthreshold
obs_to_exp_ratio_threshold=args.oethreshold

if file:
    sequences = read_seqs(file)
else:
    sequences.append(fetch_genbank(genbank_id, args.email))

for seq_r in sequences:

    islands = find_cpg_islands(seq_r.seq, island_length, cg_perc_threshold, obs_to_exp_ratio_threshold)
    num_islands.append(len(islands))
    print("Sequence ID:", seq_r.id)
    print("Islands found:")
    for island in islands:
        print(island[0], "to", island[1])

    print("")

print("Number of islands:")
for i in range(0, len(sequences)):
    print(sequences[i].id, ":", num_islands[i])



