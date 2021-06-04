#!/usr/bin/env python3
# Author: Jan Siess
# Date: 8 March 2019
# Update 1: 9 February 2021
# Purpose: To take in a FASTA sequence and return a comma-delimeted FASTA sequence as output for further use in protMutator. The output sequence will be automatically copied to the clipboard.


import os
import argparse

# ------------------------------------------------------------------#
#                               Functions                           #
# ------------------------------------------------------------------#


def read_aligned_fasta(input_file):
    global header, seq
    header = ""
    seq = ""
    fh = open(input_file, "r")
    for line in fh:
        if line.startswith(">"):
            header = line
        else:
            line = str.strip(line)
            seq = seq + line
    return seq


def add_commas(seq):
	placeholder = []
	for character in seq:
		if character == "H":
			placeholder.append("He")
		else:
			placeholder.append(character)
	
	desired = ",".join(placeholder)
	return desired


# ------------------------------------------------------------------#
#                               Main                                #
# ------------------------------------------------------------------#

if __name__ == "__main__":

    parser = argparse.ArgumentParser(prog="protMutator_template_maker.py")
    parser.add_argument(
        "-i",
        "--input",
        dest="input",
        help="The input file, in fasta format, to be processed.",
    )
    parser.add_argument(
        "-pdb",
        dest="pdb",
        help="Name of the pdb file that the input sequence will be threaded onto.",
    )
    parser.add_argument(
        "-ch",
        "--active-chains",
        dest="chains",
        help="Number of active chains in the sequences.",
    )

    options = parser.parse_args()

    # Reading in and adding commas to the input fasta
    desired_seq = read_aligned_fasta(f"{options.input}")
    seq_commas = add_commas(desired_seq)

    # Getting number of active positions (default is all of the positions)
    active_pos = [i for i in range(0, len(desired_seq))]
    active_pos = list(map(lambda x: str(x), active_pos))
    active_pos = ",".join(active_pos)

	# Creating the input file with everything provided foir use in protMutator
    f = open("mut.in", "w+")
    f.write(f"Input PDB File,{options.pdb},\n")
    f.write(f"Active Chains,{options.chains},\n")
    f.write(f"Active Positions,{active_pos},\n")
    f.write(f"{seq_commas},\n")
    f.close()
	
	
	