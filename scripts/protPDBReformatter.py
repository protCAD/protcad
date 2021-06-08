#!/usr/bin/env python3

"""
..module:: protPDBReformatter
..moduleauthor:: Jan Siess <jsiess93@gmail.com>, Douglass Pike <doughp11@gmail.com>

"""


import argparse
import numpy as np
import pandas as pd


def read_pdb_file(file):
    """
    This function is responsible for reading the input file. Each element on the
    line is split and placed into a list. The result is then filtered such that
    no empty spaces are present.
    
    :param file: Name of the PDB to be read in and parsed. 
    :param type: Namespace -- argparse type
    :return: A list of lists containing all of the data once held within the PDB.
    :rtype: list.

    """

    pdb_lines = list()
    with open(file, "r") as fh:
        lines = fh.read().splitlines()
        for l in lines:
            split_l = l.split(" ")
            pdb_lines.append(list(filter(None, split_l)))
    return pdb_lines


def renumber_atom_nums(pdb_dataframe):
    """
    This function simply reads in the PDB dataframe, sets all of the column 1
    elements to 0, and is quickly changed to reflect it's new index.

    :param pdb_dataframe: The Dataframe containing all of the PDB information.
    :param type: Pandas dataframe.
    :return: A new Pandas dataframe with renumbered atoms.
    :rtype: Pandas dataframe.

    """

    pdb_dataframe[1] = 0
    pdb_dataframe[1] = [i + 1 for i in range(len(pdb_dataframe))]
    return pdb_dataframe


def reassign_chain_id(pdb_dataframe):
    """
    This function is responsible for changing the chain IDs for each individual
    chain. A mask is constructed on pdb_dataframe to determine what values fall
    between the 'TER' tags in PDB. Atoms/residues falling in between these values
    represent one chain, and are renamed appropriately to reflect that particular
    membership.

    :param pdb_dataframe: The dataframe containing all of the PDB information.
    :param type: Pandas dataframe.
    :return: A new Pandas dataframe with chains relabeled.
    :rtype: Pandas dataframe.
    
    """

    ch = "A"
    pdb_dataframe[0] = pdb_dataframe[0].astype(str)
    pdb_dataframe[4] = pdb_dataframe[0].between("TER", "TER", inclusive=True)
    for i in range(len(pdb_dataframe)):
        if pdb_dataframe.iloc[i, 4] == False:
            pdb_dataframe.iloc[i, 4] = ch
        else:
            ch = chr(ord(ch) + 1)
            pdb_dataframe.iloc[i, 4] = 0
    return pdb_dataframe


def renumber_residues(pdb_dataframe):
    """
    This function is responsible for renumbering the residues within the PDB
    dataframe provided.


    :param pdb_dataframe: The dataframe containing all of the PDB information.
    :param type: Pandas dataframe.
    :return: A new Pandas dataframe with residues renumbered.
    :rtype: Pandas dataframe.

    """

    c = 0
    pdb_dataframe[2] = pdb_dataframe[2].astype(str)
    pdb_dataframe[5] = pdb_dataframe[2].between("N", "N", inclusive=True)
    for i in range(len(pdb_dataframe)):
        if pdb_dataframe.iloc[i, 5] == False:
            pdb_dataframe.iloc[i, 5] = c
        elif pdb_dataframe.iloc[i, 5] == True:
            c = c + 1
            pdb_dataframe.iloc[i, 5] = c
        if pdb_dataframe.iloc[i, 2] == "S1":
            c = c + 1
            pdb_dataframe.iloc[i, 5] = c
    return pdb_dataframe


if __name__ == "__main__":

    # Defining the arguments that will be able to be passed on the command line.
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--input",
        dest="inputf",
        help="Choosing the Input PDB File to be Processed and Reformatted.",
    )
    parser.add_argument(
        "-o",
        "--output",
        dest="outputf",
        help="Choosing the Name of the Ouput PDB File to be Output After Processing is Complete.",
    )
    option = parser.parse_args()
    pd.set_option("display.max_rows", None, "display.max_columns", None)


    # Feeding the dataframe into each function to renumber atoms, chains, and residues.
    pdb_li = read_pdb_file(f"{option.inputf}")
    pdb_df = pd.DataFrame(pdb_li).fillna(0)
    pdb_df = renumber_atom_nums(pdb_df)
    pdb_df = reassign_chain_id(pdb_df)
    pdb_df = renumber_residues(pdb_df).values.tolist()



    # Creating a new file and writing out the contents of the newly formed PDB list.
    f = open(f"{option.outputf}", "w+")
    for i in pdb_df:
        if i[0] == "ATOM":
            f.write(
                f"{i[0]:<6s}{i[1]:>5g} {i[2]:^4s} {i[3]:<1s} {i[4]:>1s}{i[5]:>4}    {float(i[6]):>8.3f}{float(i[7]):>8.3f}{float(i[8]):>8.3f}{float(i[9]):>6.2f} {float(i[10]):<6.2f}{i[11]:>12s}\n"
            )
        elif i[0] == "HETATM":
            f.write(
                f"{i[0]:<4s}{i[1]:>5g} {i[2]:^4s} {i[3]:<1s} {i[4]:>1s}{i[5]:>4}    {float(i[6]):>8.3f}{float(i[7]):>8.3f}{float(i[8]):>8.3f}{float(i[9]):>6.2f} {float(i[10]):<6.2f}{i[11]:>12s}\n"
            )
        else:
            f.write(f"{i[0]:<6s}{i[1]:>5g}\n")
    f.close()
