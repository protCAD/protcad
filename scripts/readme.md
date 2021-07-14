# Readme


### What are these scripts?

These scripts were created in order to help streamline one's protCAD workflow. Thus far,
we have scripts that automatically make "in files" for protMutator and that completely
reformat PDBs (e.g., renumbering atoms, residues, and relabeling chain identifiers).

### Example usage

These scripts are meant to be executed via the command line. Below is an example of what 
a command would like on the command line for each program:

#### *protMutator_template_maker.py*
As the title suggests, this script creates a template file ready to use with protMutator on the command line. If you have your sequence of interest, the pdb file that you want to thread the sequence onto, you can execute the below file to obtain an in file for protMutator automatically:

```Terminal
python3 new_comma_delim.py -i seq_of_interest.fasta -pdb structure_to_thread_seq_onto.pdb -ch 0  
```

The output will be "mut.in", which you can feed into protMutator.


#### *protPDBReformatter.py*

Similarly, as the name suggests, this file is responsible for reformatting the input PDB  provided. Provided you have your PDB, you can reformat it by executing the following on the command line: 

```Terminal
python3 protPDBReformatter.py -i PDB_you_want_to_reformat.pdb -o Name_of_new_PDB.pdb 
```

The output will be a new PDB file that will take whatever name you specify after "-o". You must add the ".pdb" file ending as well!