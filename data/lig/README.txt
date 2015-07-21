# Ligand Atom types. 
#
# When reading in ligand .lig files, all lines preceded by # are ignored. 
# Protcad will scan the first column for the atom name from the pdb file.  If a line is found
# with that atom name, protcad will use the atom types from that line.  If the atom name is not
# listed, protcad will display an error message and will assign the atom types by default from the
# first letter of the atom name. ie: C8B will be 'C' and N2C will be 'N' and so on. 
#
# A line preceded with '!' marks a ligand group definition.  Each group is marked with '! NAME
# NUM_ATOMS CONNECTIVITY'. You MUST include the Num_Atoms and 'Y' or 'N' for CONNECTIVITY
# in order for ligandTemplate.cc to work properly.  Separate by tabs.
# 
# The charges in columns 7 & 8 were calculated in insightII from the esff forcefield.
#
# column 1: atom name
# column 2: AMBER all atom atomtype
# column 3: AMBER united atom atomtype
# column 4: SUMMA atom type
# column 5: SUMMA environment type
# column 6: six-atom solvation type
# column 7: Amber All-Atom Charge
# column 8: Amber United-Atom Charge
#
# The first non-# line MUST be a ! line.  From that point on, there must not be ANY blank lines.  This is 
# important b/c ligandTemplate.cc is very picky when reading in this data. 
#
# Note: you must use 1 Tab between each column. If you do not, ligandTemplate won't parse the file.
# Note: For the ! lines, also use spaces (!/tNAME/tNumAtoms/tConnectivity).
# Note: if you do not have a value for a column, use a '-' 
# Note: ligandTemplate must count hydrogens. I made the assumption that all column 1 hydrogen names start with 'H'.
#       If you have an instance where the name doesn't start with 'H' then you will have a problem.
# Note: Hydrogens should all be listed at the end of each type.  If you have any non-hydrogens after a hydrogen,
#	then the amberUnited functions will fail.  
#
# The file "ligandList.lig" contains the name of each .lig file on separate lines.  When you add a new .lig file,
# just add the filename to a new line in this list.
#