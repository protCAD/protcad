![background](https://raw.githubusercontent.com/protCAD/protcad/master/ui/images/splash.png)


Official protCAD development tree
===================================================================================================

protCAD is an implementation of the protein design software library that originated in the 
Bill Degrado Protein Design Lab.

It is currently maintained by The Vikas Nanda Lab: https://sites.google.com/site/viknanda
The source is maintained at: https://github.com/protCAD/protcad

Publications to date on protCAD's methods and implementaions are:

-Computational Methods and their Applications for de novo Functional Protein Design and Memebrane 
 Protein Solubilization, Summa CM Thesis 2002

-Empirical estimation of local dielectric constants: Toward atomistic design of collagen mimetic 
 peptides, Biopolymers - Peptide Science. Pike & Nanda 2015; 104(4): 360-70.
 
-Computational design of a sensitive, selective phase-changing sensor protein for the VX nerve agent. 
 Science Advances McCann & Pike et al. 6 Jul 2022 Vol 8, Issue 27


 Installation
===================================================================================================

=== Install dependencies

--Windows 10

First you will need to follow instructions to install the windows ubuntu sub-system as there is no
native support for windows libraries: https://docs.microsoft.com/en-us/windows/wsl/install-win10
Then follow the Ubuntu Linux install dependency instructions and install below in the sub-system terminal.

--Ubuntu Linux:

In terminal:

sudo apt-get install g++ gfortran git;

--Mac:

In terminal install dev tools, homebrew and then dependencies:

xcode-select --install

/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

brew install gcc gfortran git


=== Install

For all systems download source and compile via these terminal commands:

git clone https://github.com/protCAD/protcad 

cd protcad

make install

Then close and re-open the terminal


 Usage
===================================================================================================

General use programs are compiled via make install and the source is available in protcad/projects.

Since make install also adds the bin directory to your path, programs in projects are immediately 
available to run in any directory.

An overview of the general use programs are described below:


--protAlign

Description: Calculates best fit and RMSD of two pdbs using SVD and aligns the second or smaller
pdb onto the first.

Input command: protAlign <inFile1.pdb> <inFile2.pdb>
Output result:  RMSD, aligned pdbs


--protDielectric

Description: Calculates the local average dielectric and solvation energy for each residue.

Input command: protDielectric <inFile.pdb>
Output result: List of local dielectrics and solvation energy for each residue


--protDihedrals

Description: Calculates backbone phi psi and backbone classification type for each residue.
Backbone Classification Type: -γ -π -α -ρ -β β ρ α π γ -γi -πi -αi -ρi -βi βi ρi αi πi γi


Input command: protDihedrals <inFile.pdb>
Output result: List of phi psi and classification type for each residue


--protEnergy

Description: Calculates the total energy of the protein in kcal/mol.

Input command: protEnergy <inFile.pdb>
Output result: Total Energy of the protein, total clashes and total backbone-backbone clashes


--protEvolver

Description: Sequence Selective Machine Learning Evolution Based Algorithm in Implicit Solvent

Input command: protEvolver <inputfile>

Input file format:
Input PDB File,xyz.pdb,
Active Chains,0,1,2,
Active Positions,0,1,2,3,5,6,7,9,10,
Random Positions,0,2,5,6,10,
Frozen Positions,4,8,
Amino Acids,A,R,N,D,C,Q,E,He,I,L,K,M,F,P,S,T,W,Y,V,G,
Backbone Relaxation,false,

Output result: Evolved model pdbs, sequences and energies written to results.out file


--protFolder (in development, nearly complete)

Description: Fold Selective Machine Learning Evolution Based Algorithm in Implicit Solvent

Input command: protFolder <inputfile>

Input file format:
Input PDB File,xyz.pdb,
Active Chains,0,1,2,
Active Positions,0,1,2,3,5,6,7,9,10,
Random Positions,0,2,5,6,10,
Backbone Types,m,c,l, p, b,t,y,a,i,g,n,d,q,r,f,h,w,k,s,v,

Backbone Classification Type: -γ -π -α -ρ -β β ρ α π γ -γi -πi -αi -ρi -βi βi ρi αi πi γi
Backbone Classification Key:   m  c  l  p  b t y a i g  n   d   q   r   f  h  w  k  s  v

Output result: Folded model pdbs, backbone type sequences and energies written to results.out file


--protInverter

Description: Generates mirror image of input pdb conformation

Input command: protInverter <inFile.pdb> <outFile.pdb>
Output result: Output pdb where pdb is mirror image of input pdb including sidechain chirality


--protMin

Description: Minimizes the energy of the structure with sidechain and local backbone motion

Input command: protMin <inFile.pdb> <outFile.pdb>
Output result: Outputs energy minimized pdb and returns start and end energy


--protMover

Description: Moves a protein structure in XYZ space via rotation and translation.

Input command:
protMover <in.pdb> <translateX> <translateY> <translateZ> <rotateX> <rotateY> <rotateZ> <out.pdb>
Output result: Output pdb rotated in degrees and translated in Angstroms


--protMutator

Description: Mutates new sequence on input protein structure, minimizes and returns model

Input command: protMutator <inputfile>

Input file format:
Input PDB File,xyz.pdb,
Active Chains,0,1,2,
Active Positions,0,1,2,3,5,6,7,9,10,
A,K,D,L,K,D,R,R,R,

Output result: Minimized mutated model pdb of amino acid muations at positions in in input file


--protOligamer

Description: Creates symmetric oligamers of input pdb with same coordinates for each chain.

Input command: protOligamer <inFile.pdb>
Output result: pdbs of symmetric oligamers (sampling parameters will need to be manually adjusted)


--protSequence

Description: Reports amino acid sequence and backbone type sequence for a protein structure.

Input command: protSequence <inFile.pdb>
Output result: Returns amino acid sequence and backbone type sequence in fasta format

Backbone Classification Type: -γ -π -α -ρ -β β ρ α π γ -γi -πi -αi -ρi -βi βi ρi αi πi γi
Backbone Classification Key:   m  c  l  p  b t y a i g  n   d   q   r   f  h  w  k  s  v


 Development process
===================================================================================================

New programs need to be added to the Makefile in protcad/ for compilation and added to projects 
directory which then can be compiled via 'make programx'.

Additional examples are in the protcad/src/tests folder for reference and usage, which can be added
to the projects directory and MakeFile for compilation.  Some tests may not compile as they were 
developed on a much earlier version of protcad but most will compile without issue.

New programs written are intended to follow the same directory structure and compilation as 
in projects.  Developers work in their own trees, then submit pull requests when they think 
their feature or bug fix is ready.

The master branch is regularly built and tested. Tags are regularly created to indicate new
official, stable release versions of protCAD.

The dev branch is a more frequently updated branch with a more experimental version being developed
but generally ready for usage. To switch to the dev branch in the protad directory use the command:

git checkout dev
make install


 Issues
===================================================================================================
Bugs and issues can be submitted into the issues section of http://www.github.com/protcad/protcad repo or a
description of the issue and result can be emailed to Douglas Pike at doughp11@gmail.com

