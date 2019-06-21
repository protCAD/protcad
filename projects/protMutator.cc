//*******************************************************************************************************
//*******************************************************************************************************
//**************************************                       ******************************************
//**************************************    protMutator 1.1    ******************************************
//**************************************                       ******************************************
//*******************************************************************************************************
//*******************************************************************************************************


//--Program setup----------------------------------------------------------------------------------------
#include <iostream>
#include <string>
#include "ensemble.h"
#include "PDBInterface.h"
#include <sstream>

UInt convertAAStringtoInt(string AA, string aminoAcidString[], UInt size);
enum aminoAcid {A,R,N,D,Dh,C,Cx,Cf,Q,E,Eh,Hd,He,Hp,I,L,K,M,F,P,O,S,T,W,Y,V,G,dA,dR,dN,dD,dDh,dC,dCx,dCf,dQ,dE,dEh,dHd,dHe,dHp,dI,dL,dK,dM,dF,dP,dO,dS,dT,dW,dY,dV,Csf,Sf4,Hca,Eoc,Oec,Saf,Hem,Cyn};
string aminoAcidString[] = {"A","R","N","D","Dh","C","Cx","Cf","Q","E","Eh","Hd","He","Hp","I","L","K","M","F","P","O","S","T","W","Y","V","G","dA","dR","dN","dD","dDh","dC","dCx","dCf","dQ","dE","dEh","dHd","dHe","dHp","dI","dL","dK","dM","dF","dP","dO","dS","dT","dW","dY","dV","Csf","Sf4","Hca","Eoc","Oec","Saf","Hem","Cyn"};
UInt aaSize = sizeof(aminoAcidString)/sizeof(aminoAcidString[0]);

int main (int argc, char* argv[])
{
	//--Program setup
	if (argc !=2)
	{
		cout << "Error: Required input file not given, should be run with input file." << endl;
		cout << "Command: protMutator <inputfile>" << endl;
		cout << "Input file format as listed below:" << endl;
		cout << "Input PDB File,xyz.pdb," << endl;
		cout << "A,K,D,L,K,D,R,R," << endl;
		cout << "A,K,E,L,K,E,R,R," << endl;
		exit(1);
	}
	string inputfile = argv[1];
	string infile;
	vector<vector<UInt> > seqs;
	vector <UInt> seq;
	
	//--read input file
	ifstream file(inputfile);
	if (file){
		string item, line;
		UInt delimitercounter, linecounter = 0;
		while(getline(file,line))
		{
			delimitercounter = 0;
			stringstream stream(line);
			while(getline(stream,item,','))
			{
				if (linecounter == 0){
					if (delimitercounter > 0){
						infile = item;
					}
				}
				else{
					UInt index = convertAAStringtoInt(item, aminoAcidString, aaSize);
					seq.push_back(index);
				}
				delimitercounter++;
			}
			linecounter++;
			seqs.push_back(seq);
			seq.clear();
		}
		file.close();
	}
	else{
		fstream inf;
		inf.open ("mutator.in", fstream::in | fstream::out | fstream::app);
		inf << "Input PDB File,xyz.pdb," << endl;
		inf << "A,K,D,L,K,D,R,R," << endl;
		inf << "A,K,E,L,K,E,R,R," << endl;
		cout << "Error: Required input file doesn't exist." << endl << "Template input file has been generated, please fill it out and rerun." << endl;
		exit(1);
	}
	
	//--Mutate same sequence on each chain and generate a model for each sequence
	for (UInt h = 0; h < seqs.size(); h++)
	{
		PDBInterface* thePDB = new PDBInterface(infile);
		ensemble* theEnsemble = thePDB->getEnsemblePointer();
		molecule* pMol = theEnsemble->getMoleculePointer(0);
		protein* bundle = static_cast<protein*>(pMol);
		UInt chainNum = bundle->getNumChains();
		for (UInt i = 0; i < chainNum; i ++)
		{
			UInt resNum = bundle->getNumResidues(i);
			for (UInt j = 0; j < resNum; j++)
			{
				bundle->activateForRepacking(i, j);
				bundle->mutateWBC(i, j, seqs[h][j]);
			}
		}
		bundle->protMin(true);
		stringstream convert;
		string countStr, outFile;
		convert << h+1, countStr = convert.str();
		outFile = countStr + ".mut.pdb";
		pdbWriter(bundle, outFile);
		delete thePDB;
	}
	return 0;
}

UInt convertAAStringtoInt(string AA, string aminoAcidString[], UInt size)
{
	UInt index;
	for (UInt i = 0; i < size; i++)
	{
		if (AA.compare(aminoAcidString[i]) == 0){index = i;}
	}
	return index;
}

