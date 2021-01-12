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
enum aminoAcid {A,R,N,D,Dh,C,Cx,Cf,Q,E,Eh,Hd,He,Hp,I,L,K,M,F,P,O,S,T,W,Y,V,G,dA,dR,dN,dD,dDh,dC,dCx,dCf,dQ,dE,dEh,dHd,dHe,dHp,dI,dL,dK,dM,dF,dP,dO,dS,dT,dW,dY,dV,SF4,HEM,NI2,CLN,CO2,MG2,OH,OXY,CLD,HIS};
string aminoAcidString[] = {"A","R","N","D","Dh","C","Cx","Cf","Q","E","Eh","Hd","He","Hp","I","L","K","M","F","P","O","S","T","W","Y","V","G","dA","dR","dN","dD","dDh","dC","dCx","dCf","dQ","dE","dEh","dHd","dHe","dHp","dI","dL","dK","dM","dF","dP","dO","dS","dT","dW","dY","dV","SF4","HEM","NI2","CLN","CO2","MG2","OH-","OXY","CLD","HIS"};
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
		cout << "Active Chains,0,1,2," << endl;
		cout << "Active Positions,0,1,2,3,5,6,7,9,10," << endl;
		cout << "A,K,D,L,K,D,R,R,R," << endl;
		exit(1);
	}
	string inputfile = argv[1];
	string infile;
	vector <UInt> seq;
	UIntVec activeChains, activeResidues;
	
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
				if (linecounter < 3){
					if (delimitercounter > 0){
						if (linecounter == 0){
							infile = item;
						}
						if (linecounter == 1){
							stringstream inputString(item);
							UInt index;
							inputString >> index;
							activeChains.push_back(index);
						}
						if (linecounter == 2){
							stringstream inputString(item);
							UInt index;
							inputString >> index;
							activeResidues.push_back(index);
						}
					}
				}
				else{
					UInt index = convertAAStringtoInt(item, aminoAcidString, aaSize);
					seq.push_back(index);
				}
				delimitercounter++;
			}
			linecounter++;
		}
		file.close();
	}
	else{
		fstream inf;
		inf.open ("mutator.in", fstream::in | fstream::out | fstream::app);
		inf << "Input PDB File,xyz.pdb," << endl;
		inf << "Active Chains,0,1,2," << endl;
		inf << "Active Positions,0,1,2,3,5,6,7,9,10," << endl;
		inf << "A,K,D,L,K,D,R,R,R," << endl;
		cout << "Error: Required input file doesn't exist." << endl << "Template input file has been generated, please fill it out and rerun." << endl;
		exit(1);
	}
		PDBInterface* thePDB = new PDBInterface(infile);
		ensemble* theEnsemble = thePDB->getEnsemblePointer();
		molecule* pMol = theEnsemble->getMoleculePointer(0);
		protein* bundle = static_cast<protein*>(pMol);
		for (UInt i = 0; i < activeChains.size(); i ++)
		{
			for (UInt j = 0; j < activeResidues.size(); j++)
			{
				bundle->activateForRepacking(activeChains[i], activeResidues[j]);
				bundle->mutateWBC(activeChains[i], activeResidues[j], seq[j]);
			}
		}
		bundle->protRelax(1000);
		string outFile;
		outFile = "mut.pdb";
		pdbWriter(bundle, outFile);
		delete thePDB;
	return 0;
}

UInt convertAAStringtoInt(string AA, string aminoAcidString[], UInt size)
{
	UInt index;
	for (UInt i = 0; i < size; i++)
	{
		if (AA.compare(aminoAcidString[i]) == 0){index = i; break;}
	}
	return index;
}

