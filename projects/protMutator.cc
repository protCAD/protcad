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
string aminoAcidString[] = {"A","R","N","D","Dh","C","Cx","Cf","Q","E","Eh","Hd","H","Hp","I","L","K","M","F","P","O","S","T","W","Y","V","G","dA","dR","dN","dD","dDh","dC","dCx","dCf","dQ","dE","dEh","dHd","dHe","dHp","dI","dL","dK","dM","dF","dP","dO","dS","dT","dW","dY","dV","SF4","HEM","NI2","CLN","CO2","MG2","OH-","OXY","CLD","HIS"};
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
		cout << "Active Chain,0," << endl;
		cout << "V57A,G23R," << endl;
		exit(1);
	}
	string inputfile = argv[1];
	string infile;
	UInt activeChain;
	vector <UIntVec> activeResidues;
	UIntVec activeRes, seq;
	vector <UIntVec> seqs;
	UInt variant = 0;
	UInt resnum, resIndex;
	
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
				if (linecounter < 2){
					if (delimitercounter > 0){
						if (linecounter == 0){
							infile = item;
						}
						if (linecounter == 1){
							stringstream inputString(item);
							UInt index;
							inputString >> index;
							activeChain = index;
						}
					}
				}
				else{
					const char *itemchar = item.c_str();
					if (isdigit(itemchar[0])){resnum = atoi(itemchar);}
					else{itemchar = &itemchar[1]; resnum = atoi(itemchar);}
					activeRes.push_back(resnum);
					UInt i = 0; string aa;
					while (item[i]){if(isalpha(item[i])){aa = item[i];}i++;}
					UInt index = convertAAStringtoInt(aa, aminoAcidString, aaSize);
					seq.push_back(index);
				}
				delimitercounter++;
			}
			linecounter++;
			if (linecounter > 2){
				activeResidues.push_back(activeRes); 
				seqs.push_back(seq);
				seq.clear(); 
				seq.resize(0);
				activeRes.clear(); 
				activeRes.resize(0); 
				variant++;
			}
		}
		file.close();
	}
	else{
		fstream inf;
		inf.open ("mutator.in", fstream::in | fstream::out | fstream::app);
		inf << "Input PDB File,xyz.pdb," << endl;
		inf << "Active Chain,0," << endl;
		inf << "V57A,G23R," << endl;
		cout << "Error: Required input file doesn't exist." << endl << "Template input file has been generated, please fill it out and rerun." << endl;
		exit(1);
	}
	for (UInt i = 0; i < seqs.size(); i++)
	{
		PDBInterface* thePDB = new PDBInterface(infile);
		ensemble* theEnsemble = thePDB->getEnsemblePointer();
		molecule* pMol = theEnsemble->getMoleculePointer(0);
		protein* prot = static_cast<protein*>(pMol);
		for (UInt j = 0; j < activeResidues[i].size(); j++)
		{
			resIndex = prot->getIndexFromResNum(activeChain, activeResidues[i][j]);
			prot->activateForRepacking(activeChain, resIndex);
			prot->mutateWBC(activeChain, resIndex, seqs[i][j]);
		}
		prot->protRelax(1000, false);
		stringstream convert; string countstr;
		convert << i+1, countstr = convert.str();
		string outFile = "cand" + countstr + "+AAKL.pdb";
		pdbWriter(prot, outFile);
		delete thePDB;
	}
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

