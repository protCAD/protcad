//*******************************************************************************************************
//*******************************************************************************************************
//**************************************                       ******************************************
//**************************************     protShaper 1.0    ******************************************
//**************************************                       ******************************************
//*******************************************************************************************************
//*******************************************************************************************************


//--Program setup----------------------------------------------------------------------------------------
#include <iostream>
#include <string>
#include "ensemble.h"
#include "PDBInterface.h"
#include <unistd.h>
#include <sstream>

UInt convertBBStringtoInt(string BB, string backboneSeq[], UInt size);
enum structure {O,L,P,B,E,X,A,I,U,J,Z,F,H,W,K,S};
string backboneSeq[] =   { "O", "L", "P", "B","E","X","A","I",  "U",  "J",  "Z",  "F", "H", "W", "K", "S"};
UInt bbSize = sizeof(backboneSeq)/sizeof(backboneSeq[0]);

int main (int argc, char* argv[])
{
	//--Program setup
	if (argc !=2)
	{
		cout << "Error: Required input file not given, should be run with input file." << endl;
		cout << "Command: protShaper <inputfile>" << endl;
		cout << "Input file format as listed below:" << endl;
		cout << "Input PDB File,xyz.pdb," << endl;
		cout << "A,K,D,L,K,D,R,R,R," << endl;
		exit(1);
	}
	int seed = (int)getpid(); srand (seed);
	stringstream convert; string startstr, outFile;
	UInt name = rand() % 100000000;
	convert << name, startstr = convert.str();
	string inputfile = argv[1];
	string infile;
	vector <UInt> seq;
	vector <double> backboneAngles(2);
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
				if (delimitercounter > 0 && linecounter == 0){
					infile = item;
				}
				else if (linecounter > 0){
					UInt index = convertBBStringtoInt(item, backboneSeq, bbSize);
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
		inf.open ("shaper.in", fstream::in | fstream::out | fstream::app);
		inf << "Input PDB File,xyz.pdb," << endl;
		inf << "A,K,D,L,K,D,R,R,R," << endl;
		cout << "Error: Required input file doesn't exist." << endl << "Template input file has been generated, please fill it out and rerun." << endl;
		exit(1);
	}
	for (UInt h = 0; h < 50; h++)
	{
		PDBInterface* thePDB = new PDBInterface(infile);
		ensemble* theEnsemble = thePDB->getEnsemblePointer();
		molecule* pMol = theEnsemble->getMoleculePointer(0);
		protein* prot = static_cast<protein*>(pMol);
		for (UInt i = 0; i < seq.size(); i ++)
		{
			if (seq[i] == 6)
			{
				prot->setDihedral(0, i,-60,0,0);
				prot->setDihedral(0, i,-40,1,0);
				prot->mutateWBC(0,i,0);
			}
			else{
				backboneAngles = prot->getRandPhiPsifromBackboneSequenceType(seq[i]);
				prot->setDihedral(0, i,backboneAngles[0],0,0);
				prot->setDihedral(0, i,backboneAngles[1],1,0);
			}
		}
		prot->protMin(true);
		UInt sec = time(NULL), timeid = sec;
		stringstream convert; string countstr;
		convert << timeid, countstr = convert.str();
		outFile = countstr + "." + startstr + ".model.pdb";
		cout << outFile << " " << prot->protEnergy() << endl;
		pdbWriter(prot, outFile);
		delete thePDB;
	}
	return 0;
}

UInt convertBBStringtoInt(string BB, string backboneSeq[], UInt size)
{
	UInt index;
	for (UInt i = 0; i < size; i++)
	{
		if (BB.compare(backboneSeq[i]) == 0){index = i; break;}
	}
	return index;
}

