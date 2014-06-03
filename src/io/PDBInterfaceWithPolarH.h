#include <string>
#include <stdlib.h>
#include <vector>
#include <iostream.h>
#include <fstream.h>
#include "typedef.h"
#include "ensemble.h"
#include "molecule.h"
#include "protein.h"
//#include "ligand.h"
//#include "residue.h"
//#include "atom.h"

#include "PDBAtomRecord.h"

#ifndef PDBINTERFACEWITHPOLARH_H
#define PDBINTERFACEWITHPOLARH_H
class PDBInterfaceWithPolarH
{

	public:
		PDBInterfaceWithPolarH();
		PDBInterfaceWithPolarH(const string& _filename);
		PDBInterfaceWithPolarH(const PDBInterface& _otherPdb);
		~PDBInterfaceWithPolarH();

	public:
		ensemble* getEnsemblePointer() {return pItsEnsemble;}

	private:
		void readData(ifstream& _infile);
		void categorizeLines();
		void parseAtomLine();
		void parseHetatmLine();
		void parseBIOMT();

	private:
		string itsFilename;
		ensemble* pItsEnsemble;
		vector<string> theLines;

		vector<UInt> atomLines;
		vector<UInt> headerLines;
		vector<UInt> hetatmLines;
		vector<UInt> titleLines;
		vector<UInt> compndLines;
		vector<UInt> sourceLines;
		vector<UInt> keywdsLines;
		vector<UInt> expdtaLines;
		vector<UInt> authorLines;
		vector<UInt> revdatLines;
		vector<UInt> jrnlLines;
		vector<UInt> remarkLines;
		vector<UInt> dbrefLines;
		vector<UInt> seqresLines;
		vector<UInt> formulLines;
		vector<UInt> helixLines;
		vector<UInt> sheetLines;
		vector<UInt> siteLines;
		vector<UInt> crystLines;
		vector<UInt> origxLines;
		vector<UInt> scaleLines;
		vector<UInt> mtrixLines;
		vector<UInt> terLines;
		vector<UInt> masterLines;
		vector<UInt> conectLines;
		vector<UInt> endLines;
		vector<UInt> ftnoteLines;
		vector<UInt> seqadvLines;
		vector<UInt> turnLines;
		vector<UInt> hetLines;
		vector<UInt> hetnamLines;
		vector<UInt> ssbondLines;
		vector<UInt> linkLines;
		vector<UInt> modresLines;
		vector<UInt> sltbrgLines;
		vector<UInt> sprsdeLines;
		vector<UInt> modelLines;
		vector<UInt> anisouLines;
		vector<UInt> userLines;
		vector<UInt> energyLines;

};
#endif
