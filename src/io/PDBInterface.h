//Modified version to include hetatm parsing.  
//Last Modified 7/10/2002 by Jeff Kearns

#include <string.h>
#include <stdlib.h>
#include <vector>
#include <iostream>
#include <fstream>
#include "typedef.h"
#include "aaBaseline.h"
#include "ensemble.h"
#include "molecule.h"
#include "protein.h"
//#include "residue.h"
//#include "atom.h"

#include "PDBAtomRecord.h"

#ifndef PDBINTERFACE_H
#define PDBINTERFACE_H
class PDBInterface
{

	public:
		PDBInterface();
		PDBInterface(const string& _filename);
		PDBInterface(const string& _filename, const bool _Hflag);
		PDBInterface(const string& _filename, const bool _Hflag, const bool _HPflag, const bool _dummyBool);
		PDBInterface(const PDBInterface& _otherPdb);
		PDBInterface(const string& _filename, bool _hetflag, bool _atomflag);
		~PDBInterface();

	public:
		ensemble* getEnsemblePointer() {return pItsEnsemble;}
		void setHetatmFlag(bool& _hetflag){itsHetatmFlag=_hetflag;}
		void setAtomFlag(bool& _atomflag){itsAtomFlag=_atomflag;}		

	private:
		void readData(ifstream& _infile);
		void categorizeLines();
		void parseAtomLine();
		void parseAtomLine(const bool _Hflag);
		void parseAtomLine(const bool _Hflag, const bool _HPflag);
		void parseHetatmLine();
		void parseBIOMT();

	private:
		string itsFilename;
		ensemble* pItsEnsemble;
		vector<string> theLines;
		bool itsHetatmFlag;
		bool itsAtomFlag;
		vector<string> resNames;

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
		static aaBaseline itsAABaseline;

};
#endif
