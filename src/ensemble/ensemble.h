//	____ _  _ ____ ____ _  _ ___  _    ____  _  _ 
//	|___ |\ | [__  |___ |\/| |__] |    |___  |__| 
//	|___ | \| ___] |___ |  | |__] |___ |___ o|  | 
//	                                              
// filename: ensemble.h
// contents: class ensemble is defined

#include "assert.h"
#include <vector>
#include "ran.h"
#include "typedef.h"
#include "molecule.h"
//#include "CMath.h"

#ifndef ENSEMBLE_H
#define ENSEMBLE_H


class ensemble
{
public:
	friend class annealer;

	// Constructor and Destructor declaration
	ensemble();
	~ensemble();

	// Accessors
	void add(molecule* _pMolecule);
	void remove(molecule* _pMolecule);
	
	int modify(ran& _ran);
	int modify(ran& _ran, vector <int> _position);
	int mutate(vector <int> _position, UInt _resType);
	int mutateWithSymmetry(vector <int> _position, UInt _resType);

	int symmetryLinkMolAndChain(UInt _molecule1, UInt _chain1, UInt _molecule2, UInt _chain2);
	void printLinkageInfo();

        molecule * getMoleculePointer(UInt _index);
        UInt getNumMolecules(){return itsMolecules.size();}
	double getPositionEnergy(vector <int> _position);
	vector <int> getLastModification();
	vector <int> chooseNextTargetPosition(ran& _ran);
	UInt chooseNextMutationIdentity(ran& _ran, vector <int> _position);
	double energy();
	double getIntraEnergy(UInt _molecule) { return itsMolecules[_molecule]->intraEnergy(); }
        void setCutoffDistance(double _cutoff){itscutoffDistance=_cutoff;}
        void setAllEnsembleCutoffDistance();
        void setAllEnsembleCutoffDistance(double _cutoff){setCutoffDistance(_cutoff); setAllEnsembleCutoffDistance();}
	double getVolume(UInt _method); 
	void acceptModification();
	void rejectModification();
	void setupSystem(ran& _ran);
	void saveState(string& _filename);


private:

	// list of mol-chain pairs that are treated as independent during symmetry related
	// operations
	vector < vector < UInt > > itsIndMolAndChainList; 

	// list of symmetry linked mols and chains 
	// first index is the independent chain these are slaves of
	// second index is the number of molecules that are linked to the independent one
	// third index is the molecule (0) and the chain (1)
	vector < vector < vector < UInt > > > itsMolAndChainLinkageMap;

	int chooseMolecule(ran& _ran);
	//variable declarations
	vector<molecule*> itsMolecules;
	int itsLastModifiedMolecule;
	double itsEnergy;
	double itsVolume;
        static double itscutoffDistance;
};
#endif
