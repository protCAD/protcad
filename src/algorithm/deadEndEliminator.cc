#include "deadEndEliminator.h"
#define DEE_DEBUG

deadEndEliminator::deadEndEliminator(protein* _pProtein)
{
#ifdef DEE_DEBUG
	cout << "deadEndEliminator constructor called:"
	     << "deadEndEliminator::deadEndEliminator()" << endl;
#endif
	pItsProtein = _pProtein;
}

deadEndEliminator::~deadEndEliminator()
{
#ifdef DEE_DEBUG
	cout<< "deadEndEliminator destructor called " << endl;
#endif
}

void deadEndEliminator::run(UInt _level, double _cutoff)
{	
	itsRunLevel = _level;
	itsEnergyCutoff = _cutoff;
	itsCurrentIteration = 0;
#ifdef DEE_DEBUG
	string outfilename;
	outfilename = "orig.pdb";
	char fmodchar[3];
	int fmodint = 1;
#endif
	
	pItsProtein->stripToGlycine();
#ifdef DEE_DEBUG
	pdbWriter(static_cast<protein*>(pItsProtein),outfilename);	
#endif

	for (UInt currentChain=0; currentChain < pItsProtein->getNumChains(); currentChain++)
	{	vector <chainPosition*> theChainPositions = pItsProtein->getChainPositionVector(currentChain);
		for (UInt currentCP = 0; currentCP < theChainPositions.size(); currentCP++)
		{
			// check to see if chain position at currentCP is active	
			if (theChainPositions[currentCP])
			{	
#ifdef DEE_DEBUG
				cout << "ChainPosition: " << theChainPositions[currentCP] << endl;
#endif
				vector< vector <UIntVec> > theDatabase;
				theDatabase = theChainPositions[currentCP]->getAllowedDB();
				UInt numAllowedAA = theDatabase.size();
				for (UInt currentAA = 0; currentAA < numAllowedAA; currentAA++)
				{
#ifdef DEE_DEBUG
					cout << "current AA: " << currentAA ;
#endif
					UInt numAllowedBpt = theDatabase[currentAA].size();
					pItsProtein->mutateWBC(currentChain, currentCP, currentAA); 
#ifdef DEE_DEBUG
					cout << " numBpt: " << numAllowedBpt << endl;
#endif
					
						for (UInt currentBpt = 0; currentBpt < numAllowedBpt ; currentBpt++)
						{
							UInt numAllowedRot = theDatabase[currentAA][currentBpt].size();
							if (numAllowedRot)
							{
								for (UInt currentRot = 0; currentRot < numAllowedRot; currentRot++)
								{
									pItsProtein->setRotamerWBC(currentChain, currentCP, currentBpt, currentRot);
#ifdef DEE_DEBUG
									outfilename = "mod";
									sprintf(fmodchar,"%i",fmodint);
									string fmodstring = fmodchar;
									outfilename += fmodstring;
									outfilename += ".pdb";
									fmodint++;
									pdbWriter(static_cast<protein*>(pItsProtein),outfilename);
#endif

									double theEnergy;
									theEnergy = pItsProtein->intraEnergy();
									cout << theDatabase[currentAA][currentBpt][currentRot] << ": ";
									cout << theEnergy << ", ";
								}
							}
							else
							{	
								// no rotamers for this amino acid
#ifdef DEE_DEBUG
								outfilename = "mod";
								sprintf(fmodchar,"%i",fmodint);
								string fmodstring = fmodchar;
								outfilename += fmodstring;
								outfilename += ".pdb";
								fmodint++;
								pdbWriter(static_cast<protein*>(pItsProtein),outfilename);
#endif
								double theEnergy;
								theEnergy = pItsProtein->intraEnergy();
								cout << theEnergy << ", ";
							}
							cout << endl;
					}
				}
			}
		}
	}
}

double deadEndEliminator::min( vector< double > _theList)
{	double lowest = _theList[0];
	for (UInt i=0; i< _theList.size(); i++)
	{
		if (_theList[i] < lowest)
		{	lowest = _theList[i];
		}
	}
	return lowest;
}

double deadEndEliminator::max( vector< double > _theList)
{	double highest = _theList[0];
	for (UInt i=0; i< _theList.size(); i++)
	{
		if (_theList[i] > highest)
		{	highest = _theList[i];
		}
	}
	return highest;
}

void makeRotamerDeadEnding(UInt _chain, UInt _residueIndex,
			UInt _residueType, UInt _bpt, UInt _rotamer)
{
#ifdef DEE_DEBUG	
	cout << "Labelling as dead-ending chain: " << _chain << " residue: "
		<< _residueIndex << " type: " << _residueType << " branchpoint: "
		<< _bpt << " rotamer: " << _rotamer << endl;
	// dummy for now
#endif
	return;
}
