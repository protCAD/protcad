#include "atomIterator.h"

atomIterator::atomIterator(protein* _pProtein)
{	
        pInputProtein = _pProtein;
        hetatmFlag=false;
        initialize();
}

atomIterator::atomIterator(ligand* _pLigand)
{
	pInputLigand=_pLigand;
	hetatmFlag=true;
	initializeLigand();
}
	
atomIterator::~atomIterator()
{
}

void atomIterator::initializeLigand()
{
	heavyOnly= false;
	failFlag= 0;
	numAtomsInCurrent=0;
	itsCurrentAtomIndex = 0;

	if((numAtomsInCurrent=pInputLigand->itsAtoms.size())){
		pItsCurrentAtom =pInputLigand->itsAtoms[itsCurrentAtomIndex];
	}
	else{
		cout <<"numAtomsInCurrent != ligand size in initializeLigand() in AtomIterator.cc"<<endl;
		pItsCurrentAtom=0;
	}
	//cout <<"pItsCurrentAtom serial= " << pItsCurrentAtom->getSerialNumber() << endl;
	//cout << "Done with initializeLigand()..." << endl;
} 

void atomIterator::initialize()
{
	heavyOnly = true;
	failFlag = 0;
	itsCurrentChainIndex = 0;
	itsCurrentResidueIndex = 0;
	itsCurrentAtomIndex = 0;
	
// need to add some bounds checking here to make sure that each of these is
// fully defined upon instantiation...
// Note: Initialization in conditions of "if" loop are intended in this case....

	if( (numChainsInCurrent = pInputProtein->itsChains.size()) )
	{	pItsCurrentChain   = (pInputProtein->itsChains)[itsCurrentChainIndex];
		if( (numResiduesInCurrent = pItsCurrentChain->itsResidues.size()) )
		{	pItsCurrentResidue = (pItsCurrentChain->itsResidues)[itsCurrentResidueIndex];
			if( (numAtomsInCurrent = pItsCurrentResidue->itsAtoms.size()) )
			{	pItsCurrentAtom    = (pItsCurrentResidue->itsAtoms)[itsCurrentAtomIndex];
			}
			else
			{	pItsCurrentAtom = 0;
			}
		}
		else
		{	pItsCurrentResidue = 0;
			pItsCurrentAtom = 0;
		}
	}
	else
	{	pItsCurrentChain = 0;
		pItsCurrentResidue = 0;
		pItsCurrentAtom = 0;
	}
//	cout << "Finished Construction of Atom Iterator" << endl;
//	cout << "pItsCurrentChain = " << pItsCurrentChain << endl;
//	cout << "pItsCurrentResidue = " << pItsCurrentResidue  << endl;
//	cout << "pItsCurrentAtom = " << pItsCurrentAtom << endl;
}

atomIterator& atomIterator::operator++ (int _x)
{	increment();
	//cout << " failFlag = " << failFlag << endl;
	return *this;
}

void atomIterator::increment()
{	
        //ligands
        if(hetatmFlag){
            if ( updateCurrentAtomIndex() ){return;}
            else
            {	
		failFlag = 1;
		return;
            }
        }
	
        //proteins
	else
        {
		if ( updateCurrentAtomIndex() ){return;}
    
                else if ( updateCurrentResidueIndex() ){return;}
	
                else if ( updateCurrentChainIndex() ){return;}
	
                else
                {	
                    failFlag = 1;
                    return;
                }	
        }
	
}

bool atomIterator::updateCurrentAtomIndex()
{	
        //ligands
        if(hetatmFlag)
        {
            if (itsCurrentAtomIndex < numAtomsInCurrent-1)
            {
                    itsCurrentAtomIndex++;
                    pItsCurrentAtom = (pInputLigand->itsAtoms[itsCurrentAtomIndex]);
                    return true;
            }
            else return false;
        }
        
        //proteins
        else
        {
            if (itsCurrentAtomIndex < numAtomsInCurrent-1)
            {
		itsCurrentAtomIndex++;
		pItsCurrentAtom = (pItsCurrentResidue->itsAtoms)[itsCurrentAtomIndex];
		return true;
            }
            else return false;
        }
}

bool atomIterator::updateCurrentResidueIndex()
{	if (itsCurrentResidueIndex < numResiduesInCurrent-1)
        {       itsCurrentResidueIndex++;
                pItsCurrentResidue = (pItsCurrentChain->itsResidues)[itsCurrentResidueIndex];
                numAtomsInCurrent = pItsCurrentResidue->itsAtoms.size();
                itsCurrentAtomIndex = 0;
                pItsCurrentAtom = (pItsCurrentResidue->itsAtoms)[itsCurrentAtomIndex];
                return true;
        }
	else return false;
}

bool atomIterator::updateCurrentChainIndex()
{	if (itsCurrentChainIndex < numChainsInCurrent-1)
	{	itsCurrentChainIndex++;
		pItsCurrentChain = (pInputProtein->itsChains)[itsCurrentChainIndex];
		numResiduesInCurrent = pItsCurrentChain->itsResidues.size();
		itsCurrentResidueIndex = 0;
		pItsCurrentResidue = (pItsCurrentChain->itsResidues)[itsCurrentResidueIndex];
		numAtomsInCurrent = pItsCurrentResidue->itsAtoms.size();
		itsCurrentAtomIndex = 0;
		pItsCurrentAtom  = (pItsCurrentResidue->itsAtoms)[itsCurrentAtomIndex];
		return true;
	}
	else return false;
}
