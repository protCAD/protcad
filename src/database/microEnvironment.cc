#include "microEnvironment.h"
//#define MICROENV_DEBUG

// default value for itsMultFactor
double microEnvironment::itsScaleFactor = 1.0;

microEnvironment::microEnvironment()
{
}

microEnvironment::microEnvironment(const double _rad, const UInt _skip)
{	
	residue::setupDataBase();
	if (_rad < 2.0 ) 
	{	cout << "Error reported from microEnvironment::microEnvironment" << endl;
		cout << "Radius of " << _rad << " is nonsense value!" << endl;
		cout << "Defaulting to 4.2 Angstroms" << endl;
		itsCriticalRadius = 4.2;
	}
	else
	{
		itsCriticalRadius = _rad;
	}
	itsResidueSkippingNumber = _skip;
	pItsDB = new microEnvDB(itsCriticalRadius, itsResidueSkippingNumber);
	pItsLoResPMF = new pmf("PMF_lores_2.dat");
	itsDistDepFlag = false;
	if (!pItsDB)
	{	cout << "Error reported from microEnvironment::microEnvironment" << endl;

		exit (1);
	}

	itsEnergyPerAtom.resize(0);	
	itsEnergySum = 0.0;
	theAtomTypes.resize(0);
	environments.resize(0);
}

microEnvironment::~microEnvironment()
{
		delete pItsDB;
		delete pItsLoResPMF;
}

void microEnvironment::initialize()
{
	itsEnergyPerAtom.resize(0);
	itsEnergySum = 0.0;
	theAtomTypes.resize(0);
	environments.resize(0);
}

double microEnvironment::calculateEnergy()
{  // NOTE: this can be done in a much smarter way in the future...
   // 1) no need to run through first loop if atoms are unchanged
   // 2) don't need to recalculate ALL energies if only one change is made

	itsEnergyPerAtom.resize(0);
	itsEnergySum = 0.0;
	environments.resize(0);
	theAtomTypes.resize(0);
	atomIterator it1(static_cast<protein*>(pItsMol));
	atomIterator it2(static_cast<protein*>(pItsMol));
	atom* pAtom1;
	atom* pAtom2;
	residue* pRes1;
	vector <UInt> atomEnv;
	double atEnergy;

	//First, run through and set atomtypes 
	for (; !(it1.last()); it1++) 
	{	if ( (pAtom1 = it1.getAtomPointer()) )
		{   if ( (pRes1 = it1.getResiduePointer()) )
			{
				vector< int > tempvector;
				string aTypeName = pAtom1->getName();
				UInt rTypeIndex = pRes1->getTypeIndex();
				int index1 = residue::dataBase[rTypeIndex].getAtomIndexOf(aTypeName);
				int tempint = residue::dataBase[rTypeIndex].getAtomEnergyTypeDefinition(index1,2);
				tempvector.push_back(tempint);
				tempint = residue::dataBase[rTypeIndex].getAtomEnergyTypeDefinition(index1,3);
				tempvector.push_back(tempint);
				theAtomTypes.push_back(tempvector);
			}
		}
	}
#ifdef MICROENV_DEBUG_ATOM_TYPES
	cout << "Atom Types as assigned by microEnvironment" << endl;
	for (UInt i=0; i<theAtomTypes.size(); i++)
	{
		for (UInt j=0; j<theAtomTypes[i].size(); j++)
		{	cout << theAtomTypes[i][j] << " ";
		}
		cout << endl;
	}
#endif
	it1.initialize();
	double dist;
	UInt counter1 = 0;
	for (; !(it1.last());it1++)
	{	atomEnv.resize(0);
		UInt counter2 = 0;
		for (UInt i=0; i<4; i++)
			{	atomEnv.push_back(0);	}
		it2.initialize();
		double atDDEnergy = 0.0;
		for (; !(it2.last()); it2++)
		{	bool relevant = true;
			pAtom1 = it1.getAtomPointer();
			pAtom2 = it2.getAtomPointer();
			if (it2.getResidueIndex() < it1.getResidueIndex())
			{	if (it2.getResidueIndex() > it1.getResidueIndex() - itsResidueSkippingNumber)
				{	relevant = false;
				}
				if ( (it1.getResidueIndex() - it2.getResidueIndex() == 1) 
						&& ( pAtom1->getName() == "C")
						&& ( pAtom2->getName() == "N")
						&& itsResidueSkippingNumber == 1 )
				{	relevant = false;
				}
			}
			else
		    {	if (it2.getResidueIndex() < it1.getResidueIndex() + itsResidueSkippingNumber)
				{	relevant = false;
				}
				if ( (it2.getResidueIndex() - it1.getResidueIndex() == 1)
					&& ( pAtom2->getName() == "C")
					&& ( pAtom1->getName() == "N")
					&& itsResidueSkippingNumber == 1 )
				{   relevant = false;
				}
			}

			if (it2.getChainIndex() != it1.getChainIndex())
				relevant = true;
    
			if (itsResidueSkippingNumber == 0)
				relevant = true;

			if (relevant)
			{	dist = pAtom1->distance(pAtom2);
				if (dist <= itsCriticalRadius)
				{
					if (theAtomTypes[counter2][1] >=0 && theAtomTypes[counter1][0] >= 0)
					{
						atomEnv[theAtomTypes[counter2][1]]++;
						if (itsDistDepFlag)
						{	
							double origPMFScale = pmf::getScaleFactor();
							pmf::setScaleFactor(1.0);
							atDDEnergy += pItsLoResPMF->getEnergy(theAtomTypes[counter1][0],theAtomTypes[counter2][1],dist);
							pmf::setScaleFactor(origPMFScale);
						}
					}
				}
			}
			counter2++;
		}
#ifdef MICROENV_DEBUG_ATOM_TYPES
		cout << (it1.getResiduePointer())->getType() << " ";
		cout << pAtom1->getName() << "  : ";
		cout << theAtomTypes[counter1][0] << " : ";
		for (UInt i=0;i<4;i++)
	    	cout << atomEnv[i] << " ";	
		cout << endl;
#endif
		environments.push_back(atomEnv);
		atEnergy = 0.0;
		if (theAtomTypes[counter1][0] >=0)
		{
			atEnergy = (itsScaleFactor * pItsDB->getEnergy(theAtomTypes[counter1][0],atomEnv)) + atDDEnergy;
		}
		itsEnergyPerAtom.push_back(atEnergy);
		itsEnergySum += atEnergy;
		counter1++;
	}
#ifdef MICROENV_DEBUG
	cout << "Energy per Atom" << endl;
	printFlagStatus();
	cout << "---Start----" << endl;
	for (UInt i=0; i< itsEnergyPerAtom.size(); i++)
	{	cout << "atom " << i << " " << itsEnergyPerAtom[i] << endl;
	}
	cout << "-----End-----" << endl;
#endif
	return itsEnergySum;
}

double microEnvironment::calculateEnergy(molecule* _pMol)
{	
	pItsMol = _pMol;
	itsEnergySum = 0.0;
	//cout << itsEnergySum;
	double energy = calculateEnergy();
	//cout << " " << itsEnergySum << endl;
	return energy;
}

void microEnvironment::printFlagStatus() const
{
	cout << "itsScaleFactor = " << itsScaleFactor << endl;
	cout << "itsDistDepFlag = " << itsDistDepFlag << endl;
	cout << "microEnvDB:useSingleBodyEnergy = " << microEnvDB::useSingleBodyEnergy << endl;
}

vector<double> microEnvironment::getEnergyPerAtom() const
{
	return itsEnergyPerAtom;
}

vector<UInt> microEnvironment::getEnvAssignedAtomTypes() const
{	vector<UInt> theList;
	for (UInt i=0; i< theAtomTypes.size(); i++)
	{	theList.push_back(theAtomTypes[i][0]);
	}
	return theList;
}
