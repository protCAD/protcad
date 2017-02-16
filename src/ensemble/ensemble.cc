// filename: ensemble.cc
// contents: class ensemble implementation

#include "ensemble.h"
#include "residue.h"

double ensemble::itscutoffDistance=8.0;

ensemble::ensemble()
{	//cout<< "default ensemble constructor called" << endl;
	itsMolecules.resize(0);
	itsEnergy = 0.0;
	//cout<< "starting build" << endl;
	itsLastModifiedMolecule = -1;
	itsIndMolAndChainList.resize(0);
	itsMolAndChainLinkageMap.resize(0);
	//cout <<"Exiting default ensemble constructor..."<<endl;
}

ensemble::~ensemble()
{	//cout<< "ensemble destructor called " << endl;
	for (UInt i=0; i < itsMolecules.size(); i++)
	{	delete itsMolecules[i];
	}
}

void ensemble::add(molecule* _pMolecule)
{	//cout << "\nIn Ensemble add()...\n";
	itsMolecules.push_back(_pMolecule);
}

void ensemble::remove(molecule* _pMolecule)
{
	for (UInt i = 0; i < itsMolecules.size(); i ++)
	{
		if (_pMolecule == itsMolecules[i])
		{
			vector <molecule*> temp;
			for (UInt j = 0; j < itsMolecules.size(); j++)
			{
				if (i !=j ) temp.push_back(itsMolecules[j]);
			}
			itsMolecules = temp;
			return;
		}
	}
	cout << "ERROR in ensemble::remove(molecule* ...)\n\tMolecule not found in this ensemble." << endl;
	return;
}

int ensemble::mutate(vector <int> _position, UInt _resType)
{
	int moleculeToModify = _position[0];
	int result = -3;
	if (moleculeToModify >=0 && moleculeToModify < (int)itsMolecules.size())
	{
		cout << "MUTATING MOL:  " << moleculeToModify << " ";
		result = itsMolecules[moleculeToModify]->mutate(_position, _resType);
		if (result == 1)
		{
			itsLastModifiedMolecule = moleculeToModify;
		}
		else 
		{
			// result = -1 (abort signal)
			cout << "Ensemble level has detected abort from mutate() ..." << endl;
		}
	}
	else	// else mol id is out of range
	{
		result = -2;
		cout << "molecule position specified:  " << _position[0] << " is out of range..." << endl;
	}		
	
	return result;
}


int ensemble::symmetryLinkMolAndChain(UInt _indMol, UInt _indChain, UInt _slaveMol, UInt _slaveChain)
{
	cout << "Linking mol-chain" << _slaveMol << "," << _slaveChain << " to the independent position: " << _indMol << "," << _indChain << endl;

	bool molAndChainInList = false;
	UInt independentPosition = 0;
	for (UInt i = 0; i < itsIndMolAndChainList.size(); i++)
	{
		if (itsIndMolAndChainList[i][0] == _indMol && itsIndMolAndChainList[i][1] == _indChain)
		{
			// molecule already exists in independent list
			molAndChainInList = true;
			independentPosition = i;
		}
	}
	if (!molAndChainInList)
	{
		// add master to the independent list
		vector <UInt> tmpVec; 
		tmpVec.push_back(_indMol);
		tmpVec.push_back(_indChain);
		itsIndMolAndChainList.push_back(tmpVec);
		independentPosition = itsIndMolAndChainList.size() - 1;

		// create new entry in linkage map and add slave mol and chain to position 0
		vector < vector <UInt> > tmpSlaveVec;
		tmpSlaveVec.resize(1);
		tmpSlaveVec[0].push_back(_slaveMol);
		tmpSlaveVec[0].push_back(_slaveChain);
		itsMolAndChainLinkageMap.push_back(tmpSlaveVec);
		printLinkageInfo();
		return 1;
	}
	else
	{
		vector<UInt> tmpSlaveVec;
		tmpSlaveVec.push_back(_slaveMol);
		tmpSlaveVec.push_back(_slaveChain);
		itsMolAndChainLinkageMap[independentPosition].push_back(tmpSlaveVec);
		printLinkageInfo();
		return 1;
	}
}

void ensemble::printLinkageInfo()
{
	for (UInt i = 0; i < itsIndMolAndChainList.size(); i ++)
	{
		cout << "---------------------------" << endl;
		cout << "Links to molecule " << itsIndMolAndChainList[i][0] << " chain " << itsIndMolAndChainList[i][1] << ":" << endl;
		cout << "---------------------------" << endl;
		for (UInt j = 0; j < itsMolAndChainLinkageMap[i].size(); j++)
		{
			cout << "\tMolecule " << itsMolAndChainLinkageMap[i][j][0] << " chain " << itsMolAndChainLinkageMap[i][j][1] << endl;
		}
		cout << endl;
	}
	return;
}

int ensemble::mutateWithSymmetry(vector <int> _position, UInt _resType)
{
	int moleculeToModify = _position[0];
	int chainToModify = _position[1];
	int independentPosition = -1;
	int result;
	for (UInt i = 0; i < itsIndMolAndChainList.size(); i++)
	{
		if (moleculeToModify == (int)itsIndMolAndChainList[i][0] && chainToModify == (int)itsIndMolAndChainList[i][1])
		{
			independentPosition = i;
		}
	}
	if (independentPosition == -1)
	{
		cout << "ERROR in ensemble::mutateWithSymmetry ... position not found in independent chain list." << endl << "\tSystem may be configured incorrectly." << endl;
		return -2;
	}

	result = mutate(_position, _resType);
	if (result == -1)
	{
		cout << "ERROR in ensemble::mutateWithSymmetry ... Ensemble level dectected abort when modifying indepedendent molecule and chain: " << _position[0] << " " << _position[1] << endl;
		return result;
	}

	for (UInt i = 0; i < itsMolAndChainLinkageMap[independentPosition].size(); i++)
	{
		vector <int> tempPos;
		tempPos.resize(0);
		tempPos.push_back(itsMolAndChainLinkageMap[independentPosition][i][0]);				// add symmetry linked molecule
		tempPos.push_back(itsMolAndChainLinkageMap[independentPosition][i][1]);				// ... and chain
		tempPos.push_back(_position[2]);													// ... and target resiude

		result = mutate(tempPos, _resType);
		if (result == -1)
		{
			cout << "ERROR in ensemble::mutateWithSymmetry ... Ensemble level dectected abort when modifying symmetry linked molecule and chain: " << itsMolAndChainLinkageMap[independentPosition][i][0] << " " << itsMolAndChainLinkageMap[independentPosition][i][1] << endl;
			cout << "\tSystem may be configured incorrectly." << endl;
			return result;
		}
	}
	return 1;
}

int ensemble::modify(ran& _ran, vector <int> _position)
{
	int result = -2;
	int moleculeToModify = _position[0];
	if (moleculeToModify >=0 && moleculeToModify < (int)itsMolecules.size())
	{
		result = itsMolecules[moleculeToModify]->modify(_ran, _position);
		if (result == 1)
		{
			itsLastModifiedMolecule = moleculeToModify;
		}
	        else
	        {
	            result = -1; 
	            cout << "Ensemble level has detected abort..." << endl;
	        }
	}
	else
	{
		result = -2;
		cout << "Error in ensemble::modify() ... molecule position specified: " << _position[0] << "  is out of range." << endl;
	}
	
	return result;
}

int ensemble::modify(ran& _ran)
{	int result;
	int moleculeToModify = chooseMolecule(_ran);
	if (moleculeToModify >=0)
	{	result = itsMolecules[moleculeToModify]->modify(_ran);
		if (result == 1)
		{	itsLastModifiedMolecule = moleculeToModify;
		}
		else
		{
			result = -1;
			cout << "Ensemble level has detected abort..." << endl;
		}
	}
	else
	{	result = -2;
	}
	return result;
}

// selects a position to modify without actually making any changes to the system
vector <int> ensemble::chooseNextTargetPosition(ran& _ran)
{
	vector <int> position;
	position.resize(0);
	position.push_back(chooseMolecule(_ran));
	vector <int> chainAndRes = itsMolecules[position[0]]->chooseNextTargetPosition(_ran);
	position.push_back(chainAndRes[0]);
	position.push_back(chainAndRes[1]);

	return position;
}

UInt ensemble::chooseNextMutationIdentity(ran& _ran, vector <int> _position)
{	
	if (_position[0] >=0 && _position[0] < (int)itsMolecules.size())
	{
		return itsMolecules[_position[0]]->chooseNextMutationIdentity(_ran, _position);
	}
	else cout << "ERROR in chooseNextMutationIdentity at ensemble level ... molecule ID " << _position[0] << "is out of range." << endl;
	return 0;
}

void ensemble::setupSystem(ran& _ran)
{
	for (UInt i=0; i < itsMolecules.size(); i++)
	{
		itsMolecules[i]->setupSystem(_ran);
	}
}

void ensemble::saveState(string& _filename)
{
	for (UInt i=0; i < itsMolecules.size(); i++)
	{
		itsMolecules[i]->saveState(_filename);
	}	
}

double ensemble::getVolume(UInt _method)
{
	itsVolume = 0.0;
	
	for (UInt i = 0; i < itsMolecules.size(); i ++)
	{
		itsVolume += itsMolecules[i]->getVolume(_method);
	}
	return itsVolume;
}

vector <int> ensemble::getLastModification()
{
	vector <int> chainAndResiduePosition = itsMolecules[itsLastModifiedMolecule]->getLastModification();
	vector <int> position;
	position.resize(0);
	position.push_back(itsLastModifiedMolecule);   // add last modified molecule
	position.push_back(chainAndResiduePosition[0]); // add last modified chain
	position.push_back(chainAndResiduePosition[1]); // add last modified residue

	return position;
}

double ensemble::getPositionEnergy(vector<int> _position)
{
	if (_position.size() == 3) // if position contains molecule id, chain id and res id
	{
		return itsMolecules[_position[0]]->getPositionEnergy(_position);	

	}
	else cout << "ERROR in ensemble::getPositionEnergy ... position variable incorrectly specified." << endl;
	return 0.0;
}

double ensemble::energy()
{	itsEnergy = 0.0;
	// This is the PMF energy component + the rotamer energy component
	double pairwisePart;
	for (UInt i=0; i<itsMolecules.size(); i++)
	{	pairwisePart = (itsMolecules[i])->intraEnergy();
		itsEnergy += pairwisePart;
		cout << "PairwisePart = " << pairwisePart << " ";
	}
	return itsEnergy;
}

void ensemble::setAllEnsembleCutoffDistance()
{
    residue::setCutoffDistance(itscutoffDistance);
}


void ensemble::acceptModification()
{	
	if (itsLastModifiedMolecule >=0)
	{
		molecule* pMol = itsMolecules[itsLastModifiedMolecule];
		pMol->acceptModification();
		itsLastModifiedMolecule = -1;
	}
}

void ensemble::rejectModification()
{	if (itsLastModifiedMolecule >=0)
	{	molecule* pMol = itsMolecules[itsLastModifiedMolecule];
		pMol->rejectModification();
		itsLastModifiedMolecule = -1;
	}
}

int ensemble::chooseMolecule(ran& _ran)
{	UInt moleculeSize = itsMolecules.size();
	if (moleculeSize !=0)
	{	return int(_ran.getNext()*moleculeSize);
	}
	else
	{	cout << "Error reported from ensemble::chooseMolecule" << endl;
		cout << "No molecules in ensemble!" << endl;
		return -1;
	}
}
molecule* ensemble::getMoleculePointer (UInt _index)
{
	if (_index < itsMolecules.size())
	{	return itsMolecules[_index];
	}
	return 0;
}

