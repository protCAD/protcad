#include "rotamerLib.h"

UInt rotamerLib::howMany = 0;

rotamerLib::rotamerLib(const UInt _numOfBpt)
{	itsRotamers.resize(_numOfBpt);
	howMany++;
}

rotamerLib::rotamerLib(const rotamerLib& _otherRotamerLib)
{
	itsRotamers = _otherRotamerLib.itsRotamers;
}

rotamerLib::~rotamerLib()
{	howMany--;
}

DouVec rotamerLib::getAngles(const UInt _bpt, const UInt _rotamer) const
{	

	vector <double> angles;		
	ASSERT(_bpt < itsRotamers.size());
	ASSERT(_rotamer < itsRotamers[_bpt].size());
	return itsRotamers[_bpt][_rotamer].getAngles();
}

double rotamerLib::getEnergy(const UInt _bpt, const UInt _rotamer) const
{	
		ASSERT(_bpt < itsRotamers.size());
		//ASSERT(_rotamer < itsRotamers[_bpt].size());
		double energy = itsRotamers[_bpt][_rotamer].getEnergy();
		//cout << "rotamerLib::getEnergy(" << _bpt << ", " << _rotamer << ") = " << energy << endl;
		return energy;
}

bool rotamerLib::rotamersExist(const UInt _bpt, const UInt _rotamer) const
{	if (_bpt < itsRotamers.size())
	{	if (_rotamer < itsRotamers[_bpt].size())
		{	return true;
		}
	}
	return false;
}

void rotamerLib::addRotamer(const StrVec& _strVec)
{	

	//cout << "In rotamerLib::addRotamer" << endl;
	int nonChi = 3; // things that aren't chi are branch point, population and energy - ergo 3 nonchis (monchichi?)
	int chis;
	chis = int(_strVec.size()) - nonChi;
	if(chis <= 0)
	{	cerr << "Data Format Invalid !!!" << endl;
		exit(1);
	}

	// first column is the branch point 
	UInt locIndex = 0;
	int bpt;
	sscanf(_strVec[locIndex].c_str(), "%i", &bpt);
	bpt --;

	// starting second column is the chi definitions
	UInt angle;
	UIntVec tempAngles;
	tempAngles.resize(0);
	for(int i=0; i<chis; i++)
	{	locIndex ++;
		int anglei;
		sscanf(_strVec[locIndex].c_str(), "%i", &anglei);
		angle = (UInt)anglei;
		tempAngles.push_back(angle);
	}

	// energy is the second column after all the chis
	locIndex += 2;
	if( locIndex > _strVec.size() )
	{	cerr << "Data Format Invalid !!!" << endl;
		exit(1);
	}
	double energy;
	sscanf(_strVec[locIndex].c_str(), "%lf", &energy);

	rotamer tempRot(tempAngles,energy);
	itsRotamers[bpt].push_back(tempRot);
}

void rotamerLib::print() const
{	for(UInt i=0; i<itsRotamers.size(); i++)
	{	cout << " bpt = " << i << endl;
		for(UInt j=0; j<itsRotamers[i].size(); j++)
		{	itsRotamers[i][j].print();
		}
	}
}
