#include "rotamer.h"

double rotamer::itsScaleFactor = 1.0;

rotamer::rotamer()
{
	itsAngles.resize(0);
	itsEnergy= 0.0;
}

rotamer::rotamer(const UIntVec _angles, const double _energy)
{
	itsAngles = _angles;
	itsEnergy= _energy;
}

rotamer::rotamer(const rotamer& _otherRotamer)
{	// deep copy constructor
	itsAngles = _otherRotamer.itsAngles;
	itsEnergy = _otherRotamer.itsEnergy;
}

rotamer::~rotamer()
{
}

DouVec rotamer::getAngles() const
{	DouVec angles;
	for(UInt i=0; i<itsAngles.size(); i++)
	{	switch( itsAngles[i] )
		{	case 1:	angles.push_back(60);
					break;
			case 2: angles.push_back(180);
					break;
			case 3: angles.push_back(-60);
					break;
			default:break;
		}
	}
	return angles;
}

double rotamer::getEnergy() const
{	return itsScaleFactor * itsEnergy;
}

void rotamer::print() const
{	cout << "angles = ";
	for(UInt i=0; i<itsAngles.size(); i++)
	{	cout << itsAngles[i] << '\t';
	}
	cout << "energy = " << itsEnergy << endl;
}
