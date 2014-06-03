// filename: molecule.h
// contents: base class molecule is defined

#include "assert.h"
#include <string.h>
#include <vector>
#include <ran.h>
#include "typedef.h"
#include "point.h"
#include "enums.h"
//#include "CMath.h"

#ifndef MOLECULE_H
#define MOLECULE_H

class molecule
{
public:
	// Constructor and Destructor declaration
	molecule();
	molecule(const string& _name);
	molecule(const molecule& _rhs);
	virtual ~molecule() = 0;
private:
	void initialize();

public:
	// Accessors

	//virtual void monteCarloMinimize(UInt _maxSteps, double _threshold) = 0;

	virtual void translate(const UInt _index, const dblVec& _dblVec) = 0;
	virtual void translate(const UInt _index, const double _x,const double _y,const double _z) = 0;
	virtual void transform(const UInt _index, const dblMat& _dblMat) = 0;
	virtual void rotate(const UInt _index, const axis _axis, const double _theta) = 0;
	virtual void rotate(const UInt _index, const point& _point, const dblVec& _R_axis, const double _theta) = 0;

	virtual double getPositionEnergy(vector <int> _position) = 0;
	virtual vector<int> getLastModification() = 0;	
	virtual vector<int> chooseNextTargetPosition(ran& _ran) = 0;
	virtual UInt chooseNextMutationIdentity(ran& _ran, vector <int> _position) = 0;
//	virtual void findLowestRotamerWithSymmetry(vector <int> _position) = 0;

	virtual void setMoleculeType(UInt _type) {itsMolType = _type;}
	virtual UInt getMoleculeType() {return itsMolType;}

	virtual void acceptModification() = 0;
	virtual void rejectModification() = 0;
	virtual int modify(ran& _ran) = 0;
	virtual int modify(ran& _ran, vector <int> _position) = 0;
	virtual int mutate(vector <int> _position, UInt _resType) = 0;
	virtual void setupSystem(ran& _ran) = 0;
	virtual void saveState(string& _filename) = 0;
	virtual void resetAllBuffers() = 0;
	virtual double getVolume(UInt _method) = 0; 
	virtual double intraEnergy() = 0;

protected:
	UInt itsMolType;
	string itsName;
};
#endif
