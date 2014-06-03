// filename: ruler.h
#include "assert.h"
#include "typedef.h"
#include <iostream>
#include <vector>
#include "generalio.h"
#include "molecule.h"
#include "protein.h"
#include "chain.h"
#include "residue.h"
#include "atom.h"

#ifndef RULER_H
#define RULER_H
class ruler
{
public:
	// Constructor and Destructor declaration
	ruler();
	~ruler();

	// Accessors
	void setStationaryMolecule(molecule* _mol1) {pItsStationaryMolecule = _mol1;}
	void setMobileMolecule(molecule* _mol2) {pItsMobileMolecule = _mol2;}
	void setAtomListForStationary(vector<vector<UInt> > _atomList1)
		{atomList1 = _atomList1;}
	void setAtomListForMobile(vector<vector<UInt> > _atomList2)
		{atomList2 = _atomList2;};
	void setWeights(vector<double> _weights)
		{itsWeights = _weights;}

	void appendToList(UInt _listnum, UInt _startRes, UInt _endRes,
	   UInt _chain, string _atomNames);
	void superimposeProteins ();
	double getrmsd() const {return itsrmsd;}
	dblVec getCentroidOfStationary() const {return itsCentroid1;}
	dblVec getCentroidOfMobile() const {return itsCentroid2;}
	dblMat getRotationMatrix() const {return itsRotationMatrix;}
	dblVec getAxisOfRotation();

private:
	molecule* pItsStationaryMolecule;
	molecule* pItsMobileMolecule;
	dblMat itsRotationMatrix;
	dblVec itsCentroid1;
	dblVec itsCentroid2;
	vector<vector<UInt> > atomList1;
	vector<vector<UInt> > atomList2;
	vector<double> itsWeights;
	double itsrmsd;
};
#endif
