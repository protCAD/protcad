#include <assert.h>
#include <vector>
#include <string.h>
#include <iostream>
#include <stdio.h>
#include "typedef.h"

#ifndef PDB_ATOM_RECORD_H
#define PDB_ATOM_RECORD_H

class PDBAtomRecord
{
public: 
	PDBAtomRecord();
	PDBAtomRecord(string& _pdbAtomLine);
	PDBAtomRecord(const PDBAtomRecord& _otherRecord);
	~PDBAtomRecord();


	string getAtomName() const {return atomName;} 
	string getResName() const {return resName;} 
	string getICode() const {return iCode;} 
	string getAltLoc()  const {return altLoc;}
	string getSegID()  const {return segID;}
	int getResSeq() const {return resSeq;} 
	int getSerial() const {return serial;}
	double getOccupancy() const {return occupancy;}
	double getTempFactor() const {return tempFactor;}
	dblVec getAtomCoord() const {return atomCoord;}
	string getChainID() const {return chainID;}
	string getCharge() const { return charge;}
	string getElement() const { return element;}
	bool getHetflag() const {return hetflag;}

private:
	void convert(string& _pdbAtomLine);

private:
	string atomName;
	string resName;
	string iCode;
	string altLoc;
	string chainID;
	string segID;
	string element;
	string charge;
	int resSeq;
	int serial;
	double occupancy;
	double tempFactor;
	dblVec atomCoord;
	bool hetflag;

};
#endif
