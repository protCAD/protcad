#include <vector>
#include <string.h>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include "typedef.h"
#include "pdbEnums.h"
#include "generalio.h"

#ifndef PDB_DATA_H
#define PDB_DATA_H

#define SIZE_OF_RECORDNAME 6

class pdbData
{
	//class pdbParse;
	class pdbParse
	{
	public:
        	pdbParse()
        	{       header = new vector<string>(0);
                	parseRule = new vector<unsigned int>(0);
        	}

		int pdbHeaderCompare(const string& _recordName)
		{	for( unsigned int i=0; i<header->size() ; i++)
			{	if( _recordName == (*header)[i] )
				{	return 0; // a match found
				}
			}
			return 1; // no match
		}

        	~pdbParse()
        	{       delete header;
			delete parseRule;
        	}

        	vector<string>* header;
        	vector<unsigned int>* parseRule;
	};
public: 
	pdbData();
	virtual ~pdbData();

	//data accessors
	void pdbPutLine(ofstream& _outFile);
	virtual int  pdbGetLine(ifstream& _inFile);
	bool getDataValid() const { return dataValid;}
	string getRecordName() { return recordName; }
	
	int pdbParseAllLine();
	int pdbParseSpecificLine();
	int pdbDoTheParsing();
	virtual void conversion() {};

protected:
	int setCurrentPdbParse();
	void updateCurrentPdbParse(pdbData::pdbParse* _temp);

private:	
	static void buildPdbParseBase();
	static pdbData::pdbParse* getPdbParseBaseItem( const unsigned int _itemIndex );
	static unsigned int getPdbParseBaseSize();

protected:	
   	vector<string>* item;
	string recordName;
	bool dataValid;
	pdbData::pdbParse* currentPdbParse;
	vector<unsigned int>* currentParseRule;


	// static variables
	static string  currentLine;
	static vector<pdbData::pdbParse*>* pdbParseBase;
	static bool pdbParseBaseBuilt;
	static unsigned int howMany;
};


class pdbAtom : public pdbData
{
public:
	pdbAtom();

	// data accessors
	string getItem( pdbAtomParse _itemIndex ) const {return (*item)[_itemIndex];} 
	// pdbAtomParse is a defined enum type
	double getOccupancy() const {return occupancy;}
	double getTempFactor() const {return tempFactor;}
	dblVec getAtomCoord() const {return atomCoord;}
	unsigned int getSerial() const {return serial;}
	int getResSeq() const {return resSeq;}
	char  getChainID() const {return chainID;}
	double getCharge() const { return charge;}
	int pdbGetLine(ifstream& _inFile);	
	void conversion();

private:
	double occupancy;
	double tempFactor;
	double charge;
	char  chainID;
	dblVec atomCoord;
	unsigned int serial;
	int resSeq;

};


#endif
