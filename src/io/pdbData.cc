#include "pdbData.h"

string pdbData::currentLine="";
bool pdbData::pdbParseBaseBuilt = false;
vector<pdbData::pdbParse*>* pdbData::pdbParseBase = new vector<pdbData::pdbParse*>(0);
unsigned int pdbData::howMany = 0;

pdbData::pdbData()
{	if( !pdbParseBaseBuilt )
	{       buildPdbParseBase();
	}
	item = new vector<string>(0);
	recordName = "UNKNWN";
	currentPdbParse = 0; // a null pointer to start with
	howMany++;
}

pdbData::~pdbData()
{	delete item;
	howMany--;
}

// pdbParseBase Related Operations

void pdbData::buildPdbParseBase()
{	
	string evname = "PROTCADDIR";
	string path = getEnvironmentVariable(evname);

	path += "/data/pdb/";
	string filename = path + "PDBHEADER";
	
	ifstream pdbHeader(filename.c_str());
	if( !pdbHeader)
	{	cout << filename << " cannot be found or open " << endl;
		exit(1);
	}
	
	string oneLine;
	unsigned int tempInt;
	pdbData::pdbParse* buf;

	bool noHeader = true;
	unsigned int counter = 0; // number of group of headers

	while( pdbHeader >> oneLine )
	{	if( noHeader)
		{	buf = new pdbData::pdbParse;
		}
		if( oneLine != "******" )
		{	noHeader = false;
			counter ++;
			buf->header->push_back(oneLine);
			continue;
		}
		if( !noHeader)
		{
			filename = path + (*(buf->header))[0]; // get the first header
			ifstream parse(filename.c_str());
		
			for(; (*(buf->header))[0].size() < SIZE_OF_RECORDNAME ;)
			{	(*(buf->header))[0].append(" ");
			} // add enough space sothat it is six letters' wide
		
			if(!parse)
			{       cout << filename << " cannot be found or open " << endl;
				exit(1);
			}

			while( parse >> tempInt )
			{	buf->parseRule->push_back(tempInt);
			}

			pdbParseBase->push_back(buf);
			parse.close();
			parse.clear();
			noHeader = true;
        	}
	}

	pdbHeader.close();
	pdbHeader.clear();
	if( counter == 0)
	{	cout << "no header found in the file " << endl;
		exit(1);
	}
	pdbParseBaseBuilt = true;
}


unsigned int pdbData::getPdbParseBaseSize()
{	return pdbParseBase->size();
}

pdbData::pdbParse* pdbData::getPdbParseBaseItem( const unsigned int _itemIndex )
{	if( _itemIndex < getPdbParseBaseSize() )
	{	return (*pdbParseBase)[_itemIndex];
	}
	else
	{	cout << " _itemIndex " << _itemIndex <<" is incompatible with the pdbParseBase ";
		cout << endl;
		exit (1);
	}
}

// CurrentPdbParse Related Operations

int pdbData::setCurrentPdbParse()
{	// if do not need to be updated
	if( currentPdbParse && ! currentPdbParse->pdbHeaderCompare( recordName ) )
	{	return 0;
	}

	// if need to be updated
	for(unsigned int i=0; i< getPdbParseBaseSize() ; i++)
	{	pdbData::pdbParse* temp = (*pdbParseBase)[i];
		if( !temp->pdbHeaderCompare( recordName) )
		{	currentPdbParse = temp;
			return 0;
		}
	}
	
	// if not found any matching parsing rule
	return 1;
}

// Line Related Operations

int pdbData::pdbGetLine(ifstream& _inFile)
{	//cout << " entering pdbGetLine " << endl;

	for(unsigned int i=0 ; getline(_inFile, currentLine, '\n');i++ )
 	{	pdbParseAllLine();
		if( dataValid )
		{	
			return 1; // normal return at end_of_line
		}
		else
		{	
			continue; // keep reading

		}
	}
	return 0; // reach the end_of_file
}	

int pdbData::pdbParseAllLine()
{	return 0;
}

int pdbData::pdbParseSpecificLine()
{	unsigned int lineLength = currentLine.size();
	item->resize(0); // clean up the old data
	dataValid = false; // reset the dataValid flag

	if( SIZE_OF_RECORDNAME <= lineLength )
	{	recordName = currentLine.substr(0,SIZE_OF_RECORDNAME); // get the recordName first
		if( !currentPdbParse->pdbHeaderCompare(recordName) )
		{	dataValid = true;
			pdbDoTheParsing();
			conversion();
			return 0; // normal return
		}
		else
		{	return -1; // to skip
		}
	}
	else
	{	cout << " the line is too short and has incomplete recordName " << endl;
		cout << " not parsed and skipped " << endl;
		return -1; // to skip
	}
}

void pdbAtom::conversion()
{	double tempDouble;
	unsigned int tempU;
	
	sscanf( (*item)[sOccupancy].c_str(), "%lf", &tempDouble);
	occupancy = tempDouble;
	sscanf( (*item)[sSerial].c_str(), "%u", &tempU);
	serial = tempU;
	sscanf( (*item)[sTempFactor].c_str(), "%lf", &tempDouble);
	tempFactor = tempDouble;
	sscanf( (*item)[sChainID].c_str(),"%c", &chainID);
	sscanf( (*item)[sResSeq].c_str(), "%u", &tempU);
	resSeq = tempU;
	sscanf( (*item)[sXCoord].c_str(), "%lf", &tempDouble);
	atomCoord[0] = tempDouble;
	sscanf( (*item)[sYCoord].c_str(), "%lf", &tempDouble);
	atomCoord[1] = tempDouble;
	sscanf( (*item)[sZCoord].c_str(), "%lf", &tempDouble);
	atomCoord[2] = tempDouble;
	sscanf( (*item)[sCharge].c_str(), "%lf", &tempDouble);
	charge = tempDouble;
}


int pdbData::pdbDoTheParsing()
{	unsigned int pos = SIZE_OF_RECORDNAME; // update the current positiona
	int incompleteLine = 0;
	unsigned int lineLength = currentLine.size();
	string empty = " ";
	string rawStr;
	string strCh;
	strCh.resize(1);
	string strBuf;

	vector<unsigned int>* temp = currentPdbParse->parseRule;

	for( unsigned int i=0; i< temp->size(); i++)
	{	if( (pos+(*temp)[i]) <= lineLength )
		{	rawStr = currentLine.substr(pos, (*temp)[i]);
			strBuf.resize(0);
			// filter away the spaces and tabs
			for(unsigned int j=0; j<rawStr.size(); j++)
			{	if(rawStr[j] !=' ' && rawStr[j] !='\t' && rawStr[j] != '\n')
				{	strCh[0] = rawStr[j];
					strBuf.append(strCh);
				}
			}
			if(strBuf.size() == 0) 
			{	strBuf = empty;
			}
			item->push_back(strBuf);
			pos += (*temp)[i];
		}
		else
		{	item->push_back(empty);
			incompleteLine = 1;
		}
	}
	return incompleteLine;
}

void pdbData::pdbPutLine(ofstream& _outFile)
{	
}

pdbAtom::pdbAtom()
{	serial = 0;
	resSeq = 0;
	atomCoord.newsize(3);
	atomCoord[0] = 0.0;
	atomCoord[1] = 0.0;
	atomCoord[2] = 0.0;
	chainID = ' ';
	occupancy = 0.0;
	tempFactor = 0.0;
	charge = 0.0;
	dataValid = false;
	recordName = "ATOM  ";
	setCurrentPdbParse();
}

int pdbAtom::pdbGetLine(ifstream& _inFile)
{	for(unsigned int i=0 ; getline(_inFile, currentLine, '\n');i++ )
 	{	// ignore the blank line
		if(currentLine.size() != 0)
		{	pdbParseSpecificLine();
		}
		if( dataValid )
		{	return 1; // normal return at end_of_line
		}
		else
		{	continue; // keep reading
		}
	}
	return 0; // reach the end_of_file
}	
