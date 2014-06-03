#include "pmf.h"

double pmf::itsScaleFactor = 1.0;
bool pmf::atomTypeNameStringsBuilt = false;
vector< vector< string > > pmf::atomTypeNameStrings;

pmf::pmf(string _fname)
{	// default constructor
#ifdef PMF_DEBUG
	cout << "called pmf::pmf() " << endl;
#endif
	// read the file until the end.....
	itsFileName = _fname;
	pairwiseEnergyData.resize(0);
	distanceBinUpperLimits.resize(0);
	buildDataBase();
	if (!atomTypeNameStringsBuilt)
	{
		readAtomTypeData();
	}
	//cout << " pmf database is built " << endl;
}

pmf::pmf(const pmf& _other)
{	
#ifdef PMF_DEBUG
	cout << "called pmf deep copy constructor " << endl;
#endif
	itsFileName = _other.itsFileName;
	pairwiseEnergyData = _other.pairwiseEnergyData;
	distanceBinUpperLimits = _other.distanceBinUpperLimits;
	stepSize = _other.stepSize;
	lowerDistanceLimit = _other.lowerDistanceLimit;	
	upperDistanceLimit = _other.upperDistanceLimit;
}

pmf::~pmf()
{
}

double pmf::getEnergy(const UInt _type1, const UInt _type2, const double _distance) const
{	int dbin = findDistanceBin(_distance);
	if (dbin == -1)
	{	return 0.0;
	}
	else if (dbin == -2)
	{	// eliminate possible numerical overflow
		// bounds check
		if (_type1 < pairwiseEnergyData.size())
		{      if (_type2 < pairwiseEnergyData[_type1].size())
			{      
				double fudgeFactor = 0.0;
				double E_at_lowest_dist = pairwiseEnergyData[_type1][_type2][0];
				double slope = (100.0 - E_at_lowest_dist)/(0.0 - lowerDistanceLimit);
				//cout << "slope = " << slope << endl;
				double energy = CMath::linearInterpolate(slope,100.0,_distance);
				return itsScaleFactor * (energy - fudgeFactor);
			}
		}
	}
	else
	{
		// bounds check
		if (_type1 < pairwiseEnergyData.size())
		{   if (_type2 < pairwiseEnergyData[_type1].size())
			{	return itsScaleFactor * (pairwiseEnergyData[_type1][_type2][dbin]);
			}
		}
	}
	cout << "PMF bounds error!" << endl;
	return 0.0;
}

int pmf::findDistanceBin(const double _distance) const
{
	if (_distance >= upperDistanceLimit)
	{	return -1;
	}
	else if(_distance < lowerDistanceLimit)
	{	return -2;
	}
	else
	{	// return binarySearch(low, high, _distance);
		return linearSearch(_distance);
	}
}

UInt pmf::binarySearch(const UInt low, const UInt high, const double _distance) const
{	if (low == high)
	{	return low;
	}
	else
	{	UInt mid = low + (high - low)/2;
		if (_distance > distanceBinUpperLimits[mid])
		{	return binarySearch(mid+1, high, _distance);
		}
		else
		{	return binarySearch(low, mid, _distance);
		}
	}
}

UInt pmf::linearSearch(const double _distance) const
{
	return UInt( (_distance - lowerDistanceLimit)/stepSize);
}

void pmf::buildDataBase()
{	
	string evname = "PROTCADDIR";
	string path = getEnvironmentVariable(evname);

	path += "/data/pmf/";
	string iFile = path + itsFileName;
	ifstream inFile;
	inFile.clear();
	string currentLine;
	StrVec parsedStrings;
	parsedStrings.resize(0);
	
	// read the pmf data file
	inFile.open(iFile.c_str());
	if (!inFile)
	{	cout << "Error: unable to open input file: ";
		cout << iFile << endl;
		exit (1);
	}

	while (getline(inFile, currentLine, '\n'))
	{	// ignore the comment line
		// comment line should start with #
		if(currentLine[0] != '#')
		{	parsedStrings=Parse::parse(currentLine);
			convertToEnergyDataBase(parsedStrings);
			parsedStrings.resize(0);
		}
	}

	inFile.close();
	inFile.clear();
	// read the distance bin data file
	string fileName = "PMFDBUL.dat";
	iFile = path + fileName;
	inFile.open(iFile.c_str());
	if(!inFile)
	{	cout << "Error: unable to open input file: ";
		cout << iFile << endl;
		exit(1);
	}

	while(getline(inFile, currentLine, '\n'))
	{	if(currentLine[0] != '#')
		{	parsedStrings=Parse::parse(currentLine);
			convertToDistanceBinUpperLimits(parsedStrings);
			parsedStrings.resize(0);
		}
	}

	inFile.close();
	inFile.clear();	
	upperDistanceLimit = distanceBinUpperLimits[distanceBinUpperLimits.size()-1];
	stepSize = distanceBinUpperLimits[distanceBinUpperLimits.size()-1] - distanceBinUpperLimits[distanceBinUpperLimits.size()-2];
	lowerDistanceLimit = distanceBinUpperLimits[0]-stepSize;
}

// where specific information about parsed data is intepreted
void pmf::convertToEnergyDataBase(const StrVec& _parsedStrings)
{	UInt type;
	sscanf(_parsedStrings[0].c_str(), "%u", &type);
	// simple data format check point
	if(type < 0)
	{	cout << "Error within the pmf.dat " << endl;
		exit(1);
	}

	DouVec tmpDouVec;
	tmpDouVec.resize(0);
	double tmpDouble;
	for(UInt i=2; i<_parsedStrings.size(); i++)
	{	sscanf(_parsedStrings[i].c_str(), "%lf", &tmpDouble);
		tmpDouVec.push_back(tmpDouble);
	}
	// when new type is encountered
	if(type != pairwiseEnergyData.size() || pairwiseEnergyData.size() == 0)
	{	vector<DouVec> tmpVecDouVec;
		tmpVecDouVec.resize(0);
		pairwiseEnergyData.push_back(tmpVecDouVec);
	}
	// add the array of energy values
	pairwiseEnergyData[type].push_back(tmpDouVec);
}

void pmf::convertToDistanceBinUpperLimits(const StrVec& _parsedStrings)
{	double tempDouble;
	for(UInt i=0; i<_parsedStrings.size(); i++)
	{	sscanf(_parsedStrings[i].c_str(), "%lf", &tempDouble);
		distanceBinUpperLimits.push_back(tempDouble);
	}
}

void pmf::readAtomTypeData()
{
	string evname = "PROTCADDIR";
	string path = getEnvironmentVariable(evname);

	path += "/data/";
	string file;
	file = path + "summa.frc";
	ifstream inFile;
	inFile.open(file.c_str());
	if (!inFile)
	{
		cout << "Error: unable to open summa.frc file: "
			<< file << endl;
		exit (1);
	}
	char ch;
	string stringChar;
	stringChar.resize(1);
	string stringBuf;
	vector< string > stringVec;	
	while (inFile.get(ch))
	{	
		if (ch != ' '  && ch != '\t' && ch != '\n')
		{
			stringChar[0] = ch;
			stringBuf.append(stringChar);
		}
		else if ( (ch == ' ' || ch == '\t') && (stringBuf != "") )
		{
			// new entry on same line. add stringBuf
			// to stringVec and recycle stringChar
			stringVec.push_back(stringBuf);
			stringBuf = "";
		}
		else if (ch == '\n')
		{
			// end of line in data file. process data
			// and recycle stringVec.
			if (stringBuf.size())
				stringVec.push_back(stringBuf);
			stringBuf = "";
			atomTypeNameStrings.push_back(stringVec);
			stringVec.resize(0);
		}
	}
	inFile.close();
	inFile.clear();
	atomTypeNameStringsBuilt = true;
/*
	cout << "summa.frc contents:" << endl;
	for (UInt i=0; i<atomTypeNameStrings.size();i++)
	{	for (UInt j=0; j<atomTypeNameStrings[i].size(); j++)
		{	cout << atomTypeNameStrings[i][j] << " ";
		}
		cout <<  endl;
	}
*/
}

int pmf::getIndexFromNameString(const string _name, const UInt _searchField) const
{
/*
	cout << "pmf::atomTypeNameStrings" << endl;
	for (UInt i=0; i<atomTypeNameStrings.size(); i++)
	{	for (UInt j=0; j<atomTypeNameStrings[i].size(); j++)
		{		cout << atomTypeNameStrings[i][j] << " ";
		}
		cout << "   size: " << atomTypeNameStrings[i].size();
		cout << endl;
	}
	cout << endl;
*/

	UInt field = _searchField + 1;
	
/*
	cout << "Name passed in : " << _name << endl;
	cout << "field :" << field << endl;
*/

	for (UInt i=0; i< atomTypeNameStrings.size(); i++)
	{
		if (field < atomTypeNameStrings[i].size())
		{
			if (_name == atomTypeNameStrings[i][field])
			{
//				cout << "Returning " << i << endl;
				return i;
			}
		}
	}
//	cout << "Returning  -1" << endl;
	return -1;
}
