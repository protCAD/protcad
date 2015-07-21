#include "microEnvDB.h"
//#define MICROENVDB_DEBUG
//#define MICROENVDB_DEBUG_HIGH

vector< vector < vector < UInt > > > microEnvDB::itsMicroEnvMap;
bool microEnvDB::useSingleBodyEnergy = true;

microEnvDB::microEnvDB()
{	// default constructor
}

microEnvDB::microEnvDB(double _rad, UInt _skip)
{
	// standard constructor for this class - first argument
	// is the counting radius and the second argument is
	// the residue skipping number

	if (_rad <= 2.0)
	{	cout << "Error in building microEnvDB class." << endl;
		cout << "Radius of " << _rad << " must be in error!" << endl;
		exit(1);
	}
	itsCriticalRadius = _rad;
	itsResidueSkippingNumber = _skip;
	readClusterData();
	readEnergyData();
	readSingleBodyEnergyData();

#ifdef MICROENVDB_DEBUG
	cout << "microEnvDB database is built" << endl;
	printMicroEnvMap();
	printEnergyData();
#endif

}

microEnvDB::~microEnvDB()
{
}

double microEnvDB::getEnergy(const UInt _type1, const UIntVec& _cluster) const
{
	double jason = findEnergyInDB(_type1, _cluster);
	return jason;
}

double microEnvDB::findEnergyInDB(const UInt _type1, const UIntVec& _cluster) const
{
		
	UInt coord = _cluster[0] + _cluster[1] + _cluster[2] + _cluster[3];
	double microEnvEnergy = 0;
	if (coord < itsMicroEnvMap.size())
	{
		for (UInt i=0; i< itsMicroEnvMap[coord].size(); i++)
		{	if (_cluster == itsMicroEnvMap[coord][i])
			{
				microEnvEnergy = energyData[_type1][coord][i];
				if (microEnvEnergy == 33.0 && useSingleBodyEnergy == true)
				{	
					microEnvEnergy = itsSingleBodyEnergy[coord][_type1];
				}
				return microEnvEnergy;
			}
		}
	}
	// default value if the coordination is larger than what we 
	// have explicitly enumerated
	microEnvEnergy = energyData[_type1][itsMicroEnvMap.size()][0];
	return microEnvEnergy;
}

void microEnvDB::readClusterData()
{	if (itsMicroEnvMap.size() == 0)
	{

		string evname = "PROTCADDIR";
		string path = getEnvironmentVariable(evname);

		path += "/data/microEnv/";
		string filebase="combo";
		string fileend =".dat";

		vector< vector < UInt > > coordVector;
		for (UInt fmodint = 0;;fmodint++)
		{
			coordVector.resize(0);
			char fmodchar[3];
			sprintf(fmodchar,"%i",fmodint);
			string fmodstring = fmodchar;
			string file = path + filebase + fmodstring + fileend;
			ifstream inFile;
			inFile.open(file.c_str());
			if (!inFile)
			{      
				if (fmodint <7)
				{	cout << "Error: unable to open microEnvDB file: "
					<< file << endl;
   	          			exit (1);
				}
				return;	
		    	}
			string linebuffer;
			int number;
			// The first line contains 3 fields which can be discarded
			for (UInt i=0;i<3;i++)
			{
				inFile >> linebuffer;	
			}
	
			int counter=0;
			vector< UInt > tempVector;
			tempVector.resize(0);
			while (inFile >> number)
			{
				if (counter == 3)
				{	tempVector.push_back(number);
					coordVector.push_back(tempVector);
					tempVector.resize(0);
					counter = 0;
					continue;
				}
				tempVector.push_back(number);
				counter++;
			}
			inFile.close();
			inFile.clear();
			itsMicroEnvMap.push_back(coordVector);
		}
	}
}
	
void microEnvDB::printMicroEnvMap()
{	UInt maxCoord = itsMicroEnvMap.size();
	for (UInt i=0;i<maxCoord; i++)
	{	UInt maxEnv = itsMicroEnvMap[i].size();
		for (UInt j=0; j<maxEnv; j++)
		{	for (UInt k=0; k<4; k++)
			{
				cout << itsMicroEnvMap[i][j][k] << "  ";
			}
			cout << endl;
		}
		cout << endl << endl;
	}
	cout << endl << "-----------" << endl;
}

void microEnvDB::printEnergyData()
{	UInt numTypes = energyData.size();
	cout << "numTypes = " << numTypes << endl;
	for (UInt i=0; i<numTypes; i++)
	{	UInt maxCoord = itsMicroEnvMap.size();
		for (UInt j=0; j<maxCoord; j++)
		{	UInt maxEnv = energyData[i][j].size();
			for (UInt k=0; k<maxEnv; k++)
			{	cout << "Atom Type " << i << " : Coord " << j << " | ";
				for (UInt l=0; l<4; l++)
				{	cout << itsMicroEnvMap[j][k][l] << " ";
				}
				cout << energyData[i][j][k] << endl;
			}
			cout << endl;
		}
		cout << "Atom Type " << i << " : Coord  >" << itsMicroEnvMap.size()-1;
		cout << " | x x x x " << energyData[i][itsMicroEnvMap.size()][0];
		cout << endl << endl;
	}
	cout << endl << "-------------" << endl;
}

void microEnvDB::readEnergyData()
{
	string evname = "PROTCADDIR";
	string path = getEnvironmentVariable(evname);

	path += "/data/microEnv/";

	string filebase = "11.1_std";

	char fmoddis[5];
	sprintf(fmoddis,"%3.1f",itsCriticalRadius);
	string fmoddisstring = fmoddis;

	char fmodskip[5];
	sprintf(fmodskip,"%i",itsResidueSkippingNumber);
	string fmodskipstring = fmodskip;

	string file = path + fmoddisstring + "_" + fmodskipstring + "_" + filebase;

	//cout << "Reading energy file: " << file << endl;
	ifstream inFile;
	inFile.open(file.c_str());
	if (!inFile)
      	{      
		cout << "Error: unable to open microEnvDB file: "
			<< file << endl;
		exit (1);
      	}

	double number;
	UInt currentAtomType = 0;
	vector <double> evector;
	vector < vector <double> > coordvector;
	coordvector.resize(0);

#ifdef MICROENVDB_DEBUG
	cout << "itsMicroEnvMap.size() = " << itsMicroEnvMap.size() << endl;
	for (UInt i=0; i<itsMicroEnvMap.size(); i++)
		cout << itsMicroEnvMap[i].size() << endl;
	cout << endl;
#endif

	UInt currentCoord = 0;
	while (inFile >> number)
	{	evector.push_back(number);
		if (evector.size() == itsMicroEnvMap[currentCoord].size())
		{	coordvector.push_back(evector);
			evector.resize(0);
			currentCoord++;
		}
		if (coordvector.size() == itsMicroEnvMap.size())
		{	// coordination number greater than what we've
			// explicitly enumerated
			inFile >> number;
			evector.resize(0);
			evector.push_back(number);
			coordvector.push_back(evector);
			energyData.push_back(coordvector);
			currentCoord = 0;
			currentAtomType++;
			evector.resize(0);
			coordvector.resize(0);
		}
	}
	inFile.close();
	inFile.clear();
}

void microEnvDB::readSingleBodyEnergyData()
{
	string evname = "PROTCADDIR";
	string path = getEnvironmentVariable(evname);

	path += "/data/microEnv/";

	string filebase = "11.1_std_single_body";

	char fmoddis[5];
	sprintf(fmoddis,"%3.1f",itsCriticalRadius);
	string fmoddisstring = fmoddis;

	char fmodskip[5];
	sprintf(fmodskip,"%i",itsResidueSkippingNumber);
	string fmodskipstring = fmodskip;

	string file = path + fmoddisstring + "_" + fmodskipstring + "_" + filebase;

	//cout << "Reading energy file: " << file << endl;
	ifstream inFile;
	inFile.open(file.c_str());
	if (!inFile)
      	{      
		cout << "Error: unable to open microEnv_singlebody file: "
			<< file << endl;
		exit (1);
      	}

	double number;
	vector <double> evector;
	evector.resize(0);

	while (inFile >> number)
	{	evector.push_back(number);
		if (evector.size() == 16)
		{
			itsSingleBodyEnergy.push_back(evector);
			evector.resize(0);
			inFile >> number;
			evector.push_back(number);
		}
	}
	inFile.close();
	inFile.clear();
/*
	cout << "SingleBodyEnergy" << endl;
	for (UInt i=0; i<itsSingleBodyEnergy.size(); i++)
	{	for (UInt j=0; j<itsSingleBodyEnergy[i].size(); j++)
		{	cout << itsSingleBodyEnergy[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl;
*/

}
