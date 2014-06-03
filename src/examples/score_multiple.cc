#include "pdbReader.h"
#include "pdbWriter.h"
#include "atomIterator.h"
#include "ran.h"
#include "annealer.h"
#include "microEnvDB.h"
#include "microEnvironment.h"
#include <time.h>
#include <fstream>
#include <string>
#include "typedef.h"

bool foundInMap(int _x, int _y, vector<vector<int> > map);
void printPerAtomEnergy(microEnvironment* pBrian);
double correlate(vector<double>& x, vector<double>& y);
double rmsd_ca_all(molecule* pProt, molecule* pNative);
//int fastfind(int _x, int _y, vector<vector<int> > map, int _start, int _finish);

int main(int argc, char* argv[])
{

	if (argc != 5)
	{   cout << "Usage: score_multiple radius skip 1:calculate_all listname" << endl;
		exit (1);
	}

	//int wherefound = 1234;
	double tempDouble;
	int tempInt;
	sscanf(argv[1],"%lf", &tempDouble);
	double radius = tempDouble;
	sscanf(argv[2],"%i", &tempInt);
	int skip = tempInt;
	sscanf(argv[3],"%i", &tempInt);
	int flag = tempInt;
	string inputlistname;
	char* pArg4 = argv[4];

	for (;*pArg4 != '\0';pArg4++)
	{
		inputlistname += *pArg4;
	}

	cout << "Radius = " << radius << endl;
	cout << "Skip = " << skip << endl;
	cout << "Filename = " << inputlistname  << endl;


        molecule* pProt1 = pdbReader("1UBI.pdb");
        delete pProt1;
        microEnvironment* pBrian = new microEnvironment(radius,skip);
	pBrian->setDistanceDependenceOff();
	residue::printDataBaseData();

	time_t startTime,endTime,runTime;
	time(&startTime);

	string evname = "MISFOLD_MULTIPLE";
	string path_level1 = getEnvironmentVariable(evname);
	string listname = path_level1 + "/" + inputlistname;
	ifstream listfile;
	listfile.open(listname.c_str());
	if (!listfile)
	{	cout << "Error unable to read listfile: " << listname << endl;
		return 1;
	}
	string linebuffer;
	while (listfile >> linebuffer)
	{
		string path_level2 = path_level1 + "/" + linebuffer;
		string listname2 = path_level2 + "/list";
		ifstream listfile2;
		listfile2.open(listname2.c_str());
		if (!listfile2)
		{	cout << "Error unable to read listfile: " << listname2 << endl;
			return 1;
		}
		string linebuffer2;
		while (listfile2 >> linebuffer2)
		{	
			string path_level3 = path_level2 + "/" + linebuffer2;
			string listname3 = path_level3 + "/list";
			//cout << listname3 << endl;
			ifstream listfile3;
			listfile3.open(listname3.c_str());
			string linebuffer3;
			double native_microE_bu = 0.0;
			double native_microE_sb = 0.0;
			double native_microE_budd = 0.0;
			double native_microE_sbdd = 0.0;
			double native_ambervdw = 0.0;
			double native_amberelec = 0.0;
			double native_ambernonbond = 0.0;
			double native_pmfenergy = 0.0; 
			vector<vector<int> > native_contactmap;
			vector<double> microE_buvec;
			vector<double> microE_buddvec;
			vector<double> microE_sbvec;
			vector<double> microE_sbddvec;
			vector<double> ambervdwvec;
			vector<double> amberelecvec;
                        vector<double> ambernonbondvec;
			vector<double> pmfenergyvec;
			vector<double> percentnativevec;
// read rmsds file
			vector<double> rmsdsvec;
/*
			string filename = path_level3 + "/rmsds";
			ifstream inputfile;
			inputfile.open(filename.c_str());
			string linebuffer;
			double number;
			int counter = 1;
			while ( (inputfile >> linebuffer) )
			{ 	
				//cout << linebuffer << "  " << counter << endl;
				if (counter == 6)
				{	inputfile >> number;
				//	cout << number;
					counter = 0;
					rmsdsvec.push_back(number);
				}
				counter++;
			}
			inputfile.close();
*/
/*		
			for (UInt i=0; i<rmsdsvec.size(); i++)
			{ cout << rmsdsvec[i] << endl;
			}
*/
			//cout << "listfile3 = " << listfile3 << endl;
			int counter = 0;
			molecule* pNative = 0;
			while (listfile3 >> linebuffer3)
			{
				string filename = path_level3 + "/" + linebuffer3;
				//cout << "About to read " << filename << endl;
				molecule* pProt = pdbReader(filename);

				if (pProt == 0)
				{	cout << "Error reading " << filename << endl;
				return 1;
				}

				//pBrian->initialize();

				double microE_bu = 0.0;
				double microE_budd = 0.0;
				double microE_sb = 0.0;
				double microE_sbdd = 0.0;
				double ambervdw =0.0;
		 		double amberelec = 0.0;
				double ambernonbond = 0.0;
				double pmfenergy = 0.0;
				vector<vector<int> > contactmap;
				//double nativeContactsRecovered = 0;
				//double percentNativeContacts = 0;

                                pmf::setScaleFactor(0.0);
				microEnvironment::setScaleFactor(1.0);
                                microEnvDB::setUseSingleBodyEnergyOff();
                                pBrian->setDistanceDependenceOff();
                                microE_bu = pBrian->calculateEnergy(pProt);
				microE_buvec.push_back(microE_bu);
					
                                pmf::setScaleFactor(0.0);
				microEnvironment::setScaleFactor(1.0);
                                microEnvDB::setUseSingleBodyEnergyOn();
                                pBrian->setDistanceDependenceOff();
                                microE_sb = pBrian->calculateEnergy(pProt);
				microE_sbvec.push_back(microE_sb);
					
                                pmf::setScaleFactor(0.0);
				microEnvironment::setScaleFactor(1.0);
                                microEnvDB::setUseSingleBodyEnergyOff();
                                pBrian->setDistanceDependenceOn();
                                microE_budd = pBrian->calculateEnergy(pProt);
				microE_buddvec.push_back(microE_budd); 

                                pmf::setScaleFactor(0.0);
				microEnvironment::setScaleFactor(1.0);
                                microEnvDB::setUseSingleBodyEnergyOn();
                                pBrian->setDistanceDependenceOn();
                                microE_sbdd = pBrian->calculateEnergy(pProt);
				microE_sbddvec.push_back(microE_sbdd);
		
				if (counter == 0)
				{	rmsdsvec.push_back(0.0);
				}
				else
				{
					rmsdsvec.push_back(rmsd_ca_all(pProt,pNative));
				}
					
				//cout << "energies calculated " << endl;
				//residue::printDataBaseData();
/*
				// calculate contact map
				atomIterator ai1(static_cast<protein*>(pProt));
				atomIterator ai2(static_cast<protein*>(pProt));
				atom* pAtom1;
				atom* pAtom2;
				residue* pRes1;
				residue* pRes2;
				int atomnum1;
				int atomnum2;
				UInt atomIndex1;
				UInt atomIndex2;
				for (; !(ai1.last()); ai1++)
				{	pAtom1 = ai1.getAtomPointer();
					atomnum1 = pAtom1->getSerialNumber();
					ai2.initialize();
					for (; !(ai2.last()); ai2++)
					{	pAtom2 = ai2.getAtomPointer();
						if ((counter == 0 && pAtom2->distance(pAtom1) <= 5.0) ||
						    (counter != 0 && pAtom2->distance(pAtom2) <= 7.0))
						{	atomnum2 = pAtom2->getSerialNumber();
							pRes1 = ai1.getResiduePointer();
							pRes2 = ai2.getResiduePointer();
							atomIndex1 = ai1.getAtomIndex();
							atomIndex2 = ai2.getAtomIndex();
							bool closeNeighbors;
							if(pRes1 == pRes2)
								closeNeighbors = pRes1->isSeparatedByFewBonds(atomIndex1,atomIndex2);
							else
								closeNeighbors = pRes1->isSeparatedByFewBonds(pRes1,atomIndex1,pRes2,atomIndex2);
							if ((!closeNeighbors) && (!foundInMap(atomnum1,atomnum2,contactmap)))
							{
								vector<int> tempvec;
								tempvec.push_back(atomnum1);
								tempvec.push_back(atomnum2);
								contactmap.push_back(tempvec);
								//cout << "size of contactmap = " << contactmap.size() << endl;
							}	
						}
					}
				}
				cout << "Start Contact Map" << endl;
				for (UInt i=0; i<contactmap.size(); i++)
				{	cout << i << "  " << contactmap[i][0] << "  " << contactmap[i][1] << endl;
				}
				cout << "End Contact Map" << endl;
*/
				if (flag==1)
				{
				amberVDW::setScaleFactor(1.0);
				amberElec::setScaleFactor(0.0);
				pmf::setScaleFactor(0.0);
				rotamer::setScaleFactor(0.0);
				microEnvironment::setScaleFactor(0.0);

				ambervdw = static_cast<protein*>(pProt)->intraEnergy();
				ambervdwvec.push_back(ambervdw);

				amberVDW::setScaleFactor(0.0);
				amberElec::setScaleFactor(1.0);
				pmf::setScaleFactor(0.0);
				rotamer::setScaleFactor(0.0);
				microEnvironment::setScaleFactor(0.0);
				
				amberelec = static_cast<protein*>(pProt)->intraEnergy();
				amberelecvec.push_back(amberelec);

                                ambernonbond = amberelec + ambervdw;

				ambernonbondvec.push_back(ambernonbond);

				amberVDW::setScaleFactor(0.0);
				amberElec::setScaleFactor(0.0);
				pmf::setScaleFactor(1.0);
				rotamer::setScaleFactor(0.0);
				microEnvironment::setScaleFactor(0.0);

				pmfenergy = static_cast<protein*>(pProt)->intraEnergy();
				pmfenergyvec.push_back(pmfenergy);

				}

				if (counter==0)
				{
					native_microE_bu = microE_bu;
					native_microE_sb = microE_sb;
					native_microE_budd = microE_budd;
					native_microE_sbdd = microE_sbdd;
					native_ambervdw = ambervdw;
					native_amberelec = amberelec;
					native_ambernonbond = ambernonbond;
					native_pmfenergy = pmfenergy;
					pNative = pProt;
					pProt = 0;
				//	native_contactmap = contactmap;
				}
				cout << filename << "  ";
				if (counter == 0)
					cout << "NATIVE  ";

				cout << rmsdsvec[counter] << "  ";
				//cout << percentNativeContacts ;
				cout  << "  " << microE_bu;
				cout  << "  " << microE_sb;
				cout << "  " << microE_budd;
				cout << "  " << microE_sbdd;
				cout <<  "  " << ambervdw;
				cout  << "  " << amberelec;
				cout << "  " << ambernonbond;
				cout << "  "  << pmfenergy;
				cout  <<  endl;
				//printPerAtomEnergy(pBrian);
				counter++;
				delete pProt;
			}
			delete pNative;

// calculate correlation coefficients for all pairs
	
			double cof;
			cof = correlate(rmsdsvec, microE_buvec);
			cout << "correlation  rmsdvec to microE_bu = " << cof << endl;
			cof = correlate(rmsdsvec, microE_sbvec);
			cout << "correlation  rmsdvec to microE_sb = " << cof << endl;
			cof = correlate(rmsdsvec, microE_buddvec);
			cout << "correlation  rmsdvec to microE_budd = " << cof << endl;
			cof = correlate(rmsdsvec, microE_sbddvec);
			cout << "correlation  rmsdvec to microE_sbdd = " << cof << endl;
			cof = correlate(rmsdsvec, ambervdwvec);
			cout << "correlation  rmsdvec to ambervdw = " << cof << endl;
			cof = correlate(rmsdsvec, amberelecvec);
			cout << "correlation  rmsdvec to amberelec = " << cof << endl;
			cof = correlate(rmsdsvec, ambernonbondvec);
			cout << "correlation  rmsdvec to ambernonbond = " << cof << endl;
			cof = correlate(rmsdsvec, pmfenergyvec);
			cout << "correlation  rmsdvec to pmfenergy = " << cof << endl;
/*
			cof = correlate(percentnativevec, microE_buvec);
			cout << "correlation  percent_native to microE_bu = " << cof << endl;
			cof = correlate(percentnativevec, microE_sbvec);
			cout << "correlation  percent_native to microE_sb = " << cof << endl;
			cof = correlate(percentnativevec, microE_buddvec);
			cout << "correlation  percent_native to microE_budd = " << cof << endl;
			cof = correlate(percentnativevec, microE_sbddvec);
			cout << "correlation  percent_native to microE_sbdd = " << cof << endl;
			cof = correlate(percentnativevec, ambervdwvec);
			cout << "correlation  percent_native to ambervdw = " << cof << endl;
			cof = correlate(percentnativevec, amberelecvec);
			cout << "correlation  percent_native to amberelec = " << cof << endl;
			cof = correlate(percentnativevec, pmfenergyvec);
			cout << "correlation  percent_native to pmfenergy = " << cof << endl;
			*/

			listfile3.close();
		}
		listfile2.close(); 
	}
	listfile.close();
	time(&endTime);
	cout << "startTime = " << startTime << endl;
	cout << "endTime = " << endTime << endl;
	runTime = endTime - startTime;
	cout << "runTime = " << runTime << endl;

	return 0;
}

void printPerAtomEnergy(microEnvironment* pBrian)
{
		vector<double> perAtomEnergy = pBrian->getEnergyPerAtom();
		vector< vector<UInt> > perAtomEnv = pBrian->getEnvPerAtom();
		vector<UInt> atomTypes = pBrian->getEnvAssignedAtomTypes();
		for (UInt i=0; i< perAtomEnergy.size(); i++)
		{
			cout << i << "  - " << atomTypes[i] << " | ";
			for (UInt j=0; j<perAtomEnv[i].size(); j++)
			{	cout << perAtomEnv[i][j] << " ";
			}
			cout << perAtomEnergy[i] << endl;
		}
	return;
}

bool foundInMap(int _x, int _y, vector<vector<int> > map)
{
	UInt dim = map.size();
	for (UInt i=0; i<dim; i++)
	{	if ( (map[i][0] == _x && map[i][1] == _y) ||
		     (map[i][0] == _y && map[i][1] == _x) )
		{		return true;
		}
		if ( _x >= map[i][0] &&  _y > map[i][1] )
				return false;
	}
	return false;
}

double correlate(vector<double>& x, vector<double>& y)
{
	// check sizes
	UInt sizeofx = x.size();
	UInt sizeofy = y.size(); 

	if (sizeofx != sizeofy)
	{	cout << "Error! x and  y vectors cannot be correlated" << endl;
		return 0.0;
	}

	double xbar = 0.0;
	double ybar = 0.0;
	for (UInt i=0; i< sizeofx; i++)
	{	xbar += x[i];
		ybar += y[i];
	}

	xbar /= double(sizeofx);
	ybar /= double(sizeofx);

	vector<double> xminusxbar;
	vector<double> yminusybar;

	for (UInt i=0; i< sizeofx; i++)
	{	
		double tempdouble = x[i]- xbar;
		xminusxbar.push_back(tempdouble);
		tempdouble = y[i] -ybar;
		yminusybar.push_back(tempdouble);
	}

	double numerator = 0.0;
	for (UInt i=0; i< sizeofx; i++)
	{	numerator += xminusxbar[i]*yminusybar[i];
	}
	double sumxminusxbarsq = 0.0;
	double sumyminusybarsq = 0.0;
	for (UInt i=0; i< sizeofx; i++)
	{	 sumxminusxbarsq +=  xminusxbar[i]*  xminusxbar[i];
		 sumyminusybarsq +=  yminusybar[i]*  yminusybar[i];
	}
	double denominator = sqrt(sumxminusxbarsq * sumyminusybarsq);
	return numerator/denominator;
}

double rmsd_ca_all(molecule* pTheProtein2, molecule* pTheProtein1)
{
	if (pTheProtein1 == 0)
	{	return 999.999;
	}

	if (pTheProtein2 == 0)
	{	return 999.999;
	}

	double tempWeight = 1.0;
	vector<dblVec> coord1;
	vector<dblVec> coord2;
	vector<double> weights;

	atomIterator theIter1(static_cast<protein*>(pTheProtein1));
	atomIterator theIter2(static_cast<protein*>(pTheProtein2));

	atom* pAtom;
	dblVec pDV;

	for (;!(theIter1.last());theIter1++)
	{
		pAtom = theIter1.getAtomPointer(); 
		if ( pAtom->getName() == "CA" )
		{
			coord1.push_back(pAtom->getCoords());
			weights.push_back(tempWeight);
		}
	}

	for (;!(theIter2.last());theIter2++)
	{
		pAtom = theIter2.getAtomPointer(); 
		if ( pAtom->getName() == "CA" )
		{
			coord2.push_back(pAtom->getCoords());
		}
	}

	// Now unroll the coordinates into a single long vector
	// of size 3*Natm1
	int numatm1 = coord1.size();
	int numatm2 = coord2.size();
	int list1[numatm1];
	int list2[numatm1];
	int nat = numatm1;
	double newCoord1[numatm1*3];
	double newCoord2[numatm1*3];
	double newCoord3[numatm1*3];
	double newWeights[numatm1];
	double rotmat[9];
	double centroid1[3];
	double centroid2[3];
	double trnvec[3];
	double rmsd = 0;
	int ierr = 0;
	for (UInt i=0; i<numatm1; i++)
	{	for (UInt c1=0; c1<3;c1++)
		{	newCoord1[ (i*3) + c1] = coord1[i][c1];
			newCoord2[ (i*3) + c1] = coord2[i][c1];
		}
		list1[i] = i+1;
		list2[i] = i+1;
		newWeights[i] = weights[i];
	}

	bestfit_(newCoord1, &numatm1, newCoord2,
		&numatm2, &nat, newCoord3, list1, list2,
		&rmsd, &ierr, rotmat, centroid1, centroid2);

	return rmsd;
}
