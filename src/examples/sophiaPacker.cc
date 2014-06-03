#include <iostream>
#include <string>
#include <vector>
#include "ensemble.h"
#include "PDBInterface.h"
#include "ran.h"
#include "generalio.h"

static double PI=3.141596;

// MODIFY BACKBONE FUNCTIONS
void createBackbone(protein* _prot, double _rad, double _phase, double _square, double _offset, double _coil);
void undoBackbone(protein* _prot, double _rad, double _phase, double _square, double _offset, double _coil);


// ENERGY FUNCTIONS
double getProteinLigandEnergy(protein* _prot, vector <dblVec> _ligCoords, double _r0, double _eps, vector <UIntVec> _exclude);
	// exclude vector has four elements ->  protein chain, protein residue, protein atom, ligand atom
double getHisHemeEnergy(protein* _prot, UInt _Hc, UInt _Hr, vector <dblVec> _lig, UInt _metal, UInt _por1, UInt por2, bool _printE);
double getThrBackboneCarbonylEnergy(protein* _prot, UInt _Tc, UInt _Tr);
double getThrHisHBondEnergy(protein* _prot, UInt _Hc, UInt _Hr, UInt _Tc, UInt _Tr);

// OTHER FUNCTIONS
vector <dblVec> getLigandCoords(string _fileName);
void optimizeRotamerGeometry(protein* _prot, UInt _Hc, UInt _Hr, UInt _Tc, UInt _Tr, vector <dblVec> _ligCoords, UInt _metalIndex, UInt _por1index, UInt _por2index);
void writeStateToFile(protein* _prot, double _rad, double _phase, double _square, double _offset, double _coil,
	double _temp, double _energy, string _fileName);

int main(int argc, char* argv[])
{
	enum aminoAcid {A, R, N, D, C, Q, E, G, H, I, L, K, M, F, P, S, T, W, Y, V};
	if (argc < 2)
	{
		cout << "portpyPacker inputfile.inp" << endl;
		exit(1);
	}
	string inputFileName = argv[1];
	ifstream inFile;
	inFile.open(inputFileName.c_str());
	if (!inFile)
	{
		cout << "Unable to find or open file" << endl;
		exit(1);
	}

	string currentLine;
	vector <string> parsedStrings;
	parsedStrings.resize(0);
	getline(inFile, currentLine, '\n'); parsedStrings = Parse::parse(currentLine);
	string inputFile = parsedStrings[0];
	string outputFile = parsedStrings[1];

	getline(inFile, currentLine, '\n'); parsedStrings = Parse::parse(currentLine);
	string ligFile = parsedStrings[0];

	vector <dblVec> ligCoords = getLigandCoords(ligFile);

	// read in single helix
	PDBInterface* thePDB = new PDBInterface(inputFile);
	ensemble* theEnsemble = thePDB->getEnsemblePointer();
	molecule* theMol = theEnsemble->getMoleculePointer(0);
	protein* prot = static_cast<protein*>(theMol);

	// setup energy parameters
	residue::setCutoffDistance(10.0);
	pmf::setScaleFactor(0.0);
	rotamer::setScaleFactor(0.0);
	microEnvironment::setScaleFactor(0.0);
	amberVDW::setScaleFactor(1.0);
	amberVDW::setRadiusScaleFactor(1.0);
	amberVDW::setLinearRepulsionDampeningOff();
	amberElec::setScaleFactor(0.0);
	solvation::setItsScaleFactor(0.0);


	prot->symmetryLinkChainAtoB(1,0); // C2 pair
	prot->symmetryLinkChainAtoB(2,0); // C2 pair
	prot->symmetryLinkChainAtoB(3,0); // C2 pair
	// get random seed
	int randomSeed;
	getline(inFile, currentLine, '\n');
	sscanf(currentLine.c_str(), "%d", &randomSeed);
	cout << "Randomseed = " << randomSeed << endl;
	ran ranNumber;
	ranNumber.setSeed((UInt)randomSeed);


	// ************* READING HISTIDINE PARAMETERS ************************
	// get Histidine position and rotamer
	UInt hisres, hisrot;
	getline(inFile, currentLine, '\n'); parsedStrings = Parse::parse(currentLine);
	sscanf(parsedStrings[0].c_str(), "%u", &hisres);
	sscanf(parsedStrings[1].c_str(), "%u", &hisrot);
	prot->activateForRepacking(0,hisres);
	prot->mutate(0, hisres, H);
	prot->setRotamer(0, hisres, 0, hisrot);


	// align helices so that beta carbon of His sits on z axis and is at origin
	dblVec hisCB = prot->getCoords(0, hisres, "CB");
	dblVec tempX(3);
	tempX[0]=0.0;
	tempX[1]=-1.0;
	tempX[2]=0.0;
	double ztrans = -1* hisCB[2];
	hisCB[2]=0.0;

	double magCB = sqrt(CMath::dotProduct(hisCB,hisCB));
	double magX = sqrt(CMath::dotProduct(tempX,tempX));

	double hisCBAngle = acos(CMath::dotProduct(hisCB,tempX)/(magX*magCB));
	hisCBAngle = (180/PI)*hisCBAngle;
	cout << "his angle " << hisCBAngle << endl;
	prot->rotate(Z_axis, hisCBAngle);
	prot->translate(0,0, ztrans);

	// read in his parameters

	double offsetminhis, offsetmaxhis, phaseminhis, phasemaxhis, radminhis, radmaxhis, coilminhis, coilmaxhis, squareminhis, squaremaxhis;


	getline(inFile, currentLine, '\n'); parsedStrings = Parse::parse(currentLine);
	sscanf(parsedStrings[0].c_str(), "%lf", &offsetminhis);
	sscanf(parsedStrings[1].c_str(), "%lf", &offsetmaxhis);

	getline(inFile, currentLine, '\n'); parsedStrings = Parse::parse(currentLine);
	sscanf(parsedStrings[0].c_str(), "%lf", &phaseminhis);
	sscanf(parsedStrings[1].c_str(), "%lf", &phasemaxhis);

	getline(inFile, currentLine, '\n'); parsedStrings = Parse::parse(currentLine);
	sscanf(parsedStrings[0].c_str(), "%lf", &radminhis);
	sscanf(parsedStrings[1].c_str(), "%lf", &radmaxhis);

	getline(inFile, currentLine, '\n'); parsedStrings = Parse::parse(currentLine);
	sscanf(parsedStrings[0].c_str(), "%lf", &coilminhis);
	sscanf(parsedStrings[1].c_str(), "%lf", &coilmaxhis);

	getline(inFile, currentLine, '\n'); parsedStrings = Parse::parse(currentLine);
	sscanf(parsedStrings[0].c_str(), "%lf", &squareminhis);
	sscanf(parsedStrings[1].c_str(), "%lf", &squaremaxhis);


	double radhis = radminhis + ranNumber.getNext()*(radmaxhis-radminhis);
	double phasehis = phaseminhis + ranNumber.getNext()*(phasemaxhis-phaseminhis);
	double squarehis = squareminhis + ranNumber.getNext()*(squaremaxhis-squareminhis);
	double offsethis = offsetminhis + ranNumber.getNext()*(offsetmaxhis-offsetminhis);
	double coilhis = coilminhis + ranNumber.getNext()*(coilmaxhis-coilminhis);

	double radbesthis = radhis;
	double phasebesthis = phasehis;
	double squarebesthis = squarehis;
	double offsetbesthis = offsethis;
	double coilbesthis = coilhis;

	double intraWeight, hemeCoordWeight, proteinLigWeight;
	// read in energy scale factors
	getline(inFile, currentLine, '\n'); parsedStrings = Parse::parse(currentLine);
	sscanf(parsedStrings[0].c_str(), "%lf", &intraWeight);
	sscanf(parsedStrings[1].c_str(), "%lf", &hemeCoordWeight);
	sscanf(parsedStrings[2].c_str(), "%lf", &proteinLigWeight);

	// thr position and rotamer
	UInt thrres, thrrot;
	getline(inFile, currentLine, '\n'); parsedStrings = Parse::parse(currentLine);
	sscanf(parsedStrings[0].c_str(), "%u", &thrres);
	sscanf(parsedStrings[1].c_str(), "%u", &thrrot);
	prot->activateForRepacking(0,thrres);
	prot->mutate(0, thrres, T);
	prot->setRotamer(0, thrres, 0, thrrot);


	// read in thr-his hbond energy weight factor and thr - backbone energy weight factor

	double thrHisHBondWeight, thrBackboneWeight;
	getline(inFile, currentLine, '\n'); parsedStrings = Parse::parse(currentLine);
	sscanf(parsedStrings[0].c_str(), "%lf", &thrHisHBondWeight);
	sscanf(parsedStrings[1].c_str(), "%lf", &thrBackboneWeight);

	// read in position for heme bounding glycines
	UInt glyres;
	getline(inFile, currentLine, '\n'); parsedStrings = Parse::parse(currentLine);
	sscanf(parsedStrings[0].c_str(), "%u", &glyres);
	prot->activateForRepacking(0,glyres);
	prot->mutate(0, glyres, G);

	// READ IN TEMPERATURE and COOLING SCHEDULE

	double startTemp1, endTemp1, steps1, startTemp2, endTemp2, steps2;
	// read in simulated annealing schedule
	getline(inFile, currentLine, '\n'); parsedStrings = Parse::parse(currentLine);
	sscanf(parsedStrings[0].c_str(), "%lf", &startTemp1);
	sscanf(parsedStrings[1].c_str(), "%lf", &endTemp1);
	sscanf(parsedStrings[2].c_str(), "%lf", &steps1);
	getline(inFile, currentLine, '\n'); parsedStrings = Parse::parse(currentLine);
	sscanf(parsedStrings[0].c_str(), "%lf", &startTemp2);
	sscanf(parsedStrings[1].c_str(), "%lf", &endTemp2);
	sscanf(parsedStrings[2].c_str(), "%lf", &steps2);

	UInt metpos, porz, pory;
	getline(inFile, currentLine, '\n'); parsedStrings = Parse::parse(currentLine);
	sscanf(parsedStrings[0].c_str(), "%u", &metpos);
	sscanf(parsedStrings[1].c_str(), "%u", &porz);
	sscanf(parsedStrings[2].c_str(), "%u", &pory);


	cout << "COORDS - metal " << ligCoords[metpos][0] << " " << ligCoords[metpos][1] << " " << ligCoords[metpos][2] << endl;
	cout << "COORDS -  porz " << ligCoords[porz][0] << " " << ligCoords[porz][1] << " " << ligCoords[porz][2] << endl;
	cout << "COORDS -  pory " << ligCoords[pory][0] << " " << ligCoords[pory][1] << " " << ligCoords[pory][2] << endl;

	prot->silenceMessages();
	createBackbone(prot, radhis, phasehis, squarehis, offsethis, coilhis);
	//optimizeRotamerGeometry(prot, 0, hisres, 0, thrres, ligCoords, metpos, porz, pory);
	pdbWriter(prot, "start.pdb");
	undoBackbone(prot, radhis, phasehis, squarehis, offsethis, coilhis);
	double bestEnergy = intraWeight*prot->intraEnergy();
	double oldEnergy = bestEnergy;

	double e = 2.71828182846;
	for (UInt cycle = 1; cycle <=2; cycle ++)
	{
		double starttemp, endtemp;
		double steps;
		if (cycle == 1)
		{
			cout << "Starting Cycle 1" << endl;
			starttemp = startTemp1; endtemp = endTemp1; steps = steps1;
		}
		if (cycle == 2)
		{
			cout << "Starting Cycle 2" << endl;
			starttemp = startTemp2; endtemp = endTemp2; steps = steps2;
		}

		double delta = (endtemp - starttemp) / steps;
		cout << "start " << starttemp << " end " << endtemp << " steps " << steps << " delta " << delta << endl;
		UInt step = 0;
		for (double temp = starttemp; temp >= endtemp; temp = temp + delta)
		{
			prot->mutate(0, hisres, H);	// prevent degradation of sidechain conformation from extensive rotamer cycling
			prot->setRotamer(0, hisres, 0, hisrot);
			prot->mutate(0, thrres, T);
			prot->setRotamer(0, thrres, 0, thrrot);

			step ++;
			UInt variable = UInt(ranNumber.getNext()*5);
			//modify his backbone
			if (variable == 0 || variable >= 4 )
			{
				radhis = radminhis + (radmaxhis-radminhis)*ranNumber.getNext();
				cout << cycle << " " << step << " temp: " << temp << " changing hiscoil radius to " << radhis << endl;
			}
			if (variable == 1 || variable >= 4 )
			{
				phasehis = phaseminhis + (phasemaxhis - phaseminhis)*ranNumber.getNext();
				cout << cycle << " " << step << " temp: " << temp << " changing hiscoil phase to " << phasehis << endl;
			}
			if (variable == 2 || variable >= 4 )
			{
				squarehis = squareminhis + (squaremaxhis-squareminhis)*ranNumber.getNext();
				cout << cycle << " " << step << " temp: " << temp << " changing hiscoil square to " << squarehis << endl;
			}
			if ( variable == 3 || variable >= 4 )
			{
				offsethis = offsetminhis + (offsetmaxhis - offsetminhis)*ranNumber.getNext();
				cout << cycle << " " << step << " temp: " << temp << " changing hiscoil offset to " << offsethis << endl;
			}
			//if ( variable == 4 || variable == 5 )
		//	{
		//		coilhis = coilminhis + (coilmaxhis - coilminhis)*ranNumber.getNext();
		//		cout << cycle << " " << step << " temp: " << temp << " changing hiscoil coil to " << coilhis << endl;
		//	}

			createBackbone(prot, radhis, phasehis, squarehis, offsethis, coilhis);

		    // get his to metal distance
	   	 	dblVec HisNe2 = prot->getCoords(0, hisres, "NE2");
	    	dblVec diff = HisNe2 - ligCoords[metpos];
	    	double hismetalDist = sqrt(CMath::dotProduct(diff,diff));

			double energy = 1E10;
			if (hismetalDist < 6.0)
			{
				double intraenergy = intraWeight*prot->intraEnergy();
				// build exclude list for prot lig energy;
				vector <UIntVec> exclude;
				UIntVec tempv(2);
				for	(UInt cha = 0; cha < prot->getNumChains(); cha ++)
				{
					for (UInt res = 0; res < prot->getNumResidues(cha); res ++)
					{
						if (prot->getTypeFromResNum(cha,res) != A)
						{
							tempv[0]=cha, tempv[1]=res;
							exclude.push_back(tempv);
						}
					}
				}
				double protLigEnergy = proteinLigWeight*getProteinLigandEnergy(prot, ligCoords, 4.0, 1.0, exclude);
				double hemeCoordEnergy = 0.0;
				double thrHisHBondEnergy = 0.0;
				double thrBackboneEnergy = 0.0;
				energy = protLigEnergy + intraenergy;
				if (energy < 0) // dont waste time optimizing rotamers if energy is high!
				{
					optimizeRotamerGeometry(prot, 0, hisres, 0, thrres, ligCoords, metpos, porz, pory);
					hemeCoordEnergy = hemeCoordWeight*getHisHemeEnergy(prot, 0, hisres, ligCoords, metpos, porz, pory, true);
					thrHisHBondEnergy = thrHisHBondWeight*getThrHisHBondEnergy(prot, 0, hisres, 3, thrres); // HACK!!!!
					thrBackboneEnergy = thrBackboneWeight*getThrBackboneCarbonylEnergy(prot, 0, thrres);
				}
				energy = intraenergy + protLigEnergy + hemeCoordEnergy + thrHisHBondEnergy + thrBackboneEnergy;
				cout << "ENERGIES:  intra = " << intraenergy << " ligand sterics = " << protLigEnergy << " heme coord = " << hemeCoordEnergy
					<< " thr-his hbond = " << thrHisHBondEnergy << " thr-backbone = " << thrBackboneEnergy << endl;
			}
			else cout << "TOO LONG  " << hismetalDist << endl;

			if (energy < oldEnergy)
			{
				oldEnergy = energy;
				radbesthis = radhis; phasebesthis = phasehis; squarebesthis = squarehis; offsetbesthis = offsethis; coilbesthis = coilhis;
				cout << "ACCEPTED GOOD" << endl;
				pdbWriter(prot, "lowest.pdb");
				if (energy < bestEnergy)
		        {
					cout << "LOWEST STATE YET " << energy << " kcal/mol" << endl;
					pdbWriter(prot, "best.pdb");
					bestEnergy = energy;
					writeStateToFile(prot, radbesthis, phasebesthis, squarebesthis, offsetbesthis, coilbesthis, temp, bestEnergy, "best.out");
				}
				undoBackbone(prot, radhis, phasehis, squarehis, offsethis, coilhis);
			}
			else
			{
				double deltaEnergy = energy - oldEnergy;
				double probability = pow(e, (-1*deltaEnergy / (0.00258*temp)));
				if (ranNumber.getNext() < probability)
			{
					oldEnergy = energy;
					radbesthis = radhis; phasebesthis = phasehis; squarebesthis = squarehis; offsetbesthis = offsethis; coilbesthis = coilhis;
					cout << "ACCEPTED BOLTZMANN" << endl;
					pdbWriter(prot, "lowestboltz.pdb");
					undoBackbone(prot, radhis, phasehis, squarehis, offsethis, coilhis);
				}
				else
				{
					undoBackbone(prot, radhis, phasehis, squarehis, offsethis, coilhis);
					radhis = radbesthis; phasehis = phasebesthis; squarehis = squarebesthis; offsethis = offsetbesthis; coilhis = coilbesthis;
					cout << "REJECTED (" << energy << " kcal/mol)" << endl;
				}
			}
			writeStateToFile(prot, radbesthis, phasebesthis, squarebesthis, offsetbesthis, coilbesthis, temp, bestEnergy, "current.out");
		}
	}
	createBackbone(prot, radhis, phasehis, squarehis, offsethis, coilhis);
	optimizeRotamerGeometry(prot, 0, hisres, 0, thrres, ligCoords, metpos, porz, pory);
	pdbWriter(prot, outputFile);


	return 0;
}

void optimizeRotamerGeometry(protein* _prot, UInt _Hc, UInt _Hr, UInt _Tc, UInt _Tr, vector <dblVec> _ligCoords, UInt _metalIndex, UInt _por1index, UInt _por2index)
{
	UIntVec rotamer;

	rotamer = _prot->getCurrentRotamer(_Hc, _Hr);
	_prot->setRotamer(_Hc, _Hr, 0, rotamer[0]);

	double hischi1standard = _prot->getChi(_Hc, _Hr, 0, 0);
	double hischi2standard = _prot->getChi(_Hc, _Hr, 0, 1);

	rotamer = _prot->getCurrentRotamer(_Tc, _Tr);
	_prot->setRotamer(_Tc, _Tr, 0, rotamer[0]);
	double thrchi1standard = _prot->getChi(_Tc, _Tr, 0,0);


	// STAGE 1:  COARSE OPTIMIZATION
	double range = 35.0;
	double step = 5.0;

	double tempangle;
	double bestEnergy = 1E10;
	double besthischi1 = 0.0;
	double besthischi2 = 0.0;
	double bestthrchi1 = 0.0;
	for (double hischi1 = -1*range; hischi1 <= range; hischi1 = hischi1 + step)
	{
		tempangle = hischi1standard + hischi1;
		_prot->setChi(_Hc, _Hr, 0, 0, tempangle);
		for (double hischi2 = -1*range; hischi2 <= range; hischi2 = hischi2 + step)
		{
			tempangle = hischi2standard + hischi2;
			_prot->setChi(_Hc, _Hr, 0, 1, tempangle);
			for (double thrchi1 = -1*range; thrchi1 <= range; thrchi1 = thrchi1 + step)
			{
				tempangle = thrchi1standard + thrchi1;
				_prot->setChi(_Tc, _Tr, 0, 0, tempangle);
				double energy = 0.0;
				energy += getHisHemeEnergy(_prot, _Hc, _Hr, _ligCoords, _metalIndex, _por1index, _por2index, false);
				energy += getThrHisHBondEnergy(_prot, _Hc, _Hr, 3, _Tr); // HACK ... should have symmetry check!
				energy += getThrBackboneCarbonylEnergy(_prot, _Tc, _Tr);
				energy += _prot->getPositionEnergy(_Tc, _Tr);
				energy += _prot->getPositionEnergy(_Hc, _Hr);
				if (energy < bestEnergy)
				{
					besthischi1 = hischi1standard + hischi1;
					besthischi2 = hischi2standard + hischi2;
					bestthrchi1 = thrchi1standard + thrchi1;
					bestEnergy = energy;
				}
			}
		}
	}

	hischi1standard = besthischi1;
	hischi2standard = besthischi2;
	thrchi1standard = bestthrchi1;

	_prot->setChi(_Hc, _Hr, 0, 0, besthischi1);
	_prot->setChi(_Hc, _Hr, 0, 1, besthischi2);
	_prot->setChi(_Tc, _Tr, 0, 0, bestthrchi1);
	return;
}


double getHisHemeEnergy(protein* _prot, UInt _Hc, UInt _Hr, vector <dblVec> _lig, UInt _metal, UInt _por1, UInt _por2, bool _printE)
{
	// get his to metal distance
	dblVec HisNe2 = _prot->getCoords(_Hc, _Hr, "NE2");
	dblVec diff = HisNe2 - _lig[_metal];
	double HisPorDist = sqrt(CMath::dotProduct(diff,diff));

	// get angle his to heme atom on z axis
	dblVec HisCd2 = _prot->getCoords(_Hc, _Hr, "CD2");
	dblVec HisCe1 = _prot->getCoords(_Hc, _Hr, "CE1");
	dblVec pseudoAtom = 0.5 * (HisCe1 + HisCd2);
	dblVec BA = HisNe2 - pseudoAtom;
	dblVec BC = _lig[_metal] - _lig[_por1];
	double magBC = sqrt( CMath::dotProduct(BC,BC) );
	double magBA = sqrt( CMath::dotProduct(BA,BA) );
	double HisPorAngle1 = acos( (CMath::dotProduct(BA,BC)) / (magBA * magBC) );
	// get angle his to heme atom on x axis
	BC = _lig[_metal] - _lig[_por2];
	magBC = sqrt( CMath::dotProduct(BC,BC) );
	double HisPorAngle2 = acos( (CMath::dotProduct(BA,BC)) / (magBA * magBC) );

	// get His Por Intercept
	BA = pseudoAtom - _lig[_metal] ;
	BC = pseudoAtom - HisNe2;
	magBC = sqrt( CMath::dotProduct(BC,BC) );
	magBA = sqrt( CMath::dotProduct(BA,BA) );
	double HisPorInterceptAngle = acos( (CMath::dotProduct(BA,BC)) / (magBA * magBC) );

	double eps = 100.0;
	double r0 = 2.0;

	double distratio  = r0/HisPorDist;
	double angleFactor1 = 1.0;
	angleFactor1 *= pow(cos(HisPorAngle1-PI/2),2);
	angleFactor1 *= pow(cos(HisPorAngle2-PI/2),2);
	angleFactor1 *= pow(cos(HisPorInterceptAngle),2);
	double energyhis = eps*(5*pow(distratio,12) - 6 * pow(distratio,10));
	energyhis *= angleFactor1;

	// calculate his-heme sterics
	double hisPorSterics = 0.0;
	// skip metal
	for (UInt i = 0; i < _prot->getNumAtoms(_Hc, _Hr); i ++)
	{
		for (UInt j = 0; j < _lig.size(); j ++)
		{
			if (j != _metal)
			{
				dblVec atomCoords = _prot->getCoords(_Hc, _Hr, i);
				dblVec diff = atomCoords - _lig[j];
				double distance = sqrt(CMath::dotProduct(diff,diff));
				double distratio = 2.0/distance;
				hisPorSterics += 5.0*(pow(distratio,10)); // soften repulsion to 10
			}
		}
	}

	if (_printE)
	{
		cout << endl;
		cout << endl;
		cout << "*********************HIS HEME INTERACTION**********************" << endl;
		cout << "His " << _Hc << " " << _Hr << endl;
		cout << "	his - por distance " << HisPorDist << " hispor angle " << 180*HisPorAngle1/PI<< " " << 180*HisPorAngle2/PI << endl;
		cout << "	his - por - his angle " << 180*HisPorInterceptAngle/PI << endl;
		cout << "	heme coord energy  " << energyhis << endl;
		cout << "	his - por sterics " << hisPorSterics << endl;
		cout << "***************************************************************" << endl;
	}

	return 2*energyhis + 2*hisPorSterics;
}

void createBackbone(protein* _prot, double _rad, double _phase, double _square, double _offset, double _coil)
{
    _prot->rotate(Z_axis, _phase);
    _prot->translate(_rad, 0, _offset);
    _prot->rotate(Z_axis, _square);
    _prot->rotate(1, Z_axis, 180.0);
    _prot->rotate(2, Y_axis, 180.0);
    _prot->rotate(3, X_axis, 180.0);
    _prot->coilcoil(_coil);

    return;
}


void undoBackbone(protein* _prot, double _rad, double _phase, double _square, double _offset, double _coil)
{
    _prot->coilcoil(-1*_coil);
    _prot->rotate(1, Z_axis, 180.0);
    _prot->rotate(2, Y_axis, 180.0);
    _prot->rotate(3, X_axis, 180.0);
    _prot->rotate(Z_axis, -1*_square);
    _prot->translate(-1*_rad, 0, -1*_offset);
    _prot->rotate(Z_axis, -1*_phase);

    return;
}


vector <string> Parse::parse(string& _currentLine)
{
	StrVec parsedStrings;
	string tmpStrChi;
	string tmpStr;
	parsedStrings.resize(0);
	tmpStrChi.resize(1);
	tmpStr.resize(0);

	for (UInt i = 0; i < _currentLine.size(); i++)
	{
		if (_currentLine[i] != ' ' && _currentLine[i] != '\t' && _currentLine[i] != '\n')
		{
			tmpStrChi[0] = _currentLine[i];
			tmpStr.append(tmpStrChi);
		}
		else if (tmpStr.size() != 0)
		{
			parsedStrings.push_back(tmpStr);
			tmpStr.resize(0);
		}
		if ( i == (_currentLine.size() - 1) && tmpStr.size() != 0)
		{
			parsedStrings.push_back(tmpStr);
		}
	}

	for (UInt i = 0; i < parsedStrings.size(); i++)
	{
		cout << parsedStrings[i] << " ";
	}
	cout << endl;
	return parsedStrings;
}

vector <dblVec> getLigandCoords(string _inputFileName)
{
	ifstream inFile;
	inFile.open(_inputFileName.c_str());
	if (!inFile)
	{
		cout << "Unable to find or open file" << endl;
		exit(1);
	}

	string currentLine;
	vector <string> parsedStrings;
	parsedStrings.resize(0);
	vector <dblVec> ligCoords;
	while(getline(inFile, currentLine, '\n'))
	{
		parsedStrings = Parse::parse(currentLine);
		if (parsedStrings[0] == "ATOM")
		{
			double x,y,z;
			sscanf(parsedStrings[5].c_str(), "%lf", &x);
			sscanf(parsedStrings[6].c_str(), "%lf", &y);
			sscanf(parsedStrings[7].c_str(), "%lf", &z);
			dblVec thisAtom(3);
			thisAtom[0] = x; thisAtom[1] = y ; thisAtom[2] = z;
			cout << "xyz " << x << " " << y << " " << z << " " << endl;
			ligCoords.push_back(thisAtom);
		}
	}
	inFile.close();
	return ligCoords;
}

double getProteinLigandEnergy(protein* _prot, vector <dblVec> _ligCoords, double _r0, double _eps, vector <UIntVec> _exclude)
{
	double energy = 0.0;
	for (UInt i = 0; i < _prot->getNumChains(); i ++)
	{
		for (UInt j = 0; j < _prot->getNumResidues(i); j ++)
		{
			bool excludeRes = false;
			for (UInt n = 0; n < _exclude.size(); n ++)
			{
				if (_exclude[n][0] == i && _exclude[n][1] == j)
				{
					excludeRes = true;
				}
			}
			if (!excludeRes)
			{
				for (UInt k = 0; k < _prot->getNumAtoms(i,j); k ++)
				{
					dblVec atomCoords = _prot->getCoords(i,j,k);
					for (UInt l = 0; l < _ligCoords.size(); l ++)
					{
						dblVec diff = atomCoords - _ligCoords[l];
						double distance = sqrt(CMath::dotProduct(diff,diff));
						double distratio = _r0/distance;
						energy += _eps*(pow(distratio,12) - 2* pow(distratio,6));
					}
				}
			}
		}
	}

	return energy;
}

void writeStateToFile(protein* _prot, double _rad, double _phase, double _square, double _offset, double _coil, 
	double _temp, double _energy, string _fileName)
{
	ofstream file(_fileName.c_str(),ios::app);
	file << "his:rad " << _rad << " phase " << _phase << " square " << _square << " offset " << _offset << " coil " << _coil << " temp " << _temp << " energy " << _energy;
	cout << _fileName << " rad " << _rad << " phase " << _phase << " square " << _square << " offset " << _offset << " coil " << _coil << " temp " << _temp << " energy " << _energy;
	file << endl;
	cout << endl;
	file.close();
	return;
}

double getThrHisHBondEnergy(protein* _prot, UInt _Hc, UInt _Hr, UInt _Tc, UInt _Tr)
{
	dblVec HisNd1 = _prot->getCoords(_Hc, _Hr, "ND1");
	dblVec ThrOg = _prot->getCoords(_Tc, _Tr, "OG1");

	dblVec diff = HisNd1 - ThrOg;
	double distance = sqrt(CMath::dotProduct(diff,diff));

	dblVec HisCg = _prot->getCoords(_Hc, _Hr, "CG");
	dblVec HisCe1 = _prot->getCoords(_Hc, _Hr, "CE1");

	dblVec pseudoAtom = 0.5 * (HisCe1 + HisCg);
	dblVec BA = HisNd1 - pseudoAtom;
	dblVec BC = HisNd1 - ThrOg;

	double magBC = sqrt(CMath::dotProduct(BC,BC));
	double magBA = sqrt(CMath::dotProduct(BA,BA));

	double angle = acos( (CMath::dotProduct(BA,BC)) / (magBA * magBC) );
	double eps = 10.0;
	double r0 = 2.8;
	double distratio = r0/distance;
	double angleFactor = pow(cos(angle),2);
	double energy = angleFactor * eps*(5*pow(distratio,12) - 6 * pow(distratio,10));
	//cout << "dist " << distance << " angle " << angle << " anglefactor " << angleFactor << " energy " << energy << endl;
	return energy;
}

double getThrBackboneCarbonylEnergy(protein* _prot, UInt _Tc, UInt _Tr)
{
	UInt numRes = _prot->getNumResidues(_Tc);

	UInt min, max;
	if (( _Tr - 5) >= 0) { min = _Tr - 5; }
	else  { min = 0; }

	if (( _Tr + 5) < numRes) { max = _Tr + 5; }
	else { max = numRes - 1 ; }

	double eps = 5.0;
	double r0 = 2.8;

	dblVec ThrCb = _prot->getCoords(_Tc, _Tr, "CB");
	dblVec ThrOg = _prot->getCoords(_Tc, _Tr, "OG1");
	dblVec BC = ThrOg - ThrCb;
	double magBC = sqrt(CMath::dotProduct(BC,BC));
	double energy = 0.0;
	for (UInt i = min; i <= max ; i ++)
	{
		dblVec carbonyl = _prot->getCoords(_Tc, i, "O");
		dblVec BA = ThrOg - carbonyl;
		double magBA = sqrt(CMath::dotProduct(BA,BA));
		double angle = acos( (CMath::dotProduct(BA,BC)) / (magBA * magBC) );
		double distratio  = r0/magBA;
		energy += cos(cos(angle-1.91)) * eps*(5*pow(distratio,12) - 6 * pow(distratio,10));
	}

	return energy;
}

