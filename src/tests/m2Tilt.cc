#include <iostream>
#include <string>
#include <vector>
#include "ensemble.h"
#include "PDBInterface.h"
#include "ran.h"
#include "generalio.h"

#define PI 3.141596
#define HISTRPDIST 4.0
#define HISTRPFORCE 10
#define NMAX 5000
#define CREATEBACKBONE createHomoTetramer
#define UNDOBACKBONE undoHomoTetramer
#define SCORE getEnsembleScore
#define LARGE 1e20
#define TINY 1e-20
#define GETENERGY getEnergy
#define GETENERGYWITHRESTRAINTS getEnergyWithM2Restraints
#define VecVecUIntVec vector < vector < vector < unsigned int >  > >
#define VecStrVec vector < vector < string > > 
#define VecUIntVec vector < vector < unsigned int > >
#define VecDblVec vector < vector < double > >
#define FIXEDCHAIN 4
enum aminoAcid {A, R, N, D, C, Q, E, G, H, I, L, K, M, F, P, S, T, W, Y, V};
protein* prot;

VecVecUIntVec seqData(0);
StrVec seqEncoding(0);
bool verbose = false;
DblVec refEnergies(0);
// function declarations


void getSequence(string _infile, VecUIntVec &_symmetryData); 
void setSymmetry(VecUIntVec _symmetryData);
VecDblVec getBackbone(string _backboneFile); 
VecDblVec getMCParams(string _mcParmasFile); 
VecDblVec monteCarlo(VecDblVec _mcParams, ran _ranNumber, VecDblVec _backboneParams);
double getEnsembleScore(); 
int mapSequence(int _index); 
void report(VecDblVec _backboneParams, string _pdbFile, int _seqIndex); 
void createHomoTetramer(VecDblVec _backboneParams); 
void undoHomoTetramer(VecDblVec _backboneParams); 
VecDblVec modifyBackbone(VecDblVec _currentBackbone, ran &number); 
UInt getTypeNumFromTypeString(string _aa); 
double getEnergy(); 
double getEnergyWithM2Restraints(); 
double getHBondEnergy(dblVec _BA, dblVec _BC, double _eps, double _r0, double _theta0);
double hisAmdEnergy();
double hisTrpDistanceRestraint(double _dist, double _force); 
double expoCool(DblVec _args, UInt _iteration); 
void optimizeRotamers(VecUIntVec _active);
void alignCB(UInt _chain, UInt _res);
void printSequenceData();
void stripToGlycine();
DblVec calculateRefEnergies();

int main (int argc, char* argv[])
{
	residue::setCutoffDistance(12.0);
	pmf::setScaleFactor(0.0);
	rotamer::setScaleFactor(1.0);
	microEnvironment::setScaleFactor(0.0);
	amberVDW::setScaleFactor(1.0);
	amberVDW::setRadiusScaleFactor(0.95);
	amberVDW::setLinearRepulsionDampeningOff();
	amberElec::setScaleFactor(0.0);
	solvation::setItsScaleFactor(0.0);

	// read in PDB file
	string protFile = argv[1];
    PDBInterface* thePDB = new PDBInterface(protFile);
    ensemble* theEnsemble = thePDB->getEnsemblePointer();
    molecule* theMol = theEnsemble->getMoleculePointer(0);
    prot = static_cast<protein*>(theMol);
	// read in file containing symmetry, phenotype and sequence information
	string seqFile = argv[2];
	VecUIntVec symmetryData(0);
	getSequence(seqFile,symmetryData);
	setSymmetry(symmetryData);
	UIntVec itsIndependentChainsMap = prot->getItsIndependentChainsMap();
	if(!verbose) prot->silenceMessages();
	for (UInt i = 0; i < itsIndependentChainsMap.size() ; i++)
	{
		if (itsIndependentChainsMap[i] != FIXEDCHAIN)
		{
			prot->activateAllForRepacking(itsIndependentChainsMap[i]);
			prot->setCanonicalHelixRotamersOnly(itsIndependentChainsMap[i]);
		}
	}
	if(verbose) printSequenceData();
	// read in backbone params	
	string backboneFile = argv[3];
	VecDblVec backboneParams = getBackbone(backboneFile);

	// align histidine on axis

	UInt hisPos = 1000;
	stripToGlycine();

	dblVec oneCA = prot->getCoords(0,0,"CA");
	dblVec sevCA = prot->getCoords(0,6,"CA");
	if (sevCA[2] > oneCA[2]) 
		for (UInt i = 0; i < 4; i ++)
		{
			prot->rotate(i, Y_axis, 180.0);
		}
	for (UInt i = 1; i < seqData[0][0].size(); i ++)
	{
		if (seqData[0][0][i] == H)  hisPos = i - 1;
	}
	prot->mutateWBC(0,hisPos,H);

	if (hisPos != 1000)
	{
		alignCB(0,hisPos);
		alignCB(1,hisPos);
		alignCB(2,hisPos);
		alignCB(3,hisPos);
	}
	stripToGlycine();
	
	// silence amantadine backbone atoms
	if (prot->getNumChains() > 4)
	{
		for (UInt i = 0; i < 4; i ++)
		{
			prot->makeAtomSilent(FIXEDCHAIN,0,i);
		}
	}
	refEnergies = calculateRefEnergies();
	// read in monte carlo schema
	string mcFile = argv[4];
	VecDblVec mcParams = getMCParams(mcFile);
	UInt seed = 0;
	sscanf(argv[5], "%u", &seed);
	ran ranNumber;
	ranNumber.setSeed(seed);
	VecDblVec minBackboneParams = monteCarlo(mcParams, ranNumber, backboneParams);
	//minBackboneParams = minimize(backboneParams, ranNumber);
	report(minBackboneParams, "final.pdb", 0);


	return 0;
}

DblVec calculateRefEnergies()
{
	DblVec refE(0);  

	prot->translate(0, 500.0, 0);
	UInt numChains = prot->getNumChains();
	
	for (UInt i = 0; i < numChains; i ++)
	{
		double angle = (double)i*360.0/(double)numChains;
		prot->rotate(i, Z_axis, angle); 
	}

	if (verbose) pdbWriter(prot, "reference.pdb");

	UIntVec itsIndependentChainsMap = prot->getItsIndependentChainsMap();
	vector < vector < int > > itsChainLinkageMap = prot->getItsChainLinkageMap();

	for (UInt i = 0; i < seqData.size(); i ++)
	{
		stripToGlycine();
		mapSequence((int)i);
		prot->optimizeRotamers();
		double energy = 0.0;
		for (UInt j = 0; j < itsIndependentChainsMap.size(); j ++)
		{
			double thisE = prot->intraEnergy(itsIndependentChainsMap[j]);
			thisE = thisE * itsChainLinkageMap[j].size();
			energy += thisE;
		}
		refE.push_back(energy);
	}

	for (UInt i = 0; i < numChains; i ++)
	{
		double angle = -1*(double)i*360.0/(double)numChains;
		prot->rotate(i, Z_axis, angle); 
	}
	prot->translate(0, -500.0, 0);


	return refE;
}

void alignCB(UInt _chain, UInt _res)
{
	if (verbose) cout << " calling alignCB " << endl;
	dblVec hisCB = prot->getCoords(_chain, _res, "CB");
    dblVec tempX(3);
    tempX[0]=0.0;
    tempX[1]=-1.0;
    tempX[2]=0.0;
    hisCB[2]=0.0;

    double magCB = sqrt(CMath::dotProduct(hisCB,hisCB));
    double magX = sqrt(CMath::dotProduct(tempX,tempX));

    double hisCBAngle = acos(CMath::dotProduct(hisCB,tempX)/(magX*magCB));
    hisCBAngle = (180/PI)*hisCBAngle;
    if(verbose)cout << "his angle " << hisCBAngle << endl;
    prot->rotate(_chain, Z_axis, hisCBAngle);

	return;
}

void createHomoTetramer(VecDblVec _backboneParams)
{
	if (_backboneParams.size() != 6 )
	{
		cout << "ERROR - Backbone parameters not compatable with homotetramer" << endl;
		exit(1);
	}
	if (verbose) cout << "creating backbone ... " << endl;
	double radius = _backboneParams[0][0];
	double phase = _backboneParams[1][0];
	double cross = _backboneParams[2][0];
	double offset = _backboneParams[3][0];
	double zmove = _backboneParams[4][0];
	double zrot = _backboneParams[5][0];
	for (UInt i = 0; i < 4; i ++)
	{	
		prot->rotate(i,Z_axis, phase);
		prot->translate(i,0.0, 0.0, offset);
		prot->rotate(i,Y_axis, cross);
		prot->translate(i,0.0,0.0, -1*offset);
		prot->translate(i,0.0, radius, 0.0);
		prot->translate(i,0.0, 0.0, zmove);
		prot->rotate(i,Z_axis, zrot);
	}
	prot->rotate(0, Z_axis, 90.0);
	prot->rotate(1, Z_axis, 180.0);
	prot->rotate(2, Z_axis, 270.0);
	

	return;
}

void undoHomoTetramer(VecDblVec _backboneParams)
{
	if (_backboneParams.size() != 6 )
	{
		cout << "ERROR - Backbone parameters not compatible with homotetramer" << endl;
		exit(1);
	}
	if (verbose) cout << "undoing backbone ... " << endl;
	double radius = _backboneParams[0][0];
    double phase = _backboneParams[1][0];
    double cross = _backboneParams[2][0];
	double offset = _backboneParams[3][0];
	double zmove = _backboneParams[4][0];
	double zrot = _backboneParams[5][0];



    prot->rotate(0, Z_axis, -90.0);
    prot->rotate(1, Z_axis, -180.0);
    prot->rotate(2, Z_axis, -270.0);
	for (UInt i = 0; i < 4; i ++)
	{
		prot->rotate(i,Z_axis, -1*zrot);
		prot->translate(i,0.0, 0.0, -1*zmove);
    	prot->translate(i,0.0, -1*radius, 0.0);
		prot->translate(i,0.0,0.0,offset);
    	prot->rotate(i,Y_axis, -1*cross);
		prot->translate(i,0.0, 0.0, -1*offset);
    	prot->rotate(i,Z_axis, -1*phase);
	}
	return;
}

double expoCool(DblVec _args, UInt _i)
{
	if (verbose) cout << "calling expocool" << endl;
	double bi = 0;
	double b0 = _args[0];
	double bN = _args[1];
	double N = _args[2];

	bi = b0 * pow( bN/b0, (double)_i / N);
	return bi;
}

VecDblVec modifyBackbone(VecDblVec _currentBackbone, ran &ranNumber)
{
	if (verbose) cout << "calling modifyBackobne" << endl;
	UInt param = (UInt)(ranNumber.getNext()*(double)_currentBackbone.size());
	if (param > (_currentBackbone.size() -1)) param = _currentBackbone.size() -1;	

	double min = _currentBackbone[param][1];
	double max = _currentBackbone[param][2];
	double newVal = min + ranNumber.getNext()*(max - min);

	_currentBackbone[param][0] = newVal;
	return _currentBackbone;
}


double getEnsembleScore()
{
	if (verbose) cout << "calling getEnsembleScore" << endl;
	stripToGlycine();
	UInt wt = 0;
	for (UInt i = 0; i < seqEncoding.size(); i ++)
	{
		if (seqEncoding[i] == "w") wt = i;
	}
	mapSequence((int)wt);
	//prot->saveCurrentState();
	double wtE = GETENERGYWITHRESTRAINTS() - refEnergies[wt];

	double aD = -100.0;
	double aS =  100.0;
	double g = 3e-5;
	double bBS = -0.7;
	double bSB = -0.2;


	DblVec unfavEbs(0);
	DblVec unfavEsb(0);
	DblVec favEbs(0);
	DblVec favEsb(0);
	for (UInt i = 0; i < seqData.size(); i ++)
	{
		if (seqEncoding[i] == "ubs")
		{
			mapSequence((int)i);
			unfavEbs.push_back(GETENERGYWITHRESTRAINTS() - 4.0 * refEnergies[i]);
		}
		if (seqEncoding[i] == "usb")
		{
			mapSequence((int)i);
			unfavEsb.push_back(GETENERGYWITHRESTRAINTS() - 4.0 * refEnergies[i]);
		}
		if (seqEncoding[i] == "fbs")
		{
			mapSequence((int)i);
			favEbs.push_back(GETENERGYWITHRESTRAINTS() - 4.0 * refEnergies[i]);
		}
		if (seqEncoding[i] == "fsb")
		{
			mapSequence((int)i);
			favEsb.push_back(GETENERGYWITHRESTRAINTS() - 4.0 * refEnergies[i]);
		}		
	}
	
	double e = 2.71828182846;
	double pDbs = 0.0;
	double pDsb = 0.0;
	double pSbs = 0.0;
	double pSsb = 0.0;
	//DISRUPTIVE Big to Small Mutation
	for (UInt i = 0; i < unfavEbs.size(); i ++)
	{
		double diff = unfavEbs[i]-wt;
		double thisPD = 1.0 / ( pow(e,bBS*diff) + 1.0);
		pDbs += thisPD;
	}
	if (unfavEbs.size() !=  0) pDbs = pDbs / (double)unfavEbs.size();
	if (!(pDbs >= 0.0 && pDbs <= 1.0))
	{
		cout << "BAD PARAMS" << endl;
		return 0.0;
	}
	//DISRUPTIVE Small to Big Mutation
	for (UInt i = 0; i < unfavEsb.size(); i ++)
	{
		double diff = unfavEsb[i]-wt;
		double thisPD = 1.0 / ( pow(e,bSB*diff) + 1.0);
		pDsb += thisPD;
	}
	if (unfavEsb.size() !=  0) pDsb = pDsb / (double)unfavEsb.size();
	if (!(pDsb >= 0.0 && pDsb <= 1.0))
	{
		cout << "BAD PARAMS" << endl;
		return 0.0;
	}
	//SILENT Big to Small Mutation
	for (UInt i = 0; i < favEbs.size(); i ++)
	{
		double diff = favEbs[i]-wt;
		double thisPS = 1.0 / ( pow(e,bBS*diff) + 1.0);
		pSbs += thisPS;
	}
	if (favEbs.size() !=  0) pSbs = pSbs / (double)favEbs.size();
	if (!(pSbs >= 0.0 && pSbs <= 1.0))
	{
		cout << "BAD PARAMS" << endl;
		return 0.0;
	}
	//SILENT Small to Big Mutation
	for (UInt i = 0; i < favEsb.size(); i ++)
	{
		double diff = favEsb[i]-wt;
		double thisPS = 1.0 / ( pow(e,bSB*diff) + 1.0);
		pSsb += thisPS;
	}
	if (favEsb.size() !=  0) pSsb = pSsb / (double)favEsb.size();
	if (!(pSsb >= 0.0 && pSsb <= 1.0))
	{
		cout << "BAD PARAMS" << endl;
		return 0.0;
	}

	double disE = aD*log(pDbs + pDsb + g);
	double favE = aS*log(pSbs + pSsb + g);
	double thisScore = wtE + disE + favE;
	cout << endl;
	cout << "score: wtE " << wtE << " dis " << disE << " fav " << favE << " tot " <<thisScore << endl;
	return thisScore;
}

double getEnergy()
{
	if (verbose) cout << "calling getEnergy" << endl;
	double energy = 0.0;
	for (UInt i = 0; i < prot->getNumChains(); i ++)
	{
		energy += prot->intraEnergy(i);
		for (UInt j = i+1; j < prot->getNumChains(); j ++)
		{
			energy += prot->intraEnergy(i,j);
		}
	}
	return energy;
}

double getEnergyWithM2Restraints()
{
	if (verbose) cout << "calling getEnergyWithM2Restraints" << endl;
	double energy = 0.0;
	energy += GETENERGY();
	energy += hisTrpDistanceRestraint(HISTRPDIST, HISTRPFORCE);
	if (prot->getNumChains() > 4) energy += hisAmdEnergy();
	return energy;
}

double hisAmdEnergy()
{
	UInt hisNum = 1000;
	for (UInt i = 0; i < prot->getNumResidues(0); i ++)
	{
		if (prot->getTypeFromResNum(0,i) == H) hisNum = i;
	}	
	if (hisNum == 1000)
	{
		cout << "HIS not found in the sequence" << endl;
		return 0.0;
	}
	dblVec HisNE2 = prot->getCoords(0,hisNum,"NE2");
	dblVec AmdCB9 = prot->getCoords(4,0,"CB9");
	dblVec AmdNB0 = prot->getCoords(4,0,"NB0");

	dblVec BA = AmdNB0 - HisNE2;
	dblVec BC = AmdNB0 - AmdCB9;

	double angle = 120.0 * PI / 180.0;

	double HBenergy = getHBondEnergy(BA,BC,1.0,2.8,angle);
	cout << "HIS AMD HBond " << HBenergy << " kcals " << endl;
	return 4.0 * HBenergy;
}

double getHBondEnergy(dblVec _BA, dblVec _BC, double _eps, double _r0, double _theta0)
{
	double magBC = sqrt(CMath::dotProduct(_BC, _BC));
	double magBA = sqrt(CMath::dotProduct(_BA, _BA));

	double angle = acos( (CMath::dotProduct(_BA,_BC)) / (magBA * magBC) );

	double distratio = _r0/magBA;
	double angleFactor = pow(cos(angle - _theta0),2);
	double energy = angleFactor * _eps*(5*pow(distratio,12) - 6 * pow(distratio,10));

	return energy;
}


double hisTrpDistanceRestraint(double _dist, double _force)
{
	if (verbose) cout << " calling hipTrpDistanceRestraint" << endl;
	UInt hisNum = 1000;
	UInt trpNum = 1000;
	for (UInt i = 0; i < prot->getNumResidues(0); i ++)
	{
		if (prot->getTypeFromResNum(0,i) == H) hisNum = i;
		if (prot->getTypeFromResNum(0,i) == W) trpNum = i;
	}
	if (hisNum == 1000 || trpNum == 1000)
	{
		cout << "HIS or TRP not found in the sequence" << endl;
		return 0.0;	
	}
	dblVec HisND1 = prot->getCoords(0,hisNum,"ND1");

	double distance = 1000.0;
	for (UInt i = 1; i < 4; i ++)
	{
		dblVec TrpCG = prot->getCoords(i, trpNum, "CG");
		dblVec diff = HisND1 - TrpCG;
		double thisDist = sqrt(CMath::dotProduct(diff,diff));
		if (thisDist < distance) distance = thisDist;
	}
	double diff = distance - _dist;
	double force = _force * diff * diff;

	if (distance < _dist) force = 0.0;
	cout << "HW restraint - " << distance << " E " << force << endl;

	return force;
}

VecDblVec getBackbone(string _backboneFile)
{
	if (verbose) cout << "getBackbone reading in backbone file " << endl;
	ifstream inFile;
	inFile.open(_backboneFile.c_str());
	if (!inFile)
	{
		cout << "Unable to find backbone file " << inFile << endl;
		exit(1);
	}
	string currentLine;
	StrVec parsedStrings;
	VecDblVec backbone(0);
	while (getline(inFile, currentLine, '\n'))
	{
		parsedStrings = Parse::parse(currentLine);
		DblVec chainBackbone(0);
		for (UInt i = 0; i < parsedStrings.size(); i ++)
		{
			double param;
			sscanf(parsedStrings[i].c_str(), "%lf", &param);
			chainBackbone.push_back(param);
		}
		backbone.push_back(chainBackbone);
	}
	inFile.close();
	return backbone;
}
	
VecDblVec getMCParams(string _mcParamsFile)
{
	if (verbose) cout << "getMCParams reading in monte carlo paramter files" << endl;
	ifstream inFile;
	inFile.open(_mcParamsFile.c_str());
	if (!inFile)
	{
		cout  << "Unable to find monte carlo params file " << inFile << endl;
		exit(1);
	}
	string currentLine;
	StrVec parsedStrings;
	VecDblVec mcParams(0);
	while (getline(inFile, currentLine, '\n'))
	{
		parsedStrings = Parse::parse(currentLine);
		DblVec cycleParams(0);
		for (UInt i = 0; i < parsedStrings.size(); i ++)
		{
			double param;
			sscanf(parsedStrings[i].c_str(), "%lf", &param);
			cycleParams.push_back(param);
		}
		mcParams.push_back(cycleParams);
	}
	inFile.close();
	return mcParams;
}

void getSequence(string _infile, VecUIntVec &_symmetryData)
{
	if (verbose) cout << "getSequence reading in the sequence data" << endl;
	ifstream inFile;
	inFile.open(_infile.c_str());
	if (!inFile)
	{
		cout << "Unable to find sequence file " << _infile << endl;
		exit(1);
	}
	string currentLine;
	StrVec parsedStrings;
	VecUIntVec tempProtSequence(0);
	bool first = true;
	while(getline(inFile, currentLine, '\n'))
	{
		parsedStrings = Parse::parse(currentLine);
		cout << currentLine << endl;
		if (parsedStrings[0] == "w" || parsedStrings[0] == "fbs" || parsedStrings[0] == "ubs" 
			|| parsedStrings[0] == "fsb" || parsedStrings[0] == "usb" || parsedStrings[0] == ".")
		{
			cout << parsedStrings[0] << " ";
			if (parsedStrings[0] != ".")
			{
				if (!first) seqData.push_back(tempProtSequence);
				seqEncoding.push_back(parsedStrings[0]);
				tempProtSequence.resize(0);
				first = false;
			}
			UIntVec tempChainSequence(0);
			UInt chainId;
			sscanf(parsedStrings[1].c_str(), "%u", &chainId);
			tempChainSequence.push_back(chainId);
			for (UInt j = 2; j < parsedStrings.size(); j++)
			{
				if (verbose) cout << parsedStrings[j] << " ";
				UInt resId = getTypeNumFromTypeString(parsedStrings[j]);
				if (resId == 1000)
				{
					cout << "ERROR - bad amino acid designation ... quitting." << endl;
 					exit(1);
				}
				tempChainSequence.push_back(resId);
			}
			if (verbose) cout << endl;
			tempProtSequence.push_back(tempChainSequence);
		}
		else if (parsedStrings[0] == "sym")
		{
			UInt master, slave;
			sscanf(parsedStrings[1].c_str(), "%u", &slave);
			sscanf(parsedStrings[2].c_str(), "%u", &master);
			UIntVec AtoB(0);
			AtoB.push_back(slave);
			AtoB.push_back(master);
			_symmetryData.push_back(AtoB);
		}
		else cout << currentLine << " not identified " << endl;		
	}
	inFile.close();
	return;
}


int mapSequence(int _index)
{
	if (verbose) cout << "calling mapSequence" << endl;
	if (_index >=0 && _index < (int)seqData.size())
	{
		cout << "Mapping Sequence " << _index << endl;
		VecUIntVec protSequence = seqData[_index];
		VecUIntVec activePositions(0);
		for (UInt i = 0; i < protSequence.size(); i ++)
		{
			UIntVec chainSequence = protSequence[i];
			UInt chainnum = chainSequence[0];
			for (UInt j = 1; j < chainSequence.size(); j ++)
			{
				UInt resnum = j - 1;
				if (prot->getTypeFromResNum(chainnum, resnum) != chainSequence[j] || _index == 0)
				{
					prot->mutateWBC(chainnum, resnum, chainSequence[j]);
					if (chainSequence[j] != A && chainSequence[j] != G && chainSequence[j] != V)
					{ 
						UIntVec temp(0);
						temp.push_back(chainnum);
						temp.push_back(resnum);
						activePositions.push_back(temp);
					}
				}
			}
			if (activePositions.size() > 0) optimizeRotamers(activePositions);
		}
	}
 	else 
	{

		cout << "Index " << _index << " in mapSequence out of range" << endl;
		return 1;
	}
	return 0;
}



void report(VecDblVec _backbone, string _pdbFile, int _index)
{
	if (verbose) cout << "calling report" << endl;
	cout << "params report" << endl;
	for (UInt i = 0; i < _backbone.size(); i ++)
	{
		cout << "backbone " << i << " " << endl;
		for (UInt j = 0; j < _backbone[i].size(); j ++)
		{
			cout << _backbone[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl;
	if (_index != -1)
	{
		CREATEBACKBONE(_backbone);
		//prot->undoState();
		//prot->saveCurrentState();
		mapSequence(_index);
		pdbWriter(prot, _pdbFile);
		UNDOBACKBONE(_backbone);
		string logfile = _pdbFile + ".out";
		ofstream fout(logfile.c_str());
		for (UInt i = 0; i < _backbone.size(); i ++)
		{
			fout << "backbone " << i << " " << endl;
			for (UInt j = 0; j < _backbone[i].size(); j ++)
			{
				fout << _backbone[i][j] << " ";
			}
			fout << endl;
		}
		fout.close();
	}
}


void setSymmetry(VecUIntVec _symmetryData)
{
	if (verbose) cout << "calling setSymmetry" << endl;
	for (UInt i = 0; i < _symmetryData.size(); i ++)
	{
		if (_symmetryData[i][0] < prot->getNumChains() && _symmetryData[i][1] < prot->getNumChains())
			prot->symmetryLinkChainAtoB(_symmetryData[i][0], _symmetryData[i][1]);
		else
		{
			cout << "symmetry data not compatible with structure" << endl;
			exit(1);
		}
	}
	return;
}


VecDblVec monteCarlo(VecDblVec _mcParams, ran _ranNumber, VecDblVec _backboneParams)
{
	if (verbose) cout << "calling monteCarlo" << endl;
    UInt wt = 0;
    for (UInt i = 0; i < seqEncoding.size(); i ++)
    {
        if (seqEncoding[i] == "w") wt = i;
    }
	double e =  2.71828182846;
	VecDblVec bestBackbone = _backboneParams;
	VecDblVec currentBackbone = _backboneParams;
	for (UInt cycle = 0; cycle < _mcParams.size(); cycle ++)
	{
		cout << "cycle " << cycle << endl;
		double bestScore = LARGE;
		double currentScore = bestScore;
		// _mParams is start, end, iterations
		DblVec coolfuncargs = _mcParams[cycle];
		double iterations = _mcParams[cycle][2];

		for ( double i = 0.0; i < iterations; i ++)
		{
			double beta = expoCool(coolfuncargs, (UInt)floor(i));
			cout << "Iteration:  " << i << " temp " << beta << endl;
			VecDblVec tmpBackbone = modifyBackbone(currentBackbone,_ranNumber);
			CREATEBACKBONE(tmpBackbone);
			double score;
			if (cycle > 0  &&  cycle == (_mcParams.size() - 1)) 
			{
				stripToGlycine();
				mapSequence((int)wt);
				score = GETENERGYWITHRESTRAINTS() - refEnergies[wt];
			}
			else score = SCORE();
			UNDOBACKBONE(tmpBackbone);
			if (score < bestScore)
			{
				bestScore = score;
				bestBackbone = tmpBackbone;
				cout << "best" << endl;
				report(bestBackbone, "best.pdb", 0);
			}
			if (score < currentScore)
			{
				currentScore = score;
				currentBackbone = tmpBackbone;
				cout << "ACCEPTED GOOD " << score <<  endl;
			}
			else
			{	
				double prob = pow(e, -1.0 * beta * (score - currentScore));
				if (_ranNumber.getNext() < prob)
				{
					currentScore = score;
					currentBackbone = tmpBackbone;
					cout << "ACCEPTED METROPOLIS " << score << endl;
				}
				else
				{
					cout << "REJECTED " << score << endl;
				}
			}
			cout << "current " << currentScore << endl;
			report(currentBackbone, "current.pdb", 0);
		}
	}
	return bestBackbone;
}	
/*
void get_psum(dblMat &p, DblVec &psum)
{
    int i,j;
    double sum;

    int mpts=p.num_rows();
    int ndim=p.num_cols();
    for (j=0; j<ndim; j++)
    {
        for (sum= 0.0, i = 0 ;i < mpts; i ++)
        {
            sum += p[i][j];
            psum[j] = sum;
        }
    }
    return;
}

void swap(double &x, double &y)
{
    double temp = x;
    x = y;
    y = temp;
    return;
}

void ameoba(dblMat &p, DblVec &y, const double ftol, int &nfunk)
{
    int i, ihi, ilo, inhi, j;
    double rtol, ysave, ytry;
    int mpts=p.num_rows();
    int ndim=p.num_cols();
    DblVec psum(ndim);
    nfunk=0;
    get_psum(p,psum);
    for (;;)
    {
        // first determine highest point
        ilo = 0;
        if (y[0] > y[1])
        {
            ihi = 0; inhi = 1;
        }
        else
        {
            inhi = 0; ihi = 1;
        }
        for(i=0; i<mpts; i++)
        {
            if (y[i] <= y[ilo]) ilo = i;
            if (y[i] > y[ihi])
            {
                inhi = ihi;
                ihi = i;
            }
            else
            {
                if (y[i] > y[inhi] && i != ihi) inhi = i;
            }
        }
        rtol = 2.0*fabs(y[ihi]-y[ilo])/(fabs(y[ihi])+fabs(y[ilo])+TINY);
        if (rtol < ftol)
        {
            swap(y[0], y[ilo]);
            for (i=0;i<ndim;i++) swap (p[0][i],p[ilo][i]);
            break;
        }

        if (nfunk >= NMAX)
        {
            cout << "NMAX exceeded" << endl;
            return;
        }
        nfunk += 2;

        ytry = amotry(p,y,psum,ihi,-1.0);
        if (ytry < y[ilo]) ytry=amotry(p,y,psum,ihi,2.0);
        else if (ytry >= y[inhi])
        {
            ysave=y[ihi];
            ytry=amotry(p,y,psum,ihi,0.5);
            if (ytry >= ysave)
            {
                for (i=0; i < mpts; i++)
                {
                    if (i != ilo)
                    {
                        for (j=0; j < ndim; j++) p[i][j]=psum[j]=0.5*(p[i][j]+p[ilo][j]);
						VecDblVec backbone(0);
						backbone.push_back(psum);
						CREATEBACKBONE(backbone);
						stripToGlycine();
						mapSequence(0);
                        y[i]=SCORE();
						UNDOBACKBONE(backbone);
                    }
                }
                nfunk += ndim;
                get_psum(p,psum);
            }
        }
        else --nfunk;
    }
    return;
}


double amotry(dblMat &p, DblVec &y, DblVec &psum, const int ihi, const double fac)
{
    int j;
    double fac1, fac2, ytry;

    int ndim=p.num_cols();
    DblVec ptry(ndim);
    fac1=(1.0-fac)/ndim;
    fac2=fac1-fac;
    for (j=0; j < ndim; j++) ptry[j]=psum[j]*fac1-p[ihi][j]*fac2;
	VecDblVec backbone(0);
	backbone.push_back(ptry);
	CREATEBACKBONE(backbone);
	stripToGlycine();
	mapSequence(0);
    ytry=SCORE();
	UNDOBACKBONE(backbone);
    if (ytry < y[ihi])
    {
        y[ihi]=ytry;
        for (j=0; j<ndim; j++)
        {
            psum[j] += ptry[j] - p[ihi][j];
            p[ihi][j]=ptry[j];
        }
    }
    return ytry;
}
*/
void optimizeRotamers(VecUIntVec _active)
{
	if (verbose) cout << "calling optimizeRotamers" << endl;
	VecUIntVec newActive(0);
	VecUIntVec notActive(0);
	double cutoff = residue::getCutoffDistance();
	residue::setCutoffDistance(4.0);
	for (UInt i = 0; i < _active.size(); i ++)
	{
		UInt activeChain = _active[i][0];
		UInt activeRes = _active[i][1];

		bool touches = false;
		UInt chain = 0;

		while (!touches && chain < prot->getNumChains())
		{
			if (chain == activeChain) chain ++;
			else
			{
				for (UInt j = 0; j < prot->getNumResidues(chain); j ++)  	
				{
					double energy = prot->getResPairEnergy(activeChain, activeRes, chain, j);
					if (energy > 1e-3 || energy < -1e-3) touches = true;
				}
				chain ++;
			}
		} 
		if (touches) newActive.push_back(_active[i]);
		else notActive.push_back(_active[i]);
	}
	residue::setCutoffDistance(cutoff);
	prot->optimizeRotamers(newActive);
	prot->optimizeRotamers(notActive);
	cout << "\nbefore contact fix " << _active.size() << " after " << newActive.size() << endl;
	return;
}

void printSequenceData()
{
	if (verbose) cout << "calling printSequenceData" << endl;
	cout << "************printing sequence data*****************" << endl;
	cout << "total number of sequences " << seqData.size() << endl;
	cout << "total number of encodings " << seqEncoding.size() << endl;
	for (UInt i = 0; i < seqData.size(); i ++)
	{
		cout <<"seq " << i << " enc " << seqEncoding[i] << " num chains " 
			<< seqData[i].size() << " numRes " << seqData[i][0].size() << endl;
		for (UInt j = 0; j < seqData[i].size(); j ++)
		{
			UIntVec sequence = seqData[i][j];
			cout << " chain " << sequence[0] << ": ";
			for (UInt k = 1; k < sequence.size(); k ++)
			{
				cout << sequence[k] << " ";
			}
			cout << endl;
		}
	}
	cout << "***************************************************" << endl;
	return;
}


UInt getTypeNumFromTypeString(string _aa)
{
	if (verbose) cout << "calling getTypeNumFromTypeString" << endl;
	UInt resId = 1000;
	if (_aa == "A" || _aa == "ALA") resId = A;
	if (_aa == "R" || _aa == "ARG") resId = R;
	if (_aa == "N" || _aa == "ASN") resId = N;
	if (_aa == "D" || _aa == "ASP") resId = D;
	if (_aa == "C" || _aa == "CYS") resId = C;
	if (_aa == "Q" || _aa == "GLN") resId = Q;
	if (_aa == "E" || _aa == "GLU") resId = E;
	if (_aa == "G" || _aa == "GLY") resId = G;
	if (_aa == "H" || _aa == "HIS") resId = H;
	if (_aa == "I" || _aa == "ILE") resId = I;
	if (_aa == "L" || _aa == "LEU") resId = L;
	if (_aa == "K" || _aa == "LYS") resId = K;
	if (_aa == "M" || _aa == "MET") resId = M;
	if (_aa == "F" || _aa == "PHE") resId = F;
	if (_aa == "P" || _aa == "PRO") resId = P;
	if (_aa == "S" || _aa == "SER") resId = S;
	if (_aa == "T" || _aa == "THR") resId = T;
	if (_aa == "W" || _aa == "TRP") resId = W;
	if (_aa == "Y" || _aa == "TYR") resId = Y;
	if (_aa == "V" || _aa == "VAL") resId = V;

	return resId;
}

void stripToGlycine()
{
	if (verbose) cout << "calling stripToGlycine" << endl;
	if (verbose) pdbWriter(prot, "preStriptoGly.pdb");
	UIntVec itsIndependentChainsMap = prot->getItsIndependentChainsMap();
	for (UInt i = 0; i < itsIndependentChainsMap.size()  ; i ++)
	{
		UInt chain = itsIndependentChainsMap[i];
		if (chain != FIXEDCHAIN)
		{
			for (UInt j = 0; j < prot->getNumResidues(chain); j++)
			{
				UInt type = prot->getTypeFromResNum(chain,j);
				if (type < 20 && type != G) 
				{
					if (verbose) cout << "stripping " << chain << " " << j << " " << prot->getTypeFromResNum(chain,j) << endl;
					prot->mutateWBC(chain,j,G);
				}
			}
		}
	}
	return;
}
