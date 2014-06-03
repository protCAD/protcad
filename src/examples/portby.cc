#include <iostream>
#include <string>
#include <vector>
#include "ensemble.h"
#include "PDBInterface.h"
#include "ran.h"
#include "generalio.h"

vector <string> parse(string& _currentLine);
void createBackbone(protein* _prot, double _offset, double _square);
void undoBackbone(protein* _prot, double _offset, double _square);
double getEnergy(protein* _prot, UInt _Hc1, UInt _Hr1, UInt _Hc2, UInt _Hr2, UInt _Tc1, UInt _Tr1, UInt _Tc2, UInt _Tr2);
double getHisPorDist(protein* _prot, UInt _Hc, UInt _Hr);
double getHisPorAngle(protein* _prot, UInt _Hc, UInt _Hr);
void optimizeRotamerGeometry(protein* _prot, UInt _Hc1, UInt _Hr1, UInt _Hc2, UInt _Hr2, UInt _Tc1, UInt _Tr1, UInt _Tc2, UInt _Tr2);
int main(int argc, char* argv[])
{
	enum aminoAcid {A, R, N, D, C, Q, E, G, H, I, L, K, M, F, P, S, T, W, Y, V};
	if (argc < 2)
	{
		cout << "portby inputfile.inp" << endl;
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
	getline(inFile, currentLine, '\n'); parsedStrings = parse(currentLine);
	string inputFile = parsedStrings[0];
	string outputFile = parsedStrings[1];

	// read in prot structure
	PDBInterface* thePDB = new PDBInterface(inputFile);
	ensemble* theEnsemble = thePDB->getEnsemblePointer();
	molecule* theMol = theEnsemble->getMoleculePointer(0);
	protein* prot = static_cast<protein*>(theMol);

	residue::setCutoffDistance(10.0);
	pmf::setScaleFactor(0.0);
	rotamer::setScaleFactor(0.0);
	microEnvironment::setScaleFactor(0.0);
	amberVDW::setScaleFactor(1.0);
	amberVDW::setRadiusScaleFactor(1.0);
	amberVDW::setLinearRepulsionDampeningOff();
	amberElec::setScaleFactor(0.0);
	solvation::setItsScaleFactor(0.0);

	double squaremin, squaremax, squarestep;
	double offsetmin, offsetmax, offsetstep;

	getline(inFile, currentLine, '\n'); parsedStrings = parse(currentLine);
	sscanf(parsedStrings[0].c_str(), "%lf", &offsetmin);
	sscanf(parsedStrings[1].c_str(), "%lf", &offsetmax);
	sscanf(parsedStrings[2].c_str(), "%lf", &offsetstep);

	getline(inFile, currentLine, '\n'); parsedStrings = parse(currentLine);
	sscanf(parsedStrings[0].c_str(), "%lf", &squaremin);
	sscanf(parsedStrings[1].c_str(), "%lf", &squaremax);
	sscanf(parsedStrings[2].c_str(), "%lf", &squarestep);

	UInt hischain1, hisres1, hisrot1, hischain2, hisres2, hisrot2;
	UInt thrchain1, thrres1, thrrot1, thrchain2, thrres2, thrrot2;
	getline(inFile, currentLine, '\n'); parsedStrings = parse(currentLine);
	sscanf(parsedStrings[0].c_str(), "%u", &hischain1);
	sscanf(parsedStrings[1].c_str(), "%u", &hisres1);
	sscanf(parsedStrings[2].c_str(), "%u", &hisrot1);

	getline(inFile, currentLine, '\n'); parsedStrings = parse(currentLine);
	sscanf(parsedStrings[0].c_str(), "%u", &hischain2);
	sscanf(parsedStrings[1].c_str(), "%u", &hisres2);
	sscanf(parsedStrings[2].c_str(), "%u", &hisrot2);

	getline(inFile, currentLine, '\n'); parsedStrings = parse(currentLine);
	sscanf(parsedStrings[0].c_str(), "%u", &thrchain1);
	sscanf(parsedStrings[1].c_str(), "%u", &thrres1);
	sscanf(parsedStrings[2].c_str(), "%u", &thrrot1);

	getline(inFile, currentLine, '\n'); parsedStrings = parse(currentLine);
	sscanf(parsedStrings[0].c_str(), "%u", &thrchain2);
	sscanf(parsedStrings[1].c_str(), "%u", &thrres2);
	sscanf(parsedStrings[2].c_str(), "%u", &thrrot2);

	prot->activateForRepacking(hischain1,hisres1);
	prot->mutate(hischain1, hisres1, H);
	prot->setRotamer(hischain1, hisres1, 0, hisrot1);

	prot->activateForRepacking(hischain2,hisres2);
	prot->mutate(hischain2,hisres2, H);
	prot->setRotamer(hischain2, hisres2, 0, hisrot2);

	prot->activateForRepacking(thrchain1, thrres1);
	prot->mutate(thrchain1, thrres1, T);
	prot->setRotamer(thrchain1, thrres1, 0, thrrot1);

	prot->activateForRepacking(thrchain2, thrres2);
	prot->mutate(thrchain2, thrres2, T);
	prot->setRotamer(thrchain2, thrres2, 0, thrrot2);

	pdbWriter(prot, "start.pdb");
	double offsetbest, squarebest;
	double bestEnergy = 1E10;

	for (double offset = offsetmin; offset <= offsetmax; offset = offset + offsetstep)
	{
		for (double square = squaremin; square <= squaremax; square = square + squarestep)
		{
			createBackbone(prot, offset, square);
			optimizeRotamerGeometry(prot,hischain1, hisres1, hischain2, hisres2, thrchain1, thrres1, thrchain2, thrres2);
			double energy = getEnergy(prot, hischain1, hisres1, hischain2, hisres2, thrchain1, thrres1, thrchain2, thrres2);
			if (energy < bestEnergy)
			{
				bestEnergy = energy;
				squarebest = square;
				offsetbest = offset;
				cout << "NEW LOW - ";
				pdbWriter(prot, "current.pdb");
			}
			cout << "offset  " << offset << " square " << square << " energy " << energy << endl;
			undoBackbone(prot, offset, square);
		}
	}
	createBackbone(prot, offsetbest, squarebest);
	optimizeRotamerGeometry(prot, hischain1, hisres1, hischain2, hisres2, thrchain1, thrres1, thrchain2, thrres2);
	pdbWriter(prot, outputFile);


	return 0;
}

void optimizeRotamerGeometry(protein* _prot, UInt _Hc1, UInt _Hr1, UInt _Hc2, UInt _Hr2, UInt _Tc1, UInt _Tr1, UInt _Tc2, UInt _Tr2)
{
	UIntVec rotamer;

	rotamer = _prot->getCurrentRotamer(_Hc1, _Hr1);
	_prot->setRotamer(_Hc1, _Hr1, 0, rotamer[0]);
	double hischi1standard = _prot->getChi(_Hc1, _Hr1, 0, 0);
	double hischi2standard = _prot->getChi(_Hc1, _Hr1, 0, 1);

	rotamer = _prot->getCurrentRotamer(_Hc2, _Hr2);
	_prot->setRotamer(_Hc2, _Hr2, 0, rotamer[0]);

	rotamer = _prot->getCurrentRotamer(_Tc1, _Tr1);
	_prot->setRotamer(_Tc1, _Tr1, 0, rotamer[0]);
	double thrchi1standard = _prot->getChi(_Tc1, _Tr1, 0,0);

	rotamer = _prot->getCurrentRotamer(_Tc2, _Tr2);
	_prot->setRotamer(_Tc2, _Tr2, 0, rotamer[0]);

	double range = 30.0;
	double step = 5.0;

	double tempangle;
	double bestEnergy = 1E10;
	double besthischi1 = 0.0;
	double besthischi2 = 0.0;
	double bestthrchi1 = 0.0;
	for (double hischi1 = -1*range; hischi1 <= range; hischi1 = hischi1 + step)
	{
		tempangle = hischi1standard + hischi1;
		_prot->setChi(_Hc1, _Hr1, 0, 0, tempangle);
		_prot->setChi(_Hc2, _Hr2, 0, 0, tempangle);

		for (double hischi2 = -1*range; hischi2 <= range; hischi2 = hischi2 + step)
		{
			tempangle = hischi2standard + hischi2;
			_prot->setChi(_Hc1, _Hr1, 0, 1, tempangle);
			_prot->setChi(_Hc2, _Hr2, 0, 1, tempangle);

			for (double thrchi1 = -1*range; thrchi1 <= range; thrchi1 = thrchi1 + step)
			{
				tempangle = thrchi1standard + thrchi1;
				_prot->setChi(_Tc1, _Tr1, 0, 0, tempangle);
				_prot->setChi(_Tc2, _Tr2, 0, 0, tempangle);
				double energy = getEnergy(_prot, _Hc1, _Hr1, _Hc2, _Hr2, _Tc1, _Tr1, _Tc2, _Tr1);
				if (energy < bestEnergy)
				{
					besthischi1 = hischi1;
					besthischi2 = hischi2;
					bestthrchi1 = thrchi1;
					bestEnergy = energy;
				}
			}
		}
	}
	_prot->setChi(_Hc1, _Hr1, 0, 0, besthischi1);
	_prot->setChi(_Hc2, _Hr2, 0, 0, besthischi1);
	_prot->setChi(_Hc1, _Hr1, 0, 1, besthischi2);
	_prot->setChi(_Hc2, _Hr2, 0, 1, besthischi2);
	_prot->setChi(_Tc1, _Tr1, 0, 0, bestthrchi1);
	_prot->setChi(_Tc2, _Tr2, 0, 0, bestthrchi1);

	return;
}


double getEnergy(protein* _prot, UInt _Hc1, UInt _Hr1, UInt _Hc2, UInt _Hr2, UInt _Tc1, UInt _Tr1, UInt _Tc2, UInt _Tr2)
{
	double HisPorDist1 = getHisPorDist(_prot, _Hc1, _Hr1);
	double HisPorDist2 = getHisPorDist(_prot, _Hc2, _Hr2);
	//double HisPorHisAngle = getHisPorHisAngle(_prot, _Hc1, _Hr1, _Hc2, _Hr2);
	double HisPorAngle1 = getHisPorAngle(_prot, _Hc1, _Hr1);
	double HisPorAngle2 = getHisPorAngle(_prot, _Hc2, _Hr2);

	double eps = 10.0;
	double r0 = 2.0;

	double distratio  = r0/HisPorDist1;
	double angleFactor1 = 1.0;
	angleFactor1 *= cos(cos(HisPorAngle1));
	//angleFactor1 *= cos(cos(HisPorHisAngle));
	double energy1 = eps*(5*pow(distratio,12) - 6 * pow(distratio,10));
	energy1 *= angleFactor1;

	distratio  = r0/HisPorDist2;
	double angleFactor2 = 1.0;
	angleFactor2 *= cos(cos(HisPorAngle2));
	//angleFactor2 *= cos(cos(HisPorHisAngle));
	double energy2 = eps*(5*pow(distratio,12) - 6 * pow(distratio,10));
	energy2 *= angleFactor2;

	double energy3 = _prot->intraEnergy();

	return energy1 + energy2 + energy3;
}

double getHisPorDist(protein* _prot, UInt _Hc, UInt _Hr)
{
	dblVec HisNe2 = _prot->getCoords(_Hc, _Hr, "NE2");
	double x = HisNe2[0];
	double y = HisNe2[1];

	double distance = sqrt(x*x + y*y);
	return distance;
}

double getHisPorAngle(protein* _prot, UInt _Hc, UInt _Hr)
{
	dblVec HisCd2 = _prot->getCoords(_Hc, _Hr, "CD2");
	dblVec HisCe1 = _prot->getCoords(_Hc, _Hr, "CE1");
	dblVec HisNe2 = _prot->getCoords(_Hc, _Hr, "NE2");

	dblVec metal(3);
	metal[0]=0.0;
	metal[1]=0.0;
	metal[2]=HisNe2[2];

	dblVec pseudoAtom = 0.5 * (HisCe1 + HisCd2);
	dblVec BA = HisNe2 - pseudoAtom;
	dblVec BC = HisNe2 - metal;

	double magBC = sqrt(CMath::dotProduct(BC,BC));
	double magBA = sqrt(CMath::dotProduct(BA,BA));

	double angle = acos( (CMath::dotProduct(BA,BC)) / (magBA * magBC) );
	return angle;
}

void createBackbone(protein* _prot, double _offset, double _square)
{
	_prot->translate(0, 0,0,_offset/2);
	_prot->translate(1, 0,0,_offset/2);
	_prot->translate(2, 0,0,-1*_offset/2);
	_prot->translate(3, 0,0,-1*_offset/2);
	_prot->rotate(0, Z_axis, _square);
	_prot->rotate(1, Z_axis, _square);
	return;
}

void undoBackbone(protein* _prot, double _offset, double _square)
{
	_prot->rotate(0, Z_axis,-1* _square);
	_prot->rotate(1, Z_axis,-1* _square);
	_prot->translate(0, 0,0,-1* _offset/2);
	_prot->translate(1, 0,0,-1* _offset/2);
	_prot->translate(2, 0,0,_offset/2);
	_prot->translate(3, 0,0,_offset/2);
	return;
}

vector <string> parse(string& _currentLine)
{
	StrVec parsedStrings;
	string tmpStrChi;
	string tmpStr;
	parsedStrings.resize(0);
	tmpStrChi.resize(1);
	tmpStr.resize(0);

	for (UInt i = 0; i < _currentLine.size(); i++)
	{
		if (_currentLine[i] != ' ' && _currentLine[i] != '\t' && _currentLine[i] != '\n' && i != (_currentLine.size()-1))
		{
			tmpStrChi[0] = _currentLine[i];
			tmpStr.append(tmpStrChi);
		}
		else if (tmpStr.size() != 0)
		{
			parsedStrings.push_back(tmpStr);
			tmpStr.resize(0);
		}
	}

	for (UInt i = 0; i < parsedStrings.size(); i++)
	{
		cout << parsedStrings[i] << " ";
	}
	cout << endl;
	return parsedStrings;
}
