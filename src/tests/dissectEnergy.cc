#include "pdbReader.h"
#include "pdbWriter.h"
#include "atomIterator.h"
#include "ran.h"
#include "annealer.h"
int main (int argc, char* argv[])
{
	if (argc < 8)
	{
		cout << "dissectEnergy\n\t(1) filename\n\t(2) pmf scale\n\t(3) micro environment scale\n\t(4) AMBER vdW scale\n\t(5) AMBER electrostatics scale\n\t(6) electrostatics distant dependent dielectric(0-off, 1-on)\n\t(7) vdW radius scale factor\n\t(8) vdW linear repulsion dampening(0-off, 1-on)\n";
		exit(1);
	}

	string fileName = argv[1];
	molecule* pMol = pdbReader(fileName);

	protein* pProt = static_cast<protein*>(pMol);
	if (pProt == 0) return 1;
	cout << fileName << " read in." << endl;

	residue::setCutoffDistance(20.0);

	double tmpDbl; 
	sscanf(argv[2], "%lf", &tmpDbl);
	pmf::setScaleFactor(tmpDbl);
	rotamer::setScaleFactor(1.0);
	sscanf(argv[3], "%lf", &tmpDbl);
	microEnvironment::setScaleFactor(tmpDbl);
	sscanf(argv[4], "%lf", &tmpDbl);
	amberVDW::setScaleFactor(tmpDbl);
	sscanf(argv[5], "%lf", &tmpDbl);
	amberElec::setScaleFactor(tmpDbl);
	sscanf(argv[6], "%lf", &tmpDbl);
	if (tmpDbl == 0) amberElec::distanceDependanceOff();
	if (tmpDbl == 1) amberElec::distanceDependanceOn();
	sscanf(argv[7], "%lf", &tmpDbl);
	amberVDW::setRadiusScaleFactor(tmpDbl);
	sscanf(argv[8], "%lf", &tmpDbl);
	if (tmpDbl == 1) amberVDW::setLinearRepulsionDampeningOn();
	if (tmpDbl == 0) amberVDW::setLinearRepulsionDampeningOff();

	cout << "*****************************************" << endl;
	cout << "           energy parameters             " << endl;
	cout << endl;
	cout << "interaction cutoff: " << residue::getCutoffDistance() << endl;
	cout << "pmf scale: " << pmf::getScaleFactor() << endl;
	cout << "rotamer scale: " << rotamer::getScaleFactor() << endl;
	cout << "microEnvironment scale: " << microEnvironment::getScaleFactor() << endl;
	cout << "AMBER vdW scale: " << amberVDW::getScaleFactor() << endl;
	cout << "\trepulsive scale: " << amberVDW::getRepulsionScaleFactor() << endl;
	cout << "\tattractive scale: " << amberVDW::getAttractionScaleFactor() << endl;
	cout << "\tradius scale factor: " << amberVDW::getRadiusScaleFactor() << endl;
	cout << "\tlinear repulsion dampening: "; if (amberVDW::linearRepulsionDampening) cout << "on"; else cout << "off"; cout << endl;
//	cout << "\t1-4 atom vdW scale: " << residue::getOneFourVDWScaleFactor() << endl;
 	cout << "AMBER electrostatics scale: " << amberElec::getScaleFactor() << endl;
	cout << "\tdistance dependent dielectric: "; if (amberElec::isDistanceDependanceOn()) cout << "on"; else cout << "off"; cout << endl;
//	cout << "\t1-4 atom electrostatics scale: " << residue::getOneFourAmberElecScaleFactor() << endl;
	cout << endl;
	cout << "*****************************************" << endl;

	chain* tmpChain = new chain;
	residue * tmpRes = new residue;
	chain* tmpChain2 = new chain;
	residue * tmpRes2 = new residue;
	double totalEnergy = 0.0;
	for (UInt i = 0; i < pProt->getNumChains(); i ++)
	{
		tmpChain = pProt->getChain(i);
		for (UInt j = 0; j < tmpChain->getNumResidues(); j++)
		{	
			tmpRes = tmpChain->getResidue(j);		
			char chain = 'A' + i;	
			cout << "chain " << chain << " res " << tmpRes->getType() << " " << j << "\t";
			double interEnergy = 0.0;
			double intraEnergy = tmpRes->intraEnergy();
			for (UInt k = 0; k < pProt->getNumChains(); k++)
			{
				char chain2 = 'A' + k;
				tmpChain2 = pProt->getChain(k);
				for (UInt l = 0; l < tmpChain2->getNumResidues(); l++)
				{
					if ( l!=j )
					{
						 tmpRes2 = tmpChain2->getResidue(l);
						 interEnergy += tmpRes->interEnergy(tmpRes2);
						 cout << "intraenergy of " << tmpRes->getType() << " " << chain << j << " and " << tmpRes2->getType() << " " << chain2  << l << " is " << tmpRes->interEnergy(tmpRes2) <<  endl;
					}
				 	else if (i != k)
                    {
                         tmpRes2 = tmpChain2->getResidue(l);
                         interEnergy += tmpRes->interEnergy(tmpRes2);
                         cout << "intraenergy of " << tmpRes->getType() << " " << chain << j << " and " << tmpRes2->getType() << " " <<chain2 << l<< " is " << tmpRes->interEnergy(tmpRes2)  << endl;
                    } 
				}
			}
			cout << "inter \t" <<interEnergy << "\t intra \t" << intraEnergy << "\t rotamer \t" << tmpChain->rotamerEnergy(j) << endl;
			totalEnergy = totalEnergy + interEnergy + intraEnergy + tmpChain->rotamerEnergy(j);
		}
	}
//	cout << "total energy:  " <<  pProt->intraEnergy() << endl;
	

	// test chain connectivity;
	for (UInt i = 1; i < pProt->getNumChains(); i++)
	{
		cout << "testing residue connectivity across chains ... " << endl;
		tmpChain = pProt->getChain(0);
		tmpChain2 = pProt->getChain(i);
	
		for (UInt j = 1; j < tmpChain->getNumResidues() - 1; j ++)
		{
			for (UInt k = 1; k < tmpChain2->getNumResidues() -1; k++)
			{
				tmpRes = tmpChain->getResidue(j);
				tmpRes2 = tmpChain2->getResidue(k);

				if (tmpRes2->getNextRes() == tmpRes) cout << "error in sequential pointers found." << endl;
				if (tmpRes2->getPrevRes() == tmpRes) cout << "error in sequential pointers found - prev." << endl;
			}
		}
	}
 
	return 0;
}
