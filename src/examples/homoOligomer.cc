#include <iostream>
#include <string>
#include "ensemble.h"
#include "PDBInterface.h"
#include <fstream>
void buildParallelTrimer (protein* _prot, double _radius, double _phase, double _coil);
void undoParallelTrimer (protein* _prot, double _radius, double _phase, double _coil);
double calcEnergy(protein* _prot);
double getHBondEnergy(protein* _prot, UInt _chain1, UInt _res1, UInt _chain2, UInt _res2);
void mapSequence(protein* _prot);

int main (int argc, char* argv[])
{

//  STEP 1:  generate bundle
    if (argc != 4)
    {
        cout << "dimerBuilder  <infilename.pdb>  <outfilename.pdb> " << endl;
    exit(1);
    }
	enum aminoAcid {A, R, N, D, C, Q, E, G, H, I, L, K, M, F, P, S, T, W, Y, V};

    string infile = argv[1];
    PDBInterface* thePDB = new PDBInterface(infile);
    ensemble* theEnsemble = thePDB->getEnsemblePointer();
    molecule* pMol = theEnsemble->getMoleculePointer(0);
    protein* bundle = static_cast<protein*>(pMol);

	string basename = argv[2];
	string outfile = basename + ".pdb";
	string logfile = basename + ".log";

	ofstream fout;
	fout.open(logfile.c_str());
	double bestEnergy = 1e10;
	double energy = bestEnergy;

	residue::setCutoffDistance(6.0);
	pmf::setScaleFactor(0.0);
    rotamer::setScaleFactor(1.0);
    microEnvironment::setScaleFactor(0.0);
    amberVDW::setScaleFactor(1.0);
    amberVDW::setRadiusScaleFactor(0.95);
    amberVDW::setLinearRepulsionDampeningOn();
    amberElec::setScaleFactor(0.0);
    solvation::setItsScaleFactor(0.0);

	bundle->symmetryLinkChainAtoB(1,0);
	bundle->symmetryLinkChainAtoB(2,0);

    int resID[] = {A,Q,L,L,I,A,N,L,L,L,I,A,V,N,L,I,L,L,I,A,V,A,R,L,R,Y,L,V,G};

    int arraySize = sizeof(resID)/sizeof(resID[0]);

    for (int i = 0; i < arraySize; i ++)
    {
        bundle->activateForRepacking(0,i);
        bundle->setCanonicalHelixRotamersOnly(0,i);
        bundle->mutate(0, i, (UInt)resID[i]);
        cout << resID[i] << endl;
    }

	bundle->silenceMessages();
	for (double radius = 6.3; radius <= 6.7; radius = radius + 0.1)
	{
		for (double phase = 150.0; phase < 170.0; phase = phase + 2.0)
		{
			buildParallelTrimer(bundle, radius, phase, 190.0);
			mapSequence(bundle);
			energy = calcEnergy(bundle);
			if (energy < bestEnergy)
			{
				bestEnergy = energy;
				pdbWriter(bundle, outfile);
			}
			fout << " r " << radius << " p " << phase << " e " << energy << endl;
			cout << " r " << radius << " p " << phase << " e " << energy << endl;
 
			undoParallelTrimer(bundle, radius, phase, 190.0);
		}
	}
	return 0; 

}

void mapSequence (protein* _prot)
{

	enum aminoAcid {A, R, N, D, C, Q, E, G, H, I, L, K, M, F, P, S, T, W, Y, V};
    int _resID[] = {A,Q,L,L,I,A,N,L,L,L,I,A,V,N,L,I,L,L,I,A,V,A,R,L,R,Y,L,V,G};
    int arraySize = sizeof(_resID)/sizeof(_resID[0]);
	cout << arraySize << endl;

    for (int i = 0; i < arraySize; i ++)
    {
		_prot->mutate(0, i, 0);
        _prot->mutate(0, i, (UInt)_resID[i]);
    }

}

double calcEnergy(protein* _prot)
{
	_prot->optimizeRotamers();
	double itsEnergy = _prot->intraEnergy();
	for (UInt i = 0; i < _prot->getNumResidues(0); i ++)
	{
		double asn = getHBondEnergy(_prot, 0, 6, 1, i);
		asn += getHBondEnergy(_prot,0, 6, 2, i);
		itsEnergy += asn;
	}

/*

	double chi11, chi12, chi21, chi22, chi31, chi32;
	double bestEnergy = 0;
	_prot->saveCurrentState();	
	for (chi11 = -180.0; chi11 < 180.0; chi11 += 120.0)
	{
		cout << "*" << endl;
		for (chi12 = -180.0; chi12 < 180.0; chi12 += 120.0)
		{
			cout << "-" << endl;
			for (chi21 = -180.0; chi21 < 180.0; chi21 += 120.0)
			{
				for (chi22 = -180.0; chi22 < 180.0; chi22 += 120.0)
				{
					for (chi31 = -180.0; chi31 < 180.0; chi31 += 120.0)
					{
						for (chi32 = -180.0; chi32 < 180.0; chi32 += 120.0)
						{
							_prot->setChi(0,6, 0, 0, chi11);
							_prot->setChi(0,6, 0, 1, chi12);
							_prot->setChi(1,6, 0, 0, chi21);
							_prot->setChi(1,6, 0, 1, chi22);
							_prot->setChi(2,6, 0, 0, chi31);
							_prot->setChi(2,6, 0, 0, chi32);
							double energy = 0.0;
							for (UInt i = 0; i < _prot->getNumResidues(0); i ++)
							{
								double asn = getHBondEnergy(_prot, 0, 6, 1, i);
								asn += getHBondEnergy(_prot, 0, 6, 2, i);
								asn += getHBondEnergy(_prot, 1, 6, 2, i);
								energy += asn;
							}
							energy += _prot->getPositionEnergy(0,6) + _prot->getPositionEnergy(1,6) + _prot->getPositionEnergy(2,6);
							if (energy < bestEnergy)
							{
								bestEnergy = energy;
								_prot->commitState();
								_prot->saveCurrentState();
							}
						}
					}
				}
			}
		}
	}
	cout << endl;
	itsEnergy += bestEnergy;
	_prot->undoState();
*/
	return itsEnergy;
}


void buildParallelTrimer (protein* _prot, double _radius, double _phase, double _coil)
{
	_prot->rotate(Z_axis, _phase);
	_prot->translate(0.0, _radius, 0.0);
	_prot->rotate(1, Z_axis, 120);
	_prot->rotate(2, Z_axis, 240);
	_prot->coilcoil(_coil);
	return;
}

void undoParallelTrimer (protein* _prot, double _radius, double _phase, double _coil)
{
	_prot->coilcoil(-1.0 * _coil);
	_prot->rotate(1, Z_axis, -120);
	_prot->rotate(2, Z_axis, -240);
	_prot->translate(0.0, -1.0*_radius, 0.0);
	_prot->rotate(Z_axis, -1.0 * _phase);
	return;
}


double getHBondEnergy(protein* _prot, UInt _chain1, UInt _res1, UInt _chain2, UInt _res2)
{
	enum aminoAcid {A, R, N, D, C, Q, E, G, H, I, L, K, M, F, P, S, T, W, Y, V};
	double PI = 3.1415;

    UInt type1 = _prot->getTypeFromResNum(_chain1, _res1);
    UInt type2 = _prot->getTypeFromResNum(_chain2, _res2);

    vector <dblVec> donorList(0);
    vector <dblVec> donorBaseList(0);
    vector <dblVec> acceptorList(0);
    vector <dblVec> acceptorBaseList(0);
    vector <double> donorAngleList(0);
    vector <double> acceptorAngleList(0);
	

    for (UInt i = 1; i <= 2; i ++)
    {
        UInt chain, res, type;
        if (i == 1)
        {
            chain = _chain1;
            res = _res1;
            type = type1;
        }
        if (i == 2)
        {
            chain = _chain2;
            res = _res2;
            type = type2;
        }

        if (type == H)
        {
            dblVec donor = _prot->getCoords(chain,res,"NE2");
            dblVec temp1 = _prot->getCoords(chain,res,"CE1");
            dblVec temp2 = _prot->getCoords(chain,res,"CD2");
            dblVec donorBase = (temp1 + temp2)/2.0;
            dblVec acceptor = _prot->getCoords(chain,res,"ND1");
            temp2 = _prot->getCoords(chain,res,"CG");
            dblVec acceptorBase = (temp1 + temp2) / 2.0;
            donorList.push_back(donor);
            donorBaseList.push_back(donorBase);
            acceptorList.push_back(acceptor);
            acceptorBaseList.push_back(acceptorBase);
            donorAngleList.push_back(180.0);
            acceptorAngleList.push_back(180.0);
        }
      	if (type == T)
        {
            dblVec donor = _prot->getCoords(chain,res,"OG1");
            dblVec donorBase = _prot->getCoords(chain,res,"CB");
            dblVec acceptor = _prot->getCoords(chain,res,"OG1");
            dblVec acceptorBase = _prot->getCoords(chain,res,"CB");
            donorList.push_back(donor);
            donorBaseList.push_back(donorBase);
            acceptorList.push_back(acceptor);
            acceptorBaseList.push_back(acceptorBase);
            donorAngleList.push_back(109.0);
            acceptorAngleList.push_back(109.0);
        }
        if (type == S)
        {
            dblVec donor = _prot->getCoords(chain,res,"OG");
            dblVec donorBase = _prot->getCoords(chain,res,"CB");
            dblVec acceptor = _prot->getCoords(chain,res,"OG");
            dblVec acceptorBase = _prot->getCoords(chain,res,"CB");
            donorList.push_back(donor);
            donorBaseList.push_back(donorBase);
            acceptorList.push_back(acceptor);
            acceptorBaseList.push_back(acceptorBase);
            donorAngleList.push_back(109.0);
            acceptorAngleList.push_back(109.0);
        }
        if (type == Y)
        {
            dblVec donor = _prot->getCoords(chain,res,"OH");
            dblVec donorBase = _prot->getCoords(chain,res,"CZ");
            dblVec acceptor = _prot->getCoords(chain,res,"OH");
            dblVec acceptorBase = _prot->getCoords(chain,res,"CZ");
            donorList.push_back(donor);
            donorBaseList.push_back(donorBase);
            acceptorList.push_back(acceptor);
            acceptorBaseList.push_back(acceptorBase);
            donorAngleList.push_back(109.0);
            acceptorAngleList.push_back(109.0);
        }
        if (type == N)
        {
            dblVec donor = _prot->getCoords(chain,res,"ND2");
            dblVec donorBase = _prot->getCoords(chain,res,"CG");
            dblVec acceptor = _prot->getCoords(chain,res,"OD1");
            dblVec acceptorBase = _prot->getCoords(chain,res,"CG");
            donorList.push_back(donor);
            donorBaseList.push_back(donorBase);
            acceptorList.push_back(acceptor);
            acceptorBaseList.push_back(acceptorBase);
            donorAngleList.push_back(120.0);
            acceptorAngleList.push_back(180.0);
        }
        if (type == Q)
        {
            dblVec donor = _prot->getCoords(chain,res,"NE2");
            dblVec donorBase = _prot->getCoords(chain,res,"CD");
            dblVec acceptor = _prot->getCoords(chain,res,"OOE1");
            dblVec acceptorBase = _prot->getCoords(chain,res,"CD");
            donorList.push_back(donor);
            donorBaseList.push_back(donorBase);
            acceptorList.push_back(acceptor);
            acceptorBaseList.push_back(acceptorBase);
            donorAngleList.push_back(120.0);
            acceptorAngleList.push_back(180.0);
        }
        if (type == D)
        {
            dblVec acceptor1 = _prot->getCoords(chain,res,"OD1");
            dblVec acceptor2 = _prot->getCoords(chain,res,"OD2");
            dblVec acceptorBase = _prot->getCoords(chain,res,"CG");
            acceptorList.push_back(acceptor1);
            acceptorList.push_back(acceptor2);
            acceptorBaseList.push_back(acceptorBase);
            acceptorBaseList.push_back(acceptorBase);
            acceptorAngleList.push_back(180.0);
            acceptorAngleList.push_back(180.0);
        }
        if (type == E)
        {
            dblVec acceptor1 = _prot->getCoords(chain,res,"OE1");
            dblVec acceptor2 = _prot->getCoords(chain,res,"OE2");
            dblVec acceptorBase = _prot->getCoords(chain,res,"CD");
            acceptorList.push_back(acceptor1);
            acceptorList.push_back(acceptor2);
            acceptorBaseList.push_back(acceptorBase);
            acceptorBaseList.push_back(acceptorBase);
            acceptorAngleList.push_back(180.0);
            acceptorAngleList.push_back(180.0);
        }
        if (type == W)
        {
            dblVec donor = _prot->getCoords(chain,res,"NE1");
            dblVec temp1 = _prot->getCoords(chain,res,"CD1");
            dblVec temp2 = _prot->getCoords(chain,res,"CE2");
            dblVec donorBase = (temp1 + temp2)/2.0;
            donorList.push_back(donor);
            donorBaseList.push_back(donorBase);
            donorAngleList.push_back(180.0);
        }
		dblVec carbonyl = _prot->getCoords(chain,res,"O");
		dblVec carbonylC = _prot->getCoords(chain,res,"C");
		acceptorList.push_back(carbonyl);
		acceptorBaseList.push_back(carbonylC);
		acceptorAngleList.push_back(120);
    }
    double energy = 0.0;
	if (donorList.size() == 0) return 0.0;

    for (UInt i = 0; i < donorList.size(); i ++)
    {
        for (UInt j = 0; j < acceptorList.size(); j ++)
        {
            dblVec DDB = donorList[i] - donorBaseList[i];
            dblVec DA = donorList[i] - acceptorList[j];
            dblVec AAB = acceptorList[j] - acceptorBaseList[j];

            double magDDB = sqrt(CMath::dotProduct(DDB,DDB));
            double magDA = sqrt(CMath::dotProduct(DA,DA));
            double magAAB = sqrt(CMath::dotProduct(AAB,AAB));

			double thisE;
			if (magDA < 1e-50) thisE = 0.0;            
			else 
			{
				
				double donorAngle = acos( CMath::dotProduct(DA,DDB)/(magDDB*magDA) );
            	double acceptorAngle = acos( CMath::dotProduct(DA,AAB)/(magDA*magAAB) );

            	double distRatio = 2.8 / magDA;
            	thisE = 2.0 *(5.0*pow(distRatio,12.0)-6.0*pow(distRatio,10.0));
            	double angleFactor1 = cos(donorAngleList[i]*PI/180.0-donorAngle);
            	double angleFactor2 = cos(acceptorAngleList[j]*PI/180.0-acceptorAngle);
            	thisE = thisE * pow(angleFactor1, 2.0) * pow (angleFactor2, 2.0);
			}
            energy += thisE;
        }
    }
    return energy;
}

