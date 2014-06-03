#include <iostream>
#include <string>
#include "ensemble.h"
#include "PDBInterface.h"
#include <fstream>
void buildParallelDimer (protein* _prot, double _radius, double _phase, double _coil);
void undoParallelDimer (protein* _prot, double _radius, double _phase, double _coil);
double calcEnergy(protein* _prot);
double getHBondEnergy(protein* _prot, UInt _chain1, UInt _res1, UInt _chain2, UInt _res2);
void mapSequence(protein* _prot);

int main (int argc, char* argv[])
{

    if (argc != 4)
    {
        cout << "serDimer  <infilename.pdb>  <outfilename.pdb> " << endl;
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
    amberVDW::setRadiusScaleFactor(1.0);
    amberVDW::setLinearRepulsionDampeningOn();
    amberElec::setScaleFactor(0.0);
    solvation::setItsScaleFactor(0.0);

	bundle->symmetryLinkChainAtoB(1,0);

/*	dblVec aCB = bundle->getCoords(0, 5, "CB");
	dblVec dCB = bundle->getCoords(0, 8, "CB");

	aCB[2] = 0.0;
	dCB[2] = 0.0;

	dblVec center = (aCB + dCB) * 0.5;
	
	dblVec tempY(3);
	tempY[0] = 1.0;
	tempY[1] = 0.0;
	tempY[2] = 0.0;

	double magCB = sqrt(CMath::dotProduct(center,center));
	double magY = sqrt(CMath::dotProduct(tempY,tempY));

	double angle = acos(CMath::dotProduct(center,tempY)/(magY*magCB));
	bundle->rotate(Z_axis, -1.0*angle);
	pdbWriter(bundle, "start.pdb");
*/	
    int resID[] = {A,A,A,A,A, V,L,L,L,I,A,V, V,L,I,L,L,I,A, V,A,R,L,R,Y,L, A,A,A,A,A};
    int arraySize = sizeof(resID)/sizeof(resID[0]);

    for (int i = 0; i < arraySize; i ++)
    {
        bundle->activateForRepacking(0,i);
        bundle->setCanonicalHelixRotamersOnly(0,i);
        bundle->mutate(0, i, (UInt)resID[i]);
        cout << resID[i] << endl;
    }

	UInt apos;
	sscanf(argv[3], "%u", &apos);
	bundle->mutate(0, 5, apos);
	bundle->mutate(0,12, apos);
	bundle->mutate(0,19, apos);
    
	bundle->silenceMessages();
	for (double radius = 4.5; radius <= 5.5; radius = radius + 0.25)
	{
		for (double phase = 60.0; phase < 180.0; phase = phase + 5.0)
		{
			buildParallelDimer(bundle, radius, phase, 190.0);
			mapSequence(bundle);
			bundle->optimizeRotamers();
			energy = bundle->intraEnergy();
			if (energy < 0)
			{
				double besthbenergy = 1E50;
				double bestchirot1 = 0.0;
				double chirot1 = 0.0;
				double bestchirot2 = 0.0;
				double chirot2 = 0.0;
				for (UInt rot1 = 0; rot1 <= 2; rot1 ++)
				{
					for (UInt rot2 = 0; rot2 <= 2; rot2 ++)
					{
						for (double dchirot1 = -20.0; dchirot1 <- 20.0; dchirot1 += 20.0)
						{
							for (double dchirot2 = -20.0; dchirot2 <- 20.0; dchirot2 += 20.0)
							{
								double itsE = 0.0;
								for (UInt pos = 0; pos < bundle->getNumResidues(0); pos ++)
								{
									UInt id = bundle->getTypeFromResNum(0,pos);
									if (id == S || id == T)
									{
										bundle->setRotamer(0,pos,0, rot1);
										chirot1 = bundle->getChi(0,pos,0,0) + dchirot1;
										bundle->setRotamer(1,pos,0, rot2);
										chirot2 = bundle->getChi(1,pos,0,0) + dchirot2;
										bundle->setRelativeChi(0,pos,0,0, chirot1);	
										bundle->setRelativeChi(1,pos,0,0, chirot2);
						 			}
								}
								for (UInt pos = 0; pos < bundle->getNumResidues(0); pos ++)
								{
									UInt id = bundle->getTypeFromResNum(0,pos);
									if (id == S || id == T)
									{
										itsE += bundle->getPositionEnergy(0,pos);
										itsE += bundle->getPositionEnergy(1,pos);
									}
								}
								itsE += calcEnergy(bundle);
								if (itsE < besthbenergy)
								{
									besthbenergy = itsE;
									bestchirot1 = chirot1;
									bestchirot2 = chirot2;
								}	
							}
						}
					}
				}
				for (UInt pos = 0; pos < bundle->getNumResidues(0); pos ++)
				{
					UInt id = bundle->getTypeFromResNum(0,pos);
					if (id == S || id == T)
					{										
						bundle->setRelativeChi(0,pos,0,0, bestchirot1);
						bundle->setRelativeChi(1,pos,0,0, bestchirot2);
					}
				}
			
				double hbenergy = calcEnergy(bundle);
				energy = bundle->intraEnergy() +  hbenergy;
				if (energy < bestEnergy)
				{
					bestEnergy = energy;
					pdbWriter(bundle, outfile);
				}
				
				fout << " r " << radius << " p " << phase << " e " << energy << " hb " << hbenergy << endl;
				cout << " r " << radius << " p " << phase << " e " << energy << " hb " << hbenergy << endl;
 
			}
			undoParallelDimer(bundle, radius, phase, 190.0);
		}
	}
	return 0; 

}

void mapSequence (protein* _prot)
{
	UInt arraySize = _prot->getNumResidues(0);
    for (UInt i = 0; i < arraySize; i ++)
    {
	UInt resID = _prot->getTypeFromResNum(0,i);
    	_prot->mutate(0, i, 0);
        _prot->mutate(0, i, resID);
    }

}

double calcEnergy(protein* _prot)
{
	double itsEnergy = 0.0;
	for (UInt i = 0; i < _prot->getNumResidues(0); i ++)
	{
		
		for (UInt j = 0; j < _prot->getNumResidues(1); j ++)
		{
			double hbe = getHBondEnergy(_prot, 0, i, 1, j);
			itsEnergy += hbe;
		}
	}

	return itsEnergy;
}


void buildParallelDimer (protein* _prot, double _radius, double _phase, double _coil)
{
	_prot->rotate(Z_axis, _phase);
	_prot->translate(0.0, _radius, 0.0);
	_prot->rotate(1, Z_axis, 180);
	_prot->coilcoil(_coil);
	return;
}

void undoParallelDimer (protein* _prot, double _radius, double _phase, double _coil)
{
	_prot->coilcoil(-1.0 * _coil);
	_prot->rotate(1, Z_axis, 180);
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
            	thisE = 10.0 *(5.0*pow(distRatio,12.0)-6.0*pow(distRatio,10.0));
            	double angleFactor1 = cos(donorAngleList[i]*PI/180.0-donorAngle);
            	double angleFactor2 = cos(acceptorAngleList[j]*PI/180.0-acceptorAngle);
            	thisE = thisE * pow(angleFactor1, 2.0) * pow (angleFactor2, 2.0);
			}
            energy += thisE;
        }
    }
    return energy;
}

