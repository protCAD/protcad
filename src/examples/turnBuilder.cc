#include "typedef.h"
#include "ensemble.h"
#include "PDBInterface.h"
#include "ran.h"

bool exceedRange(protein* _prot, UInt _A, UInt _B, double _cutoff);

int main (int argc, char* argv[])
{
	enum aminoAcid {A, R, N, D, C, Q, E, G, H, I, L, K, M, F, P, S, T, W, Y, V};

	if (argc < 2)
	{
		cout << "turnBuilder input.pdb" << endl;
		exit(1);
	}

	string inputFile = argv[1];
	ifstream inFile;
	inFile.open(inputFile.c_str());
	if (!inFile)
	{
		cout << "Unable to find or open file" << endl;
		exit(1);
	}

	PDBInterface* thePDB = new PDBInterface(inputFile);
	ensemble* theEnsemble = thePDB->getEnsemblePointer();
	molecule* theMol = theEnsemble->getMoleculePointer(0);
	protein* prot = static_cast<protein*>(theMol);

	residue::setCutoffDistance(4.0);
	pmf::setScaleFactor(0.0);
	rotamer::setScaleFactor(0.0);
	microEnvironment::setScaleFactor(0.0);
	amberVDW::setScaleFactor(1.0);
	amberVDW::setRadiusScaleFactor(0.95);
	amberVDW::setLinearRepulsionDampeningOff();
	amberElec::setScaleFactor(0.0);
	solvation::setItsScaleFactor(0.0);

	UInt startPos = 25;

	prot->activateForRepacking(0, startPos - 4);
	prot->mutate(0, startPos - 4, L);
	prot->activateForRepacking(0, startPos - 8);
	prot->mutate(0, startPos - 8, L);
	prot->activateForRepacking(0, startPos + 8);
	prot->mutate(0, startPos + 8, L);
	prot->activateForRepacking(0, startPos + 1);
	prot->mutate(0, startPos + 1, V);
	
	prot->setCanonicalHelixRotamersOnly(0);
	prot->optimizeRotamers();

	// aL
	double phi1min = 72.0 - 16.0;
	double phi1max = 72.0 + 16.0;
	double psi1min = 21.0 - 19.0;
	double psi1max = 21.0 + 19.0;
	
	// four residues of polyproline
	double phi2min = -111.0 - 27.0;
	double phi2max = -111.0 + 27.0;
	double psi2min = 141.0 - 28.0;
	double psi2max = 141.0 + 28.0;

	double phi3min = -86.0 - 23.0;
	double phi3max = -86.0 + 23.0;
	double psi3min = 129.0 - 22.0;
	double psi3max = 129.0 + 22.0;

	double phi4min = -93.0 - 27.0;
	double phi4max = -93.0 + 27.0;
	double psi4min = 135.0 - 23.0;
	double psi4max = 135.0 + 23.0;

	double phi5min = -84.0 - 25.0;
	double phi5max = -84.0 + 25.0;
	double psi5min = 154.0 - 19.0;
	double psi5max = 154.0 + 19.0;
	
	prot->setPhi(0, startPos - 1, -75.0);  // make C1 position gamma
	prot->setPsi(0, startPos - 1, 10.0);
	prot->silenceMessages();
	double steps = 10.0;
	double lowEnergy = 1e20;
	for (double phi1 = phi1min; phi1 <= phi1max; phi1 = phi1 + (phi1max - phi1min)/steps) {
		prot->setPhi(0, startPos, phi1);
		for (double psi1 = psi1min; psi1 <= psi1max; psi1 = psi1 + (psi1max - psi1min)/steps) {
			prot->setPsi(0, startPos, psi1);
			for (double phi2 = phi2min; phi2 <= phi2max; phi2 = phi2 + (phi2max - phi2min)/steps) {
				prot->setPhi(0, startPos + 1, phi2);
				for (double psi2 = psi2min; psi2 <= psi2max; psi2 = psi2 + (psi2max - psi2min)/steps) {
					prot->setPsi(0, startPos + 1, psi2);
					for (double phi3 = phi3min; phi3 <= phi3max; phi3 = phi3 + (phi3max - phi3min)/steps) {
						prot->setPhi(0, startPos + 2, phi3);
						for (double psi3 = psi3min; psi3 <= psi3max; psi3 = psi3 + (psi3max - psi3min)/steps) {
							prot->setPsi(0, startPos + 2, psi3);
							for (double phi4 = phi4min; phi4 <= phi4max; phi4 = phi4 + (phi4max - phi4min)/steps) {
								prot->setPhi(0, startPos + 3, phi4);
								for (double psi4 = psi4min; psi4 <= psi4max; psi4 = psi4 + (psi4max - psi4min)/steps) {
									prot->setPsi(0, startPos + 3, psi4);
									for (double phi5 = phi5min; phi5 <= phi5max; phi5 = phi5 + (phi5max - phi5min)/steps) {
										prot->setPhi(0, startPos + 4, phi5);
										for (double psi5 = psi5min; psi5 <= psi5max; psi5 = psi5 + (psi5max - psi5min)/steps) {
											prot->setPsi(0, startPos + 4, psi5);
											bool inrange = true;
											//if (exceedRange(prot, startPos - 7, startPos + 3, 15.0)) inrange = false;
											//if (exceedRange(prot, startPos - 8, startPos + 8, 10.0)) inrange = false;
											if (exceedRange(prot, startPos - 7, startPos + 8, 10.0)) inrange = false;
											if (inrange) {
												double energy = prot->intraEnergy();
												if (energy < lowEnergy) {
													cout << energy << " " << phi1 << " " << psi1 << " " << phi2  << " " <<  psi2 << " " 
														<< phi3 << " " << psi3 << " " << phi4 << " " << psi4 << " " << phi5 << " " << psi5 << endl;
													lowEnergy = energy;
													pdbWriter(prot, "best.pdb");
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}

	return 0;
}

bool exceedRange(protein* _prot, UInt _A, UInt _B, double _cutoff)
{
	dblVec A = _prot->getCoords(0, _A, "CA");
	dblVec B = _prot->getCoords(0, _B, "CA");

	double distance = sqrt(CMath::dotProduct(A,B));
	if (distance > _cutoff) {
		return true;
	}
	return false;
}
