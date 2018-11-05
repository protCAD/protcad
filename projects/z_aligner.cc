#include <iostream>
#include <string>
#include "ensemble.h"
#include "PDBInterface.h"

#define PI 3.14159265359

int main (int argc, char* argv[])
{

    string infile = argv[1];
    // read in protein
    PDBInterface* thePDB = new PDBInterface(infile);
    ensemble* theEnsemble = thePDB->getEnsemblePointer();
    molecule* pMol = theEnsemble->getMoleculePointer(0);
    protein* _prot = static_cast<protein*>(pMol);

    // center on backbone centroid
    dblVec center = _prot->getBackBoneCentroid();
    center = center * -1.0;
    _prot->translate(center);
    double bestProjection = 100000;
    double bestPhi, bestTheta, bestPsi;

    double phimin = -2.0 * PI; double phimax = 2.0 * PI; double step = 0.1 * 2.0 * PI;
    double thetamin = -2.0 * PI; double thetamax = 2.0 * PI;
    double psimin = -2.0 * PI; double psimax = 2.0 * PI;
    for (UInt i = 0; i < 20; i ++)
    {
        cout << endl;
        cout << "*********" << endl;
        cout << " CYCLE " << i+1 << endl;
        cout << endl;

        for (double phi = phimin; phi <= phimax; phi += step)
        {
            for (double theta = thetamin; theta <= thetamax; theta += step)
            {
                for (double psi = psimin; psi <= psimax; psi += step)
                {
                    protein* prot = new protein(*_prot);
                    prot->eulerRotate(phi, theta, psi);
                    double projection = 0.0;
                    for (UInt i = 0; i < prot->getNumChains(); i ++)
                    {
                        chain* tempChain = prot->getChain(i);
                        dblVec chainCentroid = tempChain->getBackBoneCentroid();
                        dblVec chainSum(3);
                        chainSum[0] = 0.0; chainSum[1] = 0.0; chainSum[2] = 0.0;
                        for (UInt j = 0; j < tempChain->itsResidues.size(); j++)
                        {
                            residue* tempRes = tempChain->getResidue(j);
                            chainSum = chainSum + (chainCentroid - tempRes->getBackBoneCentroid());
                            chainSum[2] = 0.0;
                            double magnitude = sqrt(CMath::dotProduct(chainSum,chainSum));
                            projection += magnitude;
                        }
                    }
                    delete prot;
                    if (projection < bestProjection)
                    {
                        bestProjection = projection;
                        bestPhi = phi;
                        bestTheta = theta;
                        bestPsi = psi;
                    }
                }
            }
        }
        cout << "phi " << bestPhi << endl;
        cout << "psi " << bestPsi << endl;
        cout << "theta " << bestTheta << endl;
        cout << "best projection " << bestProjection << endl;
        cout << endl;
        phimin = bestPhi - step;
        phimax = bestPhi + step;
        thetamin = bestTheta - step;
        thetamax = bestTheta + step;
        psimin = bestPsi - step;
        psimax = bestPsi + step;
        step = step/5;
    }

    cout << "**************BEST***************" << endl;
    cout << "phi " << bestPhi << endl;
    cout << "psi " << bestPsi << endl;
    cout << "theta " << bestTheta << endl;
    cout << "projection " << bestProjection << endl;
    _prot->eulerRotate(bestPhi, bestTheta, bestPsi);
    pdbWriter(_prot, infile);
    return 0;
}





