#include <iostream>
#include <string>
#include "ensemble.h"
#include "PDBInterface.h"

void y_align(protein* prot);

int main (int argc, char* argv[])
{

    double PI = 3.141596;
	enum aminoAcid {A, R, N, D, C, Q, E, G, H, I, L, K, M, F, P, S, T, W, Y, V};

    string infile = argv[1];
    // read in protein
    PDBInterface* thePDB = new PDBInterface(infile);
    ensemble* theEnsemble = thePDB->getEnsemblePointer();
    molecule* pMol = theEnsemble->getMoleculePointer(0);
    protein* prot = static_cast<protein*>(pMol);


	// align on y_axis and center at origin
	y_align(prot);


	// align positon zero on -z axis

	dblVec cBeta = prot->getCoords(0,0,"CB");
	cBeta[1] = 0.0;
	
	dblVec negZ(3);
	negZ[0] =  0.0;
	negZ[1] =  0.0;
	negZ[2] = -1.0;

	double magCB = sqrt(CMath::dotProduct(cBeta,cBeta));

	double angle = (180.0/PI) * acos( CMath::dotProduct(cBeta,negZ) / magCB );
	prot->rotate(0, Y_axis,  angle);
	pdbWriter(prot, "temp1.pdb");

	double phase;
	double period;
	sscanf(argv[2], "%lf", &phase);
	sscanf(argv[3], "%lf", &period);

	phase = phase * 360.0 / period;

	prot->rotate(0, Y_axis, phase);

	prot->activateAllForRepacking(0);


	for (UInt i = 6; i <= 12; i ++)
	{
		prot->mutate(0,i,Y);
		prot->setRotamer(0,i,0,3);
		if (i == 6) pdbWriter(prot,"rot3.pdb");
		dblVec cGamma = prot->getCoords(0,i,"CG");
		prot->translate(-1.0*cGamma);
		dblVec coordsRot1 = prot->getCoords(0,i,"OH");
		prot->translate(cGamma);
		prot->setRotamer(0,i,0,6);
		if (i == 6) pdbWriter(prot,"rot6.pdb");
		cGamma = prot->getCoords(0,i,"CG");
		prot->translate(-1.0*cGamma);
		dblVec coordsRot2 = prot->getCoords(0,i,"OH");
		prot->translate(cGamma);

		
		dblVec posY(3), posZ(3);
		for (UInt n = 0; n < 3; n ++)
		{
			coordsRot1[n] = fabs(coordsRot1[n]);
			coordsRot2[n] = fabs(coordsRot2[n]);
			posY[n] = 0.0; posZ[n] = 0.0;
		}
		posY[1] = 1.0; posZ[2] = 1.0;

		cout << "Pos " << i-1 << " x " << coordsRot1[0] << " y " << coordsRot1[1] << " z " << coordsRot1[2] << endl;
		cout << "Pos " << i-1 << " x " << coordsRot2[0] << " y " << coordsRot2[1] << " z " << coordsRot2[2] << endl;

		double rot1mag = sqrt(CMath::dotProduct(coordsRot1,coordsRot1));
		double rot2mag = sqrt(CMath::dotProduct(coordsRot2,coordsRot2));

		double helixAngle1 = (180/PI)*acos(CMath::dotProduct(posY,coordsRot1)/rot1mag);
		double helixAngle2 = (180/PI)*acos(CMath::dotProduct(posY,coordsRot2)/rot2mag);

		double memnormAngle1 = (180/PI)*acos(CMath::dotProduct(posZ,coordsRot1)/rot1mag);
		double memnormAngle2 = (180/PI)*acos(CMath::dotProduct(posZ,coordsRot2)/rot2mag);
		
		cout << " t  g+ :  helix = " << helixAngle1 << " membrane normal " << memnormAngle1 ;
		cout << " g- g+ :  helix = " << helixAngle2 << " membrane normal " << memnormAngle2 << endl;

		cout << endl;
	}

	pdbWriter(prot, "end.pdb");

	return 0;

}

void y_align(protein* prot)
{
    dblVec center = prot->getBackBoneCentroid();
    center = center * -1.0;
    prot->translate(center);
    double bestProjection = 100000;
    double bestPhi, bestTheta, bestPsi;

    double phimin = -180.0; double phimax = 180.0; double step = 20.0;
    double thetamin = -180.0; double thetamax = 180.0;
    double psimin = -180.0; double psimax = 180.0;
    for (UInt i = 0; i < 5; i ++)
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
                            chainSum[1] = 0.0;
                            double magnitude = sqrt(CMath::dotProduct(chainSum,chainSum));
                            projection += magnitude;
                        }
                    }
                    prot->undoEulerRotate(phi, theta, psi);
                    //cout << phi << " " << theta << " " << psi << " " << projection << endl;
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
    prot->eulerRotate(bestPhi, bestTheta, bestPsi);
}


