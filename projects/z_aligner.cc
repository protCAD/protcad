#include <iostream>
#include <string>
#include "ensemble.h"
#include "PDBInterface.h"

#define PI 3.14159
dblVec SFCentroidChain(protein* _prot, UInt _chain);
dblVec SFCentroidRes(protein* _prot, UInt _chain, UInt _residue);

int main (int argc, char* argv[])
{

	string infile = argv[1];
	// read in protein
	PDBInterface* thePDB = new PDBInterface(infile);
	ensemble* theEnsemble = thePDB->getEnsemblePointer();
	molecule* pMol = theEnsemble->getMoleculePointer(0);
	protein* prot = static_cast<protein*>(pMol);

	// center on backbone centroid
    dblVec center = SFCentroidChain(prot, 0);
	center = center * -1.0;
	prot->translate(center);
	double bestProjection = 100000;
    double bestPhi = 0.0, bestTheta = 0.0, bestPsi = 0.0;

	double phimin = -2.0 * PI; double phimax = 2.0 * PI; double step = (20.0/180.0) * 2.0 * PI;
	double thetamin = -2.0 * PI; double thetamax = 2.0 * PI;
	double psimin = -2.0 * PI; double psimax = 2.0 * PI;
    for (UInt i = 0; i < 10; i ++)
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
                        dblVec chainCentroid = SFCentroidChain(prot,i);
						dblVec chainSum(3);
						chainSum[0] = 0.0; chainSum[1] = 0.0; chainSum[2] = 0.0;
						for (UInt j = 0; j < tempChain->itsResidues.size(); j++)
						{
                            chainSum = chainSum + (chainCentroid - SFCentroidRes(prot,i,j));
							chainSum[2] = 0.0;
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
	pdbWriter(prot, infile);
	return 0;
}
dblVec SFCentroidChain(protein* _prot, UInt _chain)
{
    //--initialize and clear variables
    double number = 0;
    UInt numRes, numAtoms;
    string atomType;
    dblVec coords, coordsSum(3), coordsAve(3);
    coordsSum[0] = 0.0, coordsSum[1] = 0.0, coordsSum[2] = 0.0;
    coordsAve[0] = 0.0, coordsAve[1] = 0.0, coordsAve[2] = 0.0;

    //--loop through all atoms of all residues in search of carbon
    numRes = _prot->getNumResidues(_chain);
    for (UInt i = 0; i < numRes; i++)
    {
        numAtoms = _prot->getNumAtoms(_chain, i);
        for (UInt j = 0; j < numAtoms; j++)
        {
            atomType = _prot->getTypeStringFromAtomNum(_chain, i, j);
            if (atomType != "C" && atomType != "CA" && atomType != "O" && atomType != "CB" && atomType != "N" && atomType != "SG")
            {
                number++;
                coords = _prot->getCoords(_chain, i, j);
                coordsSum[0] += coords[0];
                coordsSum[1] += coords[1];
                coordsSum[2] += coords[2];
            }
        }
    }

    //--get average of all carbon coordinates
    coordsAve[0] = ((coordsSum[0])/number);
    coordsAve[1] = ((coordsSum[1])/number);
    coordsAve[2] = ((coordsSum[2])/number);

    return coordsAve;
}

dblVec SFCentroidRes(protein* _prot, UInt _chain, UInt _residue)
{
    //--initialize and clear variables
    double number = 0;
    UInt numAtoms;
    string atomType;
    dblVec coords, coordsSum(3), coordsAve(3);
    coordsSum[0] = 0.0, coordsSum[1] = 0.0, coordsSum[2] = 0.0;
    coordsAve[0] = 0.0, coordsAve[1] = 0.0, coordsAve[2] = 0.0;

    //--loop through all atoms of all residues in search of carbon

    numAtoms = _prot->getNumAtoms(_chain, _residue);
    for (UInt j = 0; j < numAtoms; j++)
    {
        atomType = _prot->getTypeStringFromAtomNum(_chain, _residue, j);
        if (atomType != "C" && atomType != "CA" && atomType != "O" && atomType != "CB" && atomType != "N" && atomType != "SG")
        {
            number++;
            coords = _prot->getCoords(_chain, _residue, j);
            coordsSum[0] += coords[0];
            coordsSum[1] += coords[1];
            coordsSum[2] += coords[2];
        }
    }

    //--get average of all carbon coordinates
    coordsAve[0] = ((coordsSum[0])/number);
    coordsAve[1] = ((coordsSum[1])/number);
    coordsAve[2] = ((coordsSum[2])/number);

    return coordsAve;
}



