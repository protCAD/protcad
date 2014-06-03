#include <iostream>
#include <string>
#include "ensemble.h"
#include "PDBInterface.h"

vector <dblVec> getHelixAxis ( vector < vector < UIntVec > > , protein* _prot);

int main (int argc, char* argv[])
{

	double PI = 3.141596;

	string infile = argv[1];
	// read in protein
	PDBInterface* thePDB = new PDBInterface(infile);
	ensemble* theEnsemble = thePDB->getEnsemblePointer();
	molecule* pMol = theEnsemble->getMoleculePointer(0);
	protein* prot = static_cast<protein*>(pMol);

	// center on backbone centroid
	dblVec center = prot->getBackBoneCentroid();
	center = center * -1.0;
	prot->translate(center);

	vector < vector < UIntVec > > helixPosArray;

	//
	// DETERMINE HELIX BOUNDARIES
	//
	for (UInt i = 0; i < prot->getNumChains(); i ++)
	{
	    UIntVec tempHelixBeg(2);
	    tempHelixBeg[0] = i;  // initialize helix beginning to beginning of chain
   		tempHelixBeg[1] = 0;

	    for (UInt j = 1; j < prot->getNumResidues(i); j ++)
	    {
	        dblVec CAthis = prot->getCoords(i,j,"CA");
	        dblVec CAprev = prot->getCoords(i,j-1,"CA");
			dblVec diff = CAthis - CAprev;
	
	        double distance = sqrt(CMath::dotProduct(diff,diff));
	
			// if distance falls out of acceptable CA-CA range
	        if (!( distance > 3.65 && distance < 3.95 ) || (j + 1 == prot->getNumResidues(i)) ) 
	        {
	            if (j - tempHelixBeg[1] >= 10)   // if helix is at least 10 aa long
	            {
	                UIntVec tempHelixEnd(2);
	                tempHelixEnd[0] = i;
	                tempHelixEnd[1] = j;

	                vector < UIntVec > tempHelix;
	                tempHelix.push_back(tempHelixBeg);
	                tempHelix.push_back(tempHelixEnd);
	                helixPosArray.push_back(tempHelix);
	   	            if (j+1 < prot->getNumResidues(i))  // if we are not at end of the chain
	                {
	                    tempHelixBeg[0] = i;
	                    tempHelixBeg[1] = j+1;
	                }

	            }
	            else if ( j+1 < prot->getNumResidues(i)) // else if helix is not 10 aa long, reset beginning
	            {
	                tempHelixBeg[0] = i;
	                tempHelixBeg[1] = j+1;
	            }
	        }
	    }
	}

	dblVec origin(3);
	origin[0]=0.0;
	origin[1]=0.0;
	origin[2]=0.0;


	double bestProjection = 1;
	double bestPhi, bestTheta, bestPsi;

	double phimin = -180.0; double phimax = 180.0; double step = 20.0;
	double thetamin = -180.0; double thetamax = 180.0;
	double psimin = -180.0; double psimax = 180.0;

	// 
	// OPTIMIZE ROTATION PARAMETERS
	//
	for (UInt i = 0; i < 4; i ++)
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
					prot->eulerRotate(phi,theta,psi);
					vector <dblVec> vec = getHelixAxis(helixPosArray, prot);
					dblVec z_axis(3);
					z_axis[0] = 0.0;
					z_axis[1] = 0.0;
					z_axis[2] = 1.0;
					dblVec x_axis(3);
					x_axis[0] = 1.0;
					x_axis[1] = 0.0;
					x_axis[2] = 0.0;
					dblVec y_axis(3);
					y_axis[0] = 0.0;
					y_axis[1] = 1.0;
					y_axis[2] = 0.0;

					double projection = 0.0;
					for (UInt k = 0; k < vec.size(); k ++)
					{	
						double xproj = CMath::dotProduct(vec[k],x_axis);
						double yproj = CMath::dotProduct(vec[k],y_axis);
						double zproj = CMath::dotProduct(vec[k],z_axis);

						projection = fabs(zproj);// - fabs(xproj) - fabs(yproj);
					}
					prot->undoEulerRotate(phi, theta, psi);
					if (projection > bestProjection)
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
		phimin = bestPhi - 2*step;
		phimax = bestPhi + 2*step;
		thetamin = bestTheta - step;
		thetamax = bestTheta + step;
		psimin = bestPsi - step;
		psimax = bestPsi + step;
		step = step/10.0;
	}

	cout << "**************BEST***************" << endl;
	cout << "phi " << bestPhi << endl;
	cout << "psi " << bestPsi << endl;
	cout << "theta " << bestTheta << endl;
	prot->eulerRotate(bestPhi, bestTheta, bestPsi);
	pdbWriter(prot, "best.pdb");

	return 0;
}
//
// RETURN HELIX VECTORS BASED ON BOUNDARIES
//
vector <dblVec> getHelixAxis( vector < vector <UIntVec> > _helixPosArray, protein* _prot )
{
    dblVec origin(3);
    origin[0]=0.0;
    origin[1]=0.0;
    origin[2]=0.0;

    vector < dblVec > vecArray;
    for (UInt i = 0; i < _helixPosArray.size(); i ++)
    {
        vector <UIntVec> helix = _helixPosArray[i];
        UIntVec beg = helix[0];
        UIntVec end = helix[1];


/*		dblVec helixVec = origin;
		for (UInt i = beg[1]; i <= end[1]; i ++)
		{
			dblVec carbon = _prot->getCoords(beg[0], i, "C");
			dblVec oxygen = _prot->getCoords(beg[0], i, "O");
			
			dblVec COvec = oxygen-carbon;
			helixVec = helixVec + COvec;
		}
		
		vecArray.push_back(helixVec);
	}
*/
        dblVec begCentroid = origin;
        dblVec endCentroid = origin;
        for (UInt j = beg[1]; j < beg[1] + 7; j ++)
        {
            begCentroid = begCentroid +  _prot->getCoords(beg[0], j, "CA");
        }
        begCentroid = begCentroid/7.0;

        for (UInt j = end[1] - 6; j <= end[1]; j ++)
        {
            endCentroid = endCentroid +  _prot->getCoords(end[0], j, "CA");
        }
        endCentroid = endCentroid / 7.0;

        vecArray.push_back(endCentroid - begCentroid);
	}
	return vecArray;
}

