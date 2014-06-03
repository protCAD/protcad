// from Numerical Re/cipies in C++  p415

#include <vector>
#include <iostream>
#include "ensemble.h"
#include <string>
#define FUNK funk
#define NMAX 5000
void get_psum(dblMat &p, dblVec &psum);
void swap(double &x, double &y);
void ameoba(dblMat &p, dblVec &y, const double ftol, double FUNK(dblVec &), int &nfunk);
double amotry(dblMat &p, dblVec &y, dblVec &psum, double FUNK(dblVec &), const int ihi, const double fac);
double funk(dblVec &params);
double disruptFunk(dblVec &params);
double silentFunk(dblVec &params);
double oldFunk(dblVec &params);

static vector < double > wtE;
static vector < vector <double> > favEbs;
static vector < vector <double> > favEsb;
static vector < vector <double> > unfavEbs;
static vector < vector <double> > unfavEsb;
static UInt wtStruct;
static UInt numStructs;
int main(int argc, char* argv[])
{


	// read parameters in for downhill simplex method
    string paramFileName = argv[1];
    ifstream paramFile;
    paramFile.open(paramFileName.c_str());
    if (!paramFile)
    {
        cout << "Unable to find or open file" << endl;
        exit(1);
    }

	string currentLine;
	vector <string> parsedStrings;

	// read in initial guess
    getline(paramFile, currentLine, '\n'); parsedStrings = Parse::parse(currentLine);
	UInt size = parsedStrings.size();
	dblVec firstGuess(size);
	if (argc <= 3)
	{
		for (UInt i=0; i < size; i ++)
		{
			sscanf(parsedStrings[i].c_str(), "%lf", &firstGuess[i]);
		}
	}
	else
	for (UInt i = 0; i < size; i ++)
	{
		cout << "user param "  << argv[i + 3] << endl;
		sscanf(argv[i+3], "%lf", &firstGuess[i]);
	}


	// read in length scales for each parameter
    getline(paramFile, currentLine, '\n'); parsedStrings = Parse::parse(currentLine);
	dblVec lambda(size);
	for (UInt i=0; i < size; i ++)
	{
		sscanf(parsedStrings[i].c_str(), "%lf", &lambda[i]);
	}

	// build simplex matrix
	dblMat p(size+1, size, 0.0);
	for (UInt i = 0; i < size; i ++)
	{
		 p[0][i] = firstGuess[i];
		 p[i+1][i] = lambda[i];
	}
	double ftol;
	getline(paramFile, currentLine, '\n'); parsedStrings = Parse::parse(currentLine);
	sscanf(parsedStrings[0].c_str(), "%lf", &ftol);

	// read in energies
	string dataFileName = argv[2];
	ifstream dataFile;
	dataFile.open(dataFileName.c_str());
	if (!dataFile)
	{
		cout << "Unable to find or open file" << endl;
		exit(1);
	}

	// positions of information (1) wt (2) favorable start/end (3) unfavorable start/end
	UInt wtPos;
	getline(dataFile, currentLine, '\n'); parsedStrings = Parse::parse(currentLine);
	sscanf(parsedStrings[0].c_str(), "%u", &wtPos); 	
	UIntVec favBSpos(0);
	getline(dataFile, currentLine, '\n'); parsedStrings = Parse::parse(currentLine);
	for (UInt i = 0; i < parsedStrings.size(); i ++)
	{
		UInt pos;
		sscanf(parsedStrings[i].c_str(), "%u", &pos);
		favBSpos.push_back(pos);
	}
	UIntVec favSBpos(0);
	getline(dataFile, currentLine, '\n'); parsedStrings = Parse::parse(currentLine);
	for (UInt i = 0; i < parsedStrings.size(); i ++)
	{
		UInt pos;
		sscanf(parsedStrings[i].c_str(), "%u", &pos);
		favSBpos.push_back(pos);
	}
	UIntVec unfavBSpos(0);
	getline(dataFile, currentLine, '\n'); parsedStrings = Parse::parse(currentLine);
	for (UInt i = 0; i < parsedStrings.size(); i ++)
	{
		UInt pos;
		sscanf(parsedStrings[i].c_str(), "%u", &pos);
		unfavBSpos.push_back(pos);
	}
	UIntVec unfavSBpos(0);
	getline(dataFile, currentLine, '\n'); parsedStrings = Parse::parse(currentLine);
	for (UInt i = 0; i < parsedStrings.size(); i ++)
	{
		UInt pos;
		sscanf(parsedStrings[i].c_str(), "%u", &pos);
		unfavSBpos.push_back(pos);
	}


	UInt x = 0;
	while (getline(dataFile, currentLine, '\n'))
	{
		x ++;
		bool lineGood = true;
		double wtEnergy;
		parsedStrings = Parse::parse(currentLine);
		sscanf(parsedStrings[wtPos].c_str(), "%lf", &wtEnergy);
		if (!(wtEnergy > -1e50 && wtEnergy < 1e50)) lineGood = false; 

		double favEnergyBS;
		vector <double> favEnergyBSArray;
		for (UInt i = 0; i < favBSpos.size(); i ++)
		{
			sscanf(parsedStrings[i].c_str(), "%lf", &favEnergyBS);
			favEnergyBSArray.push_back(favEnergyBS);
			if  (!(favEnergyBS > -1e50 && favEnergyBS < 1e50)) lineGood = false;
		}

		double unfavEnergyBS;
		vector <double> unfavEnergyBSArray;
		for (UInt i = 0; i < unfavBSpos.size(); i ++)
		{
			sscanf(parsedStrings[i].c_str(), "%lf", &unfavEnergyBS);
			unfavEnergyBSArray.push_back(unfavEnergyBS);
			if  (!(unfavEnergyBS > -1e50 && unfavEnergyBS < 1e50)) lineGood = false;
		}

		double favEnergySB;
		vector <double> favEnergySBArray;
		for (UInt i = 0; i < favSBpos.size(); i ++)
		{
			sscanf(parsedStrings[i].c_str(), "%lf", &favEnergySB);
			favEnergySBArray.push_back(favEnergySB);
			if  (!(favEnergySB > -1e50 && favEnergySB < 1e50)) lineGood = false;
		}

		double unfavEnergySB;
		vector <double> unfavEnergySBArray;
		for (UInt i = 0; i < unfavSBpos.size(); i ++)
		{
			sscanf(parsedStrings[i].c_str(), "%lf", &unfavEnergySB);
			unfavEnergySBArray.push_back(unfavEnergySB);
			if  (!(unfavEnergySB > -1e50 && unfavEnergySB < 1e50)) lineGood = false;
		}

		if (lineGood)
		{
			unfavEbs.push_back(unfavEnergyBSArray);
			unfavEsb.push_Back(unfavEnergySBArray);
			wtE.push_back(wtEnergy);
			favEsb.push_back(favEnergySBArray);
			favEbs.push_back(favEnergyBSArray);

		}
		else cout << "WARNING ... LINE " << x << " CONTAINS BAD DATA." << endl;
	}

	cout << "number of conformations sampled:  " << wtE.size() << endl;

	// calculate y - vector of funk evaluations

	dblVec y(size + 1);

	for (UInt i = 0; i < (UInt)p.num_rows(); i ++)
	{
		dblVec pRow(p.num_cols());
		for (UInt j = 0; j < (UInt)p.num_cols(); j ++)
		{
			pRow[j] = p[i][j];
		}
		y[i] = FUNK(pRow);
	}

	int nfunk = 1;

	ameoba(p,y,ftol,FUNK,nfunk);

	cout << "Final Simplex Matrix" << endl;
	for (UInt i = 0; i < (UInt)p.num_rows(); i ++)
	{
		for (UInt j = 0; j < (UInt)p.num_cols(); j ++)
		{
			cout << p[i][j] << " ";
		}
		cout << " evaluates to:  " << y[i] << endl;
	}

	return 0;
}


double funk(dblVec &params)
{
	double TINY = 1e-20;

	double aD = params[0];
	double aS = params[1];
	double gamma = fabs(params[2])+TINY;
	double betabs = params[3];
	double betasb = params[4];

	double e = 2.7128;

	UInt nDbs = unfavEbs[0].size();
	UInt nDsb = unfavEsb[0].size();
	UInt nSbs = favEbs[0].size();
	UInt nSsb = favEsb[0].size();

	vector <double> scores;
	double meanScore = 0.0;
	for (long int i = 0; i < (long int)wtE.size(); i ++)
	{
		double pDbs = 0.0;
		for (UInt j = 0; j < nDbs; j ++)
		{
			double diff = sqrt ( fabs (unfavEbs[i][j]-wtE[i]) );
			double thisPD = 1.0 / ( pow(e,betabs*diff) + 1.0);
			//cout << thisPD << " " << " " << unfavE[i][j] << "; ";
			pDbs += thisPD;
		}

		double pDsb = 0.0;
		


		if (!(pD >= 0.0 && pD <= 1.0))
		{
			cout << "BAD PARAMS" << endl;
			return 0.0;
		}


		pS = pS / (double)favE[i].size();
		//cout << "pS " << pS << endl;
		if (! (pS >= 0.0 && pS <= 1.0))
		{
			cout << "BAD PARAMS" << endl;
            return 0.0;
        }


		double thisScore = wtE[i] + aD*log(pD + gammaD) + aS*log(pS + gammaS);
		scores.push_back(thisScore);
		ASSERT (thisScore > -1e20 && thisScore < 1e20);
		meanScore += thisScore;
	}
	meanScore = meanScore / (double)scores.size();
	ASSERT (meanScore > -1e20 && meanScore < 1e20);
	double tmpSum;
	for (UInt i = 0; i < scores.size(); i ++)
	{
		tmpSum += pow(scores[i] - meanScore, 2);
	}
	double standDev = TINY + sqrt( tmpSum / ((double)scores.size() - 1.0));
	ASSERT (standDev > -1e20 && standDev < 1e20);
	double zScore = 0.0;
	for (UInt i = 0; i < numStructs; i ++)
	{
		zScore +=(scores[wtStruct+i] - meanScore) / standDev;
	}
	zScore /= (double)numStructs;
	cout << "preScore " << zScore << " ";
	
	if (aD > -400.0) zScore += (aD + 400.0)/10.0;
	if (aD < -600.0) zScore += (-600.0 - aD)/100.0;  
	if (aS < 80) zScore += (80.0 - aS)/10.0;
	if (aS >120) zScore += (aS-120.0)/100.0;
	if (fabs(gammaS) > 1e-5 ) zScore += (fabs(gammaS)-1e-5); 
	if (fabs(gammaD) > 1e-5 ) zScore += (fabs(gammaD)-1e-5);
	if (betaS > -5e-3 ) zScore += (betaS + 5e-3)*100;
	if (betaD > -5e-3 ) zScore += (betaD + 5e-3)*100;
	if (betaS < -10.0) zScore += (-10.0 - betaS)/10.0;
	if (betaD < -10.0) zScore += (-10.0 - betaD)/10.0;
/*
	// keep A < 10000
	if ( fabs(aD) > 10000.0) zScore += fabs(aD)/5000.0 - 2.0;
	if ( fabs(aS) > 10000.0) zScore += fabs(aS)/5000.0 - 2.0;
	//	if ( aS > 10000.0 || aS <= 0.0 ) zScore += fabs(aS)/5000.0 - 2.0;
	// keep gamma smaller than 1e-3
	if ( fabs(gammaD) > 5.0e-2) zScore += fabs(gammaD)/5.0e-2 - 1.0;
	if ( fabs(gammaS) > 5.0e-2) zScore += fabs(gammaS)/5.0e-2 - 1.0;
	// keep T approximately between 1 and 1000
	if ( betaS > -1.0) zScore += fabs(1.0 - betaS);
	if ( betaD > -1.0) zScore += fabs(1.0 - betaD);
	if ( betaS < -1000.0) zScore += fabs(betaS)/2000.0 - 0.5;
	if ( betaD < -1000.0) zScore += fabs(betaD)/2000.0 - 0.5;
*/
	cout << "current z-score: " << zScore << " thisScore " << scores[wtStruct] << endl;
	return zScore;
}

void get_psum(dblMat &p, dblVec &psum)
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

void ameoba(dblMat &p, dblVec &y, const double ftol, double FUNK(dblVec &), int &nfunk)
{
	const double TINY=1.0e-10;
	int i, ihi, ilo, inhi, j;
	double rtol, ysave, ytry;
	int mpts=p.num_rows();
	int ndim=p.num_cols();
	dblVec psum(ndim);
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

		ytry = amotry(p,y,psum,FUNK,ihi,-1.0);
		if (ytry < y[ilo]) ytry=amotry(p,y,psum,FUNK,ihi,2.0);	
		else if (ytry >= y[inhi])
		{
			ysave=y[ihi];
			ytry=amotry(p,y,psum,FUNK,ihi,0.5);
			if (ytry >= ysave)
			{
				for (i=0; i < mpts; i++)
				{
					if (i != ilo) 
					{
						for (j=0; j < ndim; j++) p[i][j]=psum[j]=0.5*(p[i][j]+p[ilo][j]);
						y[i]=FUNK(psum);
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

double amotry(dblMat &p, dblVec &y, dblVec &psum, double FUNK(dblVec &), const int ihi, const double fac)
{
	int j;
	double fac1, fac2, ytry;

	int ndim=p.num_cols();
	dblVec ptry(ndim);
	fac1=(1.0-fac)/ndim;
	fac2=fac1-fac;
	for (j=0; j < ndim; j++) ptry[j]=psum[j]*fac1-p[ihi][j]*fac2;
	ytry=FUNK(ptry);
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
