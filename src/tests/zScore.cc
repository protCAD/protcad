// from Numerical Re/cipies in C++  p415

#include <vector>
#include <iostream>
#include "ensemble.h"
#include <string>
#define FUNK oldFunk
#define NMAX 5000
double funk(dblVec &params);
double disruptFunk(dblVec &params);
double silentFunk(dblVec &params);
double oldFunk(dblVec &params);

static vector < double > wtE;
static vector < vector <double> > favE;
static vector < vector <double> > unfavE;
static vector < double > radius;
static vector < double > offset;
static vector < double > cross;
static vector < double 

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
	UInt wtPos,fStart,fEnd,uStart,uEnd;
	getline(dataFile, currentLine, '\n'); parsedStrings = Parse::parse(currentLine);
	sscanf(parsedStrings[0].c_str(), "%u", &wtPos);
	sscanf(parsedStrings[1].c_str(), "%u", &fStart);
	sscanf(parsedStrings[2].c_str(), "%u", &fEnd);
	sscanf(parsedStrings[3].c_str(), "%u", &uStart);
	sscanf(parsedStrings[4].c_str(), "%u", &uEnd);
	sscanf(parsedStrings[5].c_str(), "%u", &wtStruct);
	

	UInt x = 0;
	while (getline(dataFile, currentLine, '\n'))
	{
		x ++;
		bool lineGood = true;
		double wtEnergy;
		parsedStrings = Parse::parse(currentLine);
		sscanf(parsedStrings[wtPos].c_str(), "%lf", &wtEnergy);
		if (!(wtEnergy > -1e50 && wtEnergy < 1e50)) lineGood = false; 

		double favEnergy;
		vector <double> favEnergyArray;
		for (UInt i = fStart; i <= fEnd; i ++)
		{
			sscanf(parsedStrings[i].c_str(), "%lf", &favEnergy);
			favEnergyArray.push_back(favEnergy);
			if  (!(favEnergy > -1e50 && favEnergy < 1e50)) lineGood = false;
		}
		double unfavEnergy;
		vector <double> unfavEnergyArray;
		for (UInt i = uStart; i <= uEnd; i ++)
		{
			sscanf(parsedStrings[i].c_str(), "%lf", &unfavEnergy);
			unfavEnergyArray.push_back(unfavEnergy);
			if (!(unfavEnergy > -1e50 && unfavEnergy < 1e50)) lineGood = false;
		}
		if (lineGood)
		{
			unfavE.push_back(unfavEnergyArray);
			wtE.push_back(wtEnergy);
			favE.push_back(favEnergyArray);

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

double silentFunk(dblVec &params)
{
	double TINY = 1e-20;

	double aS = fabs(params[0]) * -1.0 - TINY;
	double gammaS = fabs(params[1])+TINY;
	double betaS = fabs(params[2]) * -1.0 - TINY;
	
	double e = 2.7128;

	cout << aS <<  " " << gammaS << " " << betaS <<  endl; 
	vector <double> scores;
	double meanScore = 0.0;
	for (long int i = 0; i < (long int)wtE.size(); i ++)
	{
		double pS = 0.0;
		double wtTemp = pow(e,betaS*wtE[i]);
		for (long int j = 0; j < (long int)unfavE[i].size(); j ++)
		{
			double thisPS = wtTemp / ( pow(e,betaS*unfavE[i][j]) + wtTemp + TINY);
			//cout << thisPD << " " << " " << unfavE[i][j] << "; ";
			pS += thisPS;
		}
		pS = pS / (double)unfavE[i].size();
		//cout << "pS "<< pS << endl;
		if (!(pS >= 0.0 && pS <= 1.0))
		{
			cout << "BAD PARAMS" << endl;
			return 0.0;
		}

		double thisScore = wtE[i] + aS*log(pS+gammaS);
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
	for (UInt i = 0; i < 10; i ++)
	{
		zScore +=(scores[wtStruct] - meanScore) / standDev;
	}
	zScore = zScore / 10.0;

	cout << "current z-score: " << zScore << " thisScore " << scores[wtStruct] << endl;
	return zScore;
}
double disruptFunk(dblVec &params)
{
	double TINY = 1e-20;

	double aD = fabs(params[0]) * -1.0 - TINY;
	double gammaD = fabs(params[1])+TINY;
	double betaD = fabs(params[2]) * -1.0 - TINY;
	
	double e = 2.7128;

	cout << aD <<  " " << gammaD << " " << betaD <<  endl; 
	vector <double> scores;
	double meanScore = 0.0;
	for (long int i = 0; i < (long int)wtE.size(); i ++)
	{
		double pD = 0.0;
		double wtTemp = pow(e,betaD*wtE[i]);
		for (long int j = 0; j < (long int)unfavE[i].size(); j ++)
		{
			double thisPD = wtTemp / ( pow(e,betaD*unfavE[i][j]) + wtTemp + TINY);
			//cout << thisPD << " " << " " << unfavE[i][j] << "; ";
			pD += thisPD;
		}
		pD = pD / (double)unfavE[i].size();
		//cout << "pD " << pD << endl;
		if (!(pD >= 0.0 && pD <= 1.0))
		{
			cout << "BAD PARAMS" << endl;
			return 0.0;
		}

		double thisScore = wtE[i] + aD*log(pD+gammaD);
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
	for (UInt i = 0; i < 10; i ++)
	{
		zScore +=(scores[wtStruct] - meanScore) / standDev;
	}
	zScore = zScore / 10.0;

	cout << "current z-score: " << zScore << " thisScore " << scores[wtStruct] << endl;
	return zScore;
}

double oldFunk(dblVec &params)
{
	double TINY = 1e-20;

	double aD = fabs(params[0]) * -1.0 - TINY;
	double aS = fabs(params[1])+TINY;
	double gammaD = fabs(params[2])+TINY;
	double gammaS = fabs(params[3])+TINY;
	double betaD = fabs(params[4]) * -1.0 - TINY;
	double betaS = fabs(params[5]) * -1.0 - TINY;
	
	double e = 2.7128;

	cout << aD << " " << aS << " " << gammaD << " " << gammaS << " " << betaD << " " << betaS << endl; 
	vector <double> scores;
	double meanScore = 0.0;
	for (long int i = 0; i < (long int)wtE.size(); i ++)
	{
		double pD = 0.0;
		double wtTemp = pow(e,betaD*wtE[i]);
		for (long int j = 0; j < (long int)unfavE[i].size(); j ++)
		{
			double thisPD =  pow(e,betaD*unfavE[i][j]);
			//cout << thisPD << " " << " " << unfavE[i][j] << "; ";
			pD += thisPD;
		}
		pD = wtTemp / (wtTemp + pD); 
		//cout << "pD " << pD << endl;
		if (!(pD >= 0.0 && pD <= 1.0))
		{
			cout << "BAD PARAMS" << endl;
			return 0.0;
		}
		double pS = 0.0;
		wtTemp = pow(e,betaS*wtE[i]);
		for (long int j = 0; j < (long int)favE[i].size(); j ++)
		{
			double thisPS = pow(e,betaS*favE[i][j]);
			//cout << thisPS << " " << " " << favE[i][j] << "; ";
			pS += thisPS;
		}
		pS = wtTemp / (wtTemp + pS);
		//cout << "pS " << pS << endl;
		if (! (pS >= 0.0 && pS <= 1.0))
		{
			cout << "BAD PARAMS" << endl;
            return 0.0;
        }


		double thisScore = wtE[i] + aD*log(pD+gammaD) + aS*log(pS + gammaS);
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
	for (UInt i = 0; i < 10; i ++)
	{
		zScore +=(scores[wtStruct] - meanScore) / standDev;
	}
	zScore = zScore / 10.0;

	cout << "current z-score: " << zScore << " thisScore " << scores[wtStruct] << endl;
	return zScore;
}
double funk(dblVec &params)
{
	double TINY = 1e-20;

	double aD = fabs(params[0]) * -1.0 - TINY;
	double aS = fabs(params[1])+TINY;
	double gammaD = fabs(params[2])+TINY;
	double gammaS = fabs(params[3])+TINY;
	double betaD = fabs(params[4]) * -1.0 - TINY;
	double betaS = fabs(params[5]) * -1.0 - TINY;
	
	double e = 2.7128;

	cout << aD << " " << aS << " " << gammaD << " " << gammaS << " " << betaD << " " << betaS << endl; 
	vector <double> scores;
	double meanScore = 0.0;
	for (long int i = 0; i < (long int)wtE.size(); i ++)
	{
		double pD = 0.0;
		double wtTemp = pow(e,betaD*wtE[i]);
		for (long int j = 0; j < (long int)unfavE[i].size(); j ++)
		{
			double thisPD = wtTemp / ( pow(e,betaD*unfavE[i][j]) + wtTemp + TINY);
			//cout << thisPD << " " << " " << unfavE[i][j] << "; ";
			pD += thisPD;
		}
		pD = pD / (double)unfavE[i].size();
		//cout << "pD " << pD << endl;
		if (!(pD >= 0.0 && pD <= 1.0))
		{
			cout << "BAD PARAMS" << endl;
			return 4.0;
		}
		double pS = 0.0;
		wtTemp = pow(e,betaS*wtE[i]);
		for (long int j = 0; j < (long int)favE[i].size(); j ++)
		{
			double thisPS= wtTemp / ( pow(e,betaS*favE[i][j]) + wtTemp + TINY);
			//cout << thisPS << " " << " " << favE[i][j] << "; ";
			pS += thisPS;
		}
		pS = pS / (double)favE[i].size();
		//cout << "pS " << pS << endl;
		if (! (pS >= 0.0 && pS <= 1.0))
		{
			cout << "BAD PARAMS" << endl;
            return 4.0;
        }


		double thisScore = wtE[i] + aD*log(pD+gammaD) + aS*log(pS + gammaS);
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
	for (UInt i = 0; i < 10; i ++)
	{
		zScore +=(scores[wtStruct] - meanScore) / standDev;
	}
	zScore = zScore / 10.0;

	cout << "current z-score: " << zScore << " thisScore " << scores[wtStruct] << endl;
	return 4.0 + zScore;
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
