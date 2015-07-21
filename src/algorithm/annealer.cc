#include "annealer.h"

UInt annealer::howMany=0;

annealer::annealer(ensemble* _pEnsemble)
{	
	cout << "Annealer constructor called:"
	     << "annealer::annealer()" << endl;
	e =  2.71828182;
	pItsEnsemble = _pEnsemble;
	itsTempDecrement = -1.0;
	howMany++;
	firstrun = true;
	writeEveryStateFlag = false;
	itsLowestEnergy = 100000.0;
}

annealer::~annealer()
{	
	cout<< "Annealer destructor called " << endl;
        howMany--;
}

void annealer::run(double _tempH, double _tempL, UInt _iter, UInt _seed)
{	
	itsSeed = _seed;
	itsRan.setSeed(_seed);
	run(_tempH, _tempL, _iter);
}

void annealer::run(double _tempH, double _tempL, UInt _iter)
{
	itsTemp = _tempH;
	itsNumIterations = _iter;
	int whenToWrite = 0;
	// Decide how often to write progress report to screen
	if (_iter <= 200)
	{	whenToWrite = 10;
	}
	if (_iter <= 1000) 
	{	whenToWrite = 100;
	}
	if (_iter >= 10000 && _iter <=100000)
	{	whenToWrite = 100;
	}

	if(_iter != 0 && _tempH >= _tempL)
	{	itsTempDecrement = (_tempH - _tempL)/_iter;
	}

	if (firstrun) 
	{	itsNumNormalKeep = 0;
		itsNumMetropolisKeep = 0;
		setupSystem(itsRan);
		itsLowestEnergy = 1000000; //energy();
		//itsOldEnergy = itsLowestEnergy;
		firstrun = false;
	}
	else
	{	itsOldEnergy = energy();
	}

	for (UInt round=1; round <= _iter; round ++)
	{
		// Check whether I should decrement temperature
		itsCurrentIteration = round;
		cout << itsCurrentIteration << " ";
		if (fmod((double)itsCurrentIteration,(double)whenToWrite) == 0)
		{	cout << itsCurrentIteration << endl;
		}
		decrementTemp();

		// Make a temporary modification to the system
		vector<int> changedPosition = pItsEnsemble->chooseNextTargetPosition(itsRan);
		cout << "TARGET: "; for (UInt i = 0; i < changedPosition.size(); i++) cout << changedPosition[i] << " ";

		while (changedPosition[0] < 0 || changedPosition[1] < 0 || changedPosition[2] < 0)
		{
			cout << " ... INVALID" << endl;
			changedPosition = pItsEnsemble->chooseNextTargetPosition(itsRan);
			cout << "TARGET: "; for (UInt i = 0; i < changedPosition.size(); i++) cout << changedPosition[i] << " ";
		}
		cout << endl;

		itsOldEnergy = pItsEnsemble->getPositionEnergy(changedPosition);
		int result = pItsEnsemble->modify(itsRan, changedPosition);
		if (result == -2)
		{	cout << "Anneal cannot run because of above errors" << endl;
			cout << "Aborting anneal" << endl;
			return;
		}

		if (result != -1)
		{

			// calculate the energy of the perturbed system
			// itsNewEnergy = energy();

			//vector<int> temp = pItsEnsemble->getLastModification();
			//cout << "OBTAINED: "; for (UInt i = 0; i < temp.size(); i++) cout << changedPosition[i] << " "; cout << endl;
			itsNewEnergy = pItsEnsemble->getPositionEnergy(changedPosition);
			// calculate the energy difference between old and new
			double deltaEnergy;
			deltaEnergy = itsNewEnergy - itsOldEnergy;
			// If energy of new system is lower, keep it by
			// default
			UInt keepBool = 0;
			cout << " OldE = " << itsOldEnergy;
			cout << " NewE = " << itsNewEnergy;
			cout << " DeltaE = " << deltaEnergy;
			cout << " Temp = " << itsTemp << endl;
			if (deltaEnergy <= 0.0) 
			{
				keepBool = 1;
			}
			// If energy of new system is higher, apply Metropolis
			// criterion to determine if we should keep it...
			if (deltaEnergy > 0.0)
			{
				if (itsTemp <= 0.0)
				{
				cout << "annealer: itsTemp has reached zero!" << endl;
				itsTemp = 1.0;
				}
				itsProbAccept = pow(e,(-deltaEnergy/itsTemp));
				//cout << " ProbAccept = " << itsProbAccept;
				if (itsRan.getNext() < itsProbAccept)
					keepBool = 2;
				else
				{
					keepBool = 0;
					//cout << " reject modification" << endl;
				}
			}
			if (keepBool != 0)
			{
				acceptModification();
				if (itsNewEnergy < itsLowestEnergy)
				{
					itsLowestEnergy = itsNewEnergy;
					cout << "LOWEST = " << itsLowestEnergy << endl;
					string lowest = "Lowest.pdb";
					saveState(lowest);
				}

				itsOldEnergy = itsNewEnergy;
				if (keepBool > 1)
				{
					itsNumMetropolisKeep +=1;
					//cout << " metropolis keep" << endl;
				}
				else
				{
					itsNumNormalKeep +=1;
					//cout << " normal keep" << endl;
				}
			}
			else
			{	rejectModification();
			}
			if (writeEveryStateFlag)
			{
				cout << " saving " ;
				string jason = "Debug";
				char fmodchar[6];
				sprintf(fmodchar,"%i",round);
				string fmodstring = fmodchar;
				jason += fmodstring;
				jason += ".pdb";
				saveState(jason);
				cout << "saved" << endl;
			}
		}
		else
		{
			//cout << "Abort detected at annealer level..." << endl;
		}
	}
}

int annealer::modifyEnsemble(ran& _ran)
{
	return	pItsEnsemble->modify(_ran);
}

double annealer::energy()
{
	return pItsEnsemble->energy();
}

void annealer::setupSystem(ran& _ran)
{
		pItsEnsemble->setupSystem(_ran);
}

void annealer::saveState(string& _filename)
{
		pItsEnsemble->saveState(_filename);
}

void annealer::decrementTemp()
{
	// this temperature decrementor varies the temperature
	// roughly linearly with the number of iterations.... 
	
	// Only executed the first time through...
	if (itsTempDecrement == -1.0)
	{	itsTemp = 0.1;
		itsTempDecrement = 0.0;
	} 
	itsTemp -= itsTempDecrement;
}

void annealer::acceptModification()
{
	pItsEnsemble->acceptModification();
}

void annealer::rejectModification()
{
	pItsEnsemble->rejectModification();
}

void annealer::writeEveryState()
{
	writeEveryStateFlag = true;
}

void annealer::dontWriteEveryState()
{
	writeEveryStateFlag = false;
}
