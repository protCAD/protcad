#include "generalio.h"

string getEnvironmentVariable(const string& _evname)
{
	const char* convEVName = _evname.c_str();
	char* pEVString = getenv(convEVName);
	if (pEVString == 0)
	{	cout << "Environment variable " << _evname << " undefined." << endl;
		cout << "Please set it properly and re-execute the program." << endl;             
		exit(1);
	}               
	string EVstring;    
	for (;*pEVString != '\0';pEVString++)
	{
		EVstring += *pEVString;
	}
	return EVstring;
}

bool compareToDelimiters(string _s)
{
	string space = " ";
	string comma = ",";
	string semi = ";";
	string colon = ":";
	string tab = "\t";
	if ( _s == space || _s == comma || _s == semi || 
		_s == colon || _s == tab)
	{	return true;
	}
	return false;
}
	
vector<string> parseString(string _string)
{
	vector<string> stringlist;
	vector<UInt> startlist;
	vector<UInt> endlist;
	UInt length = _string.size();
	int start = -1;
	int end = -1;
	for (UInt i=0; i<length; i++)
	{
		if ( !(compareToDelimiters(_string.substr(i,1))) )
		{	if (start == -1)
			{	//found one
				start = i;
				startlist.push_back(start);
			}
		}
		else
		{
			if (start == -1)
			{	continue;
			}
			else if (end == -1)
			{	// found one
				end = i;
				endlist.push_back(end-1);
				start = -1;
				end = -1;
			}
		}
	}
	if ( !(compareToDelimiters(_string.substr(length-1,1))) )
	{	endlist.push_back(length-1);
	}

/*
	for (UInt i=0; i< startlist.size(); i++)
	{	cout << startlist[i] << " " << endlist[i] << endl;
	}
*/
	for (UInt i=0; i< startlist.size(); i++)
	{	stringlist.push_back(_string.substr(startlist[i],(endlist[i]-startlist[i]+1)));	
	}

	return stringlist;
}
