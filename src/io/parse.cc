// parse.cc
#include "parse.h"

vector <string> Parse::parse(string& _currentLine)
{
	StrVec parsedStrings;
	string unwantedChars;
	string tmpStrChi;
	string tmpStr;
	parsedStrings.resize(0);
	tmpStrChi.resize(1);
	tmpStr.resize(0);

        unwantedChars = " \t\n";

	for (UInt i = 0; i < _currentLine.size(); i++)
	{
		bool testloop=false;
		for (UInt j = 0; j < unwantedChars.size(); j++)
		{	if(_currentLine[i]==unwantedChars[j]) {testloop = true; } }

		if (testloop==false)
		{
			tmpStrChi[0] = _currentLine[i];
			tmpStr.append(tmpStrChi);
		}
		else
		{	if (tmpStr.size() != 0)
			{
				parsedStrings.push_back(tmpStr);
				tmpStr.resize(0);
			}
		}
		if ( i == (_currentLine.size() - 1) && tmpStr.size() != 0)
		{
			parsedStrings.push_back(tmpStr);
		}
	}

/*
	for (UInt i = 0; i < parsedStrings.size(); i++)
	{
		cout << parsedStrings[i] << " ";
	}
	cout << endl;
*/
	return parsedStrings;
}

vector <string> Parse::parse(string& _currentLine, string& _unwantedChars)
{
	StrVec parsedStrings;
	string tmpStrChi;
	string tmpStr;
	parsedStrings.resize(0);
	tmpStrChi.resize(1);
	tmpStr.resize(0);

	for (UInt i = 0; i < _currentLine.size(); i++)
	{
		bool testloop=false;
		for (UInt j = 0; j < _unwantedChars.size(); j++)
		{	if(_currentLine[i]==_unwantedChars[j]) {testloop = true; } }

		if (testloop==false)
		{
			tmpStrChi[0] = _currentLine[i];
			tmpStr.append(tmpStrChi);
		}
		else
		{	if (tmpStr.size() != 0)
			{
				parsedStrings.push_back(tmpStr);
				tmpStr.resize(0);
			}
		}
		if ( i == (_currentLine.size() - 1) && tmpStr.size() != 0)
		{
			parsedStrings.push_back(tmpStr);
		}
	}

/*
	for (UInt i = 0; i < parsedStrings.size(); i++)
	{
		cout << parsedStrings[i] << " ";
	}
	cout << endl;
*/
	return parsedStrings;
}
