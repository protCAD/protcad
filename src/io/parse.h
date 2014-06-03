#ifndef PARSE_H
#define PARSE_H

#include <stdlib.h>
#include <string.h>
#include <vector>
#include <iostream>
#include <typedef.h>

namespace Parse
{
	vector <string> parse(string& _currentLine);
	vector <string> parse(string& _currentLine, string& _unwantedChars);
}

/*
Commonly unwanted characters:
' '   space
'\t'  tab
'\n'  new line

NULL ('\0') is considered explicitly by parse() and so should NOT be listed
in the unwanted characters string
*/
#endif



