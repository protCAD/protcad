#include <string>
#include <vector>
#include <iostream>
#include "typedef.h"

string getEnvironmentVariable(const string& _evname);
bool compareToDelimiters(string _s);
vector<string> parseString(string _s);
