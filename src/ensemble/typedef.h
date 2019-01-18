
#ifdef __APPLE__
    #include "/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.14.sdk/usr/include/stdio.h"
#else
    #include "/usr/include/stdio.h"
#endif
#include <iostream>
#include <string>
#include <vector>
using namespace std;
#include "tnt/tnt.h"
#include "tnt/vec.h"
#include "tnt/cmat.h"

#ifndef TYPEDEF_H
#define TYPEDEF_H

#define PI 3.1415926535
#define KB 0.0019872041 //Boltzmann constant kcal/mol
#define EU 2.7182818284 //Eulers number (base of natural log)

typedef unsigned int UInt;
typedef vector<double> DouVec;
typedef vector<double> DblVec;
typedef vector<UInt> UIntVec;
typedef vector<string> StrVec;
typedef TNT::Vector<double> dblVec;
typedef TNT::Matrix<double> dblMat;

#ifndef CMATH_H
#include "CMath.h"
#endif

#endif
