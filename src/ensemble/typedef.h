#include <stdio.h>
#include <iostream>
#include <string>
#include <vector>
using namespace std;
#include "tnt/tnt.h"
#include "tnt/vec.h"
#include "tnt/cmat.h"

#ifndef TYPEDEF_H
#define TYPEDEF_H

//physical constants
#define PI 3.1415926535 //Pi (Ratio of a circle's circumference to its diameter)
#define KB 0.0019872041 //Boltzmann constant (kcal/mol K)
#define EU 2.7182818284 //Eulers number (base of natural log)
#define KC 332.0636     //Coulombs constant (kcal/mol)

//amino acid types
#define Daa 27
#define Nterm 53
#define Cterm 106

//data types
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
