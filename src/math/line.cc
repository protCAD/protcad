// filename: line.cpp
// contents: class line implementation

#include "line.h"

// static class member initialization

UInt line::howMany = 0;

// Constructors

line::line()
{	
#ifdef __LINE_DEBUG
	cout << "Line constructor called: "
	     << "line::line() "
	     << endl;
#endif
	point point1(0.0,0.0,0.0);
	point point2(0.0,0.0,0.0);
	itsPoints.push_back(point1);
	itsPoints.push_back(point2);
	howMany++;
}

line::line(const point&  _point1, const point& _point2)
{     	
#ifdef __LINE_DEBUG
	cout << "Line constructor called: " << endl
	     << "line::line(const point&, const point&) "
	     << endl;
#endif
	itsPoints.push_back(_point1);
	itsPoints.push_back(_point2);
	howMany++;
}			

// Deep Copy constructor

line::line(const line& rhs)
{	
#ifdef __LINE_DEBUG
	cout << "Line copy constructor called " << endl;
#endif
	itsPoints.push_back(rhs.getPoint(0));
	itsPoints.push_back(rhs.getPoint(1));
}

// Destructor
line::~line() 
{	
#ifdef __LINE_DEBUG
	cout << "Line destructor called" << endl;
#endif
	howMany--;
}

// Coordinate Related Operations

void line::transform(dblMat* _pDoubleMatrix)
{	
	itsPoints[0].transform(_pDoubleMatrix);
	itsPoints[1].transform(_pDoubleMatrix);
}

void line::transform(const dblMat& _dblMat)
{	
	itsPoints[0].transform(_dblMat);
	itsPoints[1].transform(_dblMat);
}

void line::translate(dblVec* _pDoubleVector)
{	
	itsPoints[0].translate(_pDoubleVector);
	itsPoints[1].translate(_pDoubleVector);
}

void line::translate(const dblVec& _dblVec)
{	
	itsPoints[0].translate(_dblVec);
	itsPoints[1].translate(_dblVec);
}	 

//  Distance Related Operations
 
double line::distance(line* pOtherLine) const
{	
	return calculateDistanceOfClosestApproach(itsPoints[0],
	  itsPoints[1],pOtherLine->getPoint(0), pOtherLine->getPoint(1));
}

double line::distance(const line& otherLine) const
{	return calculateDistanceOfClosestApproach(itsPoints[0],
	  itsPoints[1],otherLine.getPoint(0), otherLine.getPoint(1));
}

double line::calculateDistanceOfClosestApproach(point P0, point P1,
	point Q0, point Q1) const
// Copyright 2001, softSurfer (www.softsurfer.com)
// This code may be freely used and modified for any purpose
// providing that this copyright notice is included with it.
// SoftSurfer makes no warranty for this code, and cannot be held
// liable for any real or imagined damage resulting from its use.
// Users of this code must verify correctness for their application.

// Assume that classes are already given for the objects:
//	Point and Vector with
//		coordinates {double x, y, z;}
//		operators for:
//			Point = Point ± Vector
//			Vector = Point - Point
//			Vector = Vector ± Vector
//			Vector = Scalar * Vector
//	Line and Segment with defining points {Point P0, P1;}
//	Track with initial position and velocity vector
//		{Point P0; Vector v;}
//===================================================================

#define SMALL_NUM  0.00000001 // anything that avoids division overflow
#define norm(v) sqrt(dot_prod(v,v))	// norm = length of vector
{
	dblVec u = P1 - P0;
	dblVec v = Q1 - Q0;
	dblVec w = P0 - Q0;
	double a = dot_prod(u,u);	// always >= 0
	double b = dot_prod(u,v);
	double c = dot_prod(v,v);	// always >= 0
	double d = dot_prod(u,w);
	double e = dot_prod(v,w);
	double D = a*c - b*b;		// always >= 0
	double sc, tc;

	// compute the line parameters of the two closest points
	if (D < SMALL_NUM)
	{	// the lines are almost parallel
		sc = 0.0;
		tc = (b>c ? d/b : e/c);	// use the largest denominator
	}
	else
	{
		sc = (b*e - c*d) / D;
		tc = (a*e - b*d) / D;
	}

	// get the difference of the two closest points
	dblVec dP = w + (sc * u) - (tc * v);	// = L1(sc) - L2(tc)

	return norm(dP);	// return the closest distance
}
