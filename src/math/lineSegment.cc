// filename: lineSegment.cpp
// contents: class lineSegment implementation

#include "lineSegment.h"

// static class member initialization

UInt lineSegment::howMany = 0;

// Constructors

lineSegment::lineSegment()
{	
#ifdef __LINESEGMENT_DEBUG
	cout << "LineSegment constructor called: "
	     << "lineSegment::lineSegment() "
	     << endl;
#endif
	howMany++;
}

lineSegment::lineSegment(const point&  _point1, const point& _point2)
{     	
#ifdef __LINESEGMENT_DEBUG
	cout << "LineSegment constructor called: " << endl
	     << "lineSegment::lineSegment(const point&, const point&) "
	     << endl;
#endif
	itsPoints.push_back(_point1);
	itsPoints.push_back(_point2);
	howMany++;
}			

// Deep Copy constructor

lineSegment::lineSegment(const lineSegment& rhs)
{	
#ifdef __LINESEGMENT_DEBUG
	cout << "LineSegment copy constructor called " << endl;
#endif
	itsPoints.push_back(rhs.getPoint(0));
	itsPoints.push_back(rhs.getPoint(1));
}

// Destructor
lineSegment::~lineSegment() 
{	
#ifdef __LINESEGMENT_DEBUG
	cout << "LineSegment destructor called" << endl;
#endif
	howMany--;
}

//  Distance Related Operations
 
double lineSegment::distance(lineSegment* pOtherLineSegment) const
{	
	return calculateDistanceOfClosestApproach(itsPoints[0],
	  itsPoints[1],pOtherLineSegment->getPoint(0), pOtherLineSegment->getPoint(1));
}

double lineSegment::distance(const lineSegment& otherLineSegment) const
{	return calculateDistanceOfClosestApproach(itsPoints[0],
	  itsPoints[1],otherLineSegment.getPoint(0), otherLineSegment.getPoint(1));
}

double lineSegment::calculateDistanceOfClosestApproach(point P0, point P1,
	point Q0, point Q1) const
// Copyright 2001, softSurfer (www.softsurfer.com)
// This code may be freely used and modified for any purpose
// providing that this copyright notice is included with it.
// SoftSurfer makes no warranty for this code, and cannot be held
// liable for any real or imagined damage resulting from its use.
// Users of this code must verify correctness for their application.

#define SMALL_NUM  0.00000001 // anything that avoids division overflow
#define norm(v) sqrt(dot_prod(v,v))	// norm = length of vector

//		Input:  two 3D line segments S1 and S2
//		Return: the shortest distance between S1 and S2

{
	cout << "Calling lineSegment::calculate" << endl;
	dblVec u = P1 - P0;
	dblVec v = Q1 - Q0;
	dblVec w = P0 - Q0;
	double a = dot_prod(u,u);			// always >= 0
	double b = dot_prod(u,v);
	double c = dot_prod(v,v);			// always >= 0
	double d = dot_prod(u,w);
	double e = dot_prod(v,w);
	double D = a*c - b*b;			// always >= 0
	double sc, sN, sD = D;		// sc = sN / sD, default sD = D >= 0
	double tc, tN, tD = D;		// tc = tN / tD, default tD = D >= 0

	// compute the line parameters of the two closest points
	if (D < SMALL_NUM)
	{	
	
		//cout << "Lines almost parallel" << endl;
		// the lines are almost parallel
		sN = 0.0;
		tN = e;
		tD = c;
		//cout << "sN = " << sN << " tN = " << tN << " tD= " << tD << endl;
	}
	else
	{	// get the closest points on the infinite lines
		sN = (b*e - c*d);
		tN = (a*e - b*d);
		//cout << "sN  = " << sN << " tN = " << tN << endl;
		if (sN < 0)
		{	// sc < 0 => the s=0 edge is visible
			sN = 0.0;
			tN = e;
			tD = c;
		}
		else if (sN > sD)
		{	// sc > 1 => the s=1 edge is visible
			sN = sD;
			tN = e + b;
			tD = c;
		}
	}

	if (tN <= 0)
	{	// tc < 0 => the t=0 edge is visible
		tN = 0.0;
		// recompute sc for this edge
		if (-d < 0)
			sN = 0.0;
		else if (-d > a)
			sN = sD;
		else
		{
			sN = -d;
			sD = a;
		}
	}
	else if (tN > tD)
	{	// tc > 1 => the t=1 edge is visible
		tN = tD;
		// recompute sc for this edge
		if ((-d + b) < 0)
			sN = 0;
		else if ((-d + b) > a)
			sN = sD;
		else
		{
			sN = (-d + b);
			sD = a;
		}
	}
	// finally do the division to get sc and tc
	sc = sN / sD;
	tc = tN / tD;

	// get the difference of the two closest points
	dblVec dP = w + (sc * u) - (tc * v);	// = S1(sc) - S2(tc)

	return norm(dP);	// return the closest distance
}
