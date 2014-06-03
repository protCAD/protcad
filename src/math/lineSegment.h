// filename: lineSegment.h
// contents: class lineSegment is defined
#include "assert.h"
#include <vector>
#include "typedef.h"
#include "point.h"
#include "line.h"

#ifndef LINESEGMENT_H
#define LINESEGMENT_H

class lineSegment : public line
{
public:
	// Constructor and Destructor declaration
	lineSegment();
	lineSegment(const point& _point1, const point& _point2);
	lineSegment(const lineSegment& _rhs);
	virtual ~lineSegment();

	// Distance Related Operations
	virtual double distance(const lineSegment& otherLine) const;
	virtual double distance(lineSegment* pOtherLine) const;
//	virtual double distance(const point& somePoint) const;
//	virtual double distance(point* pSomePoint) const;
	
	// Static Variable Accessors
	static UInt getHowMany() {return howMany; }

protected:
	virtual double calculateDistanceOfClosestApproach(point P0,
		point P1, point Q0, point Q1) const;
	static UInt howMany;

};
#endif
