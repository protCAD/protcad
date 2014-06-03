#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "unitSphere.h"

UInt spherePoint::itsSize = 114;

float spherePoint::gVertices[] =
{
0.0000f,0.0000f,1.0000f,  0.0000f,0.3827f,0.9239f, -0.1464f,0.3536f,0.9239f,
-0.2706f,0.2706f,0.9239f,  -0.3536f,0.1464f,0.9239f, -0.3827f,0.0000f,0.9239f,
-0.3536f,-0.1464f,0.9239f,  -0.2706f,-0.2706f,0.9239f, -0.1464f,-0.3536f,0.9239f,
0.0000f,-0.3827f,0.9239f,  0.1464f,-0.3536f,0.9239f, 0.2706f,-0.2706f,0.9239f,
0.3536f,-0.1464f,0.9239f,  0.3827f,0.0000f,0.9239f, 0.3536f,0.1464f,0.9239f,
0.2706f,0.2706f,0.9239f,  0.1464f,0.3536f,0.9239f,  0.0000f,0.7071f,0.7071f,
-0.2706f,0.6533f,0.7071f,  -0.5000f,0.5000f,0.7071f, -0.6533f,0.2706f,0.7071f,
-0.7071f,0.0000f,0.7071f,  -0.6533f,-0.2706f,0.7071f, -0.5000f,-0.5000f,0.7071f,
-0.2706f,-0.6533f,0.7071f,  0.0000f,-0.7071f,0.7071f, 0.2706f,-0.6533f,0.7071f,
0.5000f,-0.5000f,0.7071f,  0.6533f,-0.2706f,0.7071f, 0.7071f,0.0000f,0.7071f,
0.6533f,0.2706f,0.7071f,  0.5000f,0.5000f,0.7071f,  0.2706f,0.6533f,0.7071f,
0.0000f,0.9239f,0.3827f,  -0.3536f,0.8536f,0.3827f, -0.6533f,0.6533f,0.3827f,
-0.8536f,0.3536f,0.3827f,  -0.9239f,0.0000f,0.3827f, -0.8536f,-0.3536f,0.3827f,
-0.6533f,-0.6533f,0.3827f,  -0.3536f,-0.8536f,0.3827f, 0.0000f,-0.9239f,0.3827f,
0.3536f,-0.8536f,0.3827f,  0.6533f,-0.6533f,0.3827f, 0.8536f,-0.3536f,0.3827f,
0.9239f,0.0000f,0.3827f,  0.8536f,0.3536f,0.3827f,  0.6533f,0.6533f,0.3827f,
0.3536f,0.8536f,0.3827f,  0.0000f,1.0000f,0.0000f, -0.3827f,0.9239f,0.0000f,
-0.7071f,0.7071f,0.0000f,  -0.9239f,0.3827f,0.0000f, -1.0000f,0.0000f,0.0000f,
-0.9239f,-0.3827f,0.0000f,  -0.7071f,-0.7071f,0.0000f, -0.3827f,-0.9239f,0.0000f,
0.0000f,-1.0000f,0.0000f,  0.3827f,-0.9239f,0.0000f, 0.7071f,-0.7071f,0.0000f,
0.9239f,-0.3827f,0.0000f,  1.0000f,0.0000f,0.0000f, 0.9239f,0.3827f,0.0000f,
0.7071f,0.7071f,0.0000f,  0.3827f,0.9239f,0.0000f, 0.0000f,0.9239f,-0.3827f,
-0.3536f,0.8536f,-0.3827f,  -0.6533f,0.6533f,-0.3827f, -0.8536f,0.3536f,-0.3827f,
-0.9239f,0.0000f,-0.3827f,  -0.8536f,-0.3536f,-0.3827f, -0.6533f,-0.6533f,-0.3827f,
-0.3536f,-0.8536f,-0.3827f,  0.0000f,-0.9239f,-0.3827f, 0.3536f,-0.8536f,-0.3827f,
0.6533f,-0.6533f,-0.3827f,  0.8536f,-0.3536f,-0.3827f, 0.9239f,0.0000f,-0.3827f,
0.8536f,0.3536f,-0.3827f,  0.6533f,0.6533f,-0.3827f, 0.3536f,0.8536f,-0.3827f,
0.0000f,0.7071f,-0.7071f,  -0.2706f,0.6533f,-0.7071f, -0.5000f,0.5000f,-0.7071f,
-0.6533f,0.2706f,-0.7071f,  -0.7071f,0.0000f,-0.7071f, -0.6533f,-0.2706f,-0.7071f,
-0.5000f,-0.5000f,-0.7071f,  -0.2706f,-0.6533f,-0.7071f, 0.0000f,-0.7071f,-0.7071f,
0.2706f,-0.6533f,-0.7071f,  0.5000f,-0.5000f,-0.7071f, 0.6533f,-0.2706f,-0.7071f,
0.7071f,0.0000f,-0.7071f,  0.6533f,0.2706f,-0.7071f, 0.5000f,0.5000f,-0.7071f,
0.2706f,0.6533f,-0.7071f,  0.0000f,0.3827f,-0.9239f, -0.1464f,0.3536f,-0.9239f,
-0.2706f,0.2706f,-0.9239f,  -0.3536f,0.1464f,-0.9239f, -0.3827f,0.0000f,-0.9239f,
-0.3536f,-0.1464f,-0.9239f,  -0.2706f,-0.2706f,-0.9239f, -0.1464f,-0.3536f,-0.9239f,
0.0000f,-0.3827f,-0.9239f,  0.1464f,-0.3536f,-0.9239f, 0.2706f,-0.2706f,-0.9239f,
0.3536f,-0.1464f,-0.9239f,  0.3827f,0.0000f,-0.9239f, 0.3536f,0.1464f,-0.9239f,
0.2706f,0.2706f,-0.9239f,  0.1464f,0.3536f,-0.9239f, 0.0000f,0.0000f,-1.0000f
};

void spherePoint::setSphereSize(UInt _size)
{
	itsSize = _size;
	return;
}


vector <point*> spherePoint::getAllSpherePoints()
{
	vector <point*> allPoints;
	point* tempPoint = new point;
	for (UInt i = 0; i < itsSize; i++)
	{
		cout << "called" << endl;
		tempPoint->setCoords(gVertices[i*3], gVertices[i*3+1], gVertices[i*3+2]);
		allPoints.push_back(tempPoint);
	}
	cout << "size " << allPoints.size() << " " << allPoints[60]->getX() << endl;
	return allPoints;
}

point* spherePoint::getSpherePoint(UInt _index)
{
	point* tempPoint = new point;
	if (_index >=0 && _index < itsSize)
	{
		tempPoint->setCoords(gVertices[_index*3], gVertices[_index*3+1], gVertices[_index*3+2]);

		return tempPoint;
	}
	else
	{
		cout << "Point out of range ... specify a number between 0 and 114" << endl;
		exit(1);
	}
}

double spherePoint::getX(UInt _index)
{
	if (_index >=0 && _index < itsSize*3)
	{
		return gVertices[_index*3];
	}
	else
	{
		cout << "spherepoint y index out of range!" << endl;
		exit(1);
	}
}

double spherePoint::getY(UInt _index)
{
	if (_index >=0 && _index < itsSize*3)
	{
		return gVertices[_index*3 + 1];
	}
	else
	{
		cout << "spherepoint y index out of range!" << endl;
		exit(1);
	}
}

double spherePoint::getZ(UInt _index)
{
	if (_index >=0 && _index < itsSize*3)
	{
		return gVertices[_index*3 + 2];
	}
	else
	{
		cout << "spherepoint z index out of range!" << endl;
		exit(1);
	}
}
