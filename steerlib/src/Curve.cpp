// Copyright (c) 2015 Nilay Chakraborty
// Copyright (c) 2015 Mahyar Khayatkhoei
// Copyright (c) 2009-2014 Shawn Singh, Glen Berseth, Mubbasir Kapadia, Petros Faloutsos, Glenn Reinman
// See license.txt for complete license.
//

#include <algorithm>
#include <vector>
#include <util/Geometry.h>
#include <util/Curve.h>
#include <util/Color.h>
#include <util/DrawLib.h>
#include "Globals.h"

using namespace Util;

Curve::Curve(const CurvePoint& startPoint, int curveType) : type(curveType)
{
	controlPoints.push_back(startPoint);
}

Curve::Curve(const std::vector<CurvePoint>& inputPoints, int curveType) : type(curveType)
{
	controlPoints = inputPoints;
	sortControlPoints();
}

// Add one control point to the vector controlPoints
void Curve::addControlPoint(const CurvePoint& inputPoint)
{
	controlPoints.push_back(inputPoint);
	sortControlPoints();
}

// Add a vector of control points to the vector controlPoints
void Curve::addControlPoints(const std::vector<CurvePoint>& inputPoints)
{
	for (int i = 0; i < inputPoints.size(); i++)
		controlPoints.push_back(inputPoints[i]);
	sortControlPoints();
}

// Draw the curve shape on screen, usign window as step size (bigger window: less accurate shape)
void Curve::drawCurve(Color curveColor, float curveThickness, int window)
{
#ifdef ENABLE_GUI

	Point p;
	// Robustness: make sure there is at least two control point: start and end points
	if (!checkRobust())
	{
		return;
	}
	// Move on the curve from t=0 to t=finalPoint, using window as step size, and linearly interpolate the curve points
	glColor3f(curveColor.r,curveColor.g,curveColor.b);
	glLineWidth(curveThickness);
	glBegin(GL_LINES);
	for (float t = 0; t <= controlPoints[controlPoints.size() - 1].time; t+=0.001)
	{
		if (calculatePoint(p, t))
		{
			glVertex3f(p.x, p.y, p.z);
		}
	}
	glEnd();
	glFlush();
	return;
#endif
}

// Sort controlPoints vector in ascending order: min-first
void Curve::sortControlPoints()
{
	std::sort(controlPoints.begin(),controlPoints.end(), [](CurvePoint a , CurvePoint b)
	{
		return a.time<b.time;
	});
	return;
}

// Calculate the position on curve corresponding to the given time, outputPoint is the resulting position
bool Curve::calculatePoint(Point& outputPoint, float time)
{
	// Robustness: make sure there is at least two control point: start and end points
	if (!checkRobust())
		return false;

	// Define temporary parameters for calculation
	unsigned int nextPoint;
	float normalTime, intervalTime;

	// Find the current interval in time, supposing that controlPoints is sorted (sorting is done whenever control points are added)
	if (!findTimeInterval(nextPoint, time))
		return false;

	// Calculate position at t = time on curve
	if (type == hermiteCurve)
	{
		outputPoint = useHermiteCurve(nextPoint, time);
	}
	else if (type == catmullCurve)
	{
		outputPoint = useCatmullCurve(nextPoint, time);
	}

	// Return
	return true;
}

// Check Roboustness
bool Curve::checkRobust()
{
	if (controlPoints.size() > 1)
		return true;
	else
		return false;
}

// Find the current time interval (i.e. index of the next control point to follow according to current time)
bool Curve::findTimeInterval(unsigned int& nextPoint, float time)
{
	if (controlPoints[nextPoint+1].time == time)
	{
		nextPoint++;
		return true;
	}
	else
		return false;
	
}

// Implement Hermite curve
Point Curve::useHermiteCurve(const unsigned int nextPoint, const float time)
{
	Point newPosition;
	float normalTime, intervalTime;


	// Calculate time interval, and normal time required for later curve calculations
	normalTime = controlPoints[nextPoint-1].time;
	intervalTime = time;
	// Calculate position at t = time on Hermite curve
		newPosition.x = (2 * intervalTime*intervalTime*intervalTime - 3 * intervalTime*intervalTime + 1)*controlPoints[nextPoint - 1].position.x + (intervalTime*intervalTime*intervalTime - 2* intervalTime*intervalTime + intervalTime)*controlPoints[nextPoint - 1].tangent.x + (-2 * intervalTime*intervalTime*intervalTime + 3 * intervalTime*intervalTime)*controlPoints[nextPoint].position.x + (intervalTime*intervalTime*intervalTime - intervalTime*intervalTime)*controlPoints[nextPoint].tangent.x;
		newPosition.y = (2 * intervalTime*intervalTime*intervalTime - 3 * intervalTime*intervalTime + 1)*controlPoints[nextPoint - 1].position.y + (intervalTime*intervalTime*intervalTime - 2 * intervalTime*intervalTime + intervalTime)*controlPoints[nextPoint - 1].tangent.y + (-2 * intervalTime*intervalTime*intervalTime + 3 * intervalTime*intervalTime)*controlPoints[nextPoint].position.y + (intervalTime*intervalTime*intervalTime - intervalTime*intervalTime)*controlPoints[nextPoint].tangent.y;
		newPosition.z = (2 * intervalTime*intervalTime*intervalTime - 3 * intervalTime*intervalTime + 1)*controlPoints[nextPoint - 1].position.z + (intervalTime*intervalTime*intervalTime - 2 * intervalTime*intervalTime + intervalTime)*controlPoints[nextPoint - 1].tangent.z + (-2 * intervalTime*intervalTime*intervalTime + 3 * intervalTime*intervalTime)*controlPoints[nextPoint].position.z + (intervalTime*intervalTime*intervalTime - intervalTime*intervalTime)*controlPoints[nextPoint].tangent.z;

	// Return result
	return newPosition;
}

// Implement Catmull-Rom curve
Point Curve::useCatmullCurve(const unsigned int nextPoint, const float time)
{
	Point newPosition;

	

	// Calculate time interval, and normal time required for later curve calculations

	// Calculate position at t = time on Catmull-Rom curve
	
	// Return result
	return newPosition;
}