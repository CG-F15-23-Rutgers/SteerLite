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

	//================DELETE THIS PART AND THEN START CODING===================
	static bool flag = false;
	if (!flag)
	{
		std::cerr << "ERROR>>>>Member function drawCurve is not implemented!" << std::endl;
		flag = true;
	}
	//=========================================================================

	// Robustness: make sure there is at least two control point: start and end points

	// Move on the curve from t=0 to t=finalPoint, using window as step size, and linearly interpolate the curve points
	
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
	intervalTime = time - normalTime;
	// Calculate position at t = time on Hermite curve
	for (float t = normalTime; t <= normalTime + intervalTime; t+=0.001f)
	{
		newPosition.x = (2 *t*t*t - 3 *t*t + 1)*controlPoints[nextPoint - 1].position.x + (t*t*t - 2*t*t + t)*controlPoints[nextPoint - 1].tangent.x + (-2 *t*t*t + 3 * t*t )*controlPoints[nextPoint].position.x + (t*t*t - t*t)*controlPoints[nextPoint].tangent.x;
		newPosition.y = (2 * t*t*t - 3 * t*t + 1)*controlPoints[nextPoint - 1].position.y + (t*t*t - 2 * t*t + t)*controlPoints[nextPoint - 1].tangent.y + (-2 * t*t*t + 3 * t*t)*controlPoints[nextPoint].position.y + (t*t*t - t*t)*controlPoints[nextPoint].tangent.y;
		newPosition.z = (2 * t*t*t - 3 * t*t + 1)*controlPoints[nextPoint - 1].position.z + (t*t*t - 2 * t*t + t)*controlPoints[nextPoint - 1].tangent.z + (-2 * t*t*t + 3 * t*t)*controlPoints[nextPoint].position.z + (t*t*t - t*t)*controlPoints[nextPoint].tangent.z;

	}

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