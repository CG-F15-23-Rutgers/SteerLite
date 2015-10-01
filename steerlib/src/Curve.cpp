//
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

	float normalTime, intervalTime;


	// Calculate time interval, and normal time required for later curve calculations
	float previousTime, nextTime;

	float nextX, nextY, nextZ, previousX, previousY, previousZ, nextTanX, nextTanY, nextTanZ, previousTanX, previousTanY, previousTanZ;

	// Robustness: make sure there is at least two control point: start and end points
	if (!checkRobust())
	{
		return;
	}
	// Move on the curve from t=0 to t=finalPoint, using window as step size, and linearly interpolate the curve points
	glColor3f(curveColor.r, curveColor.g, curveColor.b);
	glLineWidth(curveThickness);
	glBegin(GL_POINTS);
	int i = 0;
	Point newPosition;
	glVertex3f(controlPoints[0].position.x, controlPoints[0].position.y, controlPoints[0].position.z);
	for (float t = controlPoints[0].time+0.001; t <= controlPoints[controlPoints.size() - 1].time; t += 0.001)
	{
		
		previousTime = controlPoints[i].time;
		
		nextTime = controlPoints[i+1].time;
	
		intervalTime = nextTime - previousTime;
		normalTime = (t - previousTime) / intervalTime;

		// Calculate position at t = time on Hermite curve
		newPosition.x = (1 + 2 * normalTime)*(1 - normalTime)*(1 - normalTime)*controlPoints[i].position.x + normalTime*(1 - normalTime)*(1 - normalTime)*intervalTime*controlPoints[i].tangent.x + normalTime*normalTime*(3 - 2 * normalTime)*controlPoints[i+1].position.x + normalTime*normalTime*(normalTime - 1)*intervalTime*controlPoints[i+1].tangent.x;

		newPosition.y = (1 + 2 * normalTime)*(1 - normalTime)*(1 - normalTime)*controlPoints[i].position.y + normalTime*(1 - normalTime)*(1 - normalTime)*intervalTime*controlPoints[i].tangent.y + normalTime*normalTime*(3 - 2 * normalTime)*controlPoints[i+1].position.y + normalTime*normalTime*(normalTime - 1)*intervalTime*controlPoints[i+1].tangent.y;

		newPosition.z = (1 + 2 * normalTime)*(1 - normalTime)*(1 - normalTime)*controlPoints[i].position.z + normalTime*(1 - normalTime)*(1 - normalTime)*intervalTime*controlPoints[i].tangent.z + normalTime*normalTime*(3 - 2 * normalTime)*controlPoints[i+1].position.z + normalTime*normalTime*(normalTime - 1)*intervalTime*controlPoints[i+1].tangent.z;

			glVertex3f(newPosition.x, newPosition.y, newPosition.z);
			if(t==nextTime)
			i++;
	}
	glEnd();
	glFlush();
	return;
	
#endif
}

// Sort controlPoints vector in ascending order: min-first
void Curve::sortControlPoints()
{

	std::sort(controlPoints.begin(), controlPoints.end(), [](CurvePoint a, CurvePoint b)
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
	if (controlPoints.size() > 1) {
		return true;
	}
	else {
		return false;
	}
}

// Find the current time interval (i.e. index of the next control point to follow according to current time)
bool Curve::findTimeInterval(unsigned int& nextPoint, float time)
{
	float previous;
	float next;

	do{
		if (nextPoint == 0) {
			previous = controlPoints[nextPoint].time;
		}
		else {
			previous = controlPoints[nextPoint - 1].time;
		}

		if (nextPoint <= controlPoints.size()) {
			next = controlPoints[nextPoint].time;
		}
		else {
			return false;
		}
		nextPoint++;
	} while (!(time >= previous && time < next));

		nextPoint--;
	return true;
}

// Implement Hermite curve
Point Curve::useHermiteCurve(const unsigned int nextPoint, const float time)
{
	Point newPosition;
	float normalTime, intervalTime;


	// Calculate time interval, and normal time required for later curve calculations
	float previousTime, nextTime;

	previousTime = controlPoints[nextPoint - 1].time;

	nextTime = controlPoints[nextPoint].time;


	intervalTime = nextTime - previousTime;
	normalTime = (time - previousTime) / intervalTime;

	// Calculate position at t = time on Hermite curve

	newPosition.x = (1 + 2 * normalTime)*(1 - normalTime)*(1 - normalTime)*controlPoints[nextPoint - 1].position.x + normalTime*(1 - normalTime)*(1 - normalTime)*intervalTime*controlPoints[nextPoint - 1].tangent.x + normalTime*normalTime*(3 - 2 * normalTime)*controlPoints[nextPoint].position.x + normalTime*normalTime*(normalTime - 1)*intervalTime*controlPoints[nextPoint].tangent.x;

	newPosition.y = (1 + 2 * normalTime)*(1 - normalTime)*(1 - normalTime)*controlPoints[nextPoint - 1].position.y + normalTime*(1 - normalTime)*(1 - normalTime)*intervalTime*controlPoints[nextPoint - 1].tangent.y + normalTime*normalTime*(3 - 2 * normalTime)*controlPoints[nextPoint].position.y + normalTime*normalTime*(normalTime - 1)*intervalTime*controlPoints[nextPoint].tangent.y;

	newPosition.z = (1 + 2 * normalTime)*(1 - normalTime)*(1 - normalTime)*controlPoints[nextPoint - 1].position.z + normalTime*(1 - normalTime)*(1 - normalTime)*intervalTime*controlPoints[nextPoint - 1].tangent.z + normalTime*normalTime*(3 - 2 * normalTime)*controlPoints[nextPoint].position.z + normalTime*normalTime*(normalTime - 1)*intervalTime*controlPoints[nextPoint].tangent.z;

	return newPosition;
}

// Implement Catmull-Rom curve
Point Curve::useCatmullCurve(const unsigned int nextPoint, const float time)
{
	Point newPosition;

	//================DELETE THIS PART AND THEN START CODING===================
	static bool flag = false;
	if (!flag)
	{
		std::cerr << "ERROR>>>>Member function useCatmullCurve is not implemented!" << std::endl;
		flag = true;
	}
	//=========================================================================


	// Calculate time interval, and normal time required for later curve calculations

	// Calculate position at t = time on Catmull-Rom curve
	
	// Return result
	return newPosition;
}