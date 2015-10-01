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

	std::sort(controlPoints.begin(), controlPoints.end(), [](CurvePoint a, CurvePoint b)
	{
		return a.time<b.time;
	});

	for (std::vector<CurvePoint>::size_type i = 0; i != controlPoints.size(); i++) {
		//std::cout << controlPoints[i].time;
	}

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

	/*//================DELETE THIS PART AND THEN START CODING===================
	static bool flag = false;
	if (!flag)
	{
		std::cerr << "ERROR>>>>Member function useHermiteCurve is not implemented!" << std::endl;
		flag = true;
	}
	//=========================================================================*/

	// Calculate time interval, and normal time required for later curve calculations
	float previousTime, nextTime;

	float nextX, nextY, nextZ, previousX, previousY, previousZ, nextTanX, nextTanY, nextTanZ, previousTanX, previousTanY, previousTanZ;

	previousTime = controlPoints[nextPoint - 1].time;
	previousX = controlPoints[nextPoint - 1].position.x;
	previousY = controlPoints[nextPoint - 1].position.y;
	previousZ = controlPoints[nextPoint - 1].position.z;
	nextTime = controlPoints[nextPoint].time;
	nextX = controlPoints[nextPoint].position.x;
	nextY = controlPoints[nextPoint].position.y;
	nextZ = controlPoints[nextPoint].position.z;
	nextTanX = controlPoints[nextPoint].tangent.x;
	nextTanY = controlPoints[nextPoint].tangent.y;
	nextTanZ = controlPoints[nextPoint].tangent.z;
	previousTanX = controlPoints[nextPoint - 1].tangent.x;
	previousTanY = controlPoints[nextPoint - 1].tangent.y;
	previousTanZ = controlPoints[nextPoint - 1].tangent.z;

	intervalTime = nextTime - previousTime;
	normalTime = (time - previousTime) / intervalTime;

	// Calculate position at t = time on Hermite curve
	float aX, bX, cX, dX;
	aX = -2 * (nextX - previousX) + (nextTanX + previousTanX)/(intervalTime);
	bX = 3 * (nextX - previousX) - (2 * nextTanX + previousTanX) / (intervalTime);
	cX = previousTanX;
	dX = previousX;
	newPosition.x = aX*pow(normalTime,3) + bX*pow(normalTime, 2) + cX*normalTime + dX;

	float aY, bY, cY, dY;
	aY = -2 * (nextY - previousY) + (nextTanY + previousTanY) / (intervalTime);
	bY = 3 * (nextY - previousY) - (2 * nextTanY + previousTanY) / (intervalTime);
	cY = previousTanY;
	dY = previousY;
	newPosition.y = aY*pow(normalTime, 3) + bY*pow(normalTime, 2) + cY*normalTime + dY;

	float aZ, bZ, cZ, dZ;
	aZ = -2 * (nextZ - previousZ) + (nextTanZ + previousTanZ) / (intervalTime);
	bZ = 3 * (nextZ - previousZ) - (2 * nextTanZ + previousTanZ) / (intervalTime);
	cZ = previousTanZ;
	dZ = previousZ;
	newPosition.z = aZ*pow(normalTime, 3) + bZ*pow(normalTime, 2) + cZ*normalTime + dZ;

	std::cout << (time - previousTime) << " " << intervalTime << std::endl;
	//std::cout << previousTime << " " << nextTime << std::endl;
	//std::cout << previousX << " " << newPosition.x << std::endl;
	//std::cout << previousY << " " << newPosition.y << std::endl;
	//std::cout << previousZ << " " << newPosition.z << std::endl;

	// Return result
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