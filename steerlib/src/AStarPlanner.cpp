//
// Copyright (c) 2009-2015 Glen Berseth, Mubbasir Kapadia, Shawn Singh, Petros Faloutsos, Glenn Reinman, Rahul Shome
// See license.txt for complete license.
//


#include <vector>
#include <stack>
#include <set>
#include <map>
#include <string>
#include <iostream>
#include <algorithm> 
#include <functional>
#include <queue>
#include <math.h>
#include "planning/AStarPlanner.h"


#define COLLISION_COST  1000
#define GRID_STEP  1
#define OBSTACLE_CLEARANCE 1
#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))
#define MAX(X,Y) ((X) > (Y) ? (X) : (Y))

namespace SteerLib
{
	AStarPlanner::AStarPlanner(){}

	AStarPlanner::~AStarPlanner(){}

	bool AStarPlanner::canBeTraversed ( int id ) 
	{
		double traversal_cost = 0;
		int current_id = id;
		unsigned int x,z;
		gSpatialDatabase->getGridCoordinatesFromIndex(current_id, x, z);
		int x_range_min, x_range_max, z_range_min, z_range_max;

		x_range_min = MAX(x-OBSTACLE_CLEARANCE, 0);
		x_range_max = MIN(x+OBSTACLE_CLEARANCE, gSpatialDatabase->getNumCellsX());

		z_range_min = MAX(z-OBSTACLE_CLEARANCE, 0);
		z_range_max = MIN(z+OBSTACLE_CLEARANCE, gSpatialDatabase->getNumCellsZ());


		for (int i = x_range_min; i<=x_range_max; i+=GRID_STEP)
		{
			for (int j = z_range_min; j<=z_range_max; j+=GRID_STEP)
			{
				int index = gSpatialDatabase->getCellIndexFromGridCoords( i, j );
				traversal_cost += gSpatialDatabase->getTraversalCost ( index );
				
			}
		}

		if ( traversal_cost > COLLISION_COST ) 
			return false;
		return true;
	}



	Util::Point AStarPlanner::getPointFromGridIndex(int id)
	{
		Util::Point p;
		gSpatialDatabase->getLocationFromIndex(id, p);
		return p;
	}



	bool AStarPlanner::computePath(std::vector<Util::Point>& agent_path,  Util::Point start, Util::Point goal, SteerLib::GridDatabase2D * _gSpatialDatabase, bool append_to_path)
	{
		if (computeEuclideanPath(agent_path, start, goal, _gSpatialDatabase, append_to_path)) {
			return true;
		}
		else {
			return false;
		}

		/*if (computeManhattanPath(agent_path, start, goal, _gSpatialDatabase, append_to_path)) {
			return true;
		}
		else {
			return false;
		}*/
	}

	bool AStarPlanner::computeEuclideanPath(std::vector<Util::Point>& agent_path, Util::Point start, Util::Point goal, SteerLib::GridDatabase2D * _gSpatialDatabase, bool append_to_path) {
		gSpatialDatabase = _gSpatialDatabase;
		float offset = 1;

		Util::Point startPoint = start;
		int startIndex = gSpatialDatabase->getCellIndexFromLocation(start.x, start.z);
		gSpatialDatabase->getLocationFromIndex(startIndex, startPoint);

		Util::Point goalPoint = goal;
		int goalIndex = gSpatialDatabase->getCellIndexFromLocation(goal.x, goal.z);
		gSpatialDatabase->getLocationFromIndex(goalIndex, goalPoint);

		float g = 0;

		// Weight A* multiplication factor
		//float w = 2;

		//Euclidean
		float h = /*w **/ (sqrt(pow((goalPoint.x - startPoint.x), 2) + pow((goalPoint.z - startPoint.z), 2)));

		AStarPlannerNode* startNode = new AStarPlannerNode(startPoint, g, g + h, NULL);

		std::vector<AStarPlannerNode*> closedSet;
		std::vector<AStarPlannerNode*> openSet;

		openSet.push_back(startNode);

		while (!openSet.empty()) {
			AStarPlannerNode* currentNode = openSet[0];
			int currentIndex = 0;
			for (int i = 1; i < openSet.size(); i++) {
				AStarPlannerNode * tempNode = openSet[i];
				//Part 1 with no tie breaks
				/*if (*tempNode < *currentNode) {
					currentNode = tempNode;
					currentIndex = i;
				}*/
				//Part 2 with larger or smaller g values
				if (*tempNode <= *currentNode) {
					currentNode = tempNode;
					currentIndex = i;
				}
			}
			int currentCellIndex = gSpatialDatabase->getCellIndexFromLocation(currentNode->point.x, currentNode->point.z);
			if (currentCellIndex == goalIndex) {
				agent_path.insert(agent_path.begin(), currentNode->point);
				while (currentNode->parent) {
					currentNode = currentNode->parent;
					agent_path.insert(agent_path.begin(), currentNode->point);
				}
				std::cout << "The number of expanded nodes is " << closedSet.size() << std::endl;
				return true;
			}
			openSet.erase(openSet.begin() + currentIndex);
			closedSet.push_back(currentNode);
			for (float x = -offset; x <= offset; x = x + offset) {
				for (float z = -offset; z <= offset; z = z + offset) {
					Util::Point currentPoint = currentNode->point;
					int currentIndex = gSpatialDatabase->getCellIndexFromLocation(currentPoint.x, currentPoint.z);
					gSpatialDatabase->getLocationFromIndex(currentIndex, currentPoint);
					Util::Point nextPoint(currentPoint.x + x, 0, currentPoint.z + z);
					int nextIndex = gSpatialDatabase->getCellIndexFromLocation(nextPoint.x, nextPoint.z);
					gSpatialDatabase->getLocationFromIndex(nextIndex, nextPoint);
					std::cout << nextIndex << std::endl;
					if (canBeTraversed(nextIndex)) {

						float dist = 1;

						//normal cases distance calculation
						g = currentNode->g + dist;

						//increasing diagonal distances
						/*float randomWeightForCost = 100;
						if (x == z) {
							g = currentNode->g + randomWeightForCost * dist;
						}
						else {
							g = currentNode->g + dist;
						}*/

						//Euclidean
						h = /*w **/ (sqrt(pow((goalPoint.x - nextPoint.x), 2) + pow((goalPoint.z - nextPoint.z), 2)));

						bool skipFlag = false;

						AStarPlannerNode * nextNode = new AStarPlannerNode(nextPoint, g, g + h, currentNode);
						for (int i = 1; i < closedSet.size(); i++) {
							if (*closedSet[i] == *nextNode) {
								skipFlag = true;
							}
						}
						for (int i = 1; i < openSet.size(); i++) {
							if (*openSet[i] == *nextNode) {
								//Part 1 with no tie breaks
								/*if (*nextNode < *openSet[i]) {
									openSet[i] = nextNode;
								}*/
								//Part 2 with smaller or larger g values for tie breaks
								if (*nextNode <= *openSet[i]) {
									openSet[i] = nextNode;
								}
								skipFlag = true;
							}
						}
						if (!skipFlag) {
							openSet.push_back(nextNode);
						}
					}
				}
			}
		}
		//TODO
		std::cout << "\nIn A*";
		return false;
	}

	bool AStarPlanner::computeManhattanPath(std::vector<Util::Point>& agent_path, Util::Point start, Util::Point goal, SteerLib::GridDatabase2D * _gSpatialDatabase, bool append_to_path) {
		gSpatialDatabase = _gSpatialDatabase;
		float offset = 1;

		Util::Point startPoint = start;
		int startIndex = gSpatialDatabase->getCellIndexFromLocation(start.x, start.z);
		gSpatialDatabase->getLocationFromIndex(startIndex, startPoint);

		Util::Point goalPoint = goal;
		int goalIndex = gSpatialDatabase->getCellIndexFromLocation(goal.x, goal.z);
		gSpatialDatabase->getLocationFromIndex(goalIndex, goalPoint);

		// Weight A* multiplication factor
		//float w = 8;

		float g = 0;

		//Manhatten
		float h = /*w **/ (abs(goalPoint.x - startPoint.x) + abs(goalPoint.z - startPoint.z));

		AStarPlannerNode* startNode = new AStarPlannerNode(startPoint, g, g + h, NULL);

		std::vector<AStarPlannerNode*> closedSet;
		std::vector<AStarPlannerNode*> openSet;

		openSet.push_back(startNode);

		while (!openSet.empty()) {
			AStarPlannerNode* currentNode = openSet[0];
			int currentIndex = 0;
			for (int i = 1; i < openSet.size(); i++) {
				AStarPlannerNode * tempNode = openSet[i];
				//Part 1 with no tie breaks
				/*if (*tempNode < *currentNode) {
				currentNode = tempNode;
				currentIndex = i;
				}*/
				//Part 2 with larger or smaller g values
				if (*tempNode <= *currentNode) {
					currentNode = tempNode;
					currentIndex = i;
				}
			}
			int currentCellIndex = gSpatialDatabase->getCellIndexFromLocation(currentNode->point.x, currentNode->point.z);
			if (currentCellIndex == goalIndex) {
				agent_path.insert(agent_path.begin(), currentNode->point);
				while (currentNode->parent) {
					currentNode = currentNode->parent;
					agent_path.insert(agent_path.begin(), currentNode->point);
				}
				std::cout << "The number of expanded nodes is " << closedSet.size() << std::endl;
				return true;
			}
			openSet.erase(openSet.begin() + currentIndex);
			closedSet.push_back(currentNode);
			for (float x = -offset; x <= offset; x = x + offset) {
				for (float z = -offset; z <= offset; z = z + offset) {
					//to restrict neighbors in case of manhattan distances
					if (((x == -offset) || (x == offset)) && (z != 0)) {
						continue;
					}
					Util::Point currentPoint = currentNode->point;
					int currentIndex = gSpatialDatabase->getCellIndexFromLocation(currentPoint.x, currentPoint.z);
					gSpatialDatabase->getLocationFromIndex(currentIndex, currentPoint);
					Util::Point nextPoint(currentPoint.x + x, 0, currentPoint.z + z);
					int nextIndex = gSpatialDatabase->getCellIndexFromLocation(nextPoint.x, nextPoint.z);
					gSpatialDatabase->getLocationFromIndex(nextIndex, nextPoint);
					std::cout << nextIndex << std::endl;
					if (canBeTraversed(nextIndex)) {

						float dist = 1;
						g = currentNode->g + dist;

						/*float randomWeightForCost = 100;
						if (x == z) {
							g = currentNode->g + randomWeightForCost * dist;
						}
						else {
							g = currentNode->g + dist;
						}*/

						//Manhatten
						float h = /*w **/ (abs(goalPoint.x - nextPoint.x) + abs(goalPoint.z - nextPoint.z));

						bool skipFlag = false;

						AStarPlannerNode * nextNode = new AStarPlannerNode(nextPoint, g, g + h, currentNode);
						for (int i = 1; i < closedSet.size(); i++) {
							if (*closedSet[i] == *nextNode) {
								skipFlag = true;
							}
						}
						for (int i = 1; i < openSet.size(); i++) {
							if (*openSet[i] == *nextNode) {
								//Part 1 with no tie breaks
								/*if (*nextNode < *openSet[i]) {
								openSet[i] = nextNode;
								}*/
								//Part 2 with smaller or larger g values for tie breaks
								if (*nextNode <= *openSet[i]) {
									openSet[i] = nextNode;
								}
								skipFlag = true;
							}
						}
						if (!skipFlag) {
							openSet.push_back(nextNode);
						}
					}
				}
			}
		}
		//TODO
		std::cout << "\nIn A*";
		return false;
	}
}