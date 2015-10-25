/*!
*
* \author Nilay Chakraborty(RUID:165003847)
*
*/


#include "obstacles/GJK_EPA.h"
#include "algorithm"

SteerLib::GJK_EPA::GJK_EPA()
{
}

//Look at the GJK_EPA.h header file for documentation and instructions
bool SteerLib::GJK_EPA::intersect(float& return_penetration_depth, Util::Vector& return_penetration_vector, const std::vector<Util::Vector>& _shapeA, const std::vector<Util::Vector>& _shapeB)
{
	std::vector<Util::Vector> simplex = GJK(_shapeA,_shapeB);
	if (simplex.size()!=0)
	{
		EPA(return_penetration_depth,return_penetration_vector,_shapeA,_shapeB,simplex);
		return true;
	}
	else
	{
		return false; // There is no collision
	}
}

void SteerLib::GJK_EPA::EPA(float& return_penetration_depth, Util::Vector& return_penetration_vector, const std::vector<Util::Vector>& _shapeA, const std::vector<Util::Vector>& _shapeB, std::vector<Util::Vector> simplex) {

	while (true)
	{
		SteerLib::GJK_EPA::Edge e = findClosestEdge(simplex);

		Util::Vector p = support(_shapeA,_shapeB,e.normal);
		double d = p.x*e.normal.x + p.y * e.normal.y + p.z *e.normal.z;


		if (d - e.distance < 0.00001)
		{
			return_penetration_depth = d;
			return_penetration_vector = e.normal;
			return;
		}
		else
		{
			simplex[e.index] = p;

		}

	}

}

SteerLib::GJK_EPA::Edge SteerLib::GJK_EPA::findClosestEdge(std::vector<Util::Vector> simplex)
{
	SteerLib::GJK_EPA::Edge closest;
	closest.distance = std::numeric_limits<double>::max();

	for (int i = 0; i < simplex.size(); i++)
	{
		int j = i + 1 == simplex.size() ? 0 : i + 1;

		Util::Vector a = simplex[i];
		Util::Vector b = simplex[j];

		Util::Vector e;
		e.x = b.x - a.x;
		e.y = b.y - a.y;
		e.z = b.z - a.z;

		Util::Vector oa = a;

		double dotProdEE = e.x * e.x + e.y * e.y + e.z * e.z;
		double dotProdEOA = e.x * oa.x + e.y * oa.y + e.z * oa.z;

		Util::Vector crossoa;
		crossoa.x = oa.x*dotProdEE;
		crossoa.y = oa.y*dotProdEE;
		crossoa.z = oa.z*dotProdEE;
		
		Util::Vector crosse;
		crosse.x = e.x*dotProdEOA;
		crosse.y = e.y*dotProdEOA;
		crosse.z = e.z*dotProdEOA;

		Util::Vector n;
		n.x = crossoa.x - crosse.x;
		n.y = crossoa.y - crosse.y;
		n.z = crossoa.z - crosse.z;

		n.x = n.x / sqrt(n.x *n.x +n.y*n.y + n.z*n.z);
		n.y = n.y / sqrt(n.x *n.x + n.y*n.y + n.z*n.z);
		n.z = n.z / sqrt(n.x *n.x + n.y*n.y + n.z*n.z);

		double d = n.x * a.x + n.y * a.y+ n.z * a.z;
		if (d < closest.distance)
		{
			closest.distance = d;
			closest.normal = n;
			closest.index = j;
		}

	}
	return closest;
}


Util::Vector SteerLib::GJK_EPA::getFarthestPoint(std::vector<Util::Vector> shape, Util::Vector direction) {

	int maxIndex = 0;
	double maxDotProd = shape[0].x*direction.x + shape[0].y*direction.y +shape[0].z*direction.z;
	for (int i = 1; i < shape.size(); i++)
	{
		double dotProd = shape[i].x * direction.x + shape[i].y*direction.y + shape[i].z * direction.z;
		if (dotProd > maxDotProd)
		{
			maxDotProd = dotProd;
			maxIndex = i;
		}
	}
	return shape[maxIndex];
}

Util::Vector SteerLib::GJK_EPA::support(std::vector<Util::Vector> shapeA, std::vector<Util::Vector> shapeB, Util::Vector direction) {

	Util::Vector A = getFarthestPoint(shapeA,direction);
	
	direction.x = - direction.x;
	direction.y = - direction.y;
	direction.z = - direction.z;


	Util::Vector B = getFarthestPoint(shapeB, direction);

	Util::Vector minkowDiff;
	minkowDiff.x = A.x - B.x;
	minkowDiff.y = A.y - B.y;
	minkowDiff.z = A.z - B.z;

	return minkowDiff;
}
std::vector<Util::Vector> SteerLib::GJK_EPA::GJK(const std::vector<Util::Vector>& _shapeA, const std::vector<Util::Vector>& _shapeB)
{
	Util::Vector d;

	d.x = 1;
	d.y = 0;
	d.z = -1;

	std::vector<Util::Vector> simplex;

	simplex.push_back(SteerLib::GJK_EPA::support(_shapeA, _shapeB, d));

	d.x = -1;
	d.y = 0;
	d.z = 1;
	while (true)
	{
		simplex.push_back(SteerLib::GJK_EPA::support(_shapeA, _shapeB, d));

		if ((simplex[simplex.size() - 1].x * d.x + simplex[simplex.size() - 1].y * d.y + simplex[simplex.size() - 1].z * d.z) <= 0)
		{
			return std::vector<Util::Vector>();
		}
		else
		{
			if(SteerLib::GJK_EPA::containsOrigin(simplex,d)){
				return simplex;
			}
			
		}

	}

}


bool SteerLib::GJK_EPA::containsOrigin(std::vector<Util::Vector> simplex, Util::Vector& d)
{
	Util::Vector a = simplex[simplex.size() - 1];
	Util::Vector ao;
	ao.x = -a.x;
	ao.y = -a.y;
	ao.z = -a.z;

	if (simplex.size() == 3)
	{
		Util::Vector b = simplex[0];
		Util::Vector c = simplex[1];
		
		Util::Vector ab;
		ab.x = b.x - a.x;
		ab.y = b.y - a.y;
		ab.z = b.z - a.z;

		Util::Vector ac;
		ac.x = c.x - a.x;
		ac.y = c.y - a.y;
		ac.z = c.z - a.z;

		double dotProdABAC = ab.x * ac.x + ab.y * ac.y + ab.z * ac.z;
		double dotProdABAB = ab.x * ab.x + ab.y * ab.y + ab.z * ab.z;
		Util::Vector crossAB;
		crossAB.x = ab.x * dotProdABAC;
		crossAB.y = ab.y* dotProdABAC;
		crossAB.z = ab.z* dotProdABAC;

		Util::Vector crossAC;
		crossAC.x = ac.x * dotProdABAB;
		crossAC.y = ac.y* dotProdABAB;
		crossAC.z = ac.z* dotProdABAB;

		Util::Vector ABPerp;
		ABPerp.x = crossAB.x - crossAC.x;
		ABPerp.y = crossAB.y - crossAC.y;
		ABPerp.z = crossAB.z - crossAC.z;

		ABPerp.x = ABPerp.x / sqrt(ABPerp.x*ABPerp.x + ABPerp.y*ABPerp.y + ABPerp.z*ABPerp.z);
		ABPerp.y = ABPerp.y / sqrt(ABPerp.x*ABPerp.x + ABPerp.y*ABPerp.y + ABPerp.z*ABPerp.z);
		ABPerp.z = ABPerp.z / sqrt(ABPerp.x*ABPerp.x + ABPerp.y*ABPerp.y + ABPerp.z*ABPerp.z);


		double dotProdACAC = ac.x * ac.x + ac.y * ac.y + ac.z * ac.z;
		crossAB.x = ab.x * dotProdACAC;
		crossAB.y = ab.y* dotProdACAC;
		crossAB.z = ab.z* dotProdACAC;

		crossAC.x = ac.x * dotProdABAC;
		crossAC.y = ac.y* dotProdABAC;
		crossAC.z = ac.z* dotProdABAC;

		Util::Vector ACPerp;
		ACPerp.x = crossAC.x - crossAB.x;
		ACPerp.y = crossAC.y - crossAB.y;
		ACPerp.z = crossAC.z - crossAB.z;

		ACPerp.x = ACPerp.x / sqrt(ACPerp.x*ACPerp.x + ACPerp.y*ACPerp.y + ACPerp.z*ACPerp.z);
		ACPerp.y = ACPerp.y / sqrt(ACPerp.x*ACPerp.x + ACPerp.y*ACPerp.y + ACPerp.z*ACPerp.z);
		ACPerp.z = ACPerp.z / sqrt(ACPerp.x*ACPerp.x + ACPerp.y*ACPerp.y + ACPerp.z*ACPerp.z);

		if ((ABPerp.x*ao.x+ABPerp.y*ao.y+ABPerp.z*ao.z)>0)
		{
			simplex.erase(simplex.begin()+1);
			d.x = ABPerp.x;
			d.y = ABPerp.y;
			d.z = ABPerp.z;

		}
		else
		{
			if ((ACPerp.x*ao.x + ACPerp.y*ao.y + ACPerp.z*ao.z) > 0)
			{
				simplex.erase(simplex.begin());
				d.x = ACPerp.x;
				d.y = ACPerp.y;
				d.z = ACPerp.z;

			}
			else {
				return true;
			}

		}

	}
	else
	{
		Util::Vector b = simplex[0];
		Util::Vector ab;
		ab.x = b.x - a.x;
		ab.y = b.y - a.y;
		ab.z = b.z - a.z;

		double dotProdABAO = ab.x * ao.x + ab.y * ao.y + ab.z * ao.z;
		double dotProdABAB = ab.x * ab.x + ab.y * ab.y + ab.z * ab.z;

		Util::Vector crossAB;
		crossAB.x = ab.x * dotProdABAO;
		crossAB.y = ab.y* dotProdABAO;
		crossAB.z = ab.z* dotProdABAO;

		Util::Vector crossAO;
		crossAO.x = ao.x * dotProdABAB;
		crossAO.y = ao.y* dotProdABAB;
		crossAO.z = ao.z* dotProdABAB;

		Util::Vector ABPerp;
		ABPerp.x = crossAB.x - crossAO.x;
		ABPerp.y = crossAB.y - crossAO.y;
		ABPerp.z = crossAB.z - crossAO.z;

		ABPerp.x = ABPerp.x / sqrt(ABPerp.x*ABPerp.x + ABPerp.y*ABPerp.y + ABPerp.z*ABPerp.z);
		ABPerp.y = ABPerp.y / sqrt(ABPerp.x*ABPerp.x + ABPerp.y*ABPerp.y + ABPerp.z*ABPerp.z);
		ABPerp.z = ABPerp.z / sqrt(ABPerp.x*ABPerp.x + ABPerp.y*ABPerp.y + ABPerp.z*ABPerp.z);

		d.x = ABPerp.x;
		d.y = ABPerp.y;
		d.z = ABPerp.z;
	}

	return false;
}
