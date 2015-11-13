class STEERLIB_API GJK_EPA
     {
-		
         public:
-
 			GJK_EPA();
 
-			static struct Edge {
-
-				Util::Vector normal;
-				double distance;
-				int index;
-			};
-
             /*
              *
              *  DO NOT CHANGE THE FUNCTION DEFINITION FOR intersect()
 @@ -136,12 +127,7 @@ namespace SteerLib
              *  DO NOT MODIFY polygon1.xml
              */
             static bool intersect(float& return_penetration_depth, Util::Vector& return_penetration_vector, const std::vector<Util::Vector>& _shapeA, const std::vector<Util::Vector>& _shapeB);
-			static void EPA(float& return_penetration_depth, Util::Vector& return_penetration_vector, const std::vector<Util::Vector>& _shapeA, const std::vector<Util::Vector>& _shapeB, std::vector<Util::Vector> simplex);
-			static Edge findClosestEdge(std::vector<Util::Vector> simplex);
-			static std::vector<Util::Vector> GJK(const std::vector<Util::Vector>& _shapeA, const std::vector<Util::Vector>& _shapeB);
-			static Util::Vector support(std::vector<Util::Vector> shapeA, std::vector<Util::Vector> shapeB, Util::Vector direction);
-			static bool containsOrigin(std::vector<Util::Vector> simplex, Util::Vector& d);
-			static Util::Vector getFarthestPoint(std::vector<Util::Vector> shape, Util::Vector direction);
+
         private:
 
     }; // class GJK_EPA
