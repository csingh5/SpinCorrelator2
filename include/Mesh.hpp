///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                          *** Mesh.hpp ***                                 //
//                                                                           //
// Base class responsible for holding discretized mesh points and weights    //
//                                                                           //
// created November 26, 2017                                                 //
// copyright Christopher N. Singh Binghamton University Physics              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#ifndef Mesh_hpp
#define Mesh_hpp

#include <Eigen/Dense>
#include <numeric>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>

class Mesh 
{
public:
	Mesh(std::string klist_file);

	Eigen::Vector3d point(int index);

	int weight(int index);
	int size();
	int norm();

private:
	std::vector<Eigen::Vector3d> m_points;
	std::vector<int> m_weights;
};

#endif /* Mesh_hpp */
