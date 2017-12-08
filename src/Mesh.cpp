///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                          *** Mesh.cpp ***                                 //
//                                                                           //
// Base class responsible for holding discretized mesh points and weights    //
//                                                                           //
// created November 26, 2017                                                 //
// copyright Christopher N. Singh Binghamton University Physics              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "Mesh.hpp"


Mesh::Mesh(std::string klist_file) {

	std::ifstream klist(klist_file);
	assert(klist.is_open());
	
	std::string line;
	Eigen::Vector3d invec;
	double trash, weight, div;

	while (std::getline(klist, line)) {
		if (line == "END") break;
		std::stringstream input(line);
		input >> trash;
		for (int i = 0; i < 3; ++i) input >> invec(i);
		input >> div;
		m_points.push_back(invec * 2 * M_PI / div);
		input >> weight;
		m_weights.push_back(weight);
	}	
}	

Eigen::Vector3d Mesh::point(int index) { return m_points[index]; }

int Mesh::weight(int index) { return m_weights[index]; }

int Mesh::size() { return m_points.size(); }

int Mesh::norm() { 
	return std::accumulate(m_weights.begin(), m_weights.end(), 0);
}
