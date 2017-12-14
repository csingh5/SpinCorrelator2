///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                          *** Eigmesh.cpp ***                              //
//                                                                           //
// class to hold the eigenvalues and eigenvectors on a mesh                  //
//                                                                           //
// created November 27, 2017                                                 //
// copyright Christopher N. Singh Binghamton University Physics              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "Eigmesh.hpp"
#include <iostream>

Eigmesh::Eigmesh(std::string klist, kelvin T) : Mesh(klist), m_T(T) {
	m_eig.reserve(size());	
}
	
eV Eigmesh::eigenvalue(int mesh_point_index, int band_index) {
	assert(m_eig.size() != 0 && "Run a diagonalization first!");
	return m_eig[mesh_point_index].eigenvalues()(band_index);	
}	

std::complex<double> Eigmesh::eigenvector(int mesh_point_index, int row, int col) {
	assert(m_eig.size() != 0 && "Run a diagonalization first!");
	return m_eig[mesh_point_index].eigenvectors()(row, col);
}	

double Eigmesh::abs2eigenvector(int mesh_point_index, int row, int col) {
	assert(m_eig.size() != 0 && "Run a diagonalization first!");
	return m_eig[mesh_point_index].eigenvectors().cwiseAbs2()(row, col);
}	

void Eigmesh::push_back(saes_t& saes) {
	m_eig.push_back(saes);	
}

nelec Eigmesh::density() {
	assert(m_eig.size() != 0 && "Run a diagonalization first!");
	nelec temp = 0;
	for (int meshpoint = 0; meshpoint < Mesh::size(); ++meshpoint) {
		for (int band = 0; band < m_eig[0].eigenvalues().size(); ++band) {
			temp += weight(meshpoint) * fermi(eigenvalue(meshpoint, band));  
		}
	}
	return temp / Mesh::norm();
}	

Eigen::ArrayXd Eigmesh::occs() {
	assert(m_eig.size() != 0 && "Run a diagonalization first!");
 	int dimension = m_eig[0].eigenvalues().size();
	Eigen::ArrayXd d = Eigen::ArrayXd::Zero(m_eig[0].eigenvalues().size());

	for (int meshpoint = 0; meshpoint < Mesh::size(); ++meshpoint) {
		for (int band = 0; band < dimension; ++band) {
			for (int col = 0; col < dimension; ++col) {
				d(band) += weight(meshpoint) * 
					fermi(eigenvalue(meshpoint, col)) * 
					abs2eigenvector(meshpoint, band, col);			
			}
		}
	}
	return d / Mesh::size();
}

nelec Eigmesh::occupation(int mesh_point_index, int band_index) {
	return fermi(m_eig[mesh_point_index].eigenvalues()(band_index));
}

int Eigmesh::dimension() { return m_eig[0].eigenvalues().size(); }

inline double Eigmesh::fermi(eV epsilon) {
	return 1.0 / (exp((1.0 / (kb * m_T)) * (epsilon - chemical_potential)) + 1.0);
}

inline nelec Eigmesh::flevel_occ() {
	Eigen::ArrayXd d = Eigmesh::occs();
	return d.segment<6>(2).sum() + d.segment<6>(11).sum();
}

void Eigmesh::tune_f_level(nelec desired_f_occupancy) {
	while (flevel_occ() < desired_f_occupancy) chemical_potential += 1e-2;
	while (flevel_occ() > desired_f_occupancy) chemical_potential -= 1e-2;
}

void Eigmesh::tune_total_density(nelec desired_total_density) {
	while (density() < desired_total_density) chemical_potential += 1e-2;
	while (density() > desired_total_density) chemical_potential -= 1e-2;
}




