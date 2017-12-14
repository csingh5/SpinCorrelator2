///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                          *** Eigmesh.hpp ***                              //
//                                                                           //
// class to hold the eigenvalues and eigenvectors on a mesh                  //
//                                                                           //
// created November 27, 2017                                                 //
// copyright Christopher N. Singh Binghamton University Physics              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#ifndef Eigmesh_hpp
#define Eigmesh_hpp

#include "Mesh.hpp"

using saes_t = Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd>;
using eV = double;
using nelec = double;
using eVperKelvin = double;
using kelvin = double;

class Eigmesh : public Mesh
{
public:
	Eigmesh(std::string klist, kelvin T);
	
	eV eigenvalue(int mesh_point_index, int band_index);
	
	std::complex<double> eigenvector(int mesh_point_index, int row, int col);
	
	double abs2eigenvector(int mesh_point_index, int row, int col);
	
	void push_back(saes_t& saes);	

	nelec density();

	Eigen::ArrayXd occs();

	int dimension();

	double occupation(int mesh_point_index, int band_index);

	void clear() { m_eig.clear(); }

	void tune_f_level(nelec desired_f_occupancy);

	void tune_total_density(nelec desired_total_density);

	nelec flevel_occ();

private:
	std::vector<saes_t> m_eig;
	eVperKelvin kb = 8.61733e-5;
	kelvin m_T;
	eV chemical_potential = 0;
	
	inline double fermi(eV epsilon);

};

#endif /* Eigmesh_hpp */
