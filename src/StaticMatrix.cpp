///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                      *** StaticMatrix.cpp ***                             //
//                                                                           //
// static spin susceptibility functor... state is eigemesh at k              //
//                                                                           //
// created November 26, 2017                                                 //
// copyright Christopher N. Singh Binghamton University Physics              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include "StaticMatrix.hpp"


StaticMatrix::StaticMatrix(std::string klist_file) : 
	mesh_at_k(klist_file),
	mesh_at_k_plus_q(klist_file)
{ 
	
	W90 w90("data/input/wan_hr.dat");

	for (int n = 0; n < mesh_at_k.size(); ++n) {
		saes_t s(w90.fft(mesh_at_k.point(n)));
		mesh_at_k.push_back(s);
	}	
}

Eigen::MatrixXcd StaticMatrix::operator()(const Eigen::Vector3d& q) {
	
	W90 w90("data/input/wan_hr.dat");

	if (mesh_at_k_plus_q.size() != 0) mesh_at_k_plus_q.clear();

	for (int n = 0; n < mesh_at_k_plus_q.size(); ++n) {
		saes_t s(w90.fft(mesh_at_k_plus_q.point(n) + q));
		mesh_at_k_plus_q.push_back(s);
	}	

	int dimension = mesh_at_k_plus_q.dimension();
	Eigen::MatrixXcd X = Eigen::MatrixXcd::Zero(dimension, dimension);

	for (int L = 0; L < dimension; ++L) {
		for (int M = 0; M < dimension; ++M) {
			X(L, M) = matrix_element(L, M);
		}
	}
	
	X /= mesh_at_k.norm();

	return X;	
}

	
std::complex<double> StaticMatrix::matrix_element(int L, int M) 
{

	std::complex<double> element = 0;
	int dimension = mesh_at_k.dimension();	

	for (int meshpoint = 0; meshpoint < mesh_at_k.size(); ++meshpoint) {
		for (int alpha = 0; alpha < dimension; ++alpha) {
			for (int beta = 0; beta < dimension; ++beta) {
				element += 
					projection_weight(meshpoint, alpha, beta, L, L, M, M) *
					thermal_occupation(meshpoint, alpha, beta, 0);
			}
		}
	}

	return element;	
}
	
