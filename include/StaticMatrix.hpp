///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                      *** StaticMatrix.hpp ***                             //
//                                                                           //
// static spin susceptibility functor... state is eigemesh at k              //
//                                                                           //
// created November 26, 2017                                                 //
// copyright Christopher N. Singh Binghamton University Physics              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#ifndef StaticMatrix_hpp
#define StaticMatrix_hpp

#include <iostream>
#include <string>
#include "Eigmesh.hpp"
#include "W90.hpp"

using meV = double;

class StaticMatrix
{
public:
	StaticMatrix(std::string klist_file);

	virtual Eigen::MatrixXcd operator()(const Eigen::Vector3d& q);


protected:
	Eigmesh mesh_at_k;
	Eigmesh mesh_at_k_plus_q;

	std::complex<double> matrix_element(int L, int M);


inline bool spin_similar(int index1, int index2) {
	if (index1 < 9 && index2 < 9)   return true;
	if (index1 >= 9 && index2 >= 9) return true;
	return false;
}

inline std::complex<double>
projection_weight(int kpoint, int alpha, int beta, int L1, int L2, int L3, int L4)
{
		std::complex<double> U1, U2, U3, U4;

		U1 = mesh_at_k_plus_q.eigenvector(kpoint, L1, alpha);
		U2 = mesh_at_k.eigenvector(kpoint, L2, beta);
		U3 = mesh_at_k.eigenvector(kpoint, L3, beta);
		U4 = mesh_at_k_plus_q.eigenvector(kpoint, L4, alpha);

		return conj(U1) * U2 * conj(U3) * U4;

		/*
		
		int dimension = mesh_at_k.dimension();
		int orbs = dimension / 2;

		U1 = mesh_at_k.eigenvector(kpoint, L1, alpha);
		U2 = mesh_at_k_plus_q.eigenvector(kpoint, L2, beta);
		U3 = mesh_at_k_plus_q.eigenvector(kpoint, L3, beta);
		U4 = mesh_at_k.eigenvector(kpoint, L4, alpha);
	
		if (L1 < orbs && L2 < orbs && L3 < orbs && L4 < orbs) {
			return U4 * conj(U1) * U2 * conj(U3);
		
		} else if (L1 >= orbs && L2 >= orbs && L3 >= orbs && L4 >= orbs) {
			return conj(U4) * U1 * conj(U2) * U3;	
		
		} else if (L4 >= orbs && L1 >= orbs && L2 < orbs && L3 < orbs) {
			return conj(U4) * U1 * U2 * conj(U3);
		
		} else if (L1 < orbs && L2 >= orbs && L3 >= orbs && L4 < orbs) {
			return U4 * conj(U1) * conj(U2) * U3;
		
		} else {
			return 0;
		}
		*/
}

inline std::complex<double>
thermal_occupation(int kpoint, int alpha, int beta, meV omega)
{

		double numerator, broadening_factor;
		std::complex<double> denominator;
		std::complex<double> i (0,1);

		broadening_factor = 1e-3;
		
		if (omega == 0) broadening_factor = 0;

		numerator = mesh_at_k_plus_q.occupation(kpoint, beta) - 
			mesh_at_k.occupation(kpoint, alpha);

		denominator = 
			omega/1000 + i * broadening_factor +
			mesh_at_k.eigenvalue(kpoint, alpha) - 
			mesh_at_k_plus_q.eigenvalue(kpoint, beta);

		return numerator / denominator;
}

};

#endif /* StaticMatrix_hpp */
