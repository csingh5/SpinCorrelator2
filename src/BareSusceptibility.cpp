///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                 *** BareSusceptibility.cpp ***                            //
//                                                                           //
// bare spin susceptibility functor... state is eigemesh at k                //
//                                                                           //
// created November 30, 2017                                                 //
// copyright Christopher N. Singh Binghamton University Physics              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include "BareSusceptibility.hpp"


Eigen::MatrixXcd 
BareSusceptibility::operator()(const Eigen::Vector3d& q, meV omega) {
	
	W90 w90("data/input/wan_hr.dat");

	mesh_at_k_plus_q.clear();

	for (int n = 0; n < mesh_at_k_plus_q.size(); ++n) {
		saes_t s(w90.fft(mesh_at_k_plus_q.point(n) + q));
		mesh_at_k_plus_q.push_back(s);
	}	

	mesh_at_k.tune_total_density(8 - 2.55);
	mesh_at_k_plus_q.tune_total_density(8 - 2.55);

	int dimension = mesh_at_k.dimension();
	int orbs = dimension / 2;
	int orbs2 = orbs * orbs;
	Eigen::MatrixXcd X = Eigen::MatrixXcd::Zero(orbs2, orbs2);

	int L12 = 0;
	int L34 = 0;
	
	for (int L1 = 0; L1 < orbs; ++L1) {
		for (int L2 = 0; L2 < orbs; ++L2) {
			for (int L3 = 0; L3 < orbs; ++L3) {
				for (int L4 = 0; L4 < orbs; ++L4) {

					X(L12, L34) = matrix_element(L1, L2, L3, L4, omega);

/*					
					if (L1 < orbs && L2 < orbs && L3 < orbs && L4 < orbs) {
						X(L12, L34) = matrix_element(L1, L2, L3, L4, omega);
					
					} else if (L1 >= orbs && L2 >= orbs && L3 >= orbs && L4 >= orbs) {
						X(L12, L34) = matrix_element(L1, L2, L3, L4, -omega);

					} else if (L1 >= orbs && L2 < orbs && L3 < orbs && L4 >= orbs) {
						X(L12, L34) = matrix_element(L1, L2, L3, L4, omega);
					
					} else if (L1 < orbs && L2 >= orbs && L3 >= orbs && L4 < orbs) {
						X(L12, L34) = matrix_element(L1, L2, L3, L4, -omega);
					
					} else {
						continue;
					}
*/

					++L34;
				}
			}
			L34 = 0;
			++L12;
		}	

	}

	return X;	
}



std::complex<double> 
BareSusceptibility::matrix_element(int L1, int L2, int L3, int L4, meV omega) 
{

	if (!spin_similar(L1, L4)) return 0;
	
	std::complex<double> element = 0;
	int dimension = mesh_at_k.dimension();
	int dn = dimension / 2;

	for (int meshpoint = 0; meshpoint < mesh_at_k.size(); ++meshpoint) {
		for (int alpha = 0; alpha < dimension; ++alpha) {
			for (int beta = 0; beta < dimension; ++beta) {
				element +=
				   	static_cast<double>(mesh_at_k.weight(meshpoint)) *	
					(projection_weight(meshpoint, alpha, beta, L1, L2, L3, L4) -
					 projection_weight(meshpoint, alpha, beta, L1+dn, L2+dn, L3+dn, L4+dn)) *
					thermal_occupation(meshpoint, alpha, beta, omega);
			}
		}
	}

	return element / (2.0 * mesh_at_k.norm());	
}


