///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                   *** BareSusceptibility.cpp ***                          //
//                                                                           //
// bare spin susceptibility functor... state is eigemesh at k                //
//                                                                           //
// created December 12, 2017                                                 //
// copyright Christopher N. Singh Binghamton University Physics              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "BareSusceptibility.hpp"

BareSusceptibility::BareSusceptibility(std::string klist_file, kelvin T, 
		std::string spinchannel) : 
	mesh_at_k(klist_file, T),
	mesh_at_k_plus_q(klist_file, T),
	m_spinchannel(spinchannel)
{ 
	
	W90 w90("data/input/wan_hr.dat");

	for (int n = 0; n < mesh_at_k.size(); ++n) {
		saes_t s(w90.fft(mesh_at_k.point(n)));
		mesh_at_k.push_back(s);
	}	
}



Eigen::MatrixXcd BareSusceptibility::operator()(const Eigen::Vector3d& q, 
	meV omega) 
{
	
	W90 w90("data/input/wan_hr.dat");

	mesh_at_k_plus_q.clear();

	for (int n = 0; n < mesh_at_k_plus_q.size(); ++n) {
		saes_t s(w90.fft(mesh_at_k_plus_q.point(n) + q));
		mesh_at_k_plus_q.push_back(s);
	}	

	int dimension = mesh_at_k.dimension();
	int numorbs = dimension / 2;
	int numorbs2 = numorbs * numorbs;
	int L12, L34;

	L12 = 0;
	L34 = 0;

	Eigen::MatrixXcd X = Eigen::MatrixXcd::Zero(numorbs2, numorbs2);

	for (int L1 = 0; L1 < numorbs; ++L1) {
		for (int L2 = 0; L2 < numorbs; ++L2) {
			for (int L3 = 0; L3 < numorbs; ++L3) {
				for (int L4 = 0; L4 < numorbs; ++L4) {
					X(L12, L34) = matrix_element(L1, L2, L3, L4, omega);		
					++L34;			
				}
			}
			L34 = 0;
			++L12;
		}
	}
	
	return X;	
}



inline std::complex<double> 
BareSusceptibility::matrix_element(int L1, int L2, int L3, int L4, meV omega) 
{
	std::complex<double> element = 0;
	int dimension = mesh_at_k.dimension();	

	for (int meshpoint = 0; meshpoint < mesh_at_k.size(); ++meshpoint) {
		for (int alpha = 0; alpha < dimension; ++alpha) {
			for (int beta = 0; beta < dimension; ++beta) {
				element +=
					static_cast<double>(mesh_at_k.weight(meshpoint)) *	
					projection_weight(meshpoint, alpha, beta, L1, L2, L3, L4) *
					thermal_occupation(meshpoint, alpha, beta, omega);
			}
		}
	}

	return element / static_cast<double>(mesh_at_k.norm());	
}



inline std::complex<double>
BareSusceptibility::projection_weight(int kpoint, int alpha, int beta, int L1, int L2, int L3, int L4)
{
		std::complex<double> U1, U2, U3, U4;

		U1 = mesh_at_k.eigenvector(kpoint, L1, alpha);
		U2 = mesh_at_k_plus_q.eigenvector(kpoint, L2, beta);
		U3 = mesh_at_k_plus_q.eigenvector(kpoint, L3, beta);
		U4 = mesh_at_k.eigenvector(kpoint, L4, alpha);

		if (m_spinchannel == "uu") return conj(U1) * U2 * conj(U3) * U4;
		if (m_spinchannel == "dd") return U1 * conj(U2) * U3 * conj(U4);
		if (m_spinchannel == "ud") return U1 * U2 * conj(U3) * conj(U4);

		std::cout << "ERROR IN PROJECTION WEIGHT!!!" << std::endl;
		exit(-1);
}

inline std::complex<double>
BareSusceptibility::thermal_occupation(int kpoint, int alpha, int beta, meV omega)
{

		double numerator;
		std::complex<double> denominator = 0;
		std::complex<double> ieta (0,5e-4);

		if (omega == 0) ieta = 0;

		numerator = mesh_at_k_plus_q.occupation(kpoint, beta) - 
			mesh_at_k.occupation(kpoint, alpha);

		denominator = 
			(omega/1000) + ieta +
			mesh_at_k.eigenvalue(kpoint, alpha) - 
			mesh_at_k_plus_q.eigenvalue(kpoint, beta);

		return numerator / denominator;
}


