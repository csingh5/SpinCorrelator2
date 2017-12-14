///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                   *** BareSusceptibility.hpp ***                          //
//                                                                           //
// bare spin susceptibility functor... state is eigemesh at k                //
//                                                                           //
// created December 9, 2017                                                  //
// copyright Christopher N. Singh Binghamton University Physics              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#ifndef BareSusceptibility_hpp
#define BareSusceptibility_hpp

#include <iostream>
#include <string>
#include "Eigmesh.hpp"
#include "W90.hpp"

using meV = double;

class BareSusceptibility
{
public:
	BareSusceptibility(std::string klist_file, kelvin T, 
			std::string spinchannel = "uu");
	
	Eigen::MatrixXcd operator()(const Eigen::Vector3d& q, meV omega);

protected:
	Eigmesh mesh_at_k;
	Eigmesh mesh_at_k_plus_q;
	std::string m_spinchannel;

	inline std::complex<double> 
	matrix_element(int L1, int L2, int L3, int L4, meV omega);

	inline std::complex<double>
	projection_weight(int kpoint, int alpha, int beta, int L1, int L2, int L3, int L4);

	inline std::complex<double>
	thermal_occupation(int kpoint, int alpha, int beta, meV omega);

};

#endif /* BareSusceptibility_hpp */
