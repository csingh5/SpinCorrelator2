///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                *** TransverseSusceptibility.cpp ***                       //
//                                                                           //
// RPA spin susceptibility functor... state is eigemesh at k                 //
//                                                                           //
// created December 12, 2017                                                 //
// copyright Christopher N. Singh Binghamton University Physics              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include "TransverseSusceptibility.hpp"

TransverseSusceptibility::TransverseSusceptibility(std::string klist_file, kelvin T) : 
		Xud(klist_file, T, "ud")
	{ 
		// Empty constructor body	
	}


Eigen::MatrixXcd 
TransverseSusceptibility::operator()(const Eigen::Vector3d& q, meV omega, eV U, eV J) 
{

	xud = Xud(q, omega);

	return (I + xud * (4*Ubar(U, J) - 2*Vbar(U, J))).inverse() * xud;	
}


Eigen::MatrixXd 
TransverseSusceptibility::Ubar(eV U, eV J) 
{
	Eigen::MatrixXd ubar = Eigen::MatrixXd::Zero(numorbs2, numorbs2);

	int L12 = 0, L34 = 0;

	for (int L1 = 0; L1 < numorbs; ++L1) {
		for (int L2 = 0; L2 < numorbs; ++L2) {
			for (int L3 = 0; L3 < numorbs; ++L3) {
				for (int L4 = 0; L4 < numorbs; ++L4) {

					if (L1 == L2 && L2 != L3 && L3 == L4) ubar(L12, L34) = (U - 3*J) / 4;
					if (L1 == L4 && L4 != L3 && L3 == L2) ubar(L12, L34) = (3*J - U) / 4;
							
					++L34;			
				}
			}
			L34 = 0;
			++L12;
		}
	}

	return ubar;
}



Eigen::MatrixXd 
TransverseSusceptibility::Vbar(eV U, eV J) 
{
	Eigen::MatrixXd vbar = Eigen::MatrixXd::Zero(numorbs2, numorbs2);

	int L12 = 0, L34 = 0;

	for (int L1 = 0; L1 < numorbs; ++L1) {
		for (int L2 = 0; L2 < numorbs; ++L2) {
			for (int L3 = 0; L3 < numorbs; ++L3) {
				for (int L4 = 0; L4 < numorbs; ++L4) {

					if (L1 == L2 && L2 == L3 && L3 == L4) vbar(L12, L34) = U / 2;
					if (L1 == L2 && L2 != L3 && L3 == L4) vbar(L12, L34) = (U - 2*J) / 2;
					if (L1 == L3 && L3 != L2 && L2 == L4) vbar(L12, L34) = J / 2;
					if (L1 == L4 && L4 != L3 && L3 == L2) vbar(L12, L34) = J / 2;
							
					++L34;			
				}
			}
			L34 = 0;
			++L12;
		}
	}

	return vbar;
}
