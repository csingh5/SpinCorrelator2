///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                *** LongitudinalSusceptibility.cpp ***                     //
//                                                                           //
// RPA spin susceptibility functor... state is eigemesh at k                 //
//                                                                           //
// created December 12, 2017                                                 //
// copyright Christopher N. Singh Binghamton University Physics              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include "LongitudinalSusceptibility.hpp"

LongitudinalSusceptibility::LongitudinalSusceptibility(std::string klist_file, kelvin T) : 
		Xuu(klist_file, T, "uu"), 
		Xdd(klist_file, T, "dd")
	{ 
		// Empty constructor body	
	}


Eigen::MatrixXcd 
LongitudinalSusceptibility::operator()(const Eigen::Vector3d& q, meV omega, eV U, eV J) 
{

	xuu = Xuu(q, omega);
	xdd = Xdd(q, omega);

	Eigen::MatrixXcd m11 = M11(q, omega, U, J);
	Eigen::MatrixXcd m12 = M12(q, omega, U, J);
	Eigen::MatrixXcd m21 = M21(q, omega, U, J);
	Eigen::MatrixXcd m22 = M22(q, omega, U, J);

	
	return (I + m22.inverse()*m21) * 
		   (I - (m11.inverse()*m12) * (m22.inverse()*m21)).inverse() *
		   (m11.inverse()*xuu) + (I + m11.inverse()*m12) +
		   (I - (m22.inverse()*m21)*(m11.inverse()*m12)).inverse() * 
		   (m22.inverse() * xdd);

}


Eigen::MatrixXcd 
LongitudinalSusceptibility::M11(const Eigen::Vector3d& q, meV omega, eV U, eV J) {
	return I + 4*xuu * Ubar(U, J);
}

Eigen::MatrixXcd 
LongitudinalSusceptibility::M12(const Eigen::Vector3d& q, meV omega, eV U, eV J) {
	return 2*xuu * Vbar(U, J);
}

Eigen::MatrixXcd 
LongitudinalSusceptibility::M21(const Eigen::Vector3d& q, meV omega, eV U, eV J) {
	return 2*xdd * Vbar(U, J);
}

Eigen::MatrixXcd 
LongitudinalSusceptibility::M22(const Eigen::Vector3d& q, meV omega, eV U, eV J) {
	return I + 4*xdd * Ubar(U, J);
}

Eigen::MatrixXd 
LongitudinalSusceptibility::Ubar(eV U, eV J) 
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
LongitudinalSusceptibility::Vbar(eV U, eV J) 
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
