///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//              *** LongitudinalSusceptibility.hpp ***                       //
//                                                                           //
// RPA spin susceptibility functor... state is eigemesh at k                 //
//                                                                           //
// created December 9, 2017                                                  //
// copyright Christopher N. Singh Binghamton University Physics              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#ifndef LongitudinalSusceptibility_hpp
#define LongitudinalSusceptibility_hpp

#include "BareSusceptibility.hpp"

class LongitudinalSusceptibility 
{
public:
	LongitudinalSusceptibility(std::string klist_file, kelvin T);

	Eigen::MatrixXcd operator()(const Eigen::Vector3d& q, meV omega, eV U, eV J);

private:
	BareSusceptibility Xuu;
	BareSusceptibility Xdd;
	Eigen::MatrixXcd xuu, xdd;
	int numorbs = 9;
	int numorbs2 = numorbs * numorbs;
	
	Eigen::MatrixXd I = Eigen::MatrixXd::Identity(numorbs2, numorbs2);
	
	Eigen::MatrixXcd M11(const Eigen::Vector3d& q, meV omega, eV U, eV J);

	Eigen::MatrixXcd M12(const Eigen::Vector3d& q, meV omega, eV U, eV J);

	Eigen::MatrixXcd M21(const Eigen::Vector3d& q, meV omega, eV U, eV J);

	Eigen::MatrixXcd M22(const Eigen::Vector3d& q, meV omega, eV U, eV J);

	Eigen::MatrixXd Ubar(eV U, eV J);

	Eigen::MatrixXd Vbar(eV U, eV J);

};





#endif /* LongitudinalSusceptibility_hpp */
