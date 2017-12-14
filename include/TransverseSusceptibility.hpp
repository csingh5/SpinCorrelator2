///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                *** TransverseSusceptibility.hpp ***                       //
//                                                                           //
// RPA spin susceptibility functor... state is eigemesh at k                 //
//                                                                           //
// created December 9, 2017                                                  //
// copyright Christopher N. Singh Binghamton University Physics              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#ifndef TransverseSusceptibility_hpp
#define TransverseSusceptibility_hpp

#include "BareSusceptibility.hpp"

class TransverseSusceptibility 
{
public:
	TransverseSusceptibility(std::string klist_file, kelvin T);

	Eigen::MatrixXcd operator()(const Eigen::Vector3d& q, meV omega, eV U, eV J);

private:
	BareSusceptibility Xud;
	Eigen::MatrixXcd xud;
	int numorbs = 9;
	int numorbs2 = numorbs * numorbs;
	
	Eigen::MatrixXd I = Eigen::MatrixXd::Identity(numorbs2, numorbs2);
	
	Eigen::MatrixXd Ubar(eV U, eV J);

	Eigen::MatrixXd Vbar(eV U, eV J);

};





#endif /* TransverseSusceptibility_hpp */
