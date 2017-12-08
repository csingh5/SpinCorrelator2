///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                 *** BareSusceptibility.hpp ***                            //
//                                                                           //
// bare spin susceptibility functor... state is eigemesh at k                //
//                                                                           //
// created November 29, 2017                                                 //
// copyright Christopher N. Singh Binghamton University Physics              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#ifndef BareSusceptibility_hpp
#define BareSusceptibility_hpp

#include "StaticMatrix.hpp"


class BareSusceptibility : public StaticMatrix
{
public:
	BareSusceptibility(std::string klist_file) : StaticMatrix(klist_file) { }

	
Eigen::MatrixXcd operator()(const Eigen::Vector3d& q, meV omega);

private:

std::complex<double> matrix_element(int L1, int L2, int L3, int L4, meV omega);


};


#endif /* BareSusceptibility_hpp */
