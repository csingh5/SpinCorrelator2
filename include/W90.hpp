///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                          *** w90.hpp ***                                  //
//                                                                           //
// class to manage methods relating to wannier90 hamiltonian                 //
// templatized over dimension of the mesh                                    //
//                                                                           //
// created November 27, 2017                                                 //
// copyright Christopher N. Singh Binghamton University Physics              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#ifndef W90_hpp
#define W90_hpp

#include <cassert>
#include <Eigen/Dense>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>

class W90
{
public:
	W90(std::string hr_file);

	Eigen::MatrixXcd fft(const Eigen::Vector3d& k);

private:
	Eigen::ArrayXXd hr;
	
	inline std::complex<double> tij(int nn, const Eigen::Vector3d& k);
};


#endif /* W90_hpp */
