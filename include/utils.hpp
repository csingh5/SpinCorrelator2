///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                          *** utils.hpp ***                                //
//                                                                           //
//                                                                           //
// created November 29, 2017                                                 //
// copyright Christopher N. Singh Binghamton University Physics              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#ifndef utils_hpp
#define utils_hpp

#include "Eigen/Dense"
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <cassert>


Eigen::ArrayXcd 
intra_orbital_channel(const Eigen::MatrixXcd& X) 
{
	
	int L12 = 0, L34 = 0, numorbs = 9, i = 0;

	Eigen::ArrayXcd value = Eigen::ArrayXcd::Zero(numorbs);

	for (int L1 = 0; L1 < numorbs; ++L1) {
		for (int L2 = 0; L2 < numorbs; ++L2) { 
			for (int L3 = 0; L3 < numorbs; ++L3) { 
				for (int L4 = 0; L4 < numorbs; ++L4) {
					if (L1 == L2 && L2 == L3 && L3 == L4) value(i++) = X(L12, L34);	
					++L34;
				}
			}
			L34 = 0;
			++L12;
		}
	}

	return value;
}

std::complex<double>
inter_orbital_sum(const Eigen::MatrixXcd& X) 
{
	return X.sum() - intra_orbital_channel(X).sum();
}




Eigen::MatrixXcd toEigenMatrix(std::string matrix_file) {
	std::ifstream file(matrix_file);
	assert(file.is_open());

	std::complex<double> value;
	std::vector<std::complex<double>> data;

	Eigen::MatrixXcd mat;

	while (file >> value) data.push_back(value);	

	int dim = static_cast<int> (sqrt(data.size()));

	mat = Eigen::Map<Eigen::MatrixXcd> (&data[0], dim, dim);
	
	return mat;	
}


void count_num_negatives(std::string folder) {

	Eigen::MatrixXcd m = toEigenMatrix(folder + "omega" + std::to_string(1.0));	

	int c = 1, n = 0;
	int dimension = static_cast<int> (std::sqrt(m.size()));

	for (double i = 0; i < 20; i += 0.5) {
		c++;
		m = toEigenMatrix(folder + "omega" + std::to_string(i));
		for (int j = 0; j < dimension; ++j) {
			for (int k = 0; k < dimension; ++k) {
				if (m(j, k).imag() < 0) ++n;
			}
		}
	}

	n *= 100.0 / (c*dimension*dimension);

	std::cout << n << "% negative in all elements over omega sweep" << std::endl;

}

void count_diag_negatives(std::string folder) {

	Eigen::MatrixXcd m = toEigenMatrix(folder + "omega" + std::to_string(1.0));	

	int c = 1, n = 0;
	int dimension = static_cast<int> (std::sqrt(m.size()));

	for (double i = 0; i < 20; i += 0.5) {
		c++;
		m = toEigenMatrix(folder + "omega" + std::to_string(i));
		for (int j = 0; j < dimension; ++j) {
			if (m.diagonal().imag()(j) < 0) ++n;
		}
	}

	n *= 100.0 / (c*dimension);

	std::cout << n << "% negative in diagonal elements over omega sweep" << std::endl;

}





#endif 
