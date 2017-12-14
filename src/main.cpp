///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                          *** main.cpp ***                                 //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include <iostream>
#include "Eigmesh.hpp"
#include "Mesh.hpp"
#include "W90.hpp"
#include "BareSusceptibility.hpp"
#include "TransverseSusceptibility.hpp"
#include "LongitudinalSusceptibility.hpp"
#include "utils.hpp"


int main(int argc, char *argv[])
{
	
	eV U = 1;
	eV J = U / 10;
	meV omega_min = 0;
	meV omega_max = 20;
	meV omega_stp = 0.1;
	kelvin T = 2;

	std::string klist = "data/input/reduced.klist";
	std::string outdir = "data/output/T" + std::to_string(static_cast<int>(T)) + "/";
	
	LongitudinalSusceptibility Xzz(klist, T);
	TransverseSusceptibility Xpm(klist, T);

	Eigen::MatrixXcd xzz, xpm;	

	// Figure 2 Alekseev a), b), c), d)
	Eigen::Vector3d q1  = {0.000, 0.695, 0.695};  
	Eigen::Vector3d q2  = {0.980, 0.690, 0.690};  
	Eigen::Vector3d q3  = {0.760, 0.760, 0.760};  
	Eigen::Vector3d q4  = {0.770, 0.260, 0.260};  
	
	// Figure 3 Alekseev a), b)
	Eigen::Vector3d q5 = {0.380, 0.380, 0.380};  
	Eigen::Vector3d q6 = {0.490, 0.490, 0.490};  

	Eigen::Vector3d p00 = {M_PI, 0, 0};
	Eigen::Vector3d pp0 = {M_PI, M_PI, 0};
	Eigen::Vector3d ppp = {M_PI, M_PI , M_PI};

	std::ofstream realzz, realpm;
	std::string realzz_filename = "rpa_realzz_omega_sweep.dat";
	std::string realpm_filename = "rpa_realpm_omega_sweep.dat";

	std::ofstream imagzz, imagpm;
	std::string imagzz_filename = "rpa_imagzz_omega_sweep.dat";
	std::string imagpm_filename = "rpa_imagpm_omega_sweep.dat";

	for (meV omega = omega_min; omega <= omega_max; omega += omega_stp) {

		xzz = Xzz(ppp, omega, U, J);
		xpm = Xpm(ppp, omega, U, J);

		Eigen::ArrayXcd intra_xzz = intra_orbital_channel(xzz);
		Eigen::ArrayXcd intra_xpm = intra_orbital_channel(xpm);
		
		realzz.open(outdir + "ppp/" + realzz_filename, std::ios::app);
		realpm.open(outdir + "ppp/" + realpm_filename, std::ios::app);
		imagzz.open(outdir + "ppp/" + imagzz_filename, std::ios::app);
		imagpm.open(outdir + "ppp/" + imagpm_filename, std::ios::app);
		
		realzz << std::scientific;
		realzz << omega << "\t";
		realzz << xzz.real().sum() << "\t";
		realzz << xzz.real().trace() << "\t";
		realzz << inter_orbital_sum(xzz).real() << "\t";
		realzz << intra_xzz.real().sum() << "\t";
		realzz << intra_xzz(0).real() + intra_xzz(1).real() << "\t";
		realzz << intra_xzz.real().sum() - (intra_xzz(0).real() + intra_xzz(1).real()) << "\n";
	
		realpm << std::scientific;
		realpm << omega << "\t";
		realpm << xpm.real().sum() << "\t";
		realpm << xpm.real().trace() << "\t";
		realpm << inter_orbital_sum(xpm).real() << "\t";
		realpm << intra_xpm.real().sum() << "\t";
		realpm << intra_xpm(0).real() + intra_xpm(1).real() << "\t";
		realpm << intra_xpm.real().sum() - (intra_xpm(0).real() + intra_xpm(1).real()) << "\n";

		imagzz << std::scientific;
		imagzz << omega << "\t";
		imagzz << xzz.imag().sum() << "\t";
		imagzz << xzz.imag().trace() << "\t";
		imagzz << inter_orbital_sum(xzz).imag() << "\t";
		imagzz << intra_xzz.imag().sum() << "\t";
		imagzz << intra_xzz(0).imag() + intra_xzz(1).imag() << "\t";
		imagzz << intra_xzz.imag().sum() - (intra_xzz(0).imag() + intra_xzz(1).imag()) << "\n";
	
		imagpm << std::scientific;
		imagpm << omega << "\t";
		imagpm << xpm.imag().sum() << "\t";
		imagpm << xpm.imag().trace() << "\t";
		imagpm << inter_orbital_sum(xpm).imag() << "\t";
		imagpm << intra_xpm.imag().sum() << "\t";
		imagpm << intra_xpm(0).imag() + intra_xpm(1).imag() << "\t";
		imagpm << intra_xpm.imag().sum() - (intra_xpm(0).imag() + intra_xpm(1).imag()) << "\n";

		realzz.close();
		realpm.close();
		imagzz.close();
		imagpm.close();

		std::cout << 100.0 * omega / omega_max << "\n";
	}




	return 0;
}
