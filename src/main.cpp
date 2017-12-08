///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                          *** main.cpp ***                                 //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include "mpi.h"
#include <iostream>
#include "Eigmesh.hpp"
#include "Mesh.hpp"
#include "W90.hpp"
#include "StaticMatrix.hpp"
#include "BareSusceptibility.hpp"

int main(int argc, char *argv[])
{

	int world_rank, world_size;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

	BareSusceptibility X("data/input/reduced_klist");

	Eigen::Vector3d q1 = {0.380, 0.380, 0.380}; 
	Eigen::Vector3d q2 = {0.980, 0.690, 0.690}; 
	Eigen::Vector3d q3 = {0.000, 0.695, 0.695}; 
	Eigen::Vector3d q4 = {0.770, 0.260, 0.260}; 

	std::ofstream q1file("data/output/T2/q1/omega" + std::to_string(world_rank/2.0));
	q1file << X(q1, world_rank/2.0) << std::endl;
	q1file.close();

	if (world_rank == 0) std::cout << "finished q1" << std::endl;
		
	std::ofstream q2file("data/output/T2/q2/omega" + std::to_string(world_rank/2.0));
	q2file << X(q2, world_rank/2.0) << std::endl;
	q2file.close();

	if (world_rank == 0) std::cout << "finished q2" << std::endl;
	
	std::ofstream q3file("data/output/T2/q3/omega" + std::to_string(world_rank/2.0));
	q3file << X(q4, world_rank/2.0) << std::endl;
	q3file.close();

	if (world_rank == 0) std::cout << "finished q3" << std::endl;
	
	std::ofstream q4file("data/output/T2/q4/omega" + std::to_string(world_rank/2.0));
	q4file << X(q3, world_rank/2.0) << std::endl;
	q4file.close();
	
	if (world_rank == 0) std::cout << "finished q4" << std::endl;
	
	MPI_Finalize();	

	return 0;
}
