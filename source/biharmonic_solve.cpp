#include "biharmonic_solve.h"
#include <igl/min_quad_with_fixed.h>
#include <igl/harmonic.h>
#include "../include/algorithm.h"

void biharmonic_solve(
  const igl::min_quad_with_fixed_data<double> & data,
  const Eigen::MatrixXd & bc,
  Eigen::MatrixXd & D)
{
	//code 1
	{
		Eigen::MatrixXd B(data.n, 3), Beq;
		B = Eigen::MatrixXd::Zero(data.n, 3);
		igl::min_quad_with_fixed_solve(data, B, bc, Beq, D);
	}


	//code 2
	{
		//const int size = data.n;
		//const int dims = bc.cols();

		//min_quad_with_fixed_solve(data, Eigen::MatrixXd::Zero(size, dims),
		//	bc, Eigen::MatrixXd::Zero(0, dims), D);
	}
}
