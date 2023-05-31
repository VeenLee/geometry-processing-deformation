#include "biharmonic_solve.h"
#include <igl/min_quad_with_fixed.h>
#include <igl/harmonic.h>
#include "../include/algorithm.h"

void biharmonic_solve(
  const igl::min_quad_with_fixed_data<double> & data,
  const Eigen::MatrixXd & bc,
  Eigen::MatrixXd & D)
{
#ifdef USE_ALGO_1
	Eigen::MatrixXd B(data.n, 3), Beq;
	B = Eigen::MatrixXd::Zero(data.n, 3);
	igl::min_quad_with_fixed_solve(data, B, bc, Beq, D);
#endif

#ifdef USE_ALGO_2
	const int size = data.n;
	const int dims = bc.cols();

	min_quad_with_fixed_solve(data, Eigen::MatrixXd::Zero(size, dims),
		bc, Eigen::MatrixXd::Zero(0, dims), D);
#endif

#ifdef USE_ALGO_3
	D = Eigen::MatrixXd::Zero(data.n, 3);
	Eigen::MatrixXd B(data.n, data.n);
	B.setZero();
	Eigen::MatrixXd Beq;
	igl::min_quad_with_fixed_solve(data, B, bc, Beq, D);
#endif

#ifdef USE_ALGO_4
	Eigen::MatrixXd B = Eigen::MatrixXd::Zero(data.n, 1);
	igl::min_quad_with_fixed_solve(data, B, bc, Eigen::MatrixXd(), D);
#endif

#ifdef USE_ALGO_5
	// Sparse type not working with igl.
  //Eigen::SparseMatrix<double> B, Beq;
  //B.resize(data.n, 3);
  //Beq.resize(data.n, 3);

  // Sparse isn't working. So use a usual matrix.
	Eigen::MatrixXd B, Beq;

	// B needs to be initialized for igl to work without assertion fail.
	B = Eigen::MatrixXd::Zero(data.n, 3);


	igl::min_quad_with_fixed_solve(data, B, bc, Beq, D);
#endif

#ifdef USE_ALGO_6
	Eigen::MatrixXd Beq, B;
	B.setZero(data.n, 3);
	igl::min_quad_with_fixed_solve(data, B, bc, Beq, D);
#endif
}
