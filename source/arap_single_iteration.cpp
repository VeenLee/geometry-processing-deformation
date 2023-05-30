#include "../include/arap_single_iteration.h"
#include <igl/polar_svd3x3.h>
#include <igl/min_quad_with_fixed.h>

void arap_single_iteration(
  const igl::min_quad_with_fixed_data<double> & data,
  const Eigen::SparseMatrix<double> & K,
  const Eigen::MatrixXd & bc,
  Eigen::MatrixXd & U)
{
    //code 1
    {
        int N = U.rows();

        // Construct C:
        Eigen::MatrixXd C(3 * N, 3);
        C = (U.transpose() * K).transpose();

        // Construct R:
        Eigen::MatrixXd R(3 * N, 3);
        for (int i = 0; i < N; i++) {
            Eigen::Matrix3d ck, rk;
            ck << C.row(i * 3),
                C.row(i * 3 + 1),
                C.row(i * 3 + 2);
            igl::polar_svd3x3(ck, rk);
            R.row(i * 3) = rk.row(0);
            R.row(i * 3 + 1) = rk.row(1);
            R.row(i * 3 + 2) = rk.row(2);
        }
        igl::min_quad_with_fixed_solve(data, K * R, bc,
            Eigen::MatrixXd(), U);
    }


    //code 2
    {
        //const int size = U.rows();
        //const int dims = U.cols();

        //const Eigen::MatrixXd C = K.transpose() * U;

        //assert(size * 3 == C.rows());

        //Eigen::MatrixXd R(C.rows(), C.cols());

        //// Compute R from K using closest rotation
        //for (int k = 0; k < size; k++) {
        //    Eigen::Matrix3d Ck = C.block<3, 3>(3 * k, 0);
        //    Eigen::Matrix3d Rk(3, 3);

        //    igl::polar_svd3x3(Ck, Rk);

        //    R.block<3, 3>(3 * k, 0) = Rk;
        //}

        //min_quad_with_fixed_solve(data, K * R, bc, Eigen::MatrixXd::Zero(0, dims), U);
    }
}
