#include "../include/arap_single_iteration.h"
#include <igl/polar_svd3x3.h>
#include <igl/min_quad_with_fixed.h>
#include "../include/algorithm.h"

void arap_single_iteration(
  const igl::min_quad_with_fixed_data<double> & data,
  const Eigen::SparseMatrix<double> & K,
  const Eigen::MatrixXd & bc,
  Eigen::MatrixXd & U)
{
#ifdef USE_ALGO_1
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
    igl::min_quad_with_fixed_solve(data, K * R, bc, Eigen::MatrixXd(), U);
#endif

#ifdef USE_ALGO_2
    const int size = U.rows();
    const int dims = U.cols();

    const Eigen::MatrixXd C = K.transpose() * U;

    assert(size * 3 == C.rows());

    Eigen::MatrixXd R(C.rows(), C.cols());

    // Compute R from K using closest rotation
    for (int k = 0; k < size; k++) {
        Eigen::Matrix3d Ck = C.block<3, 3>(3 * k, 0);
        Eigen::Matrix3d Rk(3, 3);

        igl::polar_svd3x3(Ck, Rk);

        R.block<3, 3>(3 * k, 0) = Rk;
    }

    min_quad_with_fixed_solve(data, K * R, bc, Eigen::MatrixXd::Zero(0, dims), U);
#endif

#ifdef USE_ALGO_3
    // Local step
    Eigen::MatrixXd R(3 * K.rows(), 3);
    const Eigen::MatrixXd C = K.transpose() * U;
    // Update each Rk
    Eigen::Matrix<double, 3, 3> Ck;
    Eigen::Matrix<double, 3, 3> Rk;
    for (int k = 0; k < K.rows(); k++) {
        Ck = C.block(3 * k, 0, 3, 3);
        igl::polar_svd3x3(Ck, Rk);
        R.block(3 * k, 0, 3, 3) = Rk;
    }

    // Global step
    Eigen::MatrixXd B = K * R;
    Eigen::MatrixXd Beq;
    igl::min_quad_with_fixed_solve(data, B, bc, Beq, U);
#endif

#ifdef USE_ALGO_4
    Eigen::MatrixXd C = K.transpose() * U;
    Eigen::MatrixXd R(3 * U.rows(), 3);

    for (int k = 0; k < U.rows(); k++) {
        Eigen::Matrix3d Ck, Rk;
        for (int i = 0; i < 3; i++) {
            Ck.row(i) = C.row(3 * k + i);
        }
        igl::polar_svd3x3(Ck, Rk);
        for (int i = 0; i < 3; i++) {
            R.row(3 * k + i) = Rk.row(i);
        }
    }
    igl::min_quad_with_fixed_solve(data, K * R, bc, Eigen::MatrixXd(), U);
#endif

#ifdef USE_ALGO_5
    // First do the local step
  // For local step, we need to minimize tr(V'KR) w.r.t R.
  // Note the matrix property: tr(A'B) = <A, B>_Frobenius! So trace is related to Frobenious dot product!
  // Therefore, min_R tr(V'KR) is same as min_R <K'V, R>_Frobenius => same as HW02 (finding closest rotation).

  // Assume the initalization of the final vertices to be same as original vertices U
  // C = K'V

    Eigen::MatrixXd C = K.transpose() * U;

    // For the rotations
    Eigen::MatrixXd R;
    R.resize(3 * U.rows(), 3);

    // Loop in steps of k
    for (int k = 0; k < U.rows(); k++) {

        // polar_svd3x3 doesnt work with dynamic matrices of 3x3 size
        Eigen::Matrix3d R_k;
        Eigen::Matrix3d C_k;
        C_k = C.block(k * 3, 0, 3, 3);

        // get the closest rotation
        igl::polar_svd3x3(C_k, R_k);

        // update R
        R.block(k * 3, 0, 3, 3) = R_k;
    }

    // global step
    // minimize tr (V'LV + V'KR) w.r.t V
    // Compare this with min_quad_with_solve:
    // B = KR; Z = V; 
    Eigen::MatrixXd B = K * R;

    // dummy
    Eigen::MatrixXd Beq;

    // global step, replace U
    igl::min_quad_with_fixed_solve(data, B, bc, Beq, U);
#endif

#ifdef USE_ALGO_6
    Eigen::MatrixXd C = (U.transpose() * K).transpose();
    Eigen::MatrixXd R(3 * U.rows(), 3);

    for (int i = 0; i < U.rows(); i++)
    {
        Eigen::Matrix3d R_block;
        Eigen::Matrix3d Ck = C.block(i * 3, 0, 3, 3);
        igl::polar_svd3x3(Ck, R_block);
        R.block(i * 3, 0, 3, 3) = R_block;
    }

    Eigen::MatrixXd Beq;
    Eigen::MatrixXd B = K * R;
    igl::min_quad_with_fixed_solve(data, B, bc, Beq, U);
#endif
}
