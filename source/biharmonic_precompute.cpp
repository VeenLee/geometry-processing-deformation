#include "biharmonic_precompute.h"
#include <igl/min_quad_with_fixed.h>
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/invert_diag.h>
#include "../include/algorithm.h"

void biharmonic_precompute(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::VectorXi & b,
  igl::min_quad_with_fixed_data<double> & data)
{
#ifdef USE_ALGO_1
    // Find L, M
    Eigen::SparseMatrix<double> L, M;
    igl::cotmatrix(V, F, L);
    igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_DEFAULT, M);

    // Compute M_inverse
    Eigen::SparseMatrix<double> M_inverse(M.rows(), M.rows());
    Eigen::VectorXd M_inverse_dia;
    M_inverse_dia = Eigen::ArrayXd::Ones(M.rows()) / M.diagonal().array();
    for (int i = 0; i < M.rows(); i++) {
        M_inverse.insert(i, i) = M_inverse_dia[i];
    }

    // Compute Q:
    Eigen::SparseMatrix<double> Q, Aeq;
    Q = L.transpose() * M_inverse * L;

    // Compute data:
    igl::min_quad_with_fixed_precompute(Q, b, Aeq, false, data);
#endif

#ifdef USE_ALGO_2
    data.n = V.rows();
    Eigen::SparseMatrix<double> L;
    igl::cotmatrix(V, F, L);

    Eigen::SparseMatrix<double> M;
    igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_DEFAULT, M);

    Eigen::SparseMatrix<double> M_inverse;
    igl::invert_diag(M, M_inverse);

    Eigen::SparseMatrix<double> Q;
    Q = L.transpose() * M_inverse * L;

    igl::min_quad_with_fixed_precompute(Q, b, Eigen::SparseMatrix<double>(), true, data);
#endif

#ifdef USE_ALGO_3
    data.n = V.rows();
    Eigen::SparseMatrix<double> L(V.rows(), V.cols());
    igl::cotmatrix(V, F, L);
    Eigen::SparseMatrix<double> M(V.rows(), V.cols());
    igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_DEFAULT, M);

    // Calculate M^(-1) * L
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    solver.compute(M);
    Eigen::SparseMatrix<double> Q1 = solver.solve(L);

    // Calculate L^T * M^(-1) * L
    Eigen::SparseMatrix<double> Q = L.transpose() * Q1;

    // Precompute
    Eigen::SparseMatrix<double> B(Q.rows(), Q.cols());
    B.setZero();
    Eigen::SparseMatrix<double> Aeq;
    igl::min_quad_with_fixed_precompute(Q, b, Aeq, true, data);
#endif

#ifdef USE_ALGO_4
    Eigen::SparseMatrix<double> L, M;
    igl::cotmatrix(V, F, L);
    igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_DEFAULT, M);

    Eigen::SparseMatrix<double> M_i(M.rows(), M.cols());
    for (int i = 0; i < M.rows(); i++) {
        M_i.coeffRef(i, i) = 1.0 / M.coeff(i, i);
    }

    Eigen::SparseMatrix<double> Q = L.transpose() * M_i * L;
    igl::min_quad_with_fixed_precompute(Q, b, Eigen::SparseMatrix<double>(), false, data);
#endif

#ifdef USE_ALGO_5
    data.n = V.rows();

    Eigen::SparseMatrix<double> L;
    igl::cotmatrix(V, F, L);

    Eigen::SparseMatrix<double> M;
    igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_BARYCENTRIC, M);

    // Find the inverse of the Mass matrix. 
    // Since it is a diagonal sparse,
    // manually compute element by element and create M_inv
    Eigen::SparseMatrix<double> M_inv;
    M_inv.resize(M.rows(), M.cols());
    std::vector<Eigen::Triplet<double> > triplets;

    for (int j = 0; j < M.rows(); j++) {
        triplets.push_back(Eigen::Triplet<double>(j, j, 1.0 / M.coeff(j, j)));
    }

    M_inv.setFromTriplets(triplets.begin(), triplets.end());

    // multiply by 2.0 since the igl library uses a 1/2 in front of the
    // quadratic term
    Eigen::SparseMatrix<double> Q = 2.0 * L.transpose() * M_inv * L;

    Eigen::SparseMatrix<double> Aeq;
    //Aeq.resize(V.rows(), V.rows());

    // Initially planned using Aeq and Beq below. But using Z_b, Z_bc is much simpler.
    // // for all the handle points keep the diagonal entries in Aeq = 1
    // std::vector<Eigen::Triplet<double> > Aeq_triplets;
    // for (int j=0; j<b.size(); j++) {
    //   Aeq_triplets.push_back(Eigen::Triplet<double>(b(j), b(j), 1.0));
    // }
    // Aeq.setFromTriplets(Aeq_triplets.begin(), Aeq_triplets.end());

    igl::min_quad_with_fixed_precompute(Q, b, Aeq, false, data);
#endif

#ifdef USE_ALGO_6
    Eigen::SparseMatrix<double> L, M, Q, M_in, Aeq, I;
    igl::cotmatrix(V, F, L);
    igl::massmatrix(V, F, igl::MassMatrixType::MASSMATRIX_TYPE_VORONOI, M);

    //compute M-1
    I.resize(M.rows(), M.cols());
    I.setIdentity();
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::UpLoType::Lower | Eigen::UpLoType::Upper> solver;
    solver.compute(M);
    M_in = solver.solve(I);

    //Q = Lt*M-1*L
    Q = L.transpose() * M_in * L;

    igl::min_quad_with_fixed_precompute(Q, b, Aeq, false, data);
#endif
}
