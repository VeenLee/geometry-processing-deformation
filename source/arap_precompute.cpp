#include "../include/arap_precompute.h"
#include <igl/min_quad_with_fixed.h>
#include <igl/arap_linear_block.h>
#include <igl/cotmatrix.h>
#include "../include/algorithm.h"

typedef Eigen::Triplet<double> T;

void arap_precompute(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::VectorXi & b,
  igl::min_quad_with_fixed_data<double> & data,
  Eigen::SparseMatrix<double> & K)
{
#ifdef USE_ALGO_1
    // Compute L:
    Eigen::SparseMatrix<double> L;
    igl::cotmatrix(V, F, L);

    // Compute K:
    std::vector<T> triplet_list;
    for (int i = 0; i < F.rows(); i++) {
        for (int j = 0; j < 3; j++) {
            // Find vi, and vj:
            int vi = F(i, j);
            int vj = F(i, (j + 1) % 3);
            Eigen::Vector3d Vi = V.row(vi);
            Eigen::Vector3d Vj = V.row(vj);

            // Compute eij:
            Eigen::Vector3d eij = L.coeff(vi, vj) * (Vi - Vj);

            // K could be three possible vertex:
            for (int k = 0; k < 3; k++) {
                int vk = F(i, k);
                // Loop over all 3 indexes:
                for (int l = 0; l < 3; l++) {
                    triplet_list.push_back(T(vi, 3 * vk + l, eij[l] / 6.0));
                    triplet_list.push_back(T(vj, 3 * vk + l, -eij[l] / 6.0));
                }
            }
        }
    }

    K.resize(V.rows(), 3 * V.rows());
    K.setFromTriplets(triplet_list.begin(), triplet_list.end());
    igl::min_quad_with_fixed_precompute(L, b, Eigen::SparseMatrix<double>(), false, data);
#endif

#ifdef USE_ALGO_2
    // Construct L
    Eigen::SparseMatrix<double> L;
    igl::cotmatrix(V, F, L);

    igl::min_quad_with_fixed_precompute(L, b, Eigen::SparseMatrix<double>(), false, data);

    /*
        * Note: When constructing K it’s easiest to iterate over all half-edges in the mesh
        * (by iterating over all faces and then each of the three edges).
        * Each half-edge ij contributes terms tying vi,vj to each of the (three)
        * rotations Rk that apply against their difference.
        */

    std::vector<Eigen::Triplet<double>> tripletList;

    // Construct K
    K.resize(V.rows(), 3 * V.rows());

    for (int f = 0; f < F.rows(); f++) {
        for (int k = 0; k < 3; k++) {
            for (int b = 0; b < 3; b++) {
                // edge 01
                tripletList.push_back(Eigen::Triplet<double>(F(f, 0), 3 * F(f, k) + b, L.coeff(F(f, 0), F(f, 1)) * (V(F(f, 0), b) - V(F(f, 1), b)) / 6));
                tripletList.push_back(Eigen::Triplet<double>(F(f, 1), 3 * F(f, k) + b, L.coeff(F(f, 0), F(f, 1)) * (-V(F(f, 0), b) + V(F(f, 1), b)) / 6));
                // edge 12
                tripletList.push_back(Eigen::Triplet<double>(F(f, 1), 3 * F(f, k) + b, L.coeff(F(f, 1), F(f, 2)) * (V(F(f, 1), b) - V(F(f, 2), b)) / 6));
                tripletList.push_back(Eigen::Triplet<double>(F(f, 2), 3 * F(f, k) + b, L.coeff(F(f, 1), F(f, 2)) * (-V(F(f, 1), b) + V(F(f, 2), b)) / 6));
                // edge 20
                tripletList.push_back(Eigen::Triplet<double>(F(f, 2), 3 * F(f, k) + b, L.coeff(F(f, 2), F(f, 0)) * (V(F(f, 2), b) - V(F(f, 0), b)) / 6));
                tripletList.push_back(Eigen::Triplet<double>(F(f, 0), 3 * F(f, k) + b, L.coeff(F(f, 2), F(f, 0)) * (-V(F(f, 2), b) + V(F(f, 0), b)) / 6));
            }
        }
    }

    K.setFromTriplets(tripletList.begin(), tripletList.end());
#endif

#ifdef USE_ALGO_3
    // Compute L
    Eigen::SparseMatrix<double> L(V.rows(), V.rows());
    igl::cotmatrix(V, F, L);
    // Precompute for the global step
    Eigen::SparseMatrix<double> Aeq;
    igl::min_quad_with_fixed_precompute(L, b, Aeq, false, data);

    // Prepare matrix K
    K = Eigen::SparseMatrix<double>(V.rows(), 3 * V.rows());
    K.setZero();
    Eigen::MatrixXd C(F.rows(), 3);
    igl::cotmatrix_entries(V, F, C);
    std::vector<Eigen::Triplet<double>> triplets;
    for (int f = 0; f < F.rows(); f++) {
        // Looping over faces
        int i, j;
        // Looping over edges of the f-th face
        i = 1; j = 2;
        for (int k = 0; k < 3; k++) {
            // Looping over vertices
            for (int beta = 0; beta < 3; beta++) {
                triplets.push_back({ F(f,i), 3 * F(f,k) + beta,
                     1. / 6. * (V(F(f,i), beta) - V(F(f,j), beta)) * L.coeff(F(f, i), F(f, j)) });
                triplets.push_back({ F(f,j), 3 * F(f,k) + beta,
                    -1. / 6. * (V(F(f,i), beta) - V(F(f,j), beta)) * L.coeff(F(f, i), F(f, j)) });
            }
        }
        i = 2; j = 0;
        for (int k = 0; k < 3; k++) {
            for (int beta = 0; beta < 3; beta++) {
                triplets.push_back({ F(f,i), 3 * F(f,k) + beta,
                     1. / 6. * (V(F(f,i), beta) - V(F(f,j), beta)) * L.coeff(F(f, i), F(f, j)) });
                triplets.push_back({ F(f,j), 3 * F(f,k) + beta,
                    -1. / 6. * (V(F(f,i), beta) - V(F(f,j), beta)) * L.coeff(F(f, i), F(f, j)) });
            }
        }
        i = 0; j = 1;
        for (int k = 0; k < 3; k++) {
            for (int beta = 0; beta < 3; beta++) {
                triplets.push_back({ F(f,i), 3 * F(f,k) + beta,
                     1. / 6. * (V(F(f,i), beta) - V(F(f,j), beta)) * L.coeff(F(f, i), F(f, j)) });
                triplets.push_back({ F(f,j), 3 * F(f,k) + beta,
                    -1. / 6. * (V(F(f,i), beta) - V(F(f,j), beta)) * L.coeff(F(f, i), F(f, j)) });
            }
        }
    }
    K.setFromTriplets(triplets.begin(), triplets.end());
#endif

#ifdef USE_ALGO_4
    Eigen::SparseMatrix<double> L;
    igl::cotmatrix(V, F, L);
    igl::min_quad_with_fixed_precompute(L, b, Eigen::SparseMatrix<double>(), false, data);

    Eigen::MatrixXd C(F.rows(), 3);
    igl::cotmatrix_entries(V, F, C);

    std::vector<Eigen::Triplet<double>> triplets;

    K.resize(V.rows(), 3 * V.rows());
    for (int f = 0; f < F.rows(); f++) {
        Eigen::RowVector3i face = F.row(f);
        for (int x = 0; x < 3; x++) {
            int i = face((x + 1) % 3);
            int j = face((x + 2) % 3);
            Eigen::RowVector3d edge = V.row(i) - V.row(j);
            Eigen::Vector3d eij = C(f, x) / 3.0 * edge;

            for (int a = 0; a < 3; a++) {
                int k = face((x + a) % 3);
                for (int b = 0; b < 3; b++) {
                    triplets.emplace_back(Eigen::Triplet<double>(i, 3 * k + b, eij(b)));
                    triplets.emplace_back(Eigen::Triplet<double>(j, 3 * k + b, -eij(b)));
                }
            }
        }
    }
    K.setFromTriplets(triplets.begin(), triplets.end());
#endif

#ifdef USE_ALGO_5
    // Construct K
    std::vector<Eigen::Triplet<double> > triplets;

    // get the cotmatrix
    Eigen::SparseMatrix<double> L;
    igl::cotmatrix(V, F, L);

    // iterate over faces
    for (int f = 0; f < F.rows(); f++) {

        // iterate over 3 edges of a face using modulo
        for (int e = 0; e < 3; e++) {
            int v_i_idx = F(f, e);
            int v_j_idx = F(f, (e + 1) % 3);

            Eigen::Vector3d v_i = V.row(v_i_idx);
            Eigen::Vector3d v_j = V.row(v_j_idx);

            // v_i -> v_j represents edge e_ij
            Eigen::Vector3d e_ij = (v_i - v_j) * L.coeff(v_i_idx, v_j_idx);

            // loop over 3 dimensions x, y, z (beta)
            for (int beta = 0; beta < 3; beta++) {
                // index k can be index v_i or v_j or the third vertex

                // case 1: k is i 
                int k_idx = v_i_idx;
                triplets.push_back(Eigen::Triplet<double>(v_i_idx, 3 * k_idx + beta, (1.0 / 6.0) * e_ij(beta)));
                triplets.push_back(Eigen::Triplet<double>(v_j_idx, 3 * k_idx + beta, (-1.0 / 6.0) * e_ij(beta)));

                // case 2: k is j
                k_idx = v_j_idx;
                triplets.push_back(Eigen::Triplet<double>(v_i_idx, 3 * k_idx + beta, (1.0 / 6.0) * e_ij(beta)));
                triplets.push_back(Eigen::Triplet<double>(v_j_idx, 3 * k_idx + beta, (-1.0 / 6.0) * e_ij(beta)));

                // case 3: k i the third vertex
                k_idx = F(f, (e + 2) % 3);
                triplets.push_back(Eigen::Triplet<double>(v_i_idx, 3 * k_idx + beta, (1.0 / 6.0) * e_ij(beta)));
                triplets.push_back(Eigen::Triplet<double>(v_j_idx, 3 * k_idx + beta, (-1.0 / 6.0) * e_ij(beta)));

            }

        }
    }

    K.resize(V.rows(), 3 * V.rows());
    K.setFromTriplets(triplets.begin(), triplets.end());

    // data matrix precomputed similar to biharmonic_precompute
    // need to minimize_D tr(D'LD + D'KR)

    // not using this 
    Eigen::SparseMatrix<double> Aeq;

    Eigen::SparseMatrix<double> L_times_2 = 2.0 * L;

    // pd should be false else numerical issue as L is not positive definite
    igl::min_quad_with_fixed_precompute(L_times_2, b, Aeq, false, data);
#endif

#ifdef USE_ALGO_6
    Eigen::SparseMatrix<double> A, Aeq;
    igl::cotmatrix(V, F, A);
    igl::min_quad_with_fixed_precompute(A, b, Aeq, false, data);

    //scale by 3 result too large
    A = A / 6.0;

    std::vector<T> Tlist;

    Tlist.reserve(F.rows() * 27 * 2);

    for (int idx = 0; idx < F.rows(); idx++)
    {
        for (int pos = 0; pos < 3; pos++)
        {
            int i = F(idx, pos % 3);
            int j = F(idx, (pos + 1) % 3);

            Eigen::RowVector3d u = (V.row(i) - V.row(j)) * A.coeff(i, j);

            for (int x = 0; x < 3; x++)
            {
                int k = F(idx, (pos + x) % 3);

                for (int y = 0; y < 3; y++)
                {
                    Tlist.push_back(T(i, 3 * k + y, u(y)));
                    Tlist.push_back(T(j, 3 * k + y, -u(y)));
                }
            }
        }
    }

    K.resize(V.rows(), 3 * V.rows());
    K.setFromTriplets(Tlist.begin(), Tlist.end());
#endif
}
