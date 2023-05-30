#include "../include/arap_precompute.h"
#include <igl/min_quad_with_fixed.h>
#include <igl/arap_linear_block.h>
#include <igl/cotmatrix.h>

typedef Eigen::Triplet<double> T;

void arap_precompute(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::VectorXi & b,
  igl::min_quad_with_fixed_data<double> & data,
  Eigen::SparseMatrix<double> & K)
{
    //code 1
    {
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
    }


    //code 2
    {
        //// Construct L
        //Eigen::SparseMatrix<double> L;
        //igl::cotmatrix(V, F, L);

        //igl::min_quad_with_fixed_precompute(L, b, Eigen::SparseMatrix<double>(), false, data);

        ///*
        // * Note: When constructing K it’s easiest to iterate over all half-edges in the mesh
        // * (by iterating over all faces and then each of the three edges).
        // * Each half-edge ij contributes terms tying vi,vj to each of the (three)
        // * rotations Rk that apply against their difference.
        // */

        //std::vector<Eigen::Triplet<double>> tripletList;

        //// Construct K
        //K.resize(V.rows(), 3 * V.rows());

        //for (int f = 0; f < F.rows(); f++) {
        //    for (int k = 0; k < 3; k++) {
        //        for (int b = 0; b < 3; b++) {
        //            // edge 01
        //            tripletList.push_back(Eigen::Triplet<double>(F(f, 0), 3 * F(f, k) + b, L.coeff(F(f, 0), F(f, 1)) * (V(F(f, 0), b) - V(F(f, 1), b)) / 6));
        //            tripletList.push_back(Eigen::Triplet<double>(F(f, 1), 3 * F(f, k) + b, L.coeff(F(f, 0), F(f, 1)) * (-V(F(f, 0), b) + V(F(f, 1), b)) / 6));
        //            // edge 12
        //            tripletList.push_back(Eigen::Triplet<double>(F(f, 1), 3 * F(f, k) + b, L.coeff(F(f, 1), F(f, 2)) * (V(F(f, 1), b) - V(F(f, 2), b)) / 6));
        //            tripletList.push_back(Eigen::Triplet<double>(F(f, 2), 3 * F(f, k) + b, L.coeff(F(f, 1), F(f, 2)) * (-V(F(f, 1), b) + V(F(f, 2), b)) / 6));
        //            // edge 20
        //            tripletList.push_back(Eigen::Triplet<double>(F(f, 2), 3 * F(f, k) + b, L.coeff(F(f, 2), F(f, 0)) * (V(F(f, 2), b) - V(F(f, 0), b)) / 6));
        //            tripletList.push_back(Eigen::Triplet<double>(F(f, 0), 3 * F(f, k) + b, L.coeff(F(f, 2), F(f, 0)) * (-V(F(f, 2), b) + V(F(f, 0), b)) / 6));
        //        }
        //    }
        //}

        //K.setFromTriplets(tripletList.begin(), tripletList.end());
    }
}
