#pragma once

// #define BFF_USE_MKL 1
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <Eigen/SparseCholesky>

#if defined(BFF_USE_MKL)
#include <Eigen/PardisoSupport>
#endif

#include <algorithm>
#include <array>
#include <iostream>

// Eigen typdefs
namespace bff  {
using Vector = Eigen::Matrix<double,3,1>;
using DenseMatrix = Eigen::MatrixXd;
using SparseMatrix = Eigen::SparseMatrix<double>;

#if defined(BFF_USE_MKL)
// using SparseSolver = Eigen::PardisoLDLT<SparseMatrix>;
using SparseSolver = Eigen::PardisoLLT<SparseMatrix>;
#else
// using SparseSolver = Eigen::SimplicialLDLT<SparseMatrix>;
using SparseSolver = Eigen::SimplicialLLT<SparseMatrix>;
#endif
}

namespace bff  {

// General/IO utility
/*-------------------------------------------------------------------------------------*/
/*-------------------------------------------------------------------------------------*/
template<typename T>
inline void print(const T& a) {
    std::cout << a << '\n';
}
/*-------------------------------------------------------------------------------------*/
/*-------------------------------------------------------------------------------------*/

/*-------------------------------------------------------------------------------------*/
/*-------------------------------------------------------------------------------------*/
inline double dot(const Vector& u, const Vector& v)
{
    // return u.x*v.x + u.y*v.y + u.z*v.z;
    return u.dot(v);
}

inline Vector cross(const Vector& u, const Vector& v)
{
    // return Vector(u.y*v.z - u.z*v.y, u.z*v.x - u.x*v.z, u.x*v.y - u.y*v.x);
    return u.cross(v);
}
/*-------------------------------------------------------------------------------------*/
/*-------------------------------------------------------------------------------------*/


// Dense matrix utility
/*-------------------------------------------------------------------------------------*/
/*-------------------------------------------------------------------------------------*/
inline
DenseMatrix submatrix(const DenseMatrix& A, const std::vector<int>& is, const std::vector<int>& js = {0}) {
    DenseMatrix out(is.size(),js.size());
    for (size_t i=0; i<is.size(); ++i) {
        for (size_t j=0; j<js.size(); ++j) {
            out(i,j) = A(is[i],js[j]);
        }
    }
    return out;
}

inline DenseMatrix hcat(const DenseMatrix& A, const DenseMatrix& B)
{
    size_t m = A.rows();
    size_t n1 = A.cols();
    size_t n2 = B.cols();
    DenseMatrix C(m, n1 + n2);

    for (size_t i = 0; i < m; i++) {
        for (size_t j = 0; j < n1; j++) {
            C(i, j) = A(i, j);
        }

        for (size_t j = 0; j < n2; j++) {
            C(i, n1 + j) = B(i, j);
        }
    }

    return C;
}

inline DenseMatrix vcat(const DenseMatrix& A, const DenseMatrix& B)
{
    size_t m1 = A.rows();
    size_t m2 = B.rows();
    size_t n = A.cols();
    DenseMatrix C(m1 + m2, n);

    for (size_t j = 0; j < n; j++) {
        for (size_t i = 0; i < m1; i++) {
            C(i, j) = A(i, j);
        }

        for (size_t i = 0; i < m2; i++) {
            C(m1 + i, j) = B(i, j);
        }
    }

    return C;
}
/*-------------------------------------------------------------------------------------*/
/*-------------------------------------------------------------------------------------*/



// Sparse matrix utilities
/*-------------------------------------------------------------------------------------*/
/*-------------------------------------------------------------------------------------*/
inline SparseMatrix diag(const DenseMatrix& d)
{
    SparseMatrix sm(d.rows(), d.rows());

//#if USE_EIGEN_TRIPLET
#if 1
    for (int i = 0; i < d.rows(); i++) {
        sm.coeffRef(i,i) = d(i,0);
    }
#else
    using Triplet = Eigen::Triplet<double>;
    std::vector<Triplet> T;
    T.reserve(d.rows() + 1);
    for (int i = 0; i < d.rows(); i++) {
        T.push_back(Triplet(i, i, d(i,0)));
    }
    sm.setFromTriplets(T.begin(), T.end());
#endif
    sm.makeCompressed();
    return sm;
}


inline SparseMatrix submatrix(  const SparseMatrix& A,
                                const std::vector<int>& r,
                                const std::vector<int>& c) {
#if 0
     using Triplet = Eigen::Triplet<double>;
     std::vector<Triplet> T;
     T.reserve(A.nonZeros());
     double tol = 1e-14;

     for (int i=0; i<(int)r.size(); ++i) {
         for (int j=0; j<(int)c.size(); ++j) {
             //if (std::binary_search(IJ.begin(), IJ.end(), std::make_pair(i, j))) {
                 double entry = A.coeff(r[i], c[j]);
                 if (std::abs(entry) > tol) {
                     T.push_back(Triplet(i, j, entry));
                 }
             //}
         }
     }

    //for (int i = 0; i < A.outerSize(); i++) {
    //    for (typename SparseMatrix::InnerIterator it(A, i); it; ++it) {
    //        if (
    //            std::binary_search(r.begin(), r.end(), it.row()) &&
    //            std::binary_search(c.begin(), c.end(), it.col())
    //            ) {
    //            T.push_back(Triplet(it.row(), it.col(), it.value()));
    //        }
    //    }
    //}

    SparseMatrix out(r.size(), c.size());
    out.setFromTriplets(T.begin(), T.end());
    out.makeCompressed();
    return out;
#endif

    // Get the list of rows to delete
    std::vector<int> rows_to_remove;
    for (int i = 0; i < A.rows(); ++i) {
        auto res = std::find(std::begin(r), std::end(r), i);
        if (res == std::end(r)) {
            rows_to_remove.push_back(i);
        }
    }

    // Get the list of columns to delete
    std::vector<int> cols_to_remove;
    for (int i = 0; i < A.cols(); ++i) {
        auto res = std::find(std::begin(c), std::end(c), i);
        if (res == std::end(c)) {
            cols_to_remove.push_back(i);
        }
    }

    // Create a temporary non-Eigen triplet vector to store i,j,v
    // to be modified later and push_back all ijvs that need to remain
    std::vector<std::pair<std::array<Eigen::Index, 2>, double>> triplet;
    triplet.reserve(A.nonZeros());
    for (int i = 0; i < A.outerSize(); i++) {
        for (typename SparseMatrix::InnerIterator it(A, i); it; ++it) {
            if (std::binary_search(r.begin(), r.end(), it.row()) &&
                std::binary_search(c.begin(), c.end(), it.col())
                )
            {
                std::array<Eigen::Index, 2> cur = { it.row(), it.col() };
                triplet.push_back(std::make_pair(cur, it.value()));
            }
        }
    }

    std::vector<std::pair<std::array<Eigen::Index, 2>, double>> triplet_copy = triplet;
    // Map rows to ground
    for (int i = 0; i < rows_to_remove.size(); ++i) {
        for (int j = 0; j < triplet.size(); ++j) {
            if (triplet_copy[j].first[0] > rows_to_remove[i]) {
                triplet[j].first[0] -= 1;
            }
        }
    }
    // Map cols to ground
    for (int i = 0; i < cols_to_remove.size(); ++i) {
        for (int j = 0; j < triplet.size(); ++j) {
            if (triplet_copy[j].first[1] > cols_to_remove[i]) {
                triplet[j].first[1] -= 1;
            }
        }
    }

    // Create the final Eigen sparse matrix
    using Triplet = Eigen::Triplet<double>;
    std::vector<Triplet> T;
    T.reserve(A.nonZeros());

    for (int i = 0; i < triplet.size(); ++i) {
        T.push_back(Triplet((int)triplet[i].first[0], (int)triplet[i].first[1], triplet[i].second));
    }

    SparseMatrix out(r.size(), c.size());
    out.setFromTriplets(T.begin(), T.end());
    out.makeCompressed();

    return out;
}

} // namespace bff
