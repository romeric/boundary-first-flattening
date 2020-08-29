#pragma once

// #define BFF_USE_MKL 1
// #define BFF_USE_CHOLMOD 1
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <Eigen/SparseCholesky>
#if defined(BFF_USE_MKL)
#include <Eigen/PardisoSupport>
#endif
#if defined(BFF_USE_CHOLMOD)
#include <Eigen/CholmodSupport>
#endif

#include <algorithm>
#include <array>
#include <iostream>


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

// Eigen typdefs
/*-------------------------------------------------------------------------------------*/
/*-------------------------------------------------------------------------------------*/
using Vector = Eigen::Matrix<double,3,1>;
using DenseMatrix = Eigen::MatrixXd;
using SparseMatrix = Eigen::SparseMatrix<double>;

#if !defined(BFF_USE_MKL) && !defined(BFF_USE_CHOLMOD)
using SparseSolver = Eigen::SimplicialLLT<SparseMatrix>;
#elif defined(BFF_USE_MKL) && !defined(BFF_USE_CHOLMOD)
using SparseSolver = Eigen::PardisoLLT<SparseMatrix>;
#elif !defined(BFF_USE_MKL) && defined(BFF_USE_CHOLMOD)
using SparseSolver = Eigen::CholmodDecomposition<SparseMatrix>;
#else
#error CONFLICTING SPARSE SOLVERS SELECTED - CHOOSE ONE
#endif
/*-------------------------------------------------------------------------------------*/
/*-------------------------------------------------------------------------------------*/

// Vector utility
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

// submatrix extraction when r==c
inline SparseMatrix submatrix(  const SparseMatrix& A,
                                const std::vector<int>& r) {

    if ((int)r.size()==A.rows() && (int)r.size()==A.cols()) return A;

#if 0
    using Triplet = Eigen::Triplet<double>;
    std::vector<Triplet> T;
    T.reserve(A.nonZeros());
    double tol = 1e-14;

    for (int i=0; i<(int)r.size(); ++i) {
        for (int j=0; j<(int)c.size(); ++j) {
            double entry = A.coeff(r[i], c[j]);
            if (std::abs(entry) > tol) {
             T.push_back(Triplet(i, j, entry));
            }
        }
    }

    SparseMatrix out(r.size(), c.size());
    out.setFromTriplets(T.begin(), T.end());
    out.makeCompressed();
    return out;
#endif

    // Get row mapping
    int sum = 0;
    std::vector<int> mappedRows;
    mappedRows.reserve(A.rows());
    for (int i = 0; i < A.rows(); ++i) {
        auto res = std::binary_search(std::begin(r), std::end(r), i);
        if (!res) {
            sum += 1;
        }
        mappedRows.push_back(i - sum);
    }

    // Create the final Eigen sparse matrix
    using Triplet = Eigen::Triplet<double>;
    std::vector<Triplet> T;
    T.reserve(A.nonZeros());

    for (int i = 0; i < A.outerSize(); i++) {
        for (typename SparseMatrix::InnerIterator it(A, i); it; ++it) {
            if (std::binary_search(r.begin(), r.end(), it.row()) &&
                std::binary_search(r.begin(), r.end(), it.col())
                )
            {
                const auto mit = mappedRows[it.row()];
                const auto mjt = mappedRows[it.col()];
                T.push_back(Triplet(mit, mjt, it.value()));
            }
        }
    }

    SparseMatrix out(r.size(), r.size());
    out.setFromTriplets(T.begin(), T.end());
    out.makeCompressed();

    return out;
}


// submatrix extraction when r!=c
inline SparseMatrix submatrix(  const SparseMatrix& A,
                                const std::vector<int>& r,
                                const std::vector<int>& c) {

    if ((int)r.size()==A.rows() && (int)c.size()==A.cols()) return A;

#if 0
    using Triplet = Eigen::Triplet<double>;
    std::vector<Triplet> T;
    T.reserve(A.nonZeros());
    double tol = 1e-14;

    for (int i=0; i<(int)r.size(); ++i) {
        for (int j=0; j<(int)c.size(); ++j) {
            double entry = A.coeff(r[i], c[j]);
            if (std::abs(entry) > tol) {
             T.push_back(Triplet(i, j, entry));
            }
        }
    }

    SparseMatrix out(r.size(), c.size());
    out.setFromTriplets(T.begin(), T.end());
    out.makeCompressed();
    return out;
#endif

    // Get row mapping
    int sum = 0;
    std::vector<int> mappedRows;
    mappedRows.reserve(A.rows());
    for (int i = 0; i < A.rows(); ++i) {
        auto res = std::binary_search(std::begin(r), std::end(r), i);
        if (!res) {
            sum += 1;
        }
        mappedRows.push_back(i - sum);
    }

    // Get col mapping
    sum = 0;
    std::vector<int> mappedCols;
    mappedCols.reserve(A.rows());
    for (int i = 0; i < A.cols(); ++i) {
        auto res = std::binary_search(std::begin(c), std::end(c), i);
        if (!res) {
            sum += 1;
        }
        mappedCols.push_back(i - sum);
    }

    // Create the final Eigen sparse matrix
    using Triplet = Eigen::Triplet<double>;
    std::vector<Triplet> T;
    T.reserve(A.nonZeros());

    for (int i = 0; i < A.outerSize(); i++) {
        for (typename SparseMatrix::InnerIterator it(A, i); it; ++it) {
            if (std::binary_search(r.begin(), r.end(), it.row()) &&
                std::binary_search(c.begin(), c.end(), it.col())
                )
            {
                const auto mit = mappedRows[it.row()];
                const auto mjt = mappedCols[it.col()];
                T.push_back(Triplet(mit, mjt, it.value()));
            }
        }
    }

    SparseMatrix out(r.size(), c.size());
    out.setFromTriplets(T.begin(), T.end());
    out.makeCompressed();

    return out;
}

} // namespace bff
