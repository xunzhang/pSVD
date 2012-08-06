/**
 * @file   gebrd.hpp
 * @author Changsheng Jiang <jiangzuoyan@gmail.com>
 * @date   Thu Oct 27 15:29:05 2011
 *
 * @brief bidiagonalization
 *
 *
 */
#ifndef FILE_8180cdd3_0cd2_4606_8640_9c0792474dde_H
#define FILE_8180cdd3_0cd2_4606_8640_9c0792474dde_H
#include "douban/linalg/householder.hpp"
#include "douban/linalg/gemv_ax.hpp"
#include "douban/matvec/vec_kahan.hpp"
#include "douban/utility/kahan_sum.hpp"
#include "douban/matvec/io.hpp"

namespace douban {
namespace linalg {

template <class Mat, class Vec>
void orth(Mat && mat, int first, int last, Vec && vec,
          bool unit=false, double alpha=.8) {
#if 1
  vec_container<double> p(last - first);
  vec_container<double> d(unit ? 0 : last - first);
  auto M = mat_cols(mat, first, last);
  if (!unit) {
    for (int i = 0; i < M.dim1(); ++i) d[i] = M.col(i).square_sum();
  }
  auto n = vec.square_sum();
  while (true) {
    if (unit)
      p = gemv(M.trans(), vec);
    else
      p = gemv(M.trans(), vec) / d;
    vec = vec - gemv(M, p);
    auto nn = n - p.square_sum();
    if (nn >= n * alpha || first + 1 >= last) break;
    n = nn;
  }
#else
  auto n = vec.square_sum();
  while (true) {
    for (int j = first; j < last; ++j) {
      auto c = vec.inner_prod(mat.col(j));
      if (!unit) c /= mat.col(j).square_sum();
      vec += -c * mat.col(j);
    }
    auto nn = vec.square_sum();
    if (nn >= n * alpha || first + 1 >= last) break;
    n = nn;
  }
#endif
}

template <class Mat, class Vec>
void orth(Mat && mat, Vec && vec, bool unit=false) {
  orth(mat, 0, mat.dim1(), vec, unit);
}

template <class Mat>
void orth(Mat && mat, int first, int last, int dest, bool unit=false) {
  orth(mat, first, last, mat.col(dest), unit);
}

template <class Mat>
void orth_unit(Mat && mat, int first, int last) {
  for (int j = first; j < last; ++j) {
    orth(mat, first, j, j, true);
    vec_unit(mat.col(j));
  }
}

/**
 * A = P B Q', where B is upper bidiagonal, P and Q are packed as set
 * of householder transform with essential parts stored in lower and
 * upper part of A, respectively. diag, and upper diag of B are stored
 * in diag and upper diag of A.
 *
 * see Matrix Computation (3ed), Algorithm 5.4.2.
 *
 * @param iA
 * @param itu
 * @param itv
 */
template <class CA, class CTP, class CTQ>
void bidiag_householder(CA && A, CTP && tu, CTQ && tv) {
  int m = A.dim0(), n = A.dim1();
  if (!n || !m) return;
  assert(m >= n);
  assert(tu.size() == n);
  assert(tv.size() == n - 1);
  for (int j = 0; j < n; ++j) {
    auto c = A.col(j);
    auto cl = make_vec(&c, m - j, vec_index_start(j));
    auto cln = cl.norm2();
    double t_u = householder_transform(cl);
    tu[j] = t_u;
    auto Asc = make_mat(
        &A, m - j, n - (j + 1),
        vec_index_start(j), vec_index_start(j + 1));
    householder_hm(t_u, cl, Asc);
    if (j + 1 < n) {
      auto r = A.row(j);
      auto rl = make_vec(&r, n - (j + 1), vec_index_start(j + 1));
      double rln = rl.norm2();
      double t_v = householder_transform(rl);
      auto Asr = make_mat(
          &A, m - (j + 1), n - (j + 1),
          vec_index_start(j + 1), vec_index_start(j + 1));
      householder_mh(t_v, rl, Asr);
      tv[j] = t_v;
      rl[0] = rln;
    }
    cl[0] = cln;
  }
}

/**
 * Golub-Kahan-Lanzcos bidiagonalization, decomposition A = PBQ', with
 * B is bidiagonal.
 *
 * gemv(A, x) and gemv(A.trans(), x) is implicity done by Ax and Atx functions.
 *
 * P and Q is orthognal in columns.
 *
 * One-Sided version.
 *
 * @param iAx Ax(x, y) for y=gemv(A, x)
 * @param iAtx Atx(x, y) for y=gemv(A', x)
 * @param iD diag of B
 * @param iE upper diag of B
 * @param iP
 * @param iQ
 */
template <class CAX, class CATX, class CD, class CE, class CP, class CQ>
void bidiag_gkl(CAX && Ax, CATX && Atx, CD && D, CE && E, CP && P, CQ && Q) {
  assert(Q.dim1() > E.size());
  assert(P.dim1() >= D.size());
  assert(Q.dim1() >= P.dim1() && Q.dim1() <= P.dim1() + 1);
  int k = P.dim1();
  if (Q.dim1() == 0) return;
  Q.col(0) = vec_random(Q.dim0()) - .5;
  vec_unit(Q.col(0));
  double beta = 0, alpha = 0;
  for (int i = 0; i < k; ++i) {
    Ax(Q.col(i), P.col(i));
    if (i > 0) P.col(i) += -beta * P.col(i - 1);
    // orth(iP, 0, i, P.col(i), true);
    alpha = vec_unit(P.col(i));
    if (i < D.size()) D[i] = alpha;
    if (i + 1 < Q.dim1()) {
      Atx(P.col(i), Q.col(i + 1));
      // Q.col(i + 1) += - alpha * Q.col(i);
      orth(Q, 0, i + 1, Q.col(i + 1), true);
      beta = vec_unit(Q.col(i + 1));
      if (i < E.size()) E[i] = beta;
    }
  }
}

/**
 * Golub-Kahan-Lanzcos bidiagonalization, decomposition A = PBQ', with
 * B is bidiagonal.
 *
 * P and Q is orthognal in columns.
 *
 * One-Sided version.
 *
 * @param iA input matrix
 * @param iD diag of B
 * @param iE upper diag of B
 * @param iP
 * @param iQ
 */
template <class CA, class CD, class CE, class CP, class CQ>
void bidiag_gkl(CA && A, CD && D, CE && E, CP && P, CQ && Q) {
  bidiag_gkl(make_gemv_ax(&A), make_gemv_atx(&A), D, E, P, Q);
}

template <class CAX, class CATX, class CD, class CE, class CRho, class CP, class CQ>
void bidiag_gkl_restart(
    int locked, int l, int n,
    CAX && Ax, CATX && Atx, CD && D, CE && E, CRho && rho, CP && P, CQ && Q) {
  // enhancements version from SLEPc
  const double eta = 1.e-10;
  
  int rank, nprocs;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  // Step 1
  // Ax(Q.col(l), P.col(l), P.dim0() > 1000);
  auto vec_container<double> tmp(Ax.dim0());
  
  Ax(Q.col(l), tmp, P.dim0() > 1000);
 // start_indx 
  for(int i = 0; i < tmp.size(); ++i) 
  	P.col(l).get(i) = tmp.get(i);
  for(int i = tmp.size(); i < P.col(l).size(); ++i) 
  	P.col(l).get(i) = 0;
  
  MPI_Gatherv();
  // Generate P.col(l)
  
  // Step 2
  if(rank == 0) {
    for (int j = locked; j < l; ++j) {
      P.col(l) += -rho(j) * P.col(j);
    }
  }
  MPI_Bcast(&P[0], P.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
  
  // Main loop
  vec_container<double> T(n);
  for (int j = l; j < n; ++j) {
    // Step 3   
    //Atx(P.col(j), Q.col(j + 1), Q.dim0() > 1000);
    auto vec_container<double> tmp2(Atx.dim0());
    Atx(P.col(j), tmp2, Q.dim0() > 1000);
    Q.col(j+1).get() = 0;
    MPI_Gatherv;
     
    if(rank == 0) {
      auto Qj = mat_cols(Q, 0, j + 1);
      auto Tj = make_vec(&T, j + 1);
      
      // Step 4
      Tj.assign(gemv(Qj.trans(), Q.col(j + 1)), j >= 3);

      // Step 5
      double r = Q.col(j + 1).norm2();
      D[j] = vec_unit(P.col(j));
      Q.col(j + 1).scale(1. / D[j]);
      Tj = Tj / D[j];
      r /= D[j];
      Q.col(j + 1).plus_assign(- gemv(Qj, Tj), Q.dim0() > 1000);

      // Step 6
      double beta = r * r - Tj.square_sum();
      if (beta < eta * r * r) {
        Tj.assign(gemv(Qj.trans(), Q.col(j + 1)), Q.dim0() > 1000);
        r = Q.col(j + 1).square_sum();
        Q.col(j + 1).plus_assign(-gemv(Qj, Tj), Q.dim0() > 1000);
        beta = r * r - Tj.square_sum();
      }
      beta = std::sqrt(beta);
      E[j] = beta;
      Q.col(j + 1).scale(1. / E[j]);
      
      // Step 7
      if (j + 1 < n) {
        Ax(Q.col(j + 1), P.col(j + 1), P.dim0() > 1000);
        P.col(j + 1).plus_assign(- E[j] * P.col(j), P.dim0() > 1000);
      }
    }
    return ;
}

template <class CA, class CD, class CE, class CRho, class CP, class CQ>
void bidiag_gkl_restart(int locked, int l, int n,
                        CA && A, CD && D, CE && E, CRho && rho, CP && P, CQ && Q) {
  bidiag_gkl_restart(
      locked, l, n,
      make_gemv_ax(&A), make_gemv_atx(&A), D, E, rho, P, Q);
}

} // namespace linalg
} // namespace douban
#endif
