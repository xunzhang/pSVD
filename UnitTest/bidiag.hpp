/**
 * @brief bidiagonalization
 **/
#ifndef FILE_8180cdd3_0cd2_4606_8640_9c0792474dde_H
#define FILE_8180cdd3_0cd2_4606_8640_9c0792474dde_H
#include "douban/linalg/householder.hpp"
#include "douban/linalg/gemv_ax.hpp"
#include "douban/matvec/vec_kahan.hpp"
#include "douban/utility/kahan_sum.hpp"
#include "douban/matvec/io.hpp"

#include <sys/time.h>

static double currenttime(void) {
	double timestamp;
	struct timeval tv;
	gettimeofday(&tv, 0);
	timestamp = (double)((double)(tv.tv_sec*1e6) + (double)tv.tv_usec);
	return timestamp;
}

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

/*
template <class CM, class CRecv>
void local_union(CM && mat, CRecv && recv_tmp, int indx, int nprocs) {
  for(size_t i = 0; i < mat.dim0(); ++i)
    for(size_t j = 1; j < (size_t)nprocs; ++j)
      recv_tmp.get(i) += recv_tmp.get(j * mat.dim0() + i);
  for(size_t i = 0; i < mat.dim0(); ++i)
    mat.col(indx).get(i) = recv_tmp.get(i);
}
*/

template <class CM, class CRecv>
void local_union(CM && mat, CRecv && recv_tmp, int indx, int nprocs) {
  for(size_t i = 0; i < mat.dim0(); ++i)
    // mat.col(indx)[i] = 0;
    for(size_t j = 1; j < (size_t)nprocs; ++j)
      mat.col(indx)[i] += recv_tmp.get(j * mat.dim0() + i);
}

/*
template<class CM, class CV, class CR>
void parallel_gemv_task(CM && mat, CV && vec, CR && res) {
  int rank, nprocs;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  double l_start, l_end;

  int m = (int)mat.dim0();// n = (int)mat.dim1();
  int load = m / nprocs;
  int remainder = m % nprocs;
  int max_load = load + remainder;
  int *rcounts = new int [nprocs];
  int *displs = new int [nprocs];
  
  // add include file??!!
  mat_container<double> pA(max_load, n);
  vec_container<double> y_tmp(max_load);
  
  l_start = currenttime();
  // copy value of mat by each procs
  for(int i = 0; i < load; ++i)
    for(int j = 0; j < n; ++j)
      pA.get(i, j) = mat.get(i + rank * load, j);
  if(remainder != 0 && rank == nprocs - 1) 
    for(int i = load; i < max_load; ++i)
      for(int j = 0; j < n; ++j)
        pA.get(i, j) = mat.get(i + rank * load, j); 
  l_end = currenttime();
  std::cout << "prapare in step 4 cost time is : " << (l_end - l_start) / 1.0e6 << std::endl;

  l_start = currenttime();
  // General matrix vector multiplication
  y_tmp = gemv(pA, vec);
  l_end = currenttime();
  std::cout << "pA.dim0 is " << pA.dim0() << std::endl;
  std::cout << "pA.dim1 is " << pA.dim1() << std::endl;
  std::cout << "mv in step 4 cost time is : " << (l_end - l_start) / 1.0e6 << std::endl;
  
  // Prepare for MPI_Gatherv
  for(int i = 0; i < nprocs; ++i) {
    rcounts[i] = load;
    displs[i] = i * load;
  }
  if(remainder != 0)
    rcounts[nprocs - 1] = max_load;
  if(rank == nprocs - 1)
    load = max_load;
  
  // MPI_Gatherv
  MPI_Gatherv(&y_tmp[0], load, MPI_DOUBLE, &res[0], &rcounts[0], &displs[0], MPI_DOUBLE, 0, MPI_COMM_WORLD);
  l_start = currenttime();
  std::cout << "gatherv in step 4 cost time is : " << (l_start - l_end) / 1.0e6 << std::endl;

  delete [] rcounts;
  deaete [] displs;
  
  return;
}
*/

template<class CM, class CV, class CR>
void parallel_gemv_task(CM && mat, CV && vec, CR && res) {
  int rank, nprocs;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  int m = (int)mat.dim0();// n = (int)mat.dim1();
  int load = m / nprocs;
  int remainder = m % nprocs;
  int max_load = load + remainder;
  int *rcounts = new int [nprocs];
  int *displs = new int [nprocs];
  int offset;
  offset = rank * load;
  
  douban::vec_container<double> y_tmp(max_load);
   
  // General matrix vector multiplication
  y_tmp = douban::gemv(mat_rows(mat, offset, offset + max_load), vec);
  
  // Prepare for MPI_Gatherv
  for(int i = 0; i < nprocs; ++i) {
    rcounts[i] = load;
    displs[i] = i * load;
  }
  if(remainder != 0)
    rcounts[nprocs - 1] = max_load;
  if(rank == nprocs - 1)
    load = max_load;

  // MPI_Gatherv
  MPI_Gatherv(&y_tmp[0], load, MPI_DOUBLE, &res[0], &rcounts[0], &displs[0], MPI_DOUBLE, 0, MPI_COMM_WORLD);
  
  delete [] rcounts;
  delete [] displs;
  
  return;
}

template <class CAX, class CATX, class CD, class CE, class CRho, class CP, class CQ>
void bidiag_gkl_restart(
    int locked, int l, int n,
    CAX && Ax, CATX && Atx, CD && D, CE && E, CRho && rho, CP && P, CQ && Q, int s_indx, int t_s_indx) {
  // enhancements version from SLEPc
  const double eta = 1.e-10;
  double t_start = 0.0, t_end = 0.0;
  double local_start = 0.0, local_end = 0.0; 
  double t_total3 = 0.0, t_total4 = 0.0, t_total5 = 0.0, t_total6 = 0.0, t_total7 = 0.0;
  
  int rank, nprocs;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  
  // Step 1
  int recv_len = (int)P.dim0() * nprocs;
  vec_container<double> tmp(Ax.dim0());
  vec_container<double> recv_tmp(recv_len);
  
  auto m_Ax = make_gemv_ax(&Ax);
  auto m_Atx = make_gemv_ax(&Atx);
  
  m_Ax(Q.col(l), tmp, P.dim0() > 1000);

  vec_container<double> send_data(P.dim0(),0);
  for(size_t i = s_indx; i < s_indx + Ax.dim0(); ++i)
    send_data[i] = tmp.get(i-s_indx);
  MPI_Gather(&send_data[0], P.dim0(), MPI_DOUBLE, &recv_tmp[0], P.dim0(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

  P.col(l) = 0;
  // Generate truly P.col(l)
  if(rank == 0) {
    local_union(P, recv_tmp, l, nprocs);
    // Step 2 & also in rank 0
    for (int j = locked; j < l; ++j) {
      P.col(l) += -rho(j) * P.col(j);
    }
  }
  
  MPI_Bcast(&(P.col(0)[0]), P.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD); 
  //MPI_Bcast(&(P.col(l)[0]), P.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD); 
  
  // Main loop
  vec_container<double> T(n);
  int recv_l = Q.dim0() * nprocs;
  vec_container<double> recv_t(recv_l);
  for (int j = l; j < n; ++j) {
    // Step 3   
    vec_container<double> tmp2(Atx.dim0());
    
    /* for print */
    if(rank == 0)
    	t_start = currenttime();
   
    local_start = currenttime();
    m_Atx(P.col(j), tmp2, Q.dim0() > 1000);
    local_end = currenttime();
    std::cout << "parallel mv time cost is " << (local_end - local_start) / 1.0e6 << std::endl;
     
    vec_container<double> s_data(Q.dim0(), 0); 
    for(size_t i = t_s_indx; i < t_s_indx + Atx.dim0(); ++i)
      s_data[i] = tmp2[i-t_s_indx];
    MPI_Gather(&s_data[0], Q.dim0(), MPI_DOUBLE, &recv_t[0], Q.dim0(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
    local_start = currenttime();
    std::cout << "parallel mv time cost2 is " << (local_start - local_end) / 1.0e6 << std::endl;
    
    //Q.col(j+1) = 0;
    if(rank == 0) {
      // Generate truly Q.col(j+1) 
      local_union(Q, recv_t, j + 1, nprocs);
      local_end = currenttime();
      t_end = currenttime();
      std::cout << "parallel mv time cost3 is " << (local_end - local_start) / 1.0e6 << std::endl;
      std::cout << "time of step 3 is : " << (t_end - t_start) / 1.0e6 << std::endl;
      t_total3 += (t_end - t_start) / 1.0e6;
    }
      
    // Step 4
    for(size_t aa = 0; aa < Q.dim0(); ++aa) // row
      MPI_Bcast(&(Q.row(aa)[0]), j + 2, MPI_DOUBLE, 0, MPI_COMM_WORLD); 
    // MPI_Bcast(&(Q.col(0)[0]), Q.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD); 
    
    if(rank == 0)
      t_end = currenttime();
    auto Qj = mat_cols(Q, 0, j + 1);
    auto Tj = make_vec(&T, j + 1);
    
    //Tj.assign(gemv(Qj.trans(), Q.col(j + 1)), j >= 3);
    parallel_gemv_task(Qj.trans(), Q.col(j+1), Tj);
    
    if(rank == 0) {
      t_start = currenttime();
      t_total4 += (t_start - t_end) / 1.0e6;
      std::cout << "time of step 4 is : " << (t_start - t_end) / 1.0e6 << std::endl;
    }

    // Step 5
    if(rank == 0) {
      double r = Q.col(j + 1).norm2();
      D[j] = vec_unit(P.col(j));
      Q.col(j + 1).scale(1. / D[j]);
      Tj = Tj / D[j];
      r /= D[j];
      Q.col(j + 1).plus_assign(- gemv(Qj, Tj), Q.dim0() > 1000);
      
      t_end = currenttime();
      t_total5 += (t_end - t_start) / 1.0e6;
      std::cout << "time of step 5 is : " << (t_end - t_start) / 1.0e6 << std::endl;

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
      
      t_start = currenttime();
      t_total6 += (t_start - t_end) / 1.0e6;
      std::cout << "time of step 6 is : " << (t_start - t_end) / 1.0e6 << std::endl;
    } 
    
    // Step 7
    // MPI_Bcast(&(Q.col(j+1)[0]), Q.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
    // MPI_Bcast(&(Q.col(0)[0]), Q.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
    for(size_t aa = 0; aa < Q.dim0(); ++aa)
      MPI_Bcast(&(Q.col(j+1)[aa]), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (j + 1 < n) {
      if(rank == 0) 
        t_start = currenttime();
      vec_container<double> tmp3(Ax.dim0());
      vec_container<double> se_data(P.dim0(), 0);
      
      m_Ax(Q.col(j + 1), tmp3, P.dim0() > 1000);
      
      for(size_t k1 = s_indx; k1 < s_indx + Ax.dim0(); ++k1)
        se_data[k1] = tmp3[k1-s_indx];
      MPI_Gather(&se_data[0], P.dim0(), MPI_DOUBLE, &recv_tmp[0], P.dim0(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
      
      // P.col(j+1) = 0;
      if(rank == 0) {
	local_union(P, recv_tmp, j + 1, nprocs);
	P.col(j + 1).plus_assign(- E[j] * P.col(j), P.dim0() > 1000);
      }
      
      /* for print */
      if(rank == 0) {
        t_end = currenttime();
        t_total7 += (t_end - t_start) / 1.0e6;
	std::cout << "time of step 7 is : " << (t_end - t_start) / 1.0e6 << std::endl;
      }

      // MPI_Bcast(&(P.col(l)[0]), P.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
      // MPI_Bcast(&(P.col(0)[0]), P.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
      for(size_t aa = 0; aa < P.dim0(); ++aa)
        MPI_Bcast(&(P.col(j+1)[aa]), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    }  // end if
  }    // end while
  
  /* for print */
  if(rank == 0) {
    std::cout << "total step 3 time is : " << t_total3 << std::endl; 
    std::cout << "total step 4 time is : " << t_total4 << std::endl; 
    std::cout << "total step 5 time is : " << t_total5 << std::endl; 
    std::cout << "total step 6 time is : " << t_total6 << std::endl; 
    std::cout << "total step 7 time is : " << t_total7 << std::endl; 
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
