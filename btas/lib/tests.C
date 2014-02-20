#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>

#include <cstdlib>
double rgen() { return (static_cast<double>(rand())/RAND_MAX-0.5)*2; }

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#define _DEFAULT_QUANTUM 1

#include <btas/TVector.h>

#include <btas/DENSE/DArray.h>
#include <btas/QSPARSE/QSDArray.h>

#include <time_stamp.h>

using namespace std;

int DENSE_TEST(int iprint = 0)
{
  using namespace btas;

  DArray<4> a(8, 2, 2, 8); a.generate(rgen);
  DArray<2> b(8, 2);       b.generate(rgen);

  if(iprint > 0) {
    cout << "====================================================================================================" << endl;
    cout << "[DENSE_TEST] print tensor [a]: " << a << endl;
    cout << "====================================================================================================" << endl;
    cout << "[DENSE_TEST] print matrix [b]: " << b << endl;
  }

  if(1)
  {
    // Dcopy
    DArray<4> c;
    Dcopy(a, c);
    if(iprint > 0) {
      cout << "====================================================================================================" << endl;
      cout << "[DENSE_TEST] [Dcopy(a, c)] print tensor [c]: " << c << endl;
    }
    // Dscal
    Dscal(2.0, c);
    if(iprint > 0) {
      cout << "====================================================================================================" << endl;
      cout << "[DENSE_TEST] [Dscal(2.0, c)] print tensor [c]: " << c << endl;
    }
    // Daxpy
    Daxpy(1.0, a, c);
    if(iprint > 0) {
      cout << "====================================================================================================" << endl;
      cout << "[DENSE_TEST] [Daxpy(1.0, a, c)] print tensor [c]: " << c << endl;
    }
    // Ddot
    double bnorm = Ddot(b, b);
    if(iprint > 0) {
      cout << "====================================================================================================" << endl;
      cout << "[DENSE_TEST] [Ddot(b, b)] print [bnorm]: " << bnorm << endl;
    }
  }

  if(1)
  {
    // Dgemv
    DArray<2> c;
    Dgemv(Trans, 1.0, a, b, 1.0, c);
    if(iprint > 0) {
      cout << "====================================================================================================" << endl;
      cout << "[DENSE_TEST] [Dgemv(Trans, 1.0, a, b, 1.0, c)] print matrix [c]: " << c << endl;
    }
    // Dger
    DArray<4> d;
    Dger(1.0, b, b, d);
    if(iprint > 0) {
      cout << "====================================================================================================" << endl;
      cout << "[DENSE_TEST] [Dger(1.0, b, b, d)] print tensor [d]: " << d << endl;
    }
    // Dgemm
    DArray<4> e;
    Dgemm(NoTrans, NoTrans, 1.0, d, a, 1.0, e);
    if(iprint > 0) {
      cout << "====================================================================================================" << endl;
      cout << "[DENSE_TEST] [Dgemm(NoTrans, NoTrans, 1.0, d, a, 1.0, e)] print tensor [e]: " << e << endl;
    }
  }

  if(1)
  {
    // Dcontract
    DArray<4> c;
    Dcontract(1.0, a, shape(3), b, shape(0), 1.0, c);
    if(iprint > 0) {
      cout << "====================================================================================================" << endl;
      cout << "[DENSE_TEST] [Dcontract(1.0, a, shape(3), b, shape(0), 1.0, c)] print tensor [c]: " << c << endl;
    }
    DArray<2> c2;
    Dcontract(1.0, a, shape(1, 3), b, shape(1, 0), 1.0, c2);
    if(iprint > 0) {
      cout << "====================================================================================================" << endl;
      cout << "[DENSE_TEST] [Dcontract(1.0, a, shape(1, 3), b, shape(1, 0), 1.0, c2)] print tensor [c2]: " << c2 << endl;
    }
    // Dindexed_contract
    DArray<4> d;
    enum { i, j, k, l, p };
    Dindexed_contract(1.0, a, shape(i,p,k,j), b, shape(l,p), 1.0, d, shape(i,l,j,k));
    if(iprint > 0) {
      cout << "====================================================================================================" << endl;
      cout << "[DENSE_TEST] [Dindexed_contract(1.0, a, shape(i,p,k,j), b, shape(l,p), 1.0, d, shape(i,l,j,k))] print tensor [d]: " << d << endl;
    }
  }

  if(1)
  {
    // LAPACK
  }

  if(1)
  {
    // Dpermute
    DArray<4> c;
    Dpermute(a, shape(2, 0, 1, 3), c);
    if(iprint > 0) {
      cout << "====================================================================================================" << endl;
      cout << "[DENSE_TEST] [Dpermute(a, shape(2, 0, 1, 3), c)] print tensor [c]: " << c << endl;
    }
    // Ddiagonal
    DArray<3> d;
    Ddiagonal(c, shape(1, 3), d);
    if(iprint > 0) {
      cout << "====================================================================================================" << endl;
      cout << "[DENSE_TEST] [Ddiagonal(c, shape(1, 3), d)] print tensor [d]: " << d << endl;
    }
  }

  return 0;
}

int QSPARSE_TEST(int iprint = 0)
{
  using namespace btas;

  Quantum qt(0);

  Qshapes<> qi;
  qi.reserve(3);
  qi.push_back(Quantum(-1));
  qi.push_back(Quantum( 0));
  qi.push_back(Quantum(+1));

  Dshapes di(qi.size(), 2);

  TVector<Qshapes<>, 4> a_qshape = { qi,-qi, qi,-qi };
  TVector<Dshapes,   4> a_dshape = { di, di, di, di };
  QSDArray<4> a(qt, a_qshape, a_dshape); a.generate(rgen);

  TVector<Qshapes<>, 2> b_qshape = {-qi, qi };
  TVector<Dshapes,   2> b_dshape = { di, di };
  QSDArray<2> b(qt, b_qshape, b_dshape); b.generate(rgen);

  if(iprint > 0) {
    cout << "====================================================================================================" << endl;
    cout << "[QSPARSE_TEST] print tensor [a]: " << a << endl;
    cout << "====================================================================================================" << endl;
    cout << "[QSPARSE_TEST] print matrix [b]: " << b << endl;
  }

  if(1)
  {
    // QSDgemv
    QSDArray<2> c;
    QSDgemv(NoTrans, 1.0, a, b, 1.0, c);
    if(iprint > 0) {
      cout << "====================================================================================================" << endl;
      cout << "[QSPARSE_TEST] [QSDgemv(NoTrans, 1.0, a, b, 1.0, c)] print matrix [c]: " << c << endl;
    }
    // QSDger
    QSDArray<4> d;
    QSDger(1.0, b, b, d);
    if(iprint > 0) {
      cout << "====================================================================================================" << endl;
      cout << "[QSPARSE_TEST] [QSDger(1.0, b, b, d)] print tensor [d]: " << d << endl;
    }
    // QSDgemm
    QSDArray<4> e;
    QSDgemm(NoTrans, NoTrans, 1.0, d, a, 1.0, e);
    if(iprint > 0) {
      cout << "====================================================================================================" << endl;
      cout << "[QSPARSE_TEST] [QSDgemm(NoTrans, NoTrans, 1.0, d, a, 1.0, e)] print tensor [e]: " << e << endl;
    }
  }

  if(1)
  {
    double norm = QSDdotc(a, a);
    QSDscal(1.0/sqrt(norm), a);

    // QSDgesvd
     SDArray<1> s;
    QSDArray<3> u;
    QSDArray<3> v;
    QSDgesvd(LeftArrow, a, s, u, v, 12);
    if(iprint > 0) {
      cout << "====================================================================================================" << endl;
      cout << "[QSPARSE_TEST] [QSDgesvd(LeftArrow, a, s, u, v, 4)] print tensor [s]: " << s << endl;
      cout << "====================================================================================================" << endl;
      cout << "[QSPARSE_TEST] [QSDgesvd(LeftArrow, a, s, u, v, 4)] print tensor [u]: " << u << endl;
      cout << "====================================================================================================" << endl;
      cout << "[QSPARSE_TEST] [QSDgesvd(LeftArrow, a, s, u, v, 4)] print tensor [v]: " << v << endl;
    }

    // QSDgesvd with null space vector
     SDArray<1> s_rm;
    QSDArray<3> u_rm;
    QSDArray<3> v_rm;
    QSDgesvd(LeftArrow, a, s, s_rm, 1, u, u_rm, 1, v, v_rm, 12);
    if(iprint > 0) {
      cout << "====================================================================================================" << endl;
      cout << "[QSPARSE_TEST] [QSDgesvd(LeftArrow, a, s, s_rm, 1, u, u_rm, 1, v, v_rm, 4)] print tensor [s]: " << s << endl;
      cout << "====================================================================================================" << endl;
      cout << "[QSPARSE_TEST] [QSDgesvd(LeftArrow, a, s, s_rm, 1, u, u_rm, 1, v, v_rm, 4)] print tensor [u]: " << u << endl;
      cout << "====================================================================================================" << endl;
      cout << "[QSPARSE_TEST] [QSDgesvd(LeftArrow, a, s, s_rm, 1, u, u_rm, 1, v, v_rm, 4)] print tensor [v]: " << v << endl;
      cout << "====================================================================================================" << endl;
      cout << "[QSPARSE_TEST] [QSDgesvd(LeftArrow, a, s, s_rm, 1, u, u_rm, 1, v, v_rm, 4)] print tensor [s_rm]: " << s_rm << endl;
      cout << "====================================================================================================" << endl;
      cout << "[QSPARSE_TEST] [QSDgesvd(LeftArrow, a, s, s_rm, 1, u, u_rm, 1, v, v_rm, 4)] print tensor [u_rm]: " << u_rm << endl;
      cout << "====================================================================================================" << endl;
      cout << "[QSPARSE_TEST] [QSDgesvd(LeftArrow, a, s, s_rm, 1, u, u_rm, 1, v, v_rm, 4)] print tensor [v_rm]: " << v_rm << endl;
    }
  }

  if(1)
  {
    // QSDcontract
    QSDArray<4> c;
    QSDcontract(1.0, a, shape(1), b, shape(1), 1.0, c);
    if(iprint > 0) {
      cout << "====================================================================================================" << endl;
      cout << "[QSPARSE_TEST] [QSDcontract(1.0, a, shape(1), b, shape(1), 1.0, c)] print matrix [c]: " << c << endl;
    }
    // QSDcontract with conjugation
    QSDArray<4> d;
    QSDcontract(1.0, a.conjugate(), shape(2), b, shape(1), 1.0, d);
    if(iprint > 0) {
      cout << "====================================================================================================" << endl;
      cout << "[QSPARSE_TEST] [QSDcontract(1.0, a.conjugate(), shape(2), b, shape(1), 1.0, d)] print matrix [d]: " << d << endl;
    }
  }

  if(1)
  {
    // Direct sum of arrays
    QSDArray<4> x;
    QSTdsum(a, a, x);
    if(iprint > 0) {
      cout << "====================================================================================================" << endl;
      cout << "[QSPARSE_TEST] [QSTdsum(a, a, x)] print matrix [x]: " << x << endl;
    }

    // Partial direct sum of arrays
    QSDArray<4> y;
    QSTdsum(a, a, shape(1, 2), y);
    if(iprint > 0) {
      cout << "====================================================================================================" << endl;
      cout << "[QSPARSE_TEST] [QSTdsum(a, a, shape(1, 2), y)] print matrix [y]: " << y << endl;
    }
  }

  if(1)
  {
    // Erasing quantum number
    QSDArray<4> x = a;
    x.erase(2, 1); // erasing m_q_shape[2][1]
    if(iprint > 0) {
      cout << "====================================================================================================" << endl;
      cout << "[QSPARSE_TEST] [x = a; x.erase(2, 1)] print matrix [x]: " << x << endl;
    }

    // Making sub-array
    Dshapes d_1;
    d_1.push_back(1);
    d_1.push_back(2);
    Dshapes d_2;
    d_2.push_back(0);
    d_2.push_back(2);
    Dshapes d_3;
    d_3.push_back(0);
    d_3.push_back(1);
    d_3.push_back(2);
    Dshapes d_4;
    d_4.push_back(0);
    d_4.push_back(1);
    TVector<Dshapes, 4> sub_index = make_array(d_1, d_2, d_3, d_4);
    QSDArray<4> y(a.subarray(sub_index));
    if(iprint > 0) {
      cout << "====================================================================================================" << endl;
      cout << "[QSPARSE_TEST] [y(a.subarray({{1,2},{0,2},{0,1,2},{0,1}})] print matrix [y]: " << y << endl;
    }

  }

  if(1)
  {
    // Merge quantum number
    TVector<Qshapes<>, 4> a_qshape = a.qshape();
    TVector<Dshapes,   4> a_dshape = a.dshape();

    TVector<Qshapes<>, 2> qrows = { a_qshape[0], a_qshape[1] };
    TVector<Dshapes,   2> drows = { a_dshape[0], a_dshape[1] };

    TVector<Qshapes<>, 2> qcols = { a_qshape[2], a_qshape[3] };
    TVector<Dshapes,   2> dcols = { a_dshape[2], a_dshape[3] };

    QSTmergeInfo<2> row_qinfo(qrows, drows);
    QSTmergeInfo<2> col_qinfo(qcols, dcols);

    QSDArray<2> c;
    QSDmerge(row_qinfo, a, col_qinfo, c);
    if(iprint > 0) {
      cout << "====================================================================================================" << endl;
      cout << "[QSPARSE_TEST] [QSDmerge(row_qinfo, a, col_qinfo, c] print matrix [c]: " << c << endl;
    }

  }

  return 0;
}

int SERIALIZE_TEST(int iprint = 0)
{
  using namespace btas;

  Quantum qt(0);

  Qshapes<> qi;
  qi.reserve(3);
  qi.push_back(Quantum(-1));
  qi.push_back(Quantum( 0));
  qi.push_back(Quantum(+1));

  Dshapes di(qi.size(), 2);

  TVector<Qshapes<>, 4> a_qshape = { qi,-qi, qi,-qi };
  TVector<Dshapes,   4> a_dshape = { di, di, di, di };
  QSDArray<4> a(qt, a_qshape, a_dshape); a.generate(rgen);

  if(iprint > 0) {
    cout << "====================================================================================================" << endl;
    cout << "[SERIALIZE_TEST] print tensor [a]: " << a << endl;
  }

  {
    ofstream fout("tests.tmp");
    boost::archive::text_oarchive oa(fout);
    oa << a;
  }

  a.clear();

  {
    ifstream finp("tests.tmp");
    boost::archive::text_iarchive ia(finp);
    ia >> a;
  }

  if(iprint > 0) {
    cout << "====================================================================================================" << endl;
    cout << "[SERIALIZE_TEST] print tensor [a] (loaded): " << a << endl;
  }

  return 0;
}

int main()
{
  time_stamp ts;

  ts.start();

  DENSE_TEST(1);

  cout << "Finished DENSE_TEST: total elapsed time = "
       << setw(8) << setprecision(6) << fixed << ts.elapsed() << " sec. " << endl;

  ts.start();

  QSPARSE_TEST(1);

  cout << "Finished QSPARSE_TEST: total elapsed time = "
       << setw(8) << setprecision(6) << fixed << ts.elapsed() << " sec. " << endl;

  ts.start();

  SERIALIZE_TEST(1);

  cout << "Finished SERIALIZE_TEST: total elapsed time = "
       << setw(8) << setprecision(6) << fixed << ts.elapsed() << " sec. " << endl;

  return 0;
}
