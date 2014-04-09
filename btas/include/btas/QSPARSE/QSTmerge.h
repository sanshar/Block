
//
/*! \file  QSTmerge.h
 *  \brief Merging and expanding array indices.
 *
 *  \par Merging function named QSTmerge
 *  - Row-merge: { Rows-info(i,j <-> r), A(i,j,k,l) } -> B(r,k,l)
 *  - Column-merge: { A(i,j,k,l), Cols-info(k,l <-> c) } -> B(i,j,c)
 *  - Row/Column-merge: { Rows-info(i,j <-> r), A(i,j,k,l), Cols-info(k,l -> c) } -> B(r,c)
 *
 *  \par Expanding function named QSTexpand
 *  - Row-expand: { Rows-info(i,j <-> r), A(r,k,l) } -> B(i,j,k,l)
 *  - Column-expand: { A(i,j,c), Cols-info(k,l <-> c) } -> B(i,j,k,l)
 *  - Row/Column-expand: { Rows-info(i,j <-> r), A(r,c), Cols-info(k,l <-> c) } -> B(i,j,k,l)
 */

#ifndef _BTAS_CXX11_QSTMERGE_H
#define _BTAS_CXX11_QSTMERGE_H 1

#include <btas/btas.h>

#include <btas/DENSE/TSubArray.h>

#include <btas/QSPARSE/QSTmergeInfo.h>
#include <btas/QSPARSE/QSTArray.h>

namespace btas {

   //! Merging row ranks
   /*! Row-merge: { Rows-info(i,j <-> r), A(i,j,k,l) } -> B(r,k,l) */
   template<typename T, size_t MR, size_t N, class Q = Quantum>
      void QSTmerge
      (const QSTmergeInfo<MR, Q>& rows_info, const QSTArray<T, N, Q>& a, QSTArray<T, 1+N-MR, Q>& b)
      {
         const size_t MC = N-MR;
         // copying a_tshape
         TVector<Qshapes<Q>, N> a_qshape(a.qshape());
         TVector<Dshapes,    N> a_dshape(a.dshape());
         // new qnum shapes
         TVector<Qshapes<Q>, 1+MC> b_qshape;
         b_qshape[0] = rows_info.qshape_merged();
         for(int i = 0; i < MC; ++i) b_qshape[1+i] = a_qshape[MR+i];
         // new dense shapes
         TVector<Dshapes, 1+MC> b_dshape;
         b_dshape[0] = rows_info.dshape_merged();
         for(int i = 0; i < MC; ++i) b_dshape[1+i] = a_dshape[MR+i];
         // resizing
         b.resize(a.q(), b_qshape, b_dshape, false);
         // strides
         int a_stride = a.stride(MR-1);
         int b_stride = b.stride(0);
         int b_n_rows = b.shape (0);
         // loop over merged blocks
         for(int i = 0; i < b_n_rows; ++i) {
            typename QSTmergeInfo<MR, Q>::const_range irow_range = rows_info.equal_range(i);
            if(irow_range.first == irow_range.second) continue;

            int ib_rows = i * b_stride;
            for(int j = 0; j < b_stride; ++j) {
               // construct merged dense-tensor
               IVector<1+MC> b_index = b.index(ib_rows + j);
               if(!b.allowed(b_index)) continue;
               TArray<T, 1+MC> block(b.dshape() & b_index); block.fill(0.0);

               // loop over dense-array of a
               IVector<1+MC> subbeg = uniform<int, 1+MC>(0);
               IVector<1+MC> subend;
               for(int i = 0; i < 1+MC; ++i) subend[i] = block.shape(i) - 1;

               bool non_zero = false;
               for(typename QSTmergeInfo<MR, Q>::const_iterator itr = irow_range.first; itr != irow_range.second; ++itr) {
                  int irow = itr->second;
                  int drow = rows_info.dshape_packed(irow);
                  subend[0] = subbeg[0] + drow - 1;
                  // merge
                  int tag = irow * a_stride + j;
                  typename QSTArray<T, N, Q>::const_iterator ita = a.find(tag);
                  if(ita != a.end()) {
                     non_zero = true;
                     block.subarray(subbeg, subend) = *ita->second;
                  }
                  subbeg[0] = subend[0] + 1;
               }
               if(non_zero)
                  b.insert(b_index, block);
            }
         }
      }

   //! Merging col ranks
   /*! Column-merge: { A(i,j,k,l), Cols-info(k,l <-> c) } -> B(i,j,c) */
   template<typename T, size_t N, size_t MC, class Q = Quantum>
      void QSTmerge
      (const QSTArray<T, N, Q>& a, const QSTmergeInfo<MC, Q>& cols_info, QSTArray<T, N-MC+1, Q>& b)
      {
         const size_t MR = N-MC;
         // copying a_tshape
         TVector<Qshapes<Q>, N> a_qshape(a.qshape());
         TVector<Dshapes,    N> a_dshape(a.dshape());
         // new qnum shapes
         TVector<Qshapes<Q>, MR+1> b_qshape;
         for(int i = 0; i < MR; ++i) b_qshape[i] = a_qshape[i];
         b_qshape[MR] = cols_info.qshape_merged();
         // new dense shapes
         TVector<Dshapes, MR+1> b_dshape;
         for(int i = 0; i < MR; ++i) b_dshape[i] = a_dshape[i];
         b_dshape[MR] = cols_info.dshape_merged();
         // resizing
         b.resize(a.q(), b_qshape, b_dshape, false);
         // strides
         int a_stride = a.stride(MR-1);
         int b_stride = b.stride(MR-1);
         int b_n_rows = b.size() / b_stride;
         // loop over merged blocks
         for(int i = 0; i < b_n_rows; ++i) {
            int ia_rows = i * a_stride;
            int ib_rows = i * b_stride;
            for(int j = 0; j < b_stride; ++j) {
               typename QSTmergeInfo<MC, Q>::const_range jcol_range = cols_info.equal_range(j);
               if(jcol_range.first == jcol_range.second) continue;
               // construct merged dense-tensor
               IVector<MR+1> b_index = b.index(ib_rows + j);
               if(!b.allowed(b_index)) continue;
               TArray<T, MR+1> block(b.dshape() & b_index); block.fill(0.0);

               // loop over dense-array of a
               IVector<MR+1> subbeg = uniform<int, MR+1>(0);
               IVector<MR+1> subend;
               for(int i = 0; i < MR+1; ++i) subend[i] = block.shape(i) - 1;

               bool non_zero = false;
               for(typename QSTmergeInfo<MC, Q>::const_iterator itc = jcol_range.first; itc != jcol_range.second; ++itc) {
                  int jcol = itc->second;
                  int dcol = cols_info.dshape_packed(jcol);
                  subend[MR] = subbeg[MR] + dcol - 1;
                  // merge
                  int tag = ia_rows + jcol;
                  typename QSTArray<T, N, Q>::const_iterator ita = a.find(tag);
                  if(ita != a.end()) {
                     non_zero = true;
                     block.subarray(subbeg, subend) = *ita->second;
                  }
                  subbeg[MR] = subend[MR] + 1;
               }
               if(non_zero)
                  b.insert(b_index, block);
            }
         }
      }

   //! Merging row and col ranks to form matrix
   /*! Row/Column-merge: { Rows-info(i,j <-> r), A(i,j,k,l), Cols-info(k,l -> c) } -> B(r,c) */
   template<typename T, size_t MR, size_t MC, class Q = Quantum>
      void QSTmerge
      (const QSTmergeInfo<MR, Q>& rows_info, const QSTArray<T, MR+MC, Q>& a, const QSTmergeInfo<MC, Q>& cols_info, QSTArray<T, 2, Q>& b)
      {
         const size_t N = MR+MC;
         // new qnum shapes
         TVector<Qshapes<Q>, 2> b_qshape;
         b_qshape[0] = rows_info.qshape_merged();
         b_qshape[1] = cols_info.qshape_merged();
         // new dense shapes
         TVector<Dshapes,    2> b_dshape;
         b_dshape[0] = rows_info.dshape_merged();
         b_dshape[1] = cols_info.dshape_merged();
         // resizing
         b.resize(a.q(), b_qshape, b_dshape, false);
         // strides
         int a_stride = a.stride(MR-1);
         int b_stride = b.stride(0);
         int b_n_rows = b.shape (0);
         // loop over merged blocks
         for(int i = 0; i < b_n_rows; ++i) {
            typename QSTmergeInfo<MR, Q>::const_range irow_range = rows_info.equal_range(i);
            if(irow_range.first == irow_range.second) continue;
            for(int j = 0; j < b_stride; ++j) {
               typename QSTmergeInfo<MC, Q>::const_range jcol_range = cols_info.equal_range(j);
               if(jcol_range.first == jcol_range.second) continue;
               // construct merged dense-tensor
               IVector<2> b_index = shape(i, j);
               if(!b.allowed(b_index)) continue;
               TArray<T, 2> block(b.dshape() & b_index); block.fill(0.0);

               // loop over dense-array of a
               IVector<2> subbeg = uniform<int, 2>(0);
               IVector<2> subend = uniform<int, 2>(0);

               bool non_zero = false;
               for(typename QSTmergeInfo<MR, Q>::const_iterator itr = irow_range.first; itr != irow_range.second; ++itr) {
                  int irow = itr->second;
                  int drow = rows_info.dshape_packed(irow);
                  subend[0] = subbeg[0] + drow - 1;
                  subbeg[1] = 0;
                  for(typename QSTmergeInfo<MC, Q>::const_iterator itc = jcol_range.first; itc != jcol_range.second; ++itc) {
                     int jcol = itc->second;
                     int dcol = cols_info.dshape_packed(jcol);
                     subend[1] = subbeg[1] + dcol - 1;
                     // merge
                     int tag = irow * a_stride + jcol;
                     typename QSTArray<T, N, Q>::const_iterator ita = a.find(tag);
                     if(ita != a.end()) {
                        non_zero = true;
                        block.subarray(subbeg, subend) = *ita->second;
                     }
                     subbeg[1] = subend[1] + 1;
                  }
                  subbeg[0] = subend[0] + 1;
               }
               if(non_zero)
                  b.insert(b_index, block);
            }
         }
      }

   //! Expanding row ranks
   /*! Row-expand: { Rows-info(i,j <-> r), A(r,k,l) } -> B(i,j,k,l) */
   template<typename T, size_t MR, size_t N, class Q = Quantum>
      void QSTexpand
      (const QSTmergeInfo<MR, Q>& rows_info, const QSTArray<T, 1+N-MR, Q>& a, QSTArray<T, N, Q>& b)
      {
         const size_t MC = N-MR;
         // copying a_tshape
         TVector<Qshapes<Q>, 1+MC> a_qshape(a.qshape());
         TVector<Dshapes,    1+MC> a_dshape(a.dshape());
         // new qnum shapes
         TVector<Qshapes<Q>, N> b_qshape;

         for(int i = 0; i < MR; ++i)
            b_qshape[i] = rows_info.qshape(i);

         for(int i = 0; i < MC; ++i) 
            b_qshape[MR+i] = a_qshape[1+i];

         // new dense shapes
         TVector<Dshapes, N> b_dshape;

         for(int i = 0; i < MR; ++i)
            b_dshape[i]    = rows_info.dshape(i);

         for(int i = 0; i < MC; ++i)
            b_dshape[MR+i] = a_dshape[1+i];

         // resizing
         b.resize(a.q(), b_qshape, b_dshape, false);

         // strides
         int a_stride = a.stride(0);
         int b_stride = b.stride(MR-1);

         // loop over merged blocks
         for(typename QSTArray<T, 1+MC, Q>::const_iterator ita = a.begin(); ita != a.end(); ++ita) {

            int i = ita->first / a_stride;
            int j = ita->first % a_stride;

            typename QSTmergeInfo<MR, Q>::const_range irow_range = rows_info.equal_range(i);

            // construct merged dense-tensor
            TArray<T, 1+MC>& block = *(ita->second);

            // loop over dense-array of a
            IVector<1+MC> subbeg = uniform<int, 1+MC>(0);
            IVector<1+MC> subend;

            for(int i = 0; i < 1+MC; ++i)
               subend[i] = block.shape(i) - 1;

            for(typename QSTmergeInfo<MR, Q>::const_iterator itr = irow_range.first; itr != irow_range.second; ++itr) {

               int irow = itr->second;
               int drow = rows_info.dshape_packed(irow);

               // skip if size of block to be created = 0
               if(drow == 0) continue;

               subend[0] = subbeg[0] + drow - 1;

               // expand
               int tag = irow * b_stride + j;

               typename QSTArray<T, N, Q>::iterator itb = b.reserve(tag);

               if(itb != b.end())
                  *itb->second = block.subarray(subbeg, subend);

               subbeg[0] = subend[0] + 1;

            }

         }

      }

   //! Expanding col ranks
   /*! Column-expand: { A(i,j,c), Cols-info(k,l <-> c) } -> B(i,j,k,l) */
   template<typename T, size_t N, size_t MC, class Q = Quantum>
      void QSTexpand
      (const QSTArray<T, N-MC+1, Q>& a, const QSTmergeInfo<MC, Q>& cols_info, QSTArray<T, N, Q>& b)
      {
         const size_t MR = N-MC;

         // copying a_tshape
         TVector<Qshapes<Q>, MR+1> a_qshape(a.qshape());
         TVector<Dshapes,    MR+1> a_dshape(a.dshape());

         // new qnum shapes
         TVector<Qshapes<Q>, N> b_qshape;

         for(int i = 0; i < MR; ++i) 
            b_qshape[i] = a_qshape[i];

         for(int i = 0; i < MC; ++i) 
            b_qshape[MR+i] = cols_info.qshape(i);

         // new dense shapes
         TVector<Dshapes, N> b_dshape;
         for(int i = 0; i < MR; ++i) 
            b_dshape[i]    = a_dshape[i];
         for(int i = 0; i < MC; ++i) 
            b_dshape[MR+i] = cols_info.dshape(i);

         // resizing
         b.resize(a.q(), b_qshape, b_dshape, false);

         // strides
         int a_stride = a.stride(MR-1);
         int b_stride = b.stride(MR-1);

         // loop over merged blocks
         for(typename QSTArray<T, MR+1, Q>::const_iterator ita = a.begin(); ita != a.end(); ++ita) {

            int i = ita->first / a_stride;
            int j = ita->first % a_stride;

            typename QSTmergeInfo<MC, Q>::const_range jcol_range = cols_info.equal_range(j);

            // construct merged dense-tensor
            TArray<T, MR+1>& block = *(ita->second);

            // loop over dense-array of a
            IVector<MR+1> subbeg = uniform<int, MR+1>(0);
            IVector<MR+1> subend;

            for(int i = 0; i < MR+1; ++i) 
               subend[i] = block.shape(i) - 1;

            for(typename QSTmergeInfo<MC, Q>::const_iterator itc = jcol_range.first; itc != jcol_range.second; ++itc) {

               int jcol = itc->second;
               int dcol = cols_info.dshape_packed(jcol);

               // skip if size of block to be created = 0
               if(dcol == 0) continue;

               subend[MR] = subbeg[MR] + dcol - 1;

               //expand
               int tag = i * b_stride + jcol;

               typename QSTArray<T, N, Q>::iterator itb = b.reserve(tag);

               if(itb != b.end())
                  *itb->second = block.subarray(subbeg, subend);

               subbeg[MR] = subend[MR] + 1;

            }
         }
      }

   //! Expanding row and col of matrix
   /*! Row/Column-expand: { Rows-info(i,j <-> r), A(r,c), Cols-info(k,l <-> c) } -> B(i,j,k,l) */
   template<typename T, size_t MR, size_t MC, class Q = Quantum>
      void QSTexpand
      (const QSTmergeInfo<MR, Q>& rows_info, const QSTArray<T, 2, Q>& a, const QSTmergeInfo<MC, Q>& cols_info, QSTArray<T, MR+MC, Q>& b)
      {
         const size_t N = MR+MC;
         // new qnum shapes
         TVector<Qshapes<Q>, N> b_qshape;
         for(int i = 0; i < MR; ++i) 
            b_qshape[i]    = rows_info.qshape(i);

         for(int i = 0; i < MC; ++i) 
            b_qshape[MR+i] = cols_info.qshape(i);
         // new dense shapes

         TVector<Dshapes,    N> b_dshape;

         for(int i = 0; i < MR; ++i) 
            b_dshape[i]    = rows_info.dshape(i);

         for(int i = 0; i < MC; ++i)
            b_dshape[MR+i] = cols_info.dshape(i);

         // resizing
         b.resize(a.q(), b_qshape, b_dshape, false);

         // strides
         int a_stride = a.stride(0);
         int b_stride = b.stride(MR-1);

         // loop over merged blocks
         for(typename QSTArray<T, 2, Q>::const_iterator ita = a.begin(); ita != a.end(); ++ita) {

            int i = ita->first / a_stride;
            int j = ita->first % a_stride;

            typename QSTmergeInfo<MR, Q>::const_range irow_range = rows_info.equal_range(i);
            typename QSTmergeInfo<MC, Q>::const_range jcol_range = cols_info.equal_range(j);

            // construct merged dense-tensor
            const TArray<T, 2>& block = *(ita->second);

            // loop over dense-array of a
            IVector<2> subbeg = uniform<int, 2>(0);
            IVector<2> subend = uniform<int, 2>(0);

            for(typename QSTmergeInfo<MR, Q>::const_iterator itr = irow_range.first; itr != irow_range.second; ++itr) {

               int irow = itr->second;
               int drow = rows_info.dshape_packed(irow);
               // skip if size of block to be created = 0
               if(drow == 0) continue;

               subend[0] = subbeg[0] + drow - 1;
               subbeg[1] = 0;

               for(typename QSTmergeInfo<MC, Q>::const_iterator itc = jcol_range.first; itc != jcol_range.second; ++itc) {

                  int jcol = itc->second;
                  int dcol = cols_info.dshape_packed(jcol);
                  // skip if size of block to be created = 0
                  if(dcol == 0) continue;

                  subend[1] = subbeg[1] + dcol - 1;

                  // expand
                  int tag = irow * b_stride + jcol;

                  typename QSTArray<T, N, Q>::iterator itb = b.reserve(tag);

                  if(itb != b.end())
                     *itb->second = block.subarray(subbeg, subend);

                  subbeg[1] = subend[1] + 1;

               }

               subbeg[0] = subend[0] + 1;

            }
         }
      }

};

#endif // _BTAS_CXX11_QSTMERGE_H
