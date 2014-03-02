#ifndef _BTAS_CXX11_QSTMERGE_INFO_H
#define _BTAS_CXX11_QSTMERGE_INFO_H 1

#include <map>

#include <btas/btas.h>
#include <btas/QSPARSE/Qshapes.h>

namespace btas {

   //! Merged index infomation by quantum number-based sparse array.
   /*!
    *  Provides reversible conversion of merged quantum number indices.
    *
    *  First, quantum number indices are packed as a vector form.
    *  Packed quantum numbers is a set of naive product of quantum number indices.
    *  e.g. product of two indices {q1, q2} * {q1, q2} = {q1*q1, q1*q2, q2*q1, q2*q2}
    *
    *  Then, packed indices are merged by quantum number
    *  e.g. {q1*q1, q1*q2, q2*q1, q2*q2} -> {q11, q12, q22}
    *
    *  Multi-map m_index_map provides map from {q11, q12, q22} to {q1*q1, {q1*q2, q2*q1}, q2*q2}
    *
    */
   template<size_t N, class Q = Quantum>
      class QSTmergeInfo {

         public:

            typedef typename std::multimap<int, int>::const_iterator const_iterator;
            typedef typename std::pair<const_iterator, const_iterator> const_range;

         public:
            //! Default constructor
            QSTmergeInfo() { }
            //! Initializer
            QSTmergeInfo (const TVector<Qshapes<Q>, N>& qshape, const TVector<Dshapes, N>& dshape) { reset(qshape, dshape); }

            //! Reset information
            void reset(const TVector<Qshapes<Q>, N>& qshape, const TVector<Dshapes, N>& dshape) {

               for(int i = 0; i < N; ++i)
                  if(qshape[i].size() != dshape[i].size())
                     BTAS_THROW(false, "btas::QSTmergeInfo::reset: qshape and dshape have different block size");

               // packing qshape and dshape into 1d-array
               Qshapes<Q> qshape_pkd(1, Q::zero());
               Dshapes    dshape_pkd(1, 1);

               for(int i = 0; i < N; ++i) {
                  qshape_pkd = qshape_pkd * qshape[i];
                  dshape_pkd = dshape_pkd * dshape[i];
               }

               // mapping packed index to merged quantum #
               std::map<Q, int> q_index_map;

               int n = 0;
               for(int i = 0; i < qshape_pkd.size(); ++i) {
                  if(q_index_map.find(qshape_pkd[i]) == q_index_map.end()) {
                     q_index_map.insert(std::make_pair(qshape_pkd[i], n++));
                  }
               }

               // copying merged quantum #
               Qshapes<Q> qshape_mgd(n, Q::zero());
               for(typename std::map<Q, int>::iterator iqmap = q_index_map.begin(); iqmap != q_index_map.end(); ++iqmap) {
                  qshape_mgd[iqmap->second] = iqmap->first;
               }

               // mapping packed index to merged index and computing merged block size
               m_index_map.clear();
               Dshapes dshape_mgd(n, 0);
               for(int i = 0; i < qshape_pkd.size(); ++i) {
                  int j = q_index_map.find(qshape_pkd[i])->second;
                  m_index_map.insert(std::make_pair(j, i));
                  dshape_mgd[j] += dshape_pkd[i];
               }
               // save as members
               m_dshape = dshape;
               m_dshape_packed = dshape_pkd;
               m_dshape_merged = dshape_mgd;

               m_qshape = qshape;
               m_qshape_merged = qshape_mgd;
            }

            const TVector<Dshapes, N>& dshape() const { return m_dshape; }
            const Dshapes& dshape(int i) const { return m_dshape[i]; }

            const Dshapes& dshape_packed() const { return m_dshape_packed; }
            const int& dshape_packed(int i) const { return m_dshape_packed[i]; }

            const Dshapes& dshape_merged() const { return m_dshape_merged; }
            const int& dshape_merged(int i) const { return m_dshape_merged[i]; }

            const TVector<Qshapes<Q>, N>& qshape() const { return m_qshape; }
            const Qshapes<Q>& qshape(int i) const { return m_qshape[i]; }

            const Qshapes<Q>& qshape_merged() const { return m_qshape_merged; }
            const Q& qshape_merged(int i) const { return m_qshape_merged[i]; }

            const_iterator begin() const { return m_index_map.begin(); }
            const_iterator end() const { return m_index_map.end(); }
            const_iterator find(int i) const { return m_index_map.find(i); }

            const_range equal_range(int i) const { return m_index_map.equal_range(i); }

         private:

            //! Original dense shapes
            TVector<Dshapes, N> m_dshape;

            //! Dense shapes for packed quantum numbers
            Dshapes m_dshape_packed;

            //! Dense shapes for merged quantum numbers
            Dshapes m_dshape_merged;

            //! Original quantum number indices
            TVector<Qshapes<Q>, N> m_qshape;

            //! Merged quantum number indices
            Qshapes<Q> m_qshape_merged;

            //! Map from merged index to packed indices
            std::multimap<int, int> m_index_map;

      };

};

#endif // _BTAS_CXX11_QSTMERGE_INFO_H
