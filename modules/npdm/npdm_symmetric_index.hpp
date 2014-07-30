#ifndef __NPDM_SYMMETRIC_INDEX_HPP
#define __NPDM_SYMMETRIC_INDEX_HPP

#include <utility>

namespace SpinAdapted {
namespace Npdm {

/// return conditional 2 bits to determine index case
/// (at most 40% faster than using if statement)
inline unsigned int __SF_get_unique_type_element (unsigned int i)
{
   // 00: i = 0
   // 01: i > 0
   // 11: i < 0
   return i ? (i | 0x7fffffff) >> 30 : 0u;
}

// ====================================================================================================
//  For 2PDM
// ====================================================================================================

/// Index type, which is used as offset to get spatial 2PDM element
enum INDEX_TYPE_2PDM
{
   T2_IJ = 1,
   T2_II = 0,

   T2_IJKL = 3, // 11
   T2_IJKK = 2, // 10
   T2_IIKL = 1, // 01
   T2_IIKK = 0  // 00
};

/// Index case, to determine parity in spin-orbital 2PDM, and non-redundancy of spatial 2PDM
enum INDEX_CASE_2PDM
{
   // I > J
   C2_IJ = 1, // 01
   C2_JI = 3, // 11
   C2_II = 0  // 00
};

/// return unique number depends on permutation pattern
inline unsigned int __SF_get_unique_type (int i, int j)
{
   return  __SF_get_unique_type_element(i-j);
}

// ====================================================================================================
// For 3PDM
// ====================================================================================================

/// Index type, which is used as offset to get spatial 3PDM element
enum INDEX_TYPE_3PDM
{
   T3_IJK = 3, // 11
   T3_IJJ = 2, // 10
   T3_IIJ = 1, // 01
   T3_III = 0, // 00

   T3_IJKLMN = 15, // 1111
   T3_IJKLMM = 14, // 1110
   T3_IJKLLM = 13, // 1101
   T3_IJKLLL = 12, // 1100
   T3_IJJLMN = 11, // 1011
   T3_IJJLMM = 10, // 1010
   T3_IJJLLM =  9, // 1001
   T3_IJJLLL =  8, // 1000
   T3_IIJLMN =  7, // 0111
   T3_IIJLMM =  6, // 0110
   T3_IIJLLM =  5, // 0101
   T3_IIJLLL =  4, // 0100
   T3_IIILMN =  3, // 0011
   T3_IIILMM =  2, // 0010
   T3_IIILLM =  1, // 0001
   T3_IIILLL =  0  // 0000
};

/// Index case, to determine parity in spin-orbital 3PDM, and non-redundancy of spatial 3PDM
enum INDEX_CASE_3PDM
{
   // I > J > K
   C3_IJK =  21, // 010101 +
   C3_IKJ =  23, // 010111 -
   C3_JIK =  53, // 110101 -
   C3_JKI =  31, // 011111 +
   C3_KIJ =  61, // 111101 +
   C3_KJI =  63, // 111111 -

   C3_IJJ =  20, // 010100 +
   C3_JIJ =  49, // 110001 -
   C3_JJI =  15, // 001111 +

   C3_IIJ =   5, // 000101 +
   C3_IJI =  19, // 010011 -
   C3_JII =  60, // 111100 +

   C3_III =   0, // 000000 +
};

/// return unique number depends on permutation pattern
inline unsigned int __SF_get_unique_type (int i, int j, int k)
{
   return  (__SF_get_unique_type_element(i-j) << 4)
         | (__SF_get_unique_type_element(i-k) << 2)
         |  __SF_get_unique_type_element(j-k);
}

// ====================================================================================================
// For 4PDM
// ====================================================================================================

/// Index type, which is used as offset to get spatial 4PDM element
enum INDEX_TYPE_4PDM
{
   T4_IJKL = 7,
   T4_IJKK = 6,
   T4_IJJK = 5,
   T4_IIJK = 4,
   T4_IJJJ = 3,
   T4_IIJJ = 2,
   T4_IIIJ = 1,
   T4_IIII = 0
};

/// Index case, to determine parity in spin-orbital 4PDM, and non-redundancy of spatial 4PDM
enum INDEX_CASE_4PDM
{
   // I > J > K > L

   // For either spin-orbital and spatial 4PDM

   C4_IJKL =1365, // +
   C4_IJLK =1367, // -
   C4_IKJL =1397, // -
   C4_IKLJ =1375, // +
   C4_ILJK =1405, // +
   C4_ILKJ =1407, // -

   C4_JIKL =3413, // -
   C4_JILK =3415, // +
   C4_JKIL =1909, // +
   C4_JKLI =1503, // -
   C4_JLIK =1917, // -
   C4_JLKI =1535, // +

   C4_KIJL =3925, // +
   C4_KILJ =3543, // -
   C4_KJIL =3957, // -
   C4_KJLI =3551, // +
   C4_KLIJ =2045, // +
   C4_KLJI =2047, // -

   C4_LIJK =4053, // -
   C4_LIKJ =4055, // +
   C4_LJIK =4085, // +
   C4_LJKI =4063, // -
   C4_LKIJ =4093, // -
   C4_LKJI =4095, // +

   // For spatial 4PDM

   C4_IJKK =1364, // +
   C4_IKJK =1393, // -
   C4_IKKJ =1359, // +

   C4_JIKK =3412, // -
   C4_JKIK =1905, // +
   C4_JKKI =1487, // -

   C4_KIJK =3861, // +
   C4_KJIK =3893, // -
   C4_KIKJ =3287, // -
   C4_KJKI =3295, // +
   C4_KKIJ =1021, // +
   C4_KKJI =1023, // -

   C4_IJJK =1349, // +
   C4_IJKJ =1363, // -
   C4_IKJJ =1404, // +

   C4_JIKJ =3351, // +
   C4_JKIJ =1853, // -
   C4_JIJK =3157, // -
   C4_JKJI =1279, // +
   C4_JJIK = 885, // +
   C4_JJKI = 479, // -

   C4_KJJI =4047, // -
   C4_KJIJ =4081, // +
   C4_KIJJ =4052, // -

   C4_IJKI =1311, // +
   C4_IKJI =1343, // -
   C4_IJIK =1141, // -
   C4_IKIJ =1149, // +
   C4_IIJK = 341, // +
   C4_IIKJ = 343, // -

   C4_JIIK =3909, // +
   C4_JIKI =3539, // -
   C4_JKII =2044, // +

   C4_KIIJ =4037, // -
   C4_KIJI =4051, // +
   C4_KJII =4092, // -

   C4_IJJJ =1344, // +
   C4_JIJJ =3092, // -
   C4_JJIJ = 817, // +
   C4_JJJI = 207, // -

   C4_IIJJ = 340, // +
   C4_IJIJ =1137, // -
   C4_IJJI =1295, // +
   C4_JIIJ =3845, // +
   C4_JIJI =3283, // -
   C4_JJII =1020, // +

   C4_IIIJ =  69, // +
   C4_IIJI = 275, // -
   C4_IJII =1084, // +
   C4_JIII =4032, // -

   C4_IIII =   0  // +
};

/// return unique number depends on permutation pattern
inline unsigned int __SF_get_unique_type (int i, int j, int k, int l)
{
   return  (__SF_get_unique_type_element(i-j) <<10)
         | (__SF_get_unique_type_element(i-k) << 8)
         | (__SF_get_unique_type_element(i-l) << 6)
         | (__SF_get_unique_type_element(j-k) << 4)
         | (__SF_get_unique_type_element(j-l) << 2)
         |  __SF_get_unique_type_element(k-l);
}

// ====================================================================================================

/// return whether parity needs to be taken in account (only for spin-orbital 1PDM)
inline bool get_parity (int i)
{
   return false;
}

/// return whether parity needs to be taken in account (only for spin-orbital 2PDM)
inline bool get_parity (int i, int j)
{
   return i < j;
}

/// return whether parity needs to be taken in account (only for spin-orbital 3PDM)
inline bool get_parity (int i, int j, int k)
{
   switch(__SF_get_unique_type(i,j,k)) // gives unique number for each case
   {
      case C3_KJI:
      case C3_JIK:
      case C3_IKJ:
         return true;
      default:
         return false;
   }
}

/// return whether parity needs to be taken in account (only for spin-orbital 4PDM)
inline bool get_parity (int i, int j, int k, int l)
{
   switch(__SF_get_unique_type(i,j,k,l)) // gives unique number for each case
   {
      case C4_IJLK:
      case C4_IKJL:
      case C4_ILKJ:
      case C4_JIKL:
      case C4_JKLI:
      case C4_JLIK:
      case C4_KILJ:
      case C4_KJIL:
      case C4_KLJI:
      case C4_LIJK:
      case C4_LJKI:
      case C4_LKIJ:
         return true;
      default:
         return false;
   }
}

/* ====================================================================================================

inline std::pair<unsigned int, unsigned int> get_spatial_type (int i, int j)
{
   unsigned int idx_type;
   unsigned int idx_offs = 0u;
   switch(__SF_get_unique_type(i,j)) // gives unique number for each case
   {
      case C_JI: ++idx_offs;
      case C_IJ: idx_type = T_IJ; break;

      case C_II: idx_type = T_II; break;

      default: abort();
   }

   return std::make_pair(idx_type, idx_offs);
}

inline std::pair<unsigned int, unsigned int> get_spatial_type (int i, int j, int k)
{
   unsigned int idx_type;
   unsigned int idx_offs = 0u;
   switch(__SF_get_unique_type(i,j,k)) // gives unique number for each case
   {
      case C_KJI: ++idx_offs;
      case C_KIJ: ++idx_offs;
      case C_JIK: ++idx_offs;
      case C_JKI: ++idx_offs;
      case C_IKJ: ++idx_offs;
      case C_IJK: idx_type = T_IJK; break;

      case C_JJI: ++idx_offs;
      case C_JIJ: ++idx_offs;
      case C_IJJ: idx_type = T_IJJ; break;

      case C_JII: ++idx_offs;
      case C_IJI: ++idx_offs;
      case C_IIJ: idx_type = T_IIJ; break;

      case C_III: idx_type = T_III; break;

      default: abort();
   }

   return std::make_pair(idx_type, idx_offs);
}

inline std::pair<unsigned int, unsigned int> get_spatial_type (int i, int j, int k, int l)
{
   unsigned int idx_type;
   unsigned int idx_offs = 0u;
   switch(__SF_get_unique_type(i,j,k,l)) // gives unique number for each case
   {
      // FIXME: not yet implemented...
      default: abort();
   }

   return std::make_pair(idx_type, idx_offs);
}
*/
} // namespace Npdm
} // namespace SpinAdapted

#endif // __NPDM_SYMMETRIC_INDEX_HPP
