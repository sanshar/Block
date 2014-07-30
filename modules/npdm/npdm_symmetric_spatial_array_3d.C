
#include <npdm_symmetric_spatial_array.h>

namespace SpinAdapted {
namespace Npdm {

const unsigned int
symmetric_spatial_array_base<3>::C_offset_ijk_lmn[6][6]
= { //*******lmn,lnm,mln,mnl,nlm,nml
    /*ijk*/ { 0u, 1u, 2u, 3u, 4u, 5u },
    /*ikj*/ { 1u, 0u, 3u, 2u, 5u, 4u },
    /*jik*/ { 2u, 4u, 0u, 5u, 1u, 3u },
    /*jki*/ { 4u, 2u, 5u, 0u, 3u, 1u },
    /*kij*/ { 3u, 5u, 1u, 4u, 0u, 2u },
    /*kji*/ { 5u, 3u, 4u, 1u, 2u, 0u } };

const unsigned int
symmetric_spatial_array_base<3>::C_offset_ijk_ijk[6][6]
= { /********ijk,ikj,jik,jki,kij,kji*/
    /*ijk*/ { 0u, 1u, 2u, 3u, 3u, 4u },
    /*ikj*/ { 1u, 0u, 3u, 2u, 4u, 3u },
    /*jik*/ { 2u, 3u, 0u, 4u, 1u, 3u },
    /*jki*/ { 3u, 2u, 4u, 0u, 3u, 1u },
    /*kij*/ { 3u, 4u, 1u, 3u, 0u, 2u },
    /*kji*/ { 4u, 3u, 3u, 1u, 2u, 0u } };

const unsigned int
symmetric_spatial_array_base<3>::C_offset_ijk_lmm[6][3]
= { /********lmm,mlm,mml*/
    /*ijk*/ { 0u, 1u, 2u },
    /*ikj*/ { 0u, 2u, 1u },
    /*jik*/ { 1u, 0u, 2u },
    /*jki*/ { 1u, 2u, 0u },
    /*kij*/ { 2u, 0u, 1u },
    /*kji*/ { 2u, 1u, 0u } };

const unsigned int
symmetric_spatial_array_base<3>::C_offset_ijk_llm[6][3]
= { /********llm,lml,mll*/
    /*ijk*/ { 0u, 1u, 2u },
    /*ikj*/ { 1u, 0u, 2u },
    /*jik*/ { 0u, 2u, 1u },
    /*jki*/ { 2u, 0u, 1u },
    /*kij*/ { 1u, 2u, 0u },
    /*kji*/ { 2u, 1u, 0u } };

//const unsigned int
//symmetric_spatial_array_base<3>::C_offset_ijk_lll[6][1]
//= { /********lll*/
//    /*ijk*/ { 0u },
//    /*ikj*/ { 0u },
//    /*jik*/ { 0u },
//    /*jki*/ { 0u },
//    /*kij*/ { 0u },
//    /*kji*/ { 0u } };

const unsigned int
symmetric_spatial_array_base<3>::C_offset_ijj_lmm[3][3]
= { /********lmm,mlm,mml*/
    /*ijj*/ { 0u, 1u, 1u },
    /*jij*/ { 1u, 0u, 1u },
    /*jji*/ { 1u, 1u, 0u } };

const unsigned int
symmetric_spatial_array_base<3>::C_offset_ijj_llm[3][3]
= { /********llm,lml,mll*/
    /*ijj*/ { 0u, 0u, 1u },
    /*jij*/ { 0u, 1u, 0u },
    /*jji*/ { 1u, 0u, 0u } };

const unsigned int
symmetric_spatial_array_base<3>::C_offset_ijj_iij[3][3]
= { /********iij,iji,jii*/
    /*ijj*/ { 0u, 0u, 1u },
    /*jij*/ { 0u, 1u, 0u },
    /*jji*/ { 1u, 0u, 0u } };

//const unsigned int
//symmetric_spatial_array_base<3>::C_offset_ijj_lll[3][1]
//= { /********lll*/
//    /*ijj*/ { 0u },
//    /*jij*/ { 0u },
//    /*jji*/ { 0u } };

const unsigned int
symmetric_spatial_array_base<3>::C_offset_iij_lmm[3][3]
= { /********lmm,mlm,mml*/
    /*iij*/ { 0u, 0u, 1u },
    /*iji*/ { 0u, 1u, 0u },
    /*jii*/ { 1u, 0u, 0u } };

const unsigned int
symmetric_spatial_array_base<3>::C_offset_iij_llm[3][3]
= { /********llm,lml,mll*/
    /*iij*/ { 0u, 1u, 1u },
    /*iji*/ { 1u, 0u, 1u },
    /*jii*/ { 1u, 1u, 0u } };

//const unsigned int
//symmetric_spatial_array_base<3>::C_offset_iij_lll[3][1]
//= { /********lll*/
//    /*iij*/ { 0u },
//    /*iji*/ { 0u },
//    /*jii*/ { 0u } };

//const unsigned int
//symmetric_spatial_array_base<3>::C_offset_iii_lll[1][1]
//= { /********lll*/
//    /*iii*/ { 0u } };

} // namespace Npdm
} // namespace SpinAdapted

