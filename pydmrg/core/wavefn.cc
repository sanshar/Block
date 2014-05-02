/*
 * File:
 */

#include "config.h"
#include <boost/serialization/vector.hpp>
#include "wavefunction.h"
#include "spinblock.h"
#include "SpinQuantum.h"
#include "MatrixBLAS.h"

using namespace SpinAdapted;

/*
 * ref wavefunction.C
 *     SpinAdapted::Wavefunction::LoadWavefunctionInfo
 *     guess_wavefunction.C
 */

//int mpirank = 0;
//char file [5000];
//sprintf(file, "%s%s%d%s%d%s%d%s%d%s", fileprefix, "/wave-", sites[0],
//        "-", *(sites.rbegin()), ".", mpirank, ".", root_id, ".tmp");
int save_wavefunction(char *filewave, Wavefunction *oldWave,
                      StateInfo *waveInfo)
{
    bool onedot = oldWave->get_onedot();

    std::ofstream ofs(filewave, std::ios::binary);
    boost::archive::binary_oarchive save_wave(ofs);
    save_wave << onedot
              << *waveInfo
              << *waveInfo->leftStateInfo
              << *(waveInfo->leftStateInfo->leftStateInfo)
              << *(waveInfo->leftStateInfo->rightStateInfo)
              << *waveInfo->rightStateInfo;
    if (!onedot) {
        save_wave << *(waveInfo->rightStateInfo->leftStateInfo)
            << *(waveInfo->rightStateInfo->rightStateInfo);
    }
    oldWave->Save(ofs);
    ofs.close();
    return 0;
}

int load_wavefunction(char *filewave, Wavefunction *oldWave,
                      StateInfo *waveInfo)
{
    bool onedot;

    //waveInfo->Allocate();

    std::ifstream ifs(filewave, std::ios::binary);
    boost::archive::binary_iarchive load_wave(ifs);
    load_wave >> onedot
              >> *waveInfo
              >> *waveInfo->leftStateInfo
              >> *(waveInfo->leftStateInfo->leftStateInfo)
              >> *(waveInfo->leftStateInfo->rightStateInfo)
              >> *waveInfo->rightStateInfo;
    if (!onedot) {
        load_wave >> *(waveInfo->rightStateInfo->leftStateInfo)
            >> *(waveInfo->rightStateInfo->rightStateInfo);
    }
    oldWave->Load(ifs);
    oldWave->set_onedot(onedot);
    ifs.close();
    return 0;
}

//void check_wavefunction(char *filewave)
//{
//    StateInfo oldStateInfo;
//    Wavefunction oldWave = load_wavefunction(filewave, oldStateInfo);
//}
