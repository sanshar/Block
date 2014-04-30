/*
 * File:
 */

#include "config.h"
#include "rotationmat.h"
#include "density.h"
#include "spinblock.h"
#include "guess_wavefunction.h"
#include "MatrixBLAS.h"
#include <boost/serialization/vector.hpp>
#include "operatorfunctions.h"

using namespace SpinAdapted;

/*
 * ref rotationmat.C
 *     SpinAdapted::LoadRotationMatrix
 *     guess_wavefunction.C
 */
int save_rotmat(char *filerotmat, std::vector<Matrix> *mat)
{
    std::ofstream ofs(filerotmat, std::ios::binary);
    boost::archive::binary_oarchive save_mat(ofs);
    save_mat << *mat;
    ofs.close();
    return 0;
}

int load_rotmat(char *filerotmat, std::vector<Matrix> *mat)
{
    std::ifstream ifs(filerotmat, std::ios::binary);
    boost::archive::binary_iarchive load_mat(ifs);
    load_mat >> *mat;
    ifs.close();
    return 0;
}

/*
 * ref density.C
 *     DensityMatrix::makedensitymatrix
 */
int update_rotmat(std::vector<Matrix> *rotateMatrix,
                  Wavefunction *wfn, SpinBlock *sys, SpinBlock *big,
                  int keptstates, int keptqstates, double noise)
{
    DensityMatrix tracedMatrix(sys->braStateInfo);
    tracedMatrix.allocate(sys->braStateInfo);
    std::vector<Wavefunction> wfns;
    wfns.push_back(*wfn);
    std::vector<double> wave_weights(1,1.l);
    double additional_noise = 0;
    bool warmup = 0;
    tracedMatrix.makedensitymatrix(wfns, *big, wave_weights, noise,
                                   additional_noise, warmup);

    /* tracedMatrix.add_onedot_noise is private function
    operatorfunctions::MultiplyProduct(*wfn, Transpose(*wfn),
                                       tracedMatrix, 1);

    // note tracedMatrix.makedensitymatrix and add_onedot_noise, dmrginp
    if (noise > 1.0e-14) {
        // In this call, add_onedot_noise only modify tracedMatrix
        std::vector<Wavefunction> wfns;
        wfns.push_back(*wfn);
        tracedMatrix.add_onedot_noise(wfns, *big, noise);
    } */

    double error = makeRotateMatrix(tracedMatrix, *rotateMatrix,
                                    keptstates, keptqstates);
    return 0;
}

/*
 * ref Sweep::Startup
 */
int guess_rotmat(std::vector<Matrix> *rotateMatrix, SpinBlock *newSystem,
                 int keep_states)
{
    int nquanta = newSystem->get_stateInfo().quanta.size();
    std::vector<DiagonalMatrix > energies(nquanta);
    rotateMatrix->resize(nquanta);
    DensityMatrix transformmatrix;
    transformmatrix.allocate(newSystem->get_stateInfo());
    SpinQuantum q(0,SpinSpace(0),IrrepSpace(0));

    double minval = 1e12;
    for (int i=0; i<nquanta; i++) {
        diagonalise(newSystem->get_op_rep(HAM,q)->operator_element(i,i), energies[i], transformmatrix(i,i));
        for (int j=0; j<energies[i].Nrows(); j++)
            if (minval > energies[i](j+1))
                minval = energies[i](j+1);
    }
    for (int i=0; i<nquanta; i++) {
        for (int j=0; j<energies[i].Nrows(); j++)
            energies[i](j+1) = 1.0/(energies[i](j+1)-minval+1);
    }

    vector<pair<int, int> > inorderwts;
    vector<vector<int> > wtsbyquanta;

    sort_weights(energies, inorderwts, wtsbyquanta);

    // make transformation matrix by various algorithms
    int keptstates = keep_states/2;
    int keptqstates = keep_states-keptstates;
    int totalstatesbydm = min(static_cast<int>(inorderwts.size()), keptstates);
    int totalstatesbyquanta = min(static_cast<int>(inorderwts.size()), keptstates + keptqstates) - totalstatesbydm;
    if (totalstatesbyquanta < 0) totalstatesbyquanta = 0;

    //pout << "\t\t\t total states using dm and quanta " << totalstatesbydm << " " << totalstatesbyquanta << endl
    double error = assign_matrix_by_dm(*rotateMatrix, energies, transformmatrix,
                                       inorderwts, wtsbyquanta, totalstatesbydm,
                                       totalstatesbyquanta, newSystem->size(),
                                       2*totalstatesbydm);
    //pout << "\t\t\t Total discarded weight "<<error<<endl;
    return 0;
}
