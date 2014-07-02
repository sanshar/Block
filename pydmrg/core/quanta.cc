/*
 * File:
 */

#include "IrrepSpace.h"
#include "SpinQuantum.h"

using namespace SpinAdapted;

int x_SpinQuantum_irrep(SpinQuantum *sq)
{
    return sq->get_symm().getirrep();
}

