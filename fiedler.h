/*
Copyright (c) 2013, Garnet K.-L. Chan

This program is integrated in Molpro with the permission of 
Garnet K.-L. Chan
*/

#ifndef FIEDLER_HEADER_H
#define FIEDLER_HEADER_H
#include <vector>

class SymmetricMatrix;

std::vector<int> fiedler_reorder(const SymmetricMatrix& m);

#endif
