/*
Copyright (c) 2013, Garnet K.-L. Chan

This program is integrated in Molpro with the permission of 
Garnet K.-L. Chan
*/

#ifndef FIEDLER_H
#define FIEDLER_H
#include <vector>
#include <newmat.h>

std::vector<int> get_fiedler(string& dumpname, ifstream& dumpFile, bool simple = false);
std::vector<int> fiedler_reorder(const SymmetricMatrix& m);

#endif
