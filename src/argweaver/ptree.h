#ifndef DLCOAL_PTREE_H
#define DLCOAL_PTREE_H

extern "C" {

void makeFtree(int nnodes, int *ptree, int ***ftree);

void freeFtree(int nnodes, int **ftree);

void printFtree(int nnodes, int **ftree);

} // extern "C"

#endif // DLCOAL_PTREE_H
