
#include <stdio.h>

extern "C" {

// creates a forward tree from a parent tree
// Note: assumes binary tree
void makeFtree(int nnodes, int *ptree, int ***ftree)
{
    *ftree = new int* [nnodes];
    int **ftree2 = *ftree;

    // initialize
    for (int i=0; i<nnodes; i++) {
        ftree2[i] = new int [2];
        ftree2[i][0] = -1;
        ftree2[i][1] = -1;
    }

    // populate
    for (int i=0; i<nnodes; i++) {
        int parent = ptree[i];

        if (parent != -1) {
            if (ftree2[parent][0] == -1)
                ftree2[parent][0] = i;
            else
                ftree2[parent][1] = i;
        }
    }
}



void freeFtree(int nnodes, int **ftree)
{
    for (int i=0; i<nnodes; i++)
        delete [] ftree[i];
    delete [] ftree;
}



void printFtree(int nnodes, int **ftree)
{
    for (int i=0; i<nnodes; i++) {
        printf("%2d: %2d %2d\n", i, ftree[i][0], ftree[i][1]);
    }
}



} // extern "C"

