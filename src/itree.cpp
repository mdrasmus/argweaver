

#include "itree.h"


namespace dlcoal {

// creates a int tree from a parent tree
// Note: assumes binary tree
intnode *make_itree(int nnodes, int *ptree)
{
    intnode *itree = new intnode [nnodes];
    
    // initialize
    for (int i=0; i<nnodes; i++) {
        itree[i].parent = ptree[i];
        itree[i].child[0] = -1;
        itree[i].child[1] = -1;
    }
    
    // populate
    for (int i=0; i<nnodes; i++) {
        int parent = ptree[i];
        
        if (parent != -1) {
            if (itree[parent].child[0] == -1)
                itree[parent].child[0] = i;
            else
                itree[parent].child[1] = i;
        }
    }

    return itree;
}


void free_itree(intnode *itree)
{
    delete [] itree;
}


} // namespace dlcoal
