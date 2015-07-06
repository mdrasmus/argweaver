#ifndef DLCOAL_ITREE_H
#define DLCOAL_ITREE_H

namespace dlcoal {

extern "C" {


struct intnode
{
    int parent;
    int child[2];
};

intnode *make_itree(int nnodes, int *ptree);
void free_itree(intnode *itree);


}

} // namespace dlcoal

#endif // DLCOAL_ITREE_H


