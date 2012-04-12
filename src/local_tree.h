//=============================================================================
// Local trees

#ifndef ARGHMM_LOCAL_TREES_H
#define ARGHMM_LOCAL_TREES_H

// c++ includes
#include <assert.h>
#include <list>
#include <vector>
#include <string.h>


namespace arghmm {

using namespace std;



// A block within a sequence alignment
class Block 
{
public:
    Block(int start, int end) : 
        start(start), end(end) {}
    int start;
    int end;
};


// A Subtree Pruning and Regrafting operation
class Spr
{
public:
    Spr() {}
    Spr(int recomb_node, int recomb_time, 
        int coal_node, int coal_time) :
        recomb_node(recomb_node), recomb_time(recomb_time),
        coal_node(coal_node), coal_time(coal_time) {}

    // sets the SPR to a null value
    void set_null()
    {
        recomb_node = -1;
        recomb_time = -1;
        coal_node = -1;
        coal_time = -1;
    }

    bool is_null()
    {
        return recomb_node == -1;
    }

    int recomb_node;
    int recomb_time;
    int coal_node;
    int coal_time;
};



// A node in a local tree
class LocalNode 
{
public:
    LocalNode() {}
    LocalNode(int parent, int left_child, int right_child, int age=-1) :
        parent(parent), age(age)
    {
        child[0] = left_child;
        child[1] = right_child;
    }

    inline bool is_leaf()
    {
        return child[0] == -1;
    }

    int parent;
    int child[2];
    int age;
};


// A local tree in a set of local trees
//
//   Leaves are always listed first in nodes array
//
class LocalTree
{
public:
    LocalTree() :
        nnodes(0),
        root(-1),
        nodes(NULL) 
    {}

    LocalTree(int *ptree, int nnodes, int *ages=NULL) :
        nodes(NULL) 
    {
        set_ptree(ptree, nnodes, ages);
    }

    ~LocalTree() {
        if (nodes)
            delete [] nodes;
    }

    // initialize a local tree by on a parent array
    void set_ptree(int *ptree, int _nnodes, int *ages=NULL) 
    {
        nnodes = _nnodes;

        // delete existing nodes if they exist
        if (nodes)
            delete [] nodes;
        nodes = new LocalNode [nnodes];

        // populate parent pointers
        for (int i=0; i<nnodes; i++) {
            nodes[i].parent = ptree[i];
            nodes[i].child[0] = -1;
            nodes[i].child[1] = -1;
        }

        // initialize age if given
        if (ages)
            for (int i=0; i<nnodes; i++)
                nodes[i].age = ages[i];

    
        // populate children pointers
        for (int i=0; i<nnodes; i++) {
            const int parent = ptree[i];
        
            if (parent != -1) {
                int *child = nodes[parent].child;
                if (child[0] == -1)
                    child[0] = i;
                else
                    child[1] = i;
            } else {
                root = i;
            }
        }
    }

    // Returns the postorder traversal of the nodes
    void get_postorder(int *order) const
    {
        char visit[nnodes];
        int i;

        // initialize array of number of visits to a node
        for (i=0; i<nnodes; i++)
            visit[i] = 0;

        // list leaves first
        for (i=0; i<nnodes; i++) {
            if (!nodes[i].is_leaf())
                break;
            order[i] = i;
        }
        
        // add the remaining nodes
        int end = i;
        for (i=0; i<nnodes; i++) {
            int parent = nodes[order[i]].parent;
            if (parent != -1) {
                visit[parent]++;
                
                // add parent to queue if both children have been seen
                if (visit[parent] == 2)
                    order[end++] = parent;
            }
        }
    }


    // Returns the number of leaves in the tree
    inline int get_num_leaves() const
    {
        return (nnodes + 1) / 2;
    }

    
    // convenience method for accessing nodes
    inline LocalNode &operator[](int name) 
    {
        return nodes[name];
    }


    int nnodes;
    int root;
    LocalNode *nodes;
};




// A tree local within a set of local trees
//
// Specifically this structure describes the block over which the
// local exists and it the SPR operation to the left of the local
// tree.
class LocalTreeSpr
{
public:
    LocalTreeSpr(int start, int end, LocalTree *tree, int *ispr) :
        block(start, end),
        tree(tree),
        spr(ispr[0], ispr[1], ispr[2], ispr[3])
    {}

    LocalTreeSpr(int start, int end, LocalTree *tree, Spr spr) :
        block(start, end),
        tree(tree),
        spr(spr)
    {}

    Block block;
    LocalTree *tree;
    Spr spr;
};



// A set of local trees that together specify an ARG
//
// This structure specifies a list of local trees, their blocks, and
// SPR operations.  Together this specifies an ancestral recombination
// graph (ARG), which because of our assumptions is an SMC-style ARG.
class LocalTrees 
{
public:
    LocalTrees(int **ptrees, int**ages, int **isprs, int *blocklens,
               int ntrees, int nnodes, int start=0) : 
        start_coord(start),
        nnodes(nnodes)
    {
        // copy data
        int pos = start;
        for (int i=0; i<ntrees; i++) {
            end_coord = pos + blocklens[i];
            trees.push_back(
                LocalTreeSpr(pos, end_coord,
                             new LocalTree(ptrees[i], nnodes, ages[i]),
                             isprs[i]));
            pos = end_coord;
        }

    }
    ~LocalTrees() {}

    // iterator for the local trees
    typedef list<LocalTreeSpr>::iterator iterator;

    // Returns iterator for first local tree
    iterator begin()
    {
        return trees.begin();
    }

    // Returns the ending iterator
    iterator end()
    {
        return trees.end();
    }


    int start_coord;
    int end_coord;
    int nnodes;
    list<LocalTreeSpr> trees;
};


void count_lineages(LocalTree *tree, int ntimes,
                    int *nbranches, int *nrecombs, int *ncoals);

// A structure that stores the number of lineages within each time segment
class LineageCounts
{
public:
    LineageCounts(int ntimes) :
        ntimes(ntimes)
    {
        nbranches = new int [ntimes];
        nrecombs = new int [ntimes];
        ncoals = new int [ntimes];
    }

    ~LineageCounts()
    {
        delete [] nbranches;
        delete [] nrecombs;
        delete [] ncoals;
    }

    // Counts the number of lineages for a tree
    inline void count(LocalTree *tree) {
        count_lineages(tree, ntimes, nbranches, nrecombs, ncoals);
    }

    int ntimes;
    int *nbranches;
    int *nrecombs;
    int *ncoals;
};


// tree methods
double get_treelen(const LocalTree *tree, const double *times, int ntimes);




} // namespace arghmm

#endif // ARGHMM_LOCAL_TREES
