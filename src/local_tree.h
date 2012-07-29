//=============================================================================
// Local trees

#ifndef ARGHMM_LOCAL_TREES_H
#define ARGHMM_LOCAL_TREES_H

// c++ includes
#include <assert.h>
#include <list>
#include <vector>
#include <string.h>
#include <stdio.h>

// arghmm includes
#include "sequences.h"


namespace arghmm {

using namespace std;



// A block within a sequence alignment
class Block 
{
public:
    Block() {}
    Block(int start, int end) : 
        start(start), end(end) {}
    int start;
    int end;

    int length() {
        return end - start;
    }
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

    bool is_null() const
    {
        return recomb_node == -1;
    }

    int recomb_node;
    int recomb_time;
    int coal_node;
    int coal_time;
};



// A node point can be either a recombination point or a coalescing point
class NodePoint
{
public:
    NodePoint(int node, int time) :
        node(node), time(time) {}

    int node;
    int time;
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
    LocalNode(const LocalNode &other) :
        parent(other.parent),
        age(other.age)
    {
        child[0] = other.child[0];
        child[1] = other.child[1];
    }

    inline bool is_leaf() const
    {
        return child[0] == -1;
    }

    inline int add_child(int child_node)
    {
        if (child[0] == -1) {
            child[0] = child_node;
            return 0;
        } else if (child[1] == -1) {
            child[1] = child_node;
            return 1;
        }

        // already have two children
        return -1;
    }

    inline void copy(const LocalNode &other)
    {
        parent = other.parent;
        age = other.age;
        child[0] = other.child[0];
        child[1] = other.child[1];        
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
        capacity(0),
        root(-1),
        nodes(NULL) 
    {}

    LocalTree(int nnodes, int capacity=-1) :
        nnodes(nnodes),
        capacity(capacity),
        root(-1)
    {
        if (capacity < nnodes)
            capacity = nnodes;
        nodes = new LocalNode [capacity];
    }


    LocalTree(int *ptree, int nnodes, int *ages=NULL, int capacity=-1) :
        nodes(NULL) 
    {
        set_ptree(ptree, nnodes, ages, capacity);
    }

    ~LocalTree() {
        if (nodes) {
            delete [] nodes;
            nodes = NULL;
        }
    }

    // initialize a local tree by on a parent array
    void set_ptree(int *ptree, int _nnodes, int *ages=NULL, int _capacity=-1) 
    {
        nnodes = _nnodes;
        if (_capacity >= 0)
            capacity = _capacity;
        if (capacity < nnodes)
            capacity = nnodes;

        // delete existing nodes if they exist
        if (nodes)
            delete [] nodes;
        nodes = new LocalNode [capacity];

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

    
    // Sets the root of the tree by finding node without a parent
    void set_root()
    {
        for (int j=0; j<nnodes; j++) {
            if (nodes[j].parent == -1) {
                root = j;
                break;
            }
        }
    }


    // Sets a new capacity for the allocated data structures
    void set_capacity(int _capacity)
    {
        if (_capacity == capacity)
            return;

        LocalNode *tmp = new LocalNode[_capacity];
        assert(tmp);

        // copy over nodes
        std::copy(nodes, nodes + capacity, tmp);   
        delete [] nodes;

        nodes = tmp;
        capacity = _capacity;
    }


    // Ensures that we have a certain capacity.  If capacity is not >=
    // than _capacity, increase capacity
    void ensure_capacity(int _capacity)
    {
        if (capacity < _capacity)
            set_capacity(_capacity);
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

    
    // Convenience method for accessing nodes
    inline LocalNode &operator[](int name) const
    {
        return nodes[name];
    }


    inline double get_dist(int node, const double *times) const 
    {
        int parent = nodes[node].parent;
        if (parent != -1)
            return times[nodes[parent].age] - times[nodes[node].age];
        else
            return 0.0;
    }


    // add a node to the tree and return its name
    inline int add_node()
    {
        nnodes++;
        if (nnodes > capacity)
            ensure_capacity(2*nnodes);
        return nnodes - 1;
    }


    // Copy tree structure from another tree
    inline void copy(const LocalTree &other)
    {
        // copy tree info
        nnodes = other.nnodes;        
        ensure_capacity(other.capacity);
        root = other.root;
        
        // copy node info
        for (int i=0; i<nnodes; i++)
            nodes[i].copy(other.nodes[i]);
    }


    // get the sibling of a node
    inline int get_sibling(int node) const {
        int parent = nodes[node].parent;
        if (parent == -1)
            return -1;
        int *c = nodes[parent].child;
        if (c[0] == node)
            return c[1];
        else
            return c[0];
    }


    // add a child to a node in the tree
    inline int add_child(int parent, int child) {
        int childi = nodes[parent].add_child(child);
        if (childi == -1)
            return -1;
        
        nodes[child].parent = parent;

        return childi;
    }

    int nnodes;        // number of nodes in tree
    int capacity;      // capacity of nodes array
    int root;          // id of root node
    LocalNode *nodes;  // nodes array
};




// A tree within a set of local trees
//
// Specifically this structure describes the block over which the
// local exists and the SPR operation to the left of the local
// tree.
class LocalTreeSpr
{
public:
     LocalTreeSpr(LocalTree *tree, int *ispr, int blocklen, int *mapping=NULL) :
        tree(tree),
        spr(ispr[0], ispr[1], ispr[2], ispr[3]),
        mapping(mapping),
        blocklen(blocklen)
    {}

     LocalTreeSpr(LocalTree *tree, Spr spr, int blocklen, int *mapping=NULL) :
        tree(tree),
        spr(spr),
        mapping(mapping),
        blocklen(blocklen)
    {}

    
    void clear() {
        if (tree) {
            delete tree;
            tree = NULL;
        }
        
        if (mapping) {
            delete [] mapping;
            mapping = NULL;
        }
    }

    void set_capacity(int _capacity)
    {
        if (tree->capacity == _capacity)
            return;

        tree->set_capacity(_capacity);

        // ensure capacity of mapping
        if (mapping) {
            int *tmp = new int [_capacity];
            assert(tmp);
            
            std::copy(mapping, mapping + _capacity, tmp);   
            delete [] mapping;
            
            mapping = tmp;
        }
    }

    void ensure_capacity(int _capacity)
    {
        if (tree->capacity < _capacity)
            set_capacity(_capacity);
    }

    
    LocalTree *tree;  // local tree
    Spr spr;          // SPR operation to the left of local tree
    int *mapping;     // node mapping between previous tree and this tree
    int blocklen;     // length of sequence block
};



// A set of local trees that together specify an ARG
//
// This structure specifies a list of local trees, their blocks, and
// SPR operations.  Together this specifies an ancestral recombination
// graph (ARG), which because of our assumptions is an SMC-style ARG.
class LocalTrees 
{
public:
    LocalTrees(int start_coord=0, int end_coord=0, int nnodes=0) :
        start_coord(start_coord),
        end_coord(end_coord),
        nnodes(nnodes) {}
    LocalTrees(int **ptrees, int**ages, int **isprs, int *blocklens,
               int ntrees, int nnodes, int capacity=-1, int start=0);
    ~LocalTrees() 
    {
        clear();
    }

    // iterator for the local trees
    typedef list<LocalTreeSpr>::iterator iterator;
    typedef list<LocalTreeSpr>::reverse_iterator reverse_iterator;
    typedef list<LocalTreeSpr>::const_iterator const_iterator;
    typedef list<LocalTreeSpr>::const_reverse_iterator const_reverse_iterator;


    // Returns iterator for first local tree
    iterator begin() { return trees.begin(); }
    reverse_iterator rbegin() { return trees.rbegin(); }
    LocalTreeSpr &front() { return trees.front(); }

    const_iterator begin() const { return trees.begin(); }
    const_reverse_iterator rbegin() const { return trees.rbegin(); }
    const LocalTreeSpr &front() const { return trees.front(); }


    // Returns the ending iterator
    iterator end() { return trees.end(); }
    reverse_iterator rend() { return trees.rend(); }
    LocalTreeSpr &back() { return trees.back(); }

    const_iterator end() const { return trees.end(); }
    const_reverse_iterator rend() const { return trees.rend(); }
    const LocalTreeSpr &back() const { return trees.back(); }


    // Returns number of leaves
    inline int get_num_leaves() const
    {
        return (nnodes + 1) / 2;
    }

    inline int length() const
    {
        return end_coord - start_coord;
    }

    inline int get_num_trees() const
    {
        return trees.size();
    }

    // Copy tree structure from another tree
    void copy(const LocalTrees &other);

    // deallocate local trees
    void clear()
    {
        for (iterator it=begin(); it!=end(); it++)
            it->clear();
        trees.clear();
    }

    // make trunk genealogy
    void make_trunk(int start, int end, int capacity=-1)
    {
        nnodes = 1;
        start_coord = start;
        end_coord = end;

        clear();

        int ptree[] = {-1};
        int ages[] = {0};
        LocalTree *tree = new LocalTree(ptree, 1, ages, capacity);
        trees.push_back(
            LocalTreeSpr(tree, Spr(-1, -1, -1, -1), end - start, NULL));
        set_default_seqids();
    }

    void set_default_seqids()
    {
        const int nleaves = get_num_leaves();
        seqids.clear();
        for (int i=0; i<nleaves; i++)
            seqids.push_back(i);
    }


    int start_coord;           // start coordinate of whole tree list
    int end_coord;             // end coordinate of whole tree list
    int nnodes;                // number of nodes in each tree
    list<LocalTreeSpr> trees;  // linked list of local trees
    vector<int> seqids;        // mapping from tree leaves to sequence ids
};


// count the lineages in a tree
void count_lineages(const LocalTree *tree, int ntimes,
                    int *nbranches, int *nrecombs, int *ncoals);
void count_lineages_internal(const LocalTree *tree, int ntimes,
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
    inline void count(const LocalTree *tree, bool internal=false) {
        if (internal)
            count_lineages_internal(tree, ntimes, nbranches, nrecombs, ncoals);
        else
            count_lineages(tree, ntimes, nbranches, nrecombs, ncoals);
    }

    int ntimes;      // number of time points
    int *nbranches;  // number of branches per time slice
    int *nrecombs;   // number of recombination points per time slice
    int *ncoals;     // number of coalescing points per time slice
};

//=============================================================================
// trees functions

void map_congruent_trees(const LocalTree *tree1, const int *seqids1,
                         const LocalTree *tree2, const int *seqids2, 
                         int *mapping);

bool remove_null_spr(LocalTrees *trees, LocalTrees::iterator it);
void remove_null_sprs(LocalTrees *trees);

LocalTrees *partition_local_trees(LocalTrees *trees, int pos,
                                  LocalTrees::iterator it, int it_start);
LocalTrees *partition_local_trees(LocalTrees *trees, int pos);
void append_local_trees(LocalTrees *trees, LocalTrees *trees2);

void uncompress_local_trees(LocalTrees *trees, 
                            const SitesMapping *sites_mapping);
void compress_local_trees(LocalTrees *trees, const SitesMapping *sites_mapping);
void assert_uncompress_local_trees(LocalTrees *trees, 
                                   const SitesMapping *sites_mapping);


//=============================================================================
// tree functions

void apply_spr(LocalTree *tree, const Spr &spr);
double get_treelen(const LocalTree *tree, const double *times, int ntimes,
                    bool use_basal=true);
double get_treelen_internal(const LocalTree *tree, const double *times, 
                            int ntimes);
double get_treelen_branch(const LocalTree *tree, const double *times, 
                          int ntimes, int node, int time, double treelen=-1.0, 
                          bool use_basal=true);
double get_basal_branch(const LocalTree *tree, const double *times, 
                        int ntimes, int node, int time);

inline void make_node_mapping(const int *ptree, int nnodes, int recomb_node, 
                              int *mapping)
{
    for (int j=0; j<nnodes; j++)
        mapping[j] = j;
            
    // parent of recomb is broken and therefore does not map
    const int parent = ptree[recomb_node];
    mapping[parent] = -1;
}

void print_local_tree(const LocalTree *tree, FILE *out=stdout);
void print_local_trees(const LocalTrees *trees, FILE *out=stdout);

//=============================================================================
// input and output

void write_newick_tree(FILE *out, const LocalTree *tree, 
                       const char *const *names,
                       const double *times, int depth, bool oneline);
bool write_newick_tree(const char *filename, const LocalTree *tree, 
                       const char *const *names, const double *times, 
                       bool oneline);
void write_local_trees(FILE *out, const LocalTrees *trees, 
                       const char *const *names, const double *times);
bool write_local_trees(const char *filename, const LocalTrees *trees, 
                       const char *const *names, const double *times);

void write_local_trees(FILE *out, const LocalTrees *trees, 
                       const Sequences &seqs, const double *times);
bool write_local_trees(const char *filename, const LocalTrees *trees, 
                       const Sequences &seqs, const double *times);


//=============================================================================
// assert functions

bool assert_tree_postorder(const LocalTree *tree, const int *order);
bool assert_tree(const LocalTree *tree);
bool assert_spr(const LocalTree *last_tree, const LocalTree *tree, 
                const Spr *spr, const int *mapping);
bool assert_trees(const LocalTrees *trees);




} // namespace arghmm

#endif // ARGHMM_LOCAL_TREES
