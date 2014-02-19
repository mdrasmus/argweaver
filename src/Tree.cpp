/*=============================================================================

  Matt Rasmussen
  Copyright 2010-2011

  Tree datastructure

=============================================================================*/

// c++ headers
#include <assert.h>
#include <stdio.h>

// spidir headers
#include "Tree.h"
#include "common.h"
#include "parsing.h"
#include "logging.h"
#include <iostream>
#include <cstring>
#include <string>
#include <map>

namespace spidir {

using namespace argweaver;

bool isNewickChar(char c) {
    static vector<int> vals(8);
    static int called=0;
    int i=(int)c;
    if (c >= (int)'a' || (c>=(int)'0' && c <=(int)'9')) return false;
    if (called==0) {
        called=1;
        vals[0]=(int)'(';
        vals[1]=(int)')';
        vals[2]=(int)',';
        vals[3]=(int)':';
        vals[4]=(int)'#';
        vals[5]=(int)'!';
        vals[6]=(int)'[';
        vals[7]=(int)']';
    }
    return (i==vals[0] || i==vals[1] || i==vals[2] || i==vals[3] ||
            i==vals[4] || i==vals[5] || i==vals[6] || i==vals[7]);
}

//create a tree from a newick string
Tree::Tree(string newick)
{
    int len = newick.length();
    Node *node = NULL;
    vector <int> stack;

    nnodes=0;
    for (int i=0; i < len; i++)
        if (newick[i]=='(') nnodes++;
    nnodes += (nnodes+1);  //add in leaves
    nodes.setCapacity(nnodes);
    for (int i=0; i < nnodes; i++) {
        nodes[i] = new Node();
    }
    root = nodes[0];
    root->name = 0;
    stack.push_back(0);
    nnodes = 1;

    for (int i=0; i < len; i++) {
        switch (newick[i]) {
        case ',':
            stack.pop_back();
        case '(':
            node = nodes[nnodes];
            if (stack.size()==0) {
                printError("bad newick: error parsing tree");
            } else {
                node->parent = nodes[stack.back()];
            }
            stack.push_back(nnodes);
            node->name = nnodes++;
            break;
        case ')': {
            stack.pop_back();
            node = nodes[stack.back()];
            break;
        }
        case ':':  { //optional dist next
            int j=i+1;
            while (j < len && !isNewickChar(newick[j]))
                j++;
            if (sscanf(&newick[i+1], "%f", &node->dist) != 1) {
                printError("bad newick: error reading distance");
                printf("&newick[%i+1]=%s\n", i, &newick[i+1]);
                printf("node->dist=%f\n", node->dist);
                exit(1);
            }
            i=j-1;
            break;
        }
        case '[': { // comment next; save in nhx field
            int count=1;
            int j=i+1;
            while (count != 0) {
                if (j==len) {
                    printError("bad newick: no closing bracket in NHX comment");
                    break;
                }
                if (newick[j]==']') count--;
                else if (newick[j]=='[') count++;
                j++;
            }
            node->nhx = newick.substr(i+1, j-i-2);
            i=j-1;
            break;
        }
        case ';':
            break;
        default:
           int j=i+1;
           while (j < len && !isNewickChar(newick[j]))
               j++;
           if (node->longname.length() > 0) {
               printError("bad newick format; got multiple names for a node");
               i=len;
               break;
           }
           node->longname = newick.substr(i, j-i);
           trim(node->longname);
           i=j-1;
           break;
        }
    }
    if (node != root) {
        printError("bad newick format: did not end with root");
    }
    //All done, now fill in children
    for (int i=0; i < nnodes; i++) {
        if (nodes[i]->parent != NULL) {
            nodes[i]->parent->addChild(nodes[i]);
        }
    }
}

void Tree::print_newick_recur(FILE *f, Node *n, bool internal_names,
                              char *branch_format_str,
                              bool show_nhx, bool oneline) {
    if (n->nchildren > 0) {
        fprintf(f, "(");
        print_newick_recur(f, n->children[0], internal_names,
                           branch_format_str, show_nhx, oneline);
        for (int i=1; i < n->nchildren; i++) {
            fprintf(f, ",");
            print_newick_recur(f, n->children[i], internal_names,
                               branch_format_str, show_nhx, oneline);
        }
        fprintf(f, ")");
        if (internal_names) fprintf(f, "%s", n->longname.c_str());
    } else {
        fprintf(f, "%s", n->longname.c_str());
    }

    if (branch_format_str != NULL && n->parent != NULL) {
        fprintf(f, ":");
        fprintf(f, branch_format_str, n->dist);
    }
    if (show_nhx && n->nhx.length() != 0)
        fprintf(f, "[%s]", n->nhx.c_str());
    if (!oneline && n->nchildren > 0) fprintf(f, "\n");
}


// return a copy of the tree
Tree *Tree::copy()
{
    Tree *tree2 = new Tree(nnodes);
    Node **nodes2 = tree2->nodes;

    for (int i=0; i<nnodes; i++) {
        nodes2[i]->setChildren(nodes[i]->nchildren);
        nodes2[i]->name = i;
        nodes2[i]->dist = nodes[i]->dist;
        nodes2[i]->longname = nodes[i]->longname;
    }

    for (int i=0; i<nnodes; i++) {
        for (int j=0; j<nodes[i]->nchildren; j++) {
            Node *child = nodes[i]->children[j];
            if (child)
                nodes2[i]->children[j] = nodes2[child->name];
            else
                nodes2[i]->children[j] = NULL;
        }
        Node *parent = nodes[i]->parent;
        if (parent)
            nodes2[i]->parent = nodes2[parent->name];
        else
            nodes2[i]->parent = NULL;
    }

    tree2->root = nodes2[root->name];

    return tree2;
}


void Tree::prune(set<string> leafs, bool allBut) {
    ExtendArray<Node*> postnodes;
    ExtendArray<Node*> newnodes;
    getTreePostOrder(this, &postnodes);
    vector<bool> is_leaf(postnodes.size());

    for (int i=0; i < postnodes.size(); i++) {
        is_leaf[i] = (postnodes[i]->nchildren == 0);
    }

    for (int i=0; i < postnodes.size(); i++) {
        if (postnodes[i]->nchildren == 0) {
            int prune;
            if (!is_leaf[i]) {
                //in this case, node was not originally a leaf, but now
                // has no children, so should be pruned
                prune = true;
            } else {
                prune = (leafs.find(postnodes[i]->longname) != leafs.end());
                if (allBut) prune=!prune;
            }
            if (prune) {
                Node *parent = postnodes[i]->parent;
                if (parent != NULL) {
                    int j, maxj=parent->nchildren;
                    for (j=0; j < maxj; j++) {
                        if (parent->children[j] == postnodes[i]) {
                            parent->children[j] = parent->children[parent->nchildren-1];
                            parent->nchildren--;
                            delete postnodes[i];
                            break;
                        }
                    }
                    if (j == maxj) {
                        fprintf(stderr, "error in tree.prune(): didn't find child in parent node\n");
                        exit(-1);
                    }
                } else {  //entire tree has been pruned!
                    for (i=0; i < nnodes; i++)
                        delete nodes[i];
                    nnodes=0;
                    root=NULL;
                }
            } else {
                newnodes.append(postnodes[i]);
            }
        } else if (postnodes[i]->nchildren == 1) {
            if (postnodes[i] == root) {
                root = postnodes[i]->children[0];
            } else {
                Node *parent = postnodes[i]->parent;
                int j, maxj=parent->nchildren;
                for (j=0; j < maxj; j++) {
                    if (parent->children[j] == postnodes[i]) {
                        parent->children[j] = postnodes[i]->children[0];
                        postnodes[i]->children[0]->dist += postnodes[i]->dist;
                        postnodes[i]->children[0]->parent = parent;
                        delete postnodes[i];
                        break;
                    }
                }
                if (j == maxj) {
                    fprintf(stderr, "error in tree.prune(): didn't find child in parent node2\n");
                    exit(-1);
                }
            }
        } else {
            newnodes.append(postnodes[i]);
        }
    }

    nodes.clear();
    for (int i=0; i < newnodes.size(); i++) {
        nodes.append(newnodes[i]);
        nodes[i]->name = i;
        nodes[i]->nhx.clear();
    }
    nnodes = nodes.size();
}



// assumes both trees have same number of nodes
// and have same leaves
void Tree::setTopology(Tree *other)
{
    assert(nnodes == other->nnodes);
    Node **onodes = other->nodes;

    for (int i=0; i<nnodes; i++) {
        Node *node = nodes[i];
        Node *onode = onodes[i];

        if (onode->parent)
            node->parent = nodes[onode->parent->name];
        else
            node->parent = NULL;


        if (node->isLeaf()) {
            assert(onode->isLeaf());
        } else {
            // copy child structure
            nodes[i]->setChildren(onodes[i]->nchildren);
            for (int j=0; j<onodes[i]->nchildren; j++) {
                node->children[j] = nodes[onode->children[j]->name];
            }
        }
    }
}


// root tree by a new branch/node
void Tree::reroot(Node *newroot, bool onBranch)
{
    // handle trivial case, newroot is root
    if (root == newroot ||
        (onBranch &&
         root->nchildren == 2 &&
         (root->children[0] == newroot ||
          root->children[1] == newroot)))
        return;


    // determine where to stop ascending
    Node *oldroot = root;
    Node *stop1=NULL, *stop2=NULL;

    if (isRooted()) {
        stop1 = root->children[0];
        stop2 = root->children[1];
    } else {
        stop1 = root;
    }

    // start the reversal
    Node *ptr1 = NULL, *ptr2 = NULL;
    float nextDist = 0;
    float rootdist;

    if (onBranch) {
        if (isRooted()) {
            // just need to stick current root somewhere else
            Node *other = newroot->parent;
            rootdist = stop1->dist + stop2->dist;

            oldroot->children[0] = newroot;
            oldroot->children[1] = other;
            newroot->parent = oldroot;
            newroot->dist /= 2.0;

            ptr1 = other;

            int oldchild = find_array(ptr1->children, ptr1->nchildren, newroot);
            assert(oldchild != -1);

            // prepare for reversing loop
            ptr1->children[oldchild] = oldroot;
            ptr2 = oldroot;
            nextDist = newroot->dist;
        } else {
            // need to add a new node to be root
            // TODO: not implemented
            assert(0);
        }
    } else {
        if (isRooted()) {
            // need to remove the root node, and make tribranch
            // TODO: not implemented
            assert(0);
        } else {
            // just need to swap node positions
            // TODO: not implemented
            assert(0);
        }
    }


    // reverse parent child relationships
    while (ptr1 != stop1 && ptr1 != stop2) {
        int oldchild = find_array(ptr1->children, ptr1->nchildren, ptr2);
        assert(oldchild != -1);

        Node *next = ptr1->parent;

        // ptr1 is now fixed
        ptr1->children[oldchild] = next;
        ptr1->parent = ptr2;

        // swap distances
        float tmpdist = ptr1->dist;
        ptr1->dist = nextDist;
        nextDist = tmpdist;

        // move pointers
        ptr2 = ptr1;
        ptr1 = next;
    }


    // handle last two nodes
    if (stop2 != NULL) {
        // make stop1 parent of stop2
        if (stop2 == ptr1) {
            Node *tmp = stop1;
            stop1 = ptr1;
            stop2 = tmp;
        }
        assert(ptr1 == stop1);

        int oldchild = find_array(stop1->children, stop1->nchildren, ptr2);
        stop1->children[oldchild] = stop2;
        stop1->parent = ptr2;
        stop1->dist = nextDist;
        stop2->parent = stop1;
        stop2->dist = rootdist;
    } else {
        assert(0);
    }


    // renumber nodes
    // - all leaves don't change numbers
    assert(root->name = nnodes-1);
}


void Tree::reroot(Node *node1, Node *node2)
{
    // determine new root
    Node *newroot;
    if (node1->parent == node2)
        newroot = node1;
    else if (node2->parent == node1)
        newroot = node2;
    else if (node1->parent == root ||
             node2->parent == root)
        // do nothing
        return;
    else
        // not a valid branch
        assert(0);

    reroot(newroot);
}





// store a hash key representing the topology into the key array
// key is a parent tree representation where the internal nodes are
// given a consistent numbering
void Tree::hashkey(int *key)
{
    // get post order of nodes
    ExtendArray<Node*> postnodes;
    getTreePostOrder(this, &postnodes);

    // order children
    ExtendArray<int> ordering(nnodes);
    for (int i=0; i<postnodes.size(); i++)
    {
        Node *node=postnodes[i];

        if (node->isLeaf()) {
            ordering[node->name] = node->name;
        } else {
            // propogate the min order to the parent
            int minorder = ordering[node->children[0]->name];
            for (int j=1; j<node->nchildren; j++) {
                int order = ordering[node->children[j]->name];
                if (order < minorder)
                    minorder = order;
            }
            ordering[node->name] = minorder;
        }
    }

    // get a sorted post ordering of nodes
    ExtendArray<Node*> sortpostnodes;
    getTreeSortedPostOrder(this, &sortpostnodes, ordering);

    // generate a unique key for this topology
    // postfix notation for a tree
    // ((A,B),C) is represented as
    // A, B, -1, C, -1
    for (int i=0; i<sortpostnodes.size(); i++) {
        Node *node = sortpostnodes[i];

        if (node->isLeaf())
            key[i] = node->name;
        else
            key[i] = -1;
    }
}


bool Tree::sameTopology(Tree *other)
{
    if (other->nnodes != nnodes)
        return false;

    typedef ExtendArray<int> TopologyKey;
    TopologyKey key1(nnodes);
    TopologyKey key2(other->nnodes);

    hashkey(key1);
    other->hashkey(key2);

    for (int i=0; i<nnodes; i++) {
        if (key1[i] != key2[i])
            return false;
    }
    return true;
}


void Tree::reorderLeaves(string *order)
{
    // count the leaves in the tree
    int nleaves = 0;
    for (int i=0; i<nnodes; i++)
        if (nodes[i]->isLeaf())
            nleaves++;

    ExtendArray<Node*> tmp(nleaves);

    // rename leaves
    for (int i=0; i<nleaves; i++) {
        bool found = false;
        for (int j=0; j<nleaves; j++) {
            if (nodes[i]->longname == order[j]) {
                found = true;
                nodes[i]->name = j;
                tmp[j] = nodes[i];
                break;
            }
        }
        assert(found);
    }

    // reorder leaves by name
    for (int i=0; i<nleaves; i++)
        nodes[i] = tmp[i];
}


// assert that the tree datastructure is self-consistent
bool Tree::assertTree()
{
    if (root == NULL) {
        fprintf(stderr, "root == NULL\n");
        return false;
    }
    if (nnodes != nodes.size()) {
        fprintf(stderr, "nnodes != nodes.size()\n");
        return false;
    }
    if (root->parent != NULL) {
        fprintf(stderr, "root->parent != NULL\n");
        return false;
    }
    /*if (root->name != nnodes - 1) {
        fprintf(stderr, "root->name != nnodes - 1\n");
        return false;
        }*/

    bool leaves = true;
    for (int i=0; i<nnodes; i++) {
        //printf("assert %d\n", i);
        if (nodes[i] == NULL) {
            fprintf(stderr, "nodes[i] == NULL\n");
            return false;
        }

        // names are correct
        if (nodes[i]->name != i) {
            fprintf(stderr, "nodes[i]->name != i\n");
            return false;
        }

        // do leaves come first
        if (nodes[i]->isLeaf()) {
            if (!leaves) {
                fprintf(stderr, "!leaves\n");
                return false;
            }
        } else
            leaves = false;

        // check parent child pointers
        for (int j=0; j<nodes[i]->nchildren; j++) {
            //printf("assert %d %d\n", i, j);
            if (nodes[i]->children[j] == NULL) {
                fprintf(stderr, "nodes[i]->children[j] == NULL\n");
                return false;
            }
            //printf("assert %d %d parent\n", i, j);
            if (nodes[i]->children[j]->parent != nodes[i]) {
                fprintf(stderr, "nodes[i]->children[j]->parent != nodes[i]\n");
                return false;
            }
        }
    }

    //printf("done\n");

    return true;
}



void getTreeSortedPostOrder(Tree *tree, ExtendArray<Node*> *nodes,
                      int *ordering, Node *node)
{
    if (!node)
        node = tree->root;

    // make a child index array
    int childperm[node->nchildren];
    int childorder[node->nchildren];
    for (int i=0; i<node->nchildren; i++) {
        childperm[i] = i;
        childorder[i] = ordering[node->children[i]->name];
    }

    // sort index array by order
    ranksort(childperm, childorder, node->nchildren);

    // recurse
    for (int i=0; i<node->nchildren; i++)
        getTreeSortedPostOrder(tree, nodes, ordering, node->children[childperm[i]]);

    // record post-process
    nodes->append(node);
}

void getTreePostOrder(Tree *tree, ExtendArray<Node*> *nodes, Node *node)
{
    if (!node)
        node = tree->root;

    // recurse
    for (int i=0; i<node->nchildren; i++)
        getTreePostOrder(tree, nodes, node->children[i]);

    // record post-process
    nodes->append(node);
}

void getTreePreOrder(Tree *tree, ExtendArray<Node*> *nodes, Node *node)
{
    if (!node)
        node = tree->root;

    // record pre-process
    nodes->append(node);

    // recurse
    for (int i=0; i<node->nchildren; i++)
        getTreePreOrder(tree, nodes, node->children[i]);
}





//=============================================================================
// primitive input/output


void printFtree(int nnodes, int **ftree)
{
    for (int i=0; i<nnodes; i++) {
        printf("%2d: %2d %2d\n", i, ftree[i][0], ftree[i][1]);
    }
}


// write out the names of internal nodes
void printTree(Tree *tree, Node *node, int depth)
{
    if (node == NULL) {
        if (tree->root != NULL) {
            printTree(tree, tree->root, 0);
            printf(";\n");
        }
    } else {
        if (node->nchildren == 0) {
            for (int i=0; i<depth; i++) printf("  ");
            printf("%d=%s:%f", node->name, node->longname.c_str(), node->dist);
        } else {
            // indent
            for (int i=0; i<depth; i++) printf("  ");
            printf("%d=(\n", node->name);

            for (int i=0; i<node->nchildren - 1; i++) {
                printTree(tree, node->children[i], depth+1);
                printf(",\n");
            }

            printTree(tree, node->children[node->nchildren-1], depth+1);
            printf("\n");

            for (int i=0; i<depth; i++) printf("  ");
            printf(")");

            if (depth > 0)
                printf(":%f", node->dist);
        }
    }
}



//some statistics
double Tree::total_branchlength() {
   ExtendArray<Node*> postnodes;
   getTreePostOrder(this, &postnodes);
   double len=0.0;

   for (int i=0; i < postnodes.size(); i++) {
       Node *node = postnodes[i];
       if (node != root) len += node->dist;
   }
   return len;
}



//Note: assumes all leaf nodes have same distance to root!
double Tree::tmrca() {
    double len=0.0;
    Node *node = root, *child;
    while (node->nchildren > 0) {
        child = node->children[0];
        len += child->dist;
        node = child;
    }
    return len;
}


double Tree::popsize() {
    ExtendArray<Node*> postnodes;
    int numleaf = (nnodes+1)/2;
    vector<double>tmpAge(nnodes);
    vector<double>ages;
    getTreePostOrder(this, &postnodes);
    for (int i=0; i < postnodes.size(); i++) {
        if (postnodes[i]->nchildren == 0) {
            tmpAge[postnodes[i]->name] = 0.0;
        } else {
            tmpAge[postnodes[i]->name] = tmpAge[postnodes[i]->children[0]->name] +
                postnodes[i]->children[0]->dist;
            ages.push_back(tmpAge[postnodes[i]->name]);
        }
    }
    std::sort(ages.begin(), ages.end());
    double lasttime=0;
    int k = numleaf;
    double popsize=0;
    for (unsigned int i=0; i < ages.size(); i++) {
        popsize += (double)k*(k-1)*(ages[i]-lasttime);
        lasttime = ages[i];
        k--;
    }
    return popsize/(4.0*numleaf-4);
}


double tmrca_half_rec(Node *node, int numnode, vector<int> numnodes,
                      double curr_tmrca) {
    if (node->nchildren != 2) {
        fprintf(stderr, "Error: tmrca_half only works for bifurcating trees\n");
    }
    if (numnodes[node->name] == numnode) return curr_tmrca;
    if (numnodes[node->children[0]->name] == numnode &&
        numnodes[node->children[1]->name] == numnode) {
        if (node->children[0]->dist >
            node->children[1]->dist)
            return curr_tmrca - node->children[0]->dist;
        else return curr_tmrca - node->children[1]->dist;
    }
    if (numnodes[node->children[0]->name] >= numnode) {
        assert(numnodes[node->children[1]->name] < numnode);
        return tmrca_half_rec(node->children[0], numnode, numnodes,
                              curr_tmrca - node->children[0]->dist);
    } else if (numnodes[node->children[1]->name] >= numnode) {
        assert(numnodes[node->children[0]->name] < numnode);
        return tmrca_half_rec(node->children[1], numnode, numnodes,
                              curr_tmrca - node->children[1]->dist);
    }
    return curr_tmrca;
}


double Tree::tmrca_half() {
    vector<int> numnodes(nnodes);
    ExtendArray<Node*> postnodes;
    getTreePostOrder(this, &postnodes);
    for (int i=0; i < postnodes.size(); i++) {
        numnodes[postnodes[i]->name] = 1;
        for (int j=0; j < postnodes[i]->nchildren; j++) {
            numnodes[postnodes[i]->name] += numnodes[postnodes[i]->children[j]->name];
        }
    }
    assert(nnodes == numnodes[root->name]);
    return tmrca_half_rec(root, (nnodes+1)/2-1, numnodes, this->tmrca());
}


double Tree::rth() {
    return this->tmrca_half()/this->tmrca();
}


//=============================================================================
// primitive tree format conversion functions

extern "C" {

    /*
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
    */


// create a tree object from a parent tree array
void ptree2tree(int nnodes, int *ptree, Tree *tree)
{
    Node **nodes = tree->nodes;

    // allocate children
    for (int i=0; i<nnodes; i++) {
        nodes[i]->allocChildren(2);
        nodes[i]->name = i;
        nodes[i]->nchildren = 0;
    }

    // store parent and child pointers
    for (int i=0; i<nnodes; i++) {
        int parent = ptree[i];

        if (parent != -1) {
            Node *parentnode = nodes[parent];
            parentnode->children[parentnode->nchildren++] = nodes[i];
            nodes[i]->parent = parentnode;
        } else {
            nodes[i]->parent = NULL;
        }
    }

    // set root
    tree->root = nodes[nnodes - 1];
    assert(tree->assertTree());
}


// create a parent tree from a tree object array
void tree2ptree(Tree *tree, int *ptree)
{
    Node **nodes = tree->nodes;
    int nnodes = tree->nnodes;

    for (int i=0; i<nnodes; i++) {
        if (nodes[i]->parent)
            ptree[i] = nodes[i]->parent->name;
        else
            ptree[i] = -1;
    }
}


Tree *makeTree(int nnodes, int *ptree)
{
    Tree *tree = new Tree(nnodes);
    ptree2tree(nnodes, ptree, tree);
    return tree;
}


void deleteTree(Tree *tree)
{
    delete tree;
}

void setTreeDists(Tree *tree, float *dists)
{
    tree->setDists(dists);
}




} // extern C



} // namespace spidir
