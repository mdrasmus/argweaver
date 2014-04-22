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
Tree::Tree(string newick, const vector<double>& times)
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
    coal_node = recomb_node = NULL;
    coal_time = recomb_time = -1;

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
            if (sscanf(&newick[i+1], "%lf", &node->dist) != 1) {
                printError("bad newick: error reading distance");
                printf("&newick[%i+1]=%s\n", i, &newick[i+1]);
                printf("node->dist=%lf\n", node->dist);
                exit(1);
            }
            i=j-1;
            break;
        }
        case '[': { // comment next; parse recomb_node or coal_node
            int count=1;
            int j=i+1;
	    double t;
	    string tmpstr;
            while (count != 0) {
                if (j==len) {
                    printError("bad newick: no closing bracket in NHX comment");
                    break;
                }
                if (newick[j]==']') count--;
                else if (newick[j]=='[') count++;
                j++;
            }
	    tmpstr = newick.substr(i+1, j-1-2);
	    if (sscanf(tmpstr.c_str(), "&&NHX:recomb_time=%lf]", &t)==1) {
	      this->recomb_time = t;
	      this->recomb_node = node;
	    } else if (sscanf(tmpstr.c_str(), "&&NHX:coal_time=%lf]", &t)==1) {
	      this->coal_time = t;
	      this->coal_node = node;
	    }
	    //            node->nhx = newick.substr(i+1, j-i-2);
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
    this->setPostNodes();
    for (int i=0; i < postnodes.size(); i++) {
	if (postnodes[i]->nchildren == 0) 
	    postnodes[i]->age = 0.0;
	else postnodes[i]->age = postnodes[i]->children[0]->age + 
	       postnodes[i]->children[0]->dist;
       // printf("node %i age %lf %s\n", postnodes[i]->name, postnodes[i]->age, postnodes[i]==root ? "ROOT" : "");
    }
    //printf("done assigning ages\n");
    if (times.size() > 0) {
      this->correct_times(times, 1);
      if (recomb_node != NULL )
	this->correct_recomb_times(times);
    }

    for (int i=0; i < nnodes; i++) {
      if (nodes[i]->longname.length() > 0)
	nodename_map[nodes[i]->longname] = i;
    }
}


//void Tree::correct_times(map<string,double> times) {
void Tree::correct_times(vector<double> times, double tol) {
    unsigned int lasttime=0, j;
    this->setPostNodes();
    for (int i=0; i < postnodes.size(); i++) {
	if (postnodes[i]->nchildren == 0) {
	    postnodes[i]->age = 0.0;
	    lasttime = 0;
	    for (j=0; j < times.size(); j++)
		if (fabs(times[j]-postnodes[i]->dist) < tol) break;
	    if (j == times.size()) {
		printError("No node has time %lf (leaf)", postnodes[i]->dist);
	    }
	    postnodes[i]->dist = times[j];
	    postnodes[i]->parent->age = times[j];
	}
	else {
	  /*	    postnodes[i]->age = postnodes[i]->children[0]->age + 
		    postnodes[i]->children[0]->dist;*/
	    double newage = postnodes[i]->age + postnodes[i]->dist;
	    for (j=lasttime; j < times.size(); j++) 
	      if (fabs(times[j]-newage) < tol) break;
	    if (j == times.size())
	      printError("No node has time %lf", newage);
	    postnodes[i]->dist = age_diff((double)times[j], postnodes[i]->age);
	    if (postnodes[i]->parent != NULL)
	      postnodes[i]->parent->age = times[j];
	    lasttime = j;
	}
    }
}


string Tree::print_newick_to_string_recur(Node *n, bool internal_names,
					  char *branch_format_str,
					  bool show_nhx, bool oneline) {
  string rv;
  char tmp[1000];
  if (n->nchildren > 0) {
    int first=0, second=1;
    if (n->children[0]->longname.size() > 0 &&
	n->children[1]->longname.size() > 0 &&
	n->children[0]->longname.compare(n->children[1]->longname) > 0) {
      first=1;
      second=0;
    }
    rv.append("(");
    rv.append(print_newick_to_string_recur(n->children[first], 
					   internal_names,
					   branch_format_str, 
					   show_nhx, oneline));
    for (int i=1; i < n->nchildren; i++) {
      rv.append(",");
      rv.append(print_newick_to_string_recur(n->children[second], 
					     internal_names,
					     branch_format_str, 
					     show_nhx, oneline));
    }
    rv.append(")");
    if (internal_names) rv.append(n->longname);
  } else {
    rv.append(n->longname);
  }
  //    fprintf(f, "(%i)", n->name);
  if (branch_format_str != NULL && n->parent != NULL) {
    rv.append(":");
    sprintf(tmp, branch_format_str, n->dist);
    rv.append(tmp);
  }
  if (show_nhx) {
    if (n == this->recomb_node) {
      sprintf(tmp,"[&&NHX:recomb_time=%.1f]", this->recomb_time);
      rv.append(tmp);
    } if (n == this->coal_node) {
      sprintf(tmp, "[&&NHX:coal_time=%.1f]", this->coal_time);
      rv.append(tmp);
    }
  }
  if (!oneline && n->nchildren > 0) rv.append("\n");
  return rv;
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
        nodes2[i]->age = nodes[i]->age;
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
    tree2->recomb_time = recomb_time;
    tree2->coal_time = coal_time;
    tree2->recomb_node = recomb_node == NULL ? NULL : nodes2[recomb_node->name];
    tree2->coal_node = coal_node == NULL ? NULL : nodes2[coal_node->name];
    tree2->nodename_map = nodename_map;
    tree2->postnodes.clear();
    return tree2;
}



//private function called by update_spr below
//given a nhx tag, looks backwards until a node name is found, 
//simultaneously counting close-parenthesis till it finds it.
// then the node is that many up from the leaf given by the node name.
int Tree::get_node_from_newick(char *newick, char *nhx) {
  int num_paren=0;
 /* printf("get_node_from_newick\n");
  printf("newick=%s\n", newick);
  printf("nhx=%s\n", nhx);*/
  while (1) {
    while (':' != nhx[0] && ')' != nhx[0]) {
      assert(nhx != newick);
      if (nhx[0]==']') {
	while (nhx[0]!='[') nhx--;
      }
      nhx--;
    }
    if (nhx[0]==':') nhx--;
    if (nhx[0]==')') {
      num_paren++;
      nhx--;
    } else {
      char *tmp=&nhx[1];
      assert(nhx[1]==':');
      nhx[1]='\0';
      while (!isNewickChar(nhx[0])) {
	nhx--;
      }
      nhx++;
      map<string,int>::iterator it = nodename_map.find(string(nhx));
      if (it ==nodename_map.end()) { //error
        printf("nodename_map size=%i\n", (int)nodename_map.size());
        printf("nhx=%s\n", nhx);
        assert(it != nodename_map.end());
      }
      int n = it->second;
      assert(nodes[n]->nchildren==0);  // should be leaf
      tmp[0]=':';
      for (int i=0; i < num_paren; i++) {
	if (nodes[n]->parent == NULL) {
	  printf("assert failed newick=%s nhx=%s node=%i num_paren=%i i=%i\n", newick, nhx, nodes[n]->name, num_paren, i);
	}
	assert(nodes[n]->parent != NULL);
	n = nodes[n]->parent->name;
      }
      //printf("returning %i\n", n);
      return n;
    }
  }
}

void Tree::correct_recomb_times(const vector<double>& times) {
    unsigned int i;
    for (i=0; i < times.size(); i++)
	if (fabs(recomb_time - times[i]) < 1) {
	    recomb_time = times[i];
	    break;
	}
    assert(i != times.size());
    for (; i < times.size(); i++) {
	if (fabs(coal_time - times[i]) < 1) {
	    coal_time = times[i];
	    return;
	}
    }
    assert(0);
}

void Tree::update_spr(char *newick, const vector<double>& times) {
  char search1[100]="[&&NHX:recomb_time=";
  char search2[100]="[&&NHX:coal_time=";
  char *x = strstr(newick, search1);
  /*fprintf(stderr, "update_spr\n");
  fprintf(stderr, "%s\n", newick);*/
  if (x == NULL) {
    recomb_node = NULL;
    coal_node = NULL;
    //    printf("no recomb node\n");
    return;
  }
  assert(1==sscanf(x, "[&&NHX:recomb_time=%lg", &recomb_time));
  recomb_node = nodes[this->get_node_from_newick(newick, x)];

  x = strstr(newick, search2);
  assert(x != NULL);
  assert(1 == sscanf(x, "[&&NHX:coal_time=%lg", &coal_time));
  coal_node = nodes[this->get_node_from_newick(newick, x)];

  if (times.size() > 0) this->correct_recomb_times(times);
  //printf("done update_spr recomb_node=%i coal_node=%i\n", recomb_node->name, coal_node->name);
}

//update the SPR on pruned tree based on node_map in big tree
void Tree::update_spr_pruned(Tree *orig_tree) {
  if (orig_tree->recomb_node == NULL) {
    recomb_node = coal_node = NULL;
    return;
  }
  int num=orig_tree->node_map.nm[orig_tree->recomb_node->name];
  if (num == -1 || nodes[num] == root) {
    recomb_node = coal_node = NULL; 
  } else {
    assert(num>=0);
    recomb_node = nodes[num];
    recomb_time = orig_tree->recomb_time;
    num = orig_tree->node_map.nm[orig_tree->coal_node->name];
    if (num == -1) {
      // coal node does not map; need to trace back until it does 
      Node *n = orig_tree->coal_node;
      while (orig_tree->node_map.nm[n->name] == -1) {
	//should never be root here; root should always map to pruned tree
	assert(n->parent != NULL);
	n = n->parent;
      }
      assert(orig_tree->coal_time-1 <= n->age);
      coal_time = n->age;
      coal_node = nodes[orig_tree->node_map.nm[n->name]];
    } else {
      assert(num >= 0);
      coal_node = nodes[num];
      coal_time = orig_tree->coal_time;
    }
    if (recomb_node == coal_node) {
      recomb_node = coal_node = NULL;
    }
  } 
  if (recomb_node != NULL) assert(coal_node != NULL);
}         


void Tree::remap_node(Node *n, int id, int *deleted_branch) {
  int old_id = node_map.nm[n->name];
  //  printf("remap_node %i %i->%i\n", n->name, old_id, id);
  if (old_id == id) return; // no change
  node_map.inv_nm[old_id].erase(n->name);
  if (old_id >= 0 && node_map.inv_nm[old_id].size() == 0) {
    assert((*deleted_branch)==-1 || (*deleted_branch)==old_id);
    *deleted_branch = old_id;
    //    printf("deleted_branch=%i\n", *deleted_branch); fflush(stdout);
  }
  node_map.nm[n->name] = id;
  node_map.inv_nm[id].insert(n->name);
}


void Tree::propogate_map(Node *n, int *deleted_branch, int count, 
			 int count_since_change, int maxcount, 
			 int maxcount_since_change) {
  Node *c0, *c1;

  if (count==maxcount) return;
  if (count_since_change == maxcount_since_change) return;
  if (n->nchildren == 0) 
    return propogate_map(n->parent, deleted_branch, count+1, 
			 count_since_change+1, maxcount, maxcount_since_change);
  c0 = n->children[0];
  c1 = n->children[1];
  //  printf("propogate_map %i (%i) %i (%i) %i (%i) count=%i\n", n->name, node_map.nm[n->name], c0->name, node_map.nm[c0->name], c1->name, node_map.nm[c1->name], count); fflush(stdout);
  int change=0;
  if (node_map.nm[c0->name] == -1 && node_map.nm[c1->name] == -1) {
    if (node_map.nm[n->name] != -1) {
      remap_node(n, -1, deleted_branch);
      change=1;
    }
  } else if (node_map.nm[c0->name] == -1 || 
	     node_map.nm[c1->name] == -1) {
    Node *c = node_map.nm[c0->name] == -1 ? c1 : c0;
    if (node_map.nm[n->name] != node_map.nm[c->name]) {
      remap_node(n, node_map.nm[c->name], deleted_branch);
      change=1;
    }
  } else {  // neither are -1
    if (node_map.nm[c0->name] == node_map.nm[c1->name]) {
      printf("here n=%i c0=%i c1=%i %i %i %i %i\n",
	     n->name, c0->name, c1->name, node_map.nm[n->name], node_map.nm[c0->name], node_map.nm[c1->name], *deleted_branch); fflush(stdout); 
      assert(0); 
    }
    if (node_map.nm[n->name] == -1 || node_map.nm[n->name]==-3 || 
	node_map.nm[n->name] == node_map.nm[c0->name] ||
	node_map.nm[n->name] == node_map.nm[c1->name]) {
      change=1;
      remap_node(n, -2, deleted_branch);
    }
  }
  if (n == root) return;
  return propogate_map(n->parent, deleted_branch, count+1, change==0 ? count+1 : 0, maxcount, maxcount_since_change);
}

double Tree::age_diff(double age1, double age2) {
    double diff = age1 - age2;
    if (diff < 0) {
	if (diff < -2) {
	    fprintf(stderr, "got age diff=%.8f (age1=%.8f, age2=%.8f)\n", diff, age1, age2);
	    fflush(stderr);
	    assert(0);
	}
	return 0.0;
    }
    return diff;
}


//by the prune function. It will be modified to map the branches from
//the full to pruned tree after the SPR operation.
// if prune_root is not NULL, it should correspond to node id in big tree which
// forms root of pruned tree, and will be updated if necessary
void Tree::apply_spr() {
  Node *recomb_parent, *recomb_grandparent, *recomb_sibling, *coal_parent=NULL;
  int x;
  /*fprintf(stderr, "apply_spr recomb_node=%i (%i) coal_node=%i (%i) root=%i\n", recomb_node->name, recomb_node->nchildren, coal_node->name, coal_node->nchildren, root->name);
  this->print_newick(stderr);
  fprintf(stderr, "\n"); fflush(stderr);*/
  if (recomb_node == NULL) return;
  if (recomb_node == root) {
    if (coal_node != root) {
      assert(0);
    }
    return;
  }
  if (recomb_node == coal_node) return;

  recomb_parent = recomb_node->parent;
  assert(recomb_parent != NULL);  //ie, recomb_node should not be root
  assert(recomb_parent->nchildren == 2);
  x = (recomb_parent->children[0] == recomb_node ? 0 : 1);
  recomb_sibling = recomb_parent->children[!x];

  //recomb_grandparent might be NULL
  recomb_grandparent = recomb_parent->parent;
  
  //coal_parent might be NULL too
  coal_parent = coal_node->parent;

  //special case; topology doesn't change; just adjust branch lengths/ages
  //(this violates SMC so shouldn't be true for an ARGweaver tree 
  //   but may be for a subtree)
  if (coal_parent == recomb_parent) {
    coal_parent->age = coal_time;
    coal_node->dist = age_diff(coal_time, coal_node->age);
    recomb_node->dist = age_diff(coal_time, recomb_node->age);
    if (recomb_grandparent != NULL)
	recomb_parent->dist = age_diff(recomb_grandparent->age, coal_time);
    //fprintf(stderr, "done trivial update SPR\n");
    return;
  }
  // similar other special case
  if (coal_node == recomb_parent) {
    coal_node->age = coal_time;
    recomb_node->dist = age_diff(coal_time, recomb_node->age);
    recomb_sibling->dist = age_diff(coal_time, recomb_sibling->age);
    if (coal_parent != NULL) 
	coal_node->dist = age_diff(coal_parent->age, coal_time);
    //fprintf(stderr, "done trivial update SPR2\n");
    return;
  }
  postnodes.clear();

  //now apply SPR
  recomb_sibling->parent = recomb_grandparent;
  if (recomb_grandparent != NULL) {
    int x1 = (recomb_grandparent->children[0]==recomb_parent ? 0 : 1);
    recomb_grandparent->children[x1] = recomb_sibling;
    recomb_sibling->dist += recomb_parent->dist;
  } else {
    root = recomb_sibling;
    recomb_sibling->parent = NULL;
  }

  //recomb_parent is extracted; re-use as new_node. one child is still recomb_node
  recomb_parent->children[!x] = coal_node;
  coal_node->dist = age_diff(coal_time, coal_node->age);
  recomb_node->dist = age_diff(coal_time, recomb_node->age);
  coal_node->parent = recomb_parent;
  recomb_parent->age = coal_time;
  if (coal_parent != NULL) {
    recomb_parent->parent = coal_parent;
    recomb_parent->dist = age_diff(coal_parent->age, coal_time);
    coal_parent->children[coal_parent->children[0]==coal_node ? 0 : 1] = recomb_parent;
  } else {
    root = recomb_parent;
    recomb_parent->parent = NULL;
  }
  //  this->print_newick(stdout, 1, NULL, 0); printf("\n"); fflush(stdout);
  if (node_map.size() > 0) {
    int deleted_branch=-1;
    //set recomb_node and recomb_parent maps to -3 = unknown
    remap_node(recomb_parent, -3, &deleted_branch);
    propogate_map(coal_node, &deleted_branch, 0, 0, 1, 1);
    propogate_map(recomb_node, &deleted_branch, 0, 0, 1, 1);
    propogate_map(recomb_sibling, &deleted_branch, 0, 0, -1, 4);
    propogate_map(recomb_parent, &deleted_branch, 0, 0, -1, 4);
    //    if (recomb_grandparent != NULL)
    //      propogate_map(recomb_grandparent, &deleted_branch, 0, 0, -1, 4);
    
    map<int,set<int> >::iterator it = node_map.inv_nm.find(-2);
    if (it != node_map.inv_nm.end() && it->second.size() > 0) {
      if (deleted_branch == -1) {
	assert(deleted_branch != -1);
      }
      set<int> rename_nodes = it->second;
      for (set<int>::iterator it2 = rename_nodes.begin(); it2 != rename_nodes.end(); ++it2) {
	remap_node(nodes[*it2], deleted_branch, &deleted_branch);
      }
    }
    for (int i=0; i < nnodes; i++) {
      if (node_map.nm[i] == -2) {
	assert(0);
      }
    }
    //    node_map.print();
  }
}


NodeMap Tree::prune(set<string> leafs, bool allBut) {
    ExtendArray<Node*> newnodes = ExtendArray<Node*>(0);
    map<int,int> node_map;  //maps original nodes to new nodes
    this->setPostNodes();
    vector<bool> is_leaf(postnodes.size());
    node_map.clear();
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
		node_map[postnodes[i]->name] = -1;
		if (recomb_node == postnodes[i]) {
		  recomb_node = NULL;
		  recomb_time = -1;
		  coal_node = NULL;
		  coal_time = -1;
		} 
		if (coal_node == postnodes[i]) {
		  coal_node = parent;
		  coal_time = parent->age;
		}
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
	      node_map[postnodes[i]->name] = newnodes.size();
	      newnodes.append(postnodes[i]);
            }
        } else if (postnodes[i]->nchildren == 1) {
	  if (postnodes[i] == root) {
	    node_map[postnodes[i]->name] = node_map[postnodes[i]->children[0]->name];
	    root = postnodes[i]->children[0];
	    root->parent = NULL;
	    if (recomb_node == postnodes[i]) {
	      recomb_node = NULL;
	      recomb_time = -1;
	      coal_node = NULL;
	      coal_time = -1;
	    }
	    if (coal_node == postnodes[i]) {
	      coal_node = root;
	    }
	    delete postnodes[i];
	  } else {
	    Node *parent = postnodes[i]->parent;
	    int j, maxj=parent->nchildren;
	    for (j=0; j < maxj; j++) {
	      if (parent->children[j] == postnodes[i]) {
		parent->children[j] = postnodes[i]->children[0];
		postnodes[i]->children[0]->dist += postnodes[i]->dist;
		postnodes[i]->children[0]->parent = parent;
		if (postnodes[i] == recomb_node)
		  recomb_node = postnodes[i]->children[0];
		if (postnodes[i] == coal_node)
		  coal_node = postnodes[i]->children[0];
		node_map[postnodes[i]->name] = node_map[postnodes[i]->children[0]->name];
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
	  node_map[postnodes[i]->name] = newnodes.size();
	  newnodes.append(postnodes[i]);
        }
    }
    if (recomb_node == root || recomb_node == coal_node) {
      recomb_node = NULL;
      coal_node = NULL;
      recomb_time = coal_time = -1;
    }

    nodes.clear();
    postnodes.clear();
    for (int i=0; i < newnodes.size(); i++) {
        nodes.append(newnodes[i]);
        nodes[i]->name = i;
	// temporary check
	int error=0;
	if (newnodes[i] == recomb_node) {
	  if (recomb_time+2 < nodes[i]->age) {
	    fprintf(stderr, "bad recomb age (too small)\n"); error=1;
	  } else if (recomb_time-2 > nodes[i]->age + nodes[i]->dist) {
	    fprintf(stderr, "bad recomb age (too big)\n"); error=1;
	  }
	}
	if (newnodes[i] == coal_node) {
	  if (coal_time+2 < nodes[i]->age) {
	    fprintf(stderr, "bad coal age (too small) coal_time=%lf age=%lf\n", coal_time, nodes[i]->age); error=1;
	  } else if (coal_node != root && coal_time-2 > nodes[i]->age + nodes[i]->dist) {
	    fprintf(stderr, "bad coal age (too big) coal_time=%lf age=%lf dist=%lf\n", coal_time, nodes[i]->age, nodes[i]->dist); error=1;
	  }
	}
	if (error) {
	  this->print_newick(stderr);
	  exit(1);
	}
    }
    nnodes = nodes.size();
    nodename_map.clear();
    for (int i=0; i < nnodes; i++) {
      if (nodes[i]->longname.length() > 0)
	nodename_map[nodes[i]->longname] = i;
    }
    return NodeMap(node_map);
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
    postnodes.clear();
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

    postnodes.clear();

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
    double nextDist = 0;
    double rootdist;

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
        double tmpdist = ptr1->dist;
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

    postnodes.clear();
    reroot(newroot);
}





// store a hash key representing the topology into the key array
// key is a parent tree representation where the internal nodes are
// given a consistent numbering
void Tree::hashkey(int *key)
{
    // get post order of nodes
    this->setPostNodes();

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
    postnodes.clear();
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


void Tree::setPostNodesRec(Node *n) {
  for (int i=0; i < n->nchildren; i++)
    setPostNodesRec(n->children[i]);
  postnodes.append(n);
}

void Tree::setPostNodes() {
  if (postnodes.size() == nnodes) return;
  postnodes.clear();
  this->setPostNodesRec(root);
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
            printf("%d=%s:%lf", node->name, node->longname.c_str(), node->dist);
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
                printf(":%lf", node->dist);
        }
    }
}



//some statistics
double Tree::total_branchlength() {
  double len=0.0;
  this->setPostNodes();
  for (int i=0; i < postnodes.size(); i++) {
    Node *node = postnodes[i];
    if (node != root) len += node->dist;
  }
  return len;
}



//Note: assumes all leaf nodes have same distance to root!
double Tree::tmrca() {
    return root->age;
}


double Tree::popsize() {
    int numleaf = (nnodes+1)/2;
    vector<double>ages;
    double lasttime=0, popsize=0;
    int k=numleaf;
    for (int i=0; i < nnodes; i++)
      if (nodes[i]->nchildren > 0)
	ages.push_back(nodes[i]->age);
    std::sort(ages.begin(), ages.end());
    for (unsigned int i=0; i < ages.size(); i++) {
        popsize += (double)k*(k-1)*(ages[i]-lasttime);
        lasttime = ages[i];
        k--;
    }
    return popsize/(4.0*numleaf-4);
}

//assume that times is sorted!
vector<double> Tree::coalCounts(vector<double> times) {
  vector<double> counts(times.size(), 0.0);
  vector<double> ages;
  unsigned int total=0;
  for (int i=0; i < nnodes; i++) {
    if (nodes[i]->nchildren > 0)
      ages.push_back(nodes[i]->age);
  }
  std::sort(ages.begin(), ages.end());
  unsigned int idx=0;
  for (unsigned int i=0; i < ages.size(); i++) {
    while (1) {
      if (fabs(ages[i]-times[idx]) < 0.00001) {
	counts[idx]++;
	total++;
	break;
      }
      idx++;
      assert(idx < times.size());
    }
  }
  assert(total == ages.size());
  return counts;
}
    

double Tree::num_zero_branches() {
    int count=0;
    for (int i=0; i < nnodes; i++) {
	if (nodes[i] != root && fabs(nodes[i]->dist) < 0.0001)
	    count++;
    }
    return count;
}
	


double tmrca_half_rec(Node *node, int numnode, vector<int> numnodes) {
    if (node->nchildren != 2) {
        fprintf(stderr, "Error: tmrca_half only works for bifurcating trees\n");
    }
    if (numnodes[node->name] == numnode) return node->age;
    if (numnodes[node->children[0]->name] == numnode &&
        numnodes[node->children[1]->name] == numnode) {
      return min(node->children[0]->age, node->children[1]->age);
    }
    if (numnodes[node->children[0]->name] >= numnode) {
        assert(numnodes[node->children[1]->name] < numnode);
        return tmrca_half_rec(node->children[0], numnode, numnodes);
    } else if (numnodes[node->children[1]->name] >= numnode) {
        assert(numnodes[node->children[0]->name] < numnode);
        return tmrca_half_rec(node->children[1], numnode, numnodes);
    }
    return node->age;
}


double Tree::tmrca_half() {
    vector<int> numnodes(nnodes);
    this->setPostNodes();
    for (int i=0; i < postnodes.size(); i++) {
        numnodes[postnodes[i]->name] = 1;
        for (int j=0; j < postnodes[i]->nchildren; j++) {
            numnodes[postnodes[i]->name] += numnodes[postnodes[i]->children[j]->name];
        }
    }
    assert(nnodes == numnodes[root->name]);
    return tmrca_half_rec(root, (nnodes+1)/2-1, numnodes);
}


double Tree::rth() {
    return this->tmrca_half()/this->tmrca();
}

double Tree::distBetweenLeaves(Node *n1, Node *n2) {
    if (n1 == n2) return 0.0;
    this->setPostNodes();
    vector<int> count(postnodes.size());
    int s=0;
    double rv=0.0;
    for (int i=0; i < postnodes.size(); i++) {
	if (postnodes[i] == n1 || postnodes[i] == n2) {
	    count[postnodes[i]->name] = 1;
	    s++;
	}
	if (postnodes[i]->nchildren == 2) {
	    count[postnodes[i]->name] = count[postnodes[i]->children[0]->name] + 
		count[postnodes[i]->children[1]->name];
	    if (count[postnodes[i]->name] == 2) break;
	}
	if (count[postnodes[i]->name]) rv += postnodes[i]->dist;
    }
    assert(s == 2);
    return rv;
}


// want to return set of nodes above which mutations happened under infinite
// sites to cause site pattern.
// assumes tree has been pruned so that only leafs with information remain (no Ns)
// this is not particularly efficient! may want to think of something faster some time.
set<Node*> Tree::lca(set<Node*> derived) {
  set<Node*>::iterator it;
  set<Node*> rv;
  
  if (derived.size() == 1) return derived;

  this->setPostNodes();
  for (int i=0; i < postnodes.size(); i++) {
    if (postnodes[i]->nchildren == 0) continue;
    if (postnodes[i] == root) {
      assert(derived.size() == 1);
      rv.insert(*(derived.begin()));
      return rv;
    }
    int count=0;
    for (int j=0; j < postnodes[i]->nchildren; j++) {
      if (derived.find(postnodes[i]->children[j]) != derived.end())
	count++;
    }
    if (count == postnodes[i]->nchildren) { // all children are derived
      for (int j=0; j < postnodes[i]->nchildren; j++)
	derived.erase(postnodes[i]->children[j]);
      derived.insert(postnodes[i]);
    }  else if (count != 0) {
      for (int j=0; j < postnodes[i]->nchildren; j++) {
	if (derived.find(postnodes[i]->children[j]) != derived.end()) {
	  rv.insert(postnodes[i]->children[j]);
	  derived.erase(postnodes[i]->children[j]);
	}
      }
    }
    if (derived.size() == 0) return rv;
  }
  fprintf(stderr, "got to end of LCA\n"); fflush(stderr);
  return rv;
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

void setTreeDists(Tree *tree, double *dists)
{
    tree->setDists(dists);
}




} // extern C



} // namespace spidir
