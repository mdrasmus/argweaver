/*=============================================================================

  Matt Rasmussen

  Newick tree reading/writing

=============================================================================*/

#ifndef ARGHMM_NEWICK_H
#define ARGHMM_NEWICK_H

#include <stdlib.h>
#include <string>

#include "Tree.h"

namespace argweaver {

using namespace spidir;

Tree *readNewickTree(FILE *infile, Tree *tree=NULL);
Tree *readNewickTree(const char *filename, Tree *tree=NULL);

void writeNewickNode(FILE *out, Node *node, int depth, bool oneline=false);
void writeNewickTree(FILE *out, Tree *tree, int depth, bool oneline=false);
bool writeNewickTree(const char *filename, Tree *tree, bool oneline=false);


} // namespace argweaver


#endif // ARGHMM_NEWICK_H

