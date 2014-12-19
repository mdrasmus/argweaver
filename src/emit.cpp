
#include "common.h"
#include "emit.h"
#include "seq.h"
#include "thread.h"

namespace argweaver {


//=============================================================================
// invariant sites
// TODO: make this mask aware


// Returns true if position 'pos' is invariant
static inline bool is_invariant_site(const char *const *seqs,
                                     const int nseqs, const int pos)
{
    const char c = seqs[0][pos];
    for (int j=1; j<nseqs; j++) {
        if (seqs[j][pos] != c) {
            return false;
        }
    }
    return true;
}


// Populates array 'invariant' with status of invariant sites
void find_invariant_sites(const char *const *seqs, int nseqs, int seqlen,
                          bool *invariant)
{
    // find invariant sites
    for (int i=0; i<seqlen; i++) {
        char c = seqs[0][i];
        bool mut = false;
        for (int j=1; j<nseqs; j++) {
            if (seqs[j][i] != c) {
                mut = true;
                break;
            }
        }
        invariant[i] = !mut;
    }
}


void find_masked_sites(const char *const *seqs, int nseqs, int seqlen,
                       bool *masked, bool *invariant)
{
    if (invariant) {
        for (int i=0; i<seqlen; i++)
            masked[i] = (seqs[0][i] == 'N' && invariant[i]);
    } else {
        for (int i=0; i<seqlen; i++)
            masked[i] = (seqs[0][i] == 'N' && is_invariant_site(seqs, nseqs, i));
    }
}


int count_alleles(const char *const *seqs,
                  const int nseqs, const int pos)
{
    int counts[4] = {0, 0, 0, 0};

    for (int j=0; j<nseqs; j++) {
        int c = dna2int[(int) seqs[j][pos]];
        if (c != -1)
            counts[c]++;
    }
    int alleles = int(counts[0] > 0) +
        int(counts[1] > 0) +
        int(counts[2] > 0) +
        int(counts[3] > 0);

    return alleles;
}


//=============================================================================
// calculate mutation probabilities


// Juke-Cantor
// probability of mutation over time 't' with mutation rate 'mu'
static inline double prob_branch(double t, double mu, bool mut)
{
    const double f = 4. / 3.;
    if (!mut)
        return .25 * (1.0 + 3. * exp(-f*mu*t));
    else
        return .25 * (1.0 - exp(-f*mu*t));
}


// get mutation probabilities
void prob_tree_mutation(const LocalTree *tree, const ArgModel *model,
                        double *muts, double *nomuts)
{
    const double *times = model->times;
    const int nnodes = tree->nnodes;
    const LocalNode *nodes = tree->nodes;
    const double mintime = model->get_mintime();

    for (int i=0; i<nnodes; i++) {
        if (i == tree->root)
            continue;
        int parent_age = nodes[nodes[i].parent].age;
        if (parent_age == model->get_removed_root_time())
            continue;

        double t = max(times[parent_age] - times[nodes[i].age], mintime);
        muts[i] = prob_branch(t, model->mu, true);
        nomuts[i] = prob_branch(t, model->mu, false);
    }
}


//============================================================================
// emissions

// table of partial likelihood values
typedef double lk_row[4];

class LikelihoodTable
{
public:
    LikelihoodTable(int seqlen, int nnodes) :
        seqlen(seqlen),
        nnodes(nnodes)
    {
        data = new lk_row* [seqlen];
        for (int i=0; i<seqlen; i++)
            data[i] = new double [nnodes][4];
    }

    ~LikelihoodTable()
    {
        for (int i=0; i<seqlen; i++)
            delete [] data[i];
        delete [] data;
    }

    int seqlen;
    int nnodes;
    lk_row **data;
};


// calculate inner partial likelihood for one node and site
inline void likelihood_site_node_inner(
    const LocalTree *tree, const int node,
    const char *const *seqs, const int pos,
    const double *muts, const double *nomuts, lk_row* inner)
{
    const LocalNode* nodes = tree->nodes;
    const int j = node;

    if (nodes[j].is_leaf()) {
        // leaf case
        const char c = seqs[j][pos];
        if (c == 'N') {
            inner[j][0] = 1.0;
            inner[j][1] = 1.0;
            inner[j][2] = 1.0;
            inner[j][3] = 1.0;
        } else {
            inner[j][0] = 0.0;
            inner[j][1] = 0.0;
            inner[j][2] = 0.0;
            inner[j][3] = 0.0;
            inner[j][dna2int[(int) c]] = 1.0;
        }
    } else {
        // internal node case
        int c1 = nodes[j].child[0];
        int c2 = nodes[j].child[1];

        for (int a=0; a<4; a++) {
            double p1 = 0.0;
            double p2 = 0.0;

            for (int b=0; b<4; b++) {
                if (a == b) {
                    p1 += inner[c1][b] * nomuts[c1];
                    p2 += inner[c2][b] * nomuts[c2];
                } else {
                    p1 += inner[c1][b] * muts[c1];
                    p2 += inner[c2][b] * muts[c2];
                }
            }

            inner[j][a] = p1 * p2;
        }
    }
}


// calculate outer partial likelihood for one node and site
inline void likelihood_site_node_outer(
    const LocalTree *tree, const int root, const int node,
    const char *const *seqs, const int pos,
    const double *muts, const double *nomuts,
    lk_row* outer, lk_row* inner)
{
    const LocalNode* nodes = tree->nodes;
    const int j = node;

    if (node == root) {
        // root case
        outer[j][0] = 1.0;
        outer[j][1] = 1.0;
        outer[j][2] = 1.0;
        outer[j][3] = 1.0;
    } else {
        // non-root case
        int sib = tree->get_sibling(j);
        int parent = nodes[j].parent;

        if (parent != root) {
            for (int a=0; a<4; a++) {
                double p1 = 0.0;
                double p2 = 0.0;
                for (int b=0; b<4; b++) {
                    if (a == b) {
                        p1 += inner[sib][b] * nomuts[sib];
                        p2 += outer[parent][b] * nomuts[parent];
                    } else {
                        p1 += inner[sib][b] * muts[sib];
                        p2 += outer[parent][b] * muts[parent];
                    }
                }
                outer[j][a] = p1 * p2;
            }
        } else {
            for (int a=0; a<4; a++) {
                double p1 = 0.0;
                for (int b=0; b<4; b++) {
                    if (a == b)
                        p1 += inner[sib][b] * nomuts[sib];
                    else
                        p1 += inner[sib][b] * muts[sib];
                }
                outer[j][a] = p1;
            }
        }
    }
}


// calculate entire inner partial likelihood table
double likelihood_site_inner(
    const LocalTree *tree, const char *const *seqs,
    const int pos, const int *order, const int norder,
    const double *muts, const double *nomuts, lk_row* inner)
{
    // iterate postorder through nodes
    for (int i=0; i<norder; i++)
        likelihood_site_node_inner(
            tree, order[i], seqs, pos, muts, nomuts, inner);

    // sum over root node
    double p = 0.0;
    int root = tree->root;
    for (int a=0; a<4; a++)
        p += inner[root][a]  * .25;

    return p;
}


// calculate eniter outer partial likelihood table
void likelihood_site_outer(
    const LocalTree *tree, const char *const *seqs, const int pos,
    const double *muts, const double *nomuts, bool internal,
    lk_row *inner, lk_row *outer)
{
    int queue[tree->nnodes];
    int top = 0;

    // process in preorder
    int maintree_root = internal ? tree->nodes[tree->root].child[1] :
        tree->root;
    queue[top++] = maintree_root;
    while (top > 0) {
        int node = queue[--top];
        likelihood_site_node_outer(tree, maintree_root, node, seqs, pos,
                                   muts, nomuts, outer, inner);

        // recurse
        if (!tree->nodes[node].is_leaf()) {
            queue[top++] = tree->nodes[node].child[0];
            queue[top++] = tree->nodes[node].child[1];
        }
    }
}


void calc_inner_outer(const LocalTree *tree, const ArgModel *model,
                      const char *const *seqs, const int seqlen,
                      const bool *invariant, bool internal,
                      lk_row **inner, lk_row **outer)
{
    // get postorder
    int norder = tree->nnodes;
    int order[tree->nnodes];
    tree->get_postorder(order);


    // get mutation probabilities and treelen
    double muts[tree->nnodes];
    double nomuts[tree->nnodes];
    prob_tree_mutation(tree, model, muts, nomuts);

    // calculate emissions for tree at each site
    for (int i=0; i<seqlen; i++) {
        if (!invariant[i]) {
            likelihood_site_inner(tree, seqs, i, order, norder,
                                  muts, nomuts, inner[i]);
            likelihood_site_outer(tree, seqs, i,
                                  muts, nomuts, internal, inner[i], outer[i]);
        }
    }
}



void likelihood_sites(const LocalTree *tree, const ArgModel *model,
                      const char *const *seqs,
                      const int seqlen, const int statei,
                      const bool *invariant,
                      double **emit, lk_row **table,
                      const int prev_node=-1, const int new_node=-1)
{
    const double *times = model->times;
    const LocalNode *nodes = tree->nodes;
    const double mintime = model->get_mintime();

    // get postorder
    int order[tree->nnodes];
    int norder;

    if (prev_node == -1) {
        tree->get_postorder(order);
        norder = tree->nnodes;
    } else {
        // find partial postorder

        // find dirty entries
        bool dirty[tree->nnodes];
        fill(dirty, dirty+tree->nnodes, false);
        for (int j=prev_node; j!=-1; j=nodes[j].parent)
            dirty[j] = true;

        // walk up root path
        norder = 0;
        for (int j=new_node; !dirty[j]; j=nodes[j].parent)
            order[norder++] = j;
        for (int j=prev_node; j!=-1; j=nodes[j].parent)
            order[norder++] = j;
    }


    // get mutation probabilities and treelen
    double muts[tree->nnodes];
    double nomuts[tree->nnodes];
    double treelen = 0.0;
    for (int i=0; i<tree->nnodes; i++) {
        if (i != tree->root) {
            double t = max(times[nodes[nodes[i].parent].age] -
                           times[nodes[i].age], mintime);
            muts[i] = prob_branch(t, model->mu, true);
            nomuts[i] = prob_branch(t, model->mu, false);
            treelen += t;
        }
    }


    // calculate invariant_lk
    double invariant_lk = .25 * exp(- model->mu * max(treelen, mintime));

    // calculate emissions for tree at each site
    for (int i=0; i<seqlen; i++) {
        if (invariant && invariant[i]) {
            // use precommuted invariant site likelihood
            if (seqs[0][i] == 'N')
                emit[i][statei] = 1.0; // masked case
            else
                emit[i][statei] = invariant_lk;
        } else {
            emit[i][statei] = likelihood_site_inner(
                tree, seqs, i, order, norder, muts, nomuts, table[i]);
        }
    }
}



double likelihood_tree(const LocalTree *tree, const ArgModel *model,
                       const char *const *seqs, const int nseqs,
                       const int start, const int end)
{
    const double *times = model->times;
    const int nnodes = tree->nnodes;
    const LocalNode *nodes = tree->nodes;
    const double mintime = model->get_mintime();
    double invariant_lk = -1;
    lk_row table[nnodes];

    // get postorder
    int order[tree->nnodes];
    tree->get_postorder(order);

    // get mutation probabilities
    double muts[tree->nnodes];
    double nomuts[tree->nnodes];
    for (int i=0; i<tree->nnodes; i++) {
        if (i != tree->root) {
            double t = max(times[nodes[nodes[i].parent].age] -
                           times[nodes[i].age], mintime);
            muts[i] = prob_branch(t, model->mu, true);
            nomuts[i] = prob_branch(t, model->mu, false);
        }
    }


    // calculate emissions for tree at each site
    double lnl = 0.0;
    for (int i=start; i<end; i++) {
        double lk;
        bool invariant = is_invariant_site(seqs, nseqs, i);

        if (invariant && invariant_lk > 0)
            // use precommuted invariant site likelihood
            lk = invariant_lk;
        else {
            lk = likelihood_site_inner(tree, seqs, i, order, tree->nnodes,
                                       muts, nomuts, table);

            // save invariant likelihood
            if (invariant)
                invariant_lk = lk;
        }

        lnl += log(lk);
    }

    return lnl;
}



//=============================================================================
// emission calculation


void get_infinite_sites_states(const States &states, const LocalTree *tree,
                               const char *const *seqs, int nseqs, int seqlen,
                               bool *invariant,
                               bool internal, bool **valid_states)
{
    const int nstates = states.size();
    const int nnodes = tree->nnodes;

    if (internal) {
        // internal branch case
        const int maintree_root = tree->nodes[tree->root].child[1];
        const int subtree_root = tree->nodes[tree->root].child[0];

        // get nodes within trees
        int nsubnodes = 0, subnodes[nnodes];
        tree->get_preorder(subtree_root, subnodes, nsubnodes);
        reverse(subnodes, subnodes + nsubnodes);
        int nmainnodes = 0, mainnodes[nnodes];
        tree->get_preorder(maintree_root, mainnodes, nmainnodes);
        reverse(mainnodes, mainnodes + nmainnodes);

        for (int i=0; i<seqlen; i++) {
            if (invariant[i])
                continue;

            // get set of all bases in the subtree
            char subset = 0;
            for (int k=0; k<nsubnodes; k++) {
                int j = subnodes[k];
                if (tree->nodes[j].is_leaf() && seqs[j][i] != 'N')
                    subset |= 1 << dna2int[(int) seqs[j][i]];
            }

            // get set of all bases in the maintree
            char mainset = 0;
            for (int k=0; k<nmainnodes; k++) {
                int j = mainnodes[k];
                if (tree->nodes[j].is_leaf() && seqs[j][i] != 'N')
                    mainset |= 1 << dna2int[(int) seqs[j][i]];
            }

            // if subtree and main have distinct bases then all states are valid
            if (!(subset & mainset)) {
                for (int j=0; j<nstates; j++)
                    valid_states[i][j] = true;
                continue;
            }


            // infer ancestral reconstruction
            char bases[nnodes];
            parsimony_ancestral_set(
                tree, seqs, i, subnodes, nsubnodes, bases);
            int cset = bases[subtree_root];
            parsimony_ancestral_set(
                tree, seqs, i, mainnodes, nmainnodes, bases);

            bool valid_nodes[nnodes];
            for (int k=0; k<nmainnodes; k++) {
                int j = mainnodes[k];
                int parent = tree->nodes[j].parent;
                valid_nodes[j] = bool(cset & bases[j]) ||
                    (j != maintree_root && (cset & bases[parent]));
            }
            for (int j=0; j<nstates; j++)
                valid_states[i][j] = valid_nodes[states[j].node];

            // check for at least one valid state
            bool valid = false;
            for (int j=0; j<nstates; j++) {
                if (valid_states[i][j]) {
                    valid = true;
                    break;
                }
            }
            if (!valid) {
                printLog(LOG_LOW, "unable to satisfy infinite sites assumption");
            }
        }

    } else {
        // external branch  case
        int newleaf = tree->get_num_leaves();
        int postorder[nnodes];
        tree->get_postorder(postorder);

        for (int i=0; i<seqlen; i++) {
            if (invariant[i])
                continue;

            // get set of all bases in the tree
            char set = 0;
            for (int j=0; j<newleaf; j++)
                if (seqs[j][i] != 'N')
                    set |= 1 << dna2int[(int) seqs[j][i]];

            char c = seqs[newleaf][i];
            char cset = ((c != 'N') ? 1 << dna2int[(int) c] : 0);
            if (!(cset & set)) {
                // newleaf has a new base, thus newleaf can go anywhere
                for (int j=0; j<nstates; j++)
                    valid_states[i][j] = true;
                continue;
            }

            // infer ancestral reconstruction
            char bases[nnodes];
            parsimony_ancestral_set(
                tree, seqs, i, postorder, nnodes, bases);

            bool valid_nodes[nnodes];
            for (int j=0; j<nnodes; j++) {
                int parent = tree->nodes[j].parent;
                valid_nodes[j] = bool(cset & bases[j]) ||
                    (parent != -1 && (cset & bases[parent]));
            }
            for (int j=0; j<nstates; j++)
                valid_states[i][j] = valid_nodes[states[j].node];

            // check for at least one valid state
            bool valid = false;
            for (int j=0; j<nstates; j++) {
                if (valid_states[i][j]) {
                    valid = true;
                    break;
                }
            }
            if (!valid) {
                printLog(LOG_LOW, "unable to satisfy infinite sites assumption");
            }
        }
    }
}


// calculate emissions for external branch resampling
void calc_emissions(const States &states, const LocalTree *tree,
                    const char *const *seqs, int nseqs, int seqlen,
                    const ArgModel *model, bool internal, double **emit)
{
    const int nstates = states.size();
    const double mintime = model->get_mintime();
    const int newleaf = tree->get_num_leaves();
    const int maintree_root = internal ? tree->nodes[tree->root].child[1] :
        tree->root;
    const int subtree_root = internal ? tree->nodes[tree->root].child[0] :
        tree->root;


    // special case: ignore fully specified local tree
    if (internal && nstates == 0) {
        for (int i=0; i<seqlen; i++)
            emit[i][0] = 1.0;
        return;
    }


    // find invariant sites
    bool *invariant = new bool [seqlen];
    bool *masked = new bool [seqlen];
    find_invariant_sites(seqs, nseqs, seqlen, invariant);
    find_masked_sites(seqs, nseqs, seqlen, masked, invariant);


    // compute inner and outer likelihood tables
    LikelihoodTable inner(seqlen, tree->nnodes);
    LikelihoodTable inner_subtree(seqlen, 1);
    LikelihoodTable outer(seqlen, tree->nnodes);
    calc_inner_outer(tree, model, seqs, seqlen, invariant, internal,
                     inner.data, outer.data);


    if (!internal) {
        // compute inner table for new leaf
        for (int i=0; i<seqlen; i++) {
            const char c = seqs[newleaf][i];
            if (c == 'N') {
                inner_subtree.data[i][0][0] = 1.0;
                inner_subtree.data[i][0][1] = 1.0;
                inner_subtree.data[i][0][2] = 1.0;
                inner_subtree.data[i][0][3] = 1.0;
            } else {
                inner_subtree.data[i][0][0] = 0.0;
                inner_subtree.data[i][0][1] = 0.0;
                inner_subtree.data[i][0][2] = 0.0;
                inner_subtree.data[i][0][3] = 0.0;
                inner_subtree.data[i][0][dna2int[(int) c]] = 1.0;
            }
        }
    }


    // calc tree lengths
    int queue[tree->nnodes];
    int top = 0;

    double maintreelen = 0.0;
    // process in preorder maintree nodes
    queue[top++] = maintree_root;
    while (top > 0) {
        int node = queue[--top];
        if (node != maintree_root)
            maintreelen += max(tree->get_dist(node, model->times), mintime);
        if (!tree->nodes[node].is_leaf()) {
            queue[top++] = tree->nodes[node].child[0];
            queue[top++] = tree->nodes[node].child[1];
        }
    }

    double subtreelen = 0.0;
    if (internal) {
        // process in preorder subtree nodes
        queue[top++] = subtree_root;
        while (top > 0) {
            int node = queue[--top];
            if (node != subtree_root)
                subtreelen += max(tree->get_dist(node, model->times), mintime);
            if (!tree->nodes[node].is_leaf()) {
                queue[top++] = tree->nodes[node].child[0];
                queue[top++] = tree->nodes[node].child[1];
            }
        }
    }


    // populate emission table
    for (int j=0; j<nstates; j++) {
        State state = states[j];

        // get nodes
        int node1 = internal ? subtree_root : 0;
        int node2 = state.node;
        int parent = tree->nodes[node2].parent;

        // get times
        double time1 = internal ? model->times[tree->nodes[node1].age] : 0.0;
        double time2 = model->times[tree->nodes[node2].age];
        double parent_time = (parent != -1) ?
            model->times[min(tree->nodes[parent].age, model->ntimes-1)] : 0.0;
        double coal_time = model->times[state.time];

        // get distances
        double dist1 = max(coal_time - time1, mintime);
        double dist2 = max(coal_time - time2, mintime);
        double dist3 = max(parent_time - coal_time, mintime);

        // get mutation probabilities
        double mut1 = prob_branch(dist1, model->mu, true);
        double mut2 = prob_branch(dist2, model->mu, true);
        double mut3 = prob_branch(dist3, model->mu, true);
        double nomut1 = prob_branch(dist1, model->mu, false);
        double nomut2 = prob_branch(dist2, model->mu, false);
        double nomut3 = prob_branch(dist3, model->mu, false);

        // get tree length
        double treelen;
        if (node2 == maintree_root)
            treelen = maintreelen + subtreelen
                + max(coal_time - time1, mintime)
                + max(coal_time - model->times[tree->nodes[maintree_root].age],
                      mintime);
        else
            treelen = maintreelen + subtreelen
                + max(coal_time - time1, mintime);

        // calculate invariant_lk
        double invariant_lk = .25 * exp(- model->mu * max(treelen, mintime));

        // fill in row of emission table
        for (int i=0; i<seqlen; i++) {
            if (masked[i]) {
                // masked site
                emit[i][j] = 1.0;
            } else if (invariant[i]) {
                // invariant site
                emit[i][j] = invariant_lk;
            } else {
                // variant site
                lk_row *in = inner.data[i];
                lk_row *out = outer.data[i];
                lk_row *in2 = internal ? in : inner_subtree.data[i];

                emit[i][j] = 0.0;
                for (int a=0; a<4; a++) {
                    double p1 = 0.0, p2 = 0.0, p3 = 0.0;
                    for (int b=0; b<4; b++) {
                        if (a == b) {
                            p1 += in2[node1][b] * nomut1;
                            p2 += in[node2][b] * nomut2;
                            p3 += out[node2][b] * nomut3;
                        } else {
                            p1 += in2[node1][b] * mut1;
                            p2 += in[node2][b] * mut2;
                            p3 += out[node2][b] * mut3;
                        }
                    }

                    if (node2 != maintree_root) {
                        emit[i][j] += p1 * p2 * p3 * .25;
                    } else {
                        emit[i][j] += p1 * p2 * .25;
                    }
                }
            }
        }
    }

    // optionally enforce infinite sites model
    if (model->infsites_penalty < 1.0) {
        bool **valid_states = new_matrix<bool>(seqlen, nstates);
        get_infinite_sites_states(states, tree, seqs, nseqs, seqlen,
                                  invariant, internal, valid_states);
        for (int i=0; i<seqlen; i++) {
            if (!invariant[i]) {
                for (int j=0; j<nstates; j++)
                    if (!valid_states[i][j])
                        emit[i][j] *= model->infsites_penalty;
            }
        }
        delete_matrix<bool>(valid_states, seqlen);
    }


    // clean up
    delete [] invariant;
    delete [] masked;
}

// calculate emissions for external branch resampling
void calc_emissions_external(const States &states, const LocalTree *tree,
                             const char *const *seqs, int nseqs, int seqlen,
                             const ArgModel *model, double **emit)
{
    calc_emissions(states, tree, seqs, nseqs, seqlen, model, false, emit);
}

// calculate emissions for internal branch resampling
void calc_emissions_internal(const States &states, const LocalTree *tree,
                             const char *const *seqs, int nseqs, int seqlen,
                             const ArgModel *model, double **emit)
{
    calc_emissions(states, tree, seqs, nseqs, seqlen, model, true, emit);
}




//=============================================================================
// counting non-compatiable sites


void parsimony_ancestral_set(const LocalTree *tree, const char * const *seqs,
                             int pos, int *postorder, int npostorder,
                             char *ancestral)
{
    const int nnodes = tree->nnodes;
    const LocalNode *nodes = tree->nodes;
    char sets[nnodes];

    // clear sets
    for (int node=0; node<nnodes; node++)
        sets[node] = 0;

    // do unweighted parsimony by postorder traversal
    int postorder2[nnodes];
    if (!postorder) {
        tree->get_postorder(postorder2);
        postorder = postorder2;
        npostorder = nnodes;
    }
    for (int i=0; i<npostorder; i++) {
        int node = postorder[i];
        if (nodes[node].is_leaf()) {
            char c = seqs[node][pos];
            if (c == 'N')
                sets[node] = 0;
            else
                sets[node] = 1 << dna2int[(int) c];
        } else {
            char lset = sets[nodes[node].child[0]];
            char rset = sets[nodes[node].child[1]];
            char intersect = lset & rset;
            if (intersect > 0)
                sets[node] = intersect;
            else
                sets[node] = lset | rset;
        }
    }

    // traceback
    int root = postorder[npostorder-1];
    ancestral[root] = sets[root];

    // traceback with preorder traversal
    for (int i=npostorder-2; i>=0; i--) {
        int node = postorder[i];
        char s = sets[node];
        char pset = ancestral[nodes[node].parent] & s;
        if (pset) {
            // use parent char if possible
            ancestral[node] = pset;
        } else {
            // otherwise do not refine set
            ancestral[node] = s;
        }
    }
}


void parsimony_ancestral_seq(const LocalTree *tree, const char * const *seqs,
                             int nseqs, int pos, char *ancestral,
                             int *postorder)
{
    const int nnodes = tree->nnodes;
    const LocalNode *nodes = tree->nodes;
    char sets[nnodes];
    int pchar;

    // clear sets
    for (int node=0; node<nnodes; node++)
        sets[node] = 0;

    // do unweighted parsimony by postorder traversal
    int postorder2[nnodes];
    if (!postorder) {
        tree->get_postorder(postorder2);
        postorder = postorder2;
    }
    for (int i=0; i<nnodes; i++) {
        int node = postorder[i];
        if (nodes[node].is_leaf()) {
            char c = seqs[node][pos];
            if (c == 'N')
                sets[node] = 1 + 2 + 4 + 8;
            else
                sets[node] = 1 << dna2int[(int) c];
        } else {
            char lset = sets[nodes[node].child[0]];
            char rset = sets[nodes[node].child[1]];
            char intersect = lset & rset;
            if (intersect > 0)
                sets[node] = intersect;
            else
                sets[node] = lset | rset;
        }
    }

    // traceback
    // arbitrary choose root base from set
    int root = postorder[nnodes-1];
    char rootset = sets[root];
    ancestral[root] = (rootset & 1) ? int2dna[0] :
        (rootset & 2) ? int2dna[1] :
        (rootset & 4) ? int2dna[2] : int2dna[3];

    // traceback with preorder traversal
    for (int i=nnodes-2; i>=0; i--) {
        int node = postorder[i];
        char s = sets[node];

        switch (s) {
        case 1: // just A
            ancestral[node] = int2dna[0];
            break;
        case 2: // just C
            ancestral[node] = int2dna[1];
            break;
        case 4: // just G
            ancestral[node] = int2dna[2];
            break;
        case 8: // just T
            ancestral[node] = int2dna[3];
            break;
        default:
            pchar = ancestral[nodes[node].parent];
            if (dna2int[pchar] & s) {
                // use parent char if possible
                ancestral[node] = pchar;
            } else {
                // use arbitrary char otherwise
                ancestral[node] = (s & 1) ? int2dna[0] :
                    (s & 2) ? int2dna[1] :
                    (s & 4) ? int2dna[2] : int2dna[3];
            }
        }
    }
}


int parsimony_cost_seq(const LocalTree *tree, const char * const *seqs,
                        int nseqs, int pos, int *postorder)
{
    const int nnodes = tree->nnodes;
    const LocalNode *nodes = tree->nodes;
    const int maxcost = 100000;
    int costs[nnodes][4];

    // do unweighted parsimony by postorder traversal
    int postorder2[nnodes];
    if (!postorder) {
        tree->get_postorder(postorder2);
        postorder = postorder2;
    }
    for (int i=0; i<nnodes; i++) {
        int node = postorder[i];
        if (tree->nodes[node].is_leaf()) {
            char c = seqs[node][pos];
            if (c == 'N') {
                for (int a=0; a<4; a++)
                    costs[node][a] = 0;
            } else {
                for (int a=0; a<4; a++)
                    costs[node][a] = maxcost;
                costs[node][dna2int[(int)c]] = 0;
            }
        } else {
            int *left_costs = costs[nodes[node].child[0]];
            int *right_costs = costs[nodes[node].child[1]];

            for (int a=0; a<4; a++) {
                int left_min = maxcost;
                int right_min = maxcost;
                for (int b=0; b<4; b++) {
                    left_min = min(left_min, int(a != b) + left_costs[b]);
                    right_min = min(right_min, int(a != b) + right_costs[b]);
                }
                costs[node][a] = left_min + right_min;
            }
        }
    }

    int root_min = maxcost;
    for (int a=0; a<4; a++)
        root_min = min(root_min, costs[tree->root][a]);

    return root_min;
}


int count_noncompat(const LocalTree *tree, const char * const *seqs,
                    int nseqs, int seqlen, int *postorder)
{
    // get postorder
    int postorder2[tree->nnodes];
    if (!postorder) {
        tree->get_postorder(postorder2);
        postorder = postorder2;
    }

    int noncompat = 0;
    for (int i=0; i<seqlen; i++)
        if (!is_invariant_site(seqs, nseqs, i)) {
            int a = count_alleles(seqs, nseqs, i);
            int c = parsimony_cost_seq(tree, seqs, nseqs, i, postorder);
            noncompat += int(c > a - 1 );
        }

    return noncompat;
}


int count_noncompat(const LocalTrees *trees, const char * const *seqs,
                    int nseqs, int seqlen)
{
    int noncompat = 0;

    int end = trees->start_coord;
    for (LocalTrees::const_iterator it=trees->begin();
         it != trees->end(); ++it)
    {
        int start = end;
        end += it->blocklen;
        LocalTree *tree = it->tree;
        int blocklen = it->blocklen;

        // get subsequence block
        char const* subseqs[nseqs];
        for (int i=0; i<nseqs; i++)
            subseqs[i] = &seqs[i][start];

        noncompat += count_noncompat(tree, subseqs, nseqs, blocklen, NULL);
    }

    return noncompat;
}



//=============================================================================
// slow literal emission calculation
// useful for testing against


void calc_emissions_external_slow(
    const States &states, const LocalTree *tree,
    const char *const *seqs, int nseqs, int seqlen,
    const ArgModel *model, double **emit)
{
    const int nstates = states.size();
    const int newleaf = tree->get_num_leaves();
    bool *invariant = new bool [seqlen];
    LikelihoodTable table(seqlen, tree->nnodes+2);

    // create local tree we can edit
    LocalTree tree2(tree->nnodes, tree->nnodes + 2);
    tree2.copy(*tree);

    // find invariant sites
    find_invariant_sites(seqs, nseqs, seqlen, invariant);

    for (int j=0; j<nstates; j++) {
        State state = states[j];
        add_tree_branch(&tree2, state.node, state.time);

        likelihood_sites(&tree2, model, seqs, seqlen, j, invariant, emit,
                         table.data, -1, -1);

        remove_tree_branch(&tree2, newleaf, NULL);
    }

    delete [] invariant;
}


void calc_emissions_internal_slow(
    const States &states, const LocalTree *tree,
    const char *const *seqs, int nseqs, int seqlen,
    const ArgModel *model, double **emit)
{
    const int nstates = states.size();
    const int subtree_root = tree->nodes[tree->root].child[0];
    const int subtree_root_age = tree->nodes[subtree_root].age;
    const int maxtime = model->ntimes + 1;

    // ignore fully specified local tree
    if (nstates == 0) {
        for (int i=0; i<seqlen; i++)
            emit[i][0] = 1.0;
        return;
    }

    bool *invariant = new bool [seqlen];
    LikelihoodTable table(seqlen, tree->nnodes+2);

    // create local tree we can edit
    LocalTree tree2(tree->nnodes, tree->nnodes + 2);
    tree2.copy(*tree);

    // find invariant sites
    find_invariant_sites(seqs, nseqs, seqlen, invariant);

    for (int j=0; j<nstates; j++) {
        State state = states[j];
        assert(subtree_root != tree2.root);

        Spr add_spr(subtree_root, subtree_root_age, state.node, state.time);
        apply_spr(&tree2, add_spr);

        likelihood_sites(&tree2, model, seqs, seqlen, j, invariant, emit,
                         table.data, -1, -1);

        Spr remove_spr(subtree_root, subtree_root_age,
                       tree2.root, maxtime);
        apply_spr(&tree2, remove_spr);
    }

    // clean up
    delete [] invariant;
}



//=============================================================================
// assert emissions

bool assert_emissions(const States &states, const LocalTree *tree,
                      const char *const *seqs, int nseqs, int seqlen,
                      const ArgModel *model)
{
    const int nstates = states.size();

    double **emit = new_matrix<double>(seqlen, nstates);
    double **emit2 = new_matrix<double>(seqlen, nstates);

    calc_emissions_external(states, tree, seqs, nseqs, seqlen, model, emit);
    calc_emissions_external_slow(states, tree, seqs, nseqs, seqlen, model, emit2);

    // compare emission tables
    for (int j=0; j<nstates; j++) {
        for (int i=0; i<seqlen; i++) {
            bool invar = is_invariant_site(seqs, nseqs, i);

            printLog(LOG_MEDIUM, ">> %d,%d: %e %e   invar=%d\n",
                     i, j, emit[i][j], emit2[i][j], invar);
            if (!fequal(emit[i][j], emit2[i][j], .0001, 1e-12))
                return false;
        }
    }

    delete_matrix<double>(emit, seqlen);
    delete_matrix<double>(emit2, seqlen);

    return true;
}


bool assert_emissions_internal(const States &states, const LocalTree *tree,
                               const char *const *seqs, int nseqs, int seqlen,
                               const ArgModel *model)
{
    const int nstates = states.size();

    double **emit = new_matrix<double>(seqlen, nstates);
    double **emit2 = new_matrix<double>(seqlen, nstates);

    calc_emissions_internal(states, tree, seqs, nseqs, seqlen,
                              model, emit);
    calc_emissions_internal_slow(states, tree, seqs, nseqs, seqlen,
                                   model, emit2);

    // compare emission tables
    for (int j=0; j<nstates; j++) {
        for (int i=0; i<seqlen; i++) {
            printLog(LOG_MEDIUM, ">> %d,%d: %e %e\n",
                     i, j, emit[i][j], emit2[i][j]);
            if (!fequal(emit[i][j], emit2[i][j], .0001, 1e-12))
                return false;
        }
    }

    delete_matrix<double>(emit, seqlen);
    delete_matrix<double>(emit2, seqlen);

    return true;
}


//=============================================================================
// C interface
extern "C" {

double **new_emissions(intstate *istates, int nstates,
                       int *ptree, int nnodes, int *ages_index,
                       char **seqs, int nseqs, int seqlen,
                       double *times, int ntimes,
                       double mu)
{
    States states;
    make_states(istates, nstates, states);
    LocalTree tree(ptree, nnodes, ages_index);
    ArgModel model(ntimes, times, NULL, 0.0, mu);

    double **emit = new_matrix<double>(seqlen, nstates);
    calc_emissions_external(states, &tree, seqs, nseqs, seqlen, &model, emit);

    return emit;
}


void delete_emissions(double **emit, int seqlen)
{
    delete_matrix<double>(emit, seqlen);
}


bool arghmm_assert_emit(
    LocalTrees *trees, int ntimes, double *times, double mu,
    char **seqs, int nseqs, int seqlen)
{
    ArgModel model(ntimes, times, NULL, 0.0, mu);
    States states;

    int end = trees->start_coord;
    for (LocalTrees::const_iterator it=trees->begin(); it!=trees->end(); ++it) {
        int start = end;
        int blocklen = it->blocklen;
        end = start + blocklen;
        LocalTree *tree = it->tree;

        get_coal_states(tree, model.ntimes, states);

        // get subsequence
        char *seqs2[nseqs];
        for (int j=0; j<nseqs; j++)
            seqs2[j] = &seqs[j][start];

        if (!assert_emissions(states, tree, seqs2, nseqs, blocklen, &model))
            return false;
    }

    return true;
}

bool arghmm_assert_emit_internal(
    LocalTrees *trees, int ntimes, double *times, double mu,
    char **seqs, int nseqs, int seqlen)
{
    ArgModel model(ntimes, times, NULL, 0.0, mu);
    const int maxtime = model.ntimes + 1;
    States states;
    int *removal_path = new int [trees->get_num_trees()];

    // randomly choose branch to remove
    LocalTrees trees2(*trees);

    // ramdomly choose a removal path
    sample_arg_removal_path_uniform(&trees2, removal_path);
    remove_arg_thread_path(&trees2, removal_path, maxtime);

    delete [] removal_path;

    int end = trees->start_coord;
    for (LocalTrees::const_iterator it=trees2.begin(); it!=trees2.end(); ++it)
    {
        int start = end;
        int blocklen = it->blocklen;
        end = start + blocklen;
        LocalTree *tree = it->tree;

        get_coal_states_internal(tree, model.ntimes, states);

        // get subsequence
        char *seqs2[nseqs];
        for (int j=0; j<nseqs; j++)
            seqs2[j] = &seqs[j][start];

        if (!assert_emissions_internal(states, tree, seqs2, nseqs, blocklen, &model))
            return false;
    }

    return true;
}



} // extern "C"



} // namespace argweaver


