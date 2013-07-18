"""
 Birth death process and reconstructed process for trees

"""

from math import *
import random
from rasmus import util, stats, treelib


def prob_birth_death1(ngenes, t, birth, death):
    """
    Returns the probability that one lineage leaves 'ngenes' genes
    after time 't'
    """

    # special cases
    if birth == death:
        if birth == 0.0:
            if ngenes == 1:
                return 1.0
            else:
                return 0.0
        ut = t / (1.0 / birth + t)
        if ngenes == 0:
            return ut
        else:
            return ((1.0 - ut)**2) * (ut**(ngenes-1))
    
    l = birth
    u = death
    r = l - u
    a = u / l

    ut = (1.0 - exp(-r*t)) / (1.0 - a * exp(-r*t))
    p0 = a*ut
    
    if ngenes == 0:
        return p0
    
    return (1.0 - p0)*(1.0 - ut) * (ut**(ngenes-1))


def prob_birth_death(genes1, genes2, t, birth, death):
    """Probability of 'genes1' genes at time 0 give rise to 'genes2' genes at
       time 't' with 'birth' and 'death' rates.
    """

    # special cases
    if birth == 0.0 and death == 0.0:
        if genes1 == genes2:
            return 1.0
        else:
            return 0.0
    
    
    l = birth
    u = death
    elut = exp((l-u)*t)
    a = u * (elut - 1.0) / (l*elut - u) # alpha
    b = l * (elut - 1.0) / (l*elut - u) # beta
    n = genes1
    i = genes2

    if genes1 < 1:
        return 0.0

    if genes2 == 0:
        return a ** n
    else:
        return sum(stats.choose(n,j) * stats.choose(n+i-j-1, n-1) *\
                   a**(n-j) * b**(i-j) * (1.0 - a - b)**j
                   for j in xrange(min(n, i)+1))
   

def birth_wait_time(t, n, T, birth, death):
    """Probability density for for next birth at time 't' given
       'n' lineages starting at time 0, evolving until time 'T' with a
       'birth' and 'death' rates for a reconstructed process.
    """

    # special case
    if birth == death:
        t2 = t - T
        nl = 1.0 / birth
        return birth * n * (-nl + t2)**n / (-nl - T)**n / (1.0 - birth * t2)
        
    l = birth
    u = death
    r = l - u
    a = u / l

    return n * r * exp(-n*r*t) * \
           ((1.0 - a * exp(-r * (T - t)))**(n-1)) / \
           ((1.0 - a * exp(-r * T))**n)


def prob_no_birth(n, T, birth, death):
    """Probability of no birth from 'n' lineages starting at time 0, 
       evolving until time 'T' with 'birth' and 'death' rates
       for a reconstructed process.
    """

    # special cases
    if birth == 0.0:
        return 1.0
    elif birth == death:
        return 1.0 / (1.0 + birth * T)**n
    
    l = birth
    u = death
    r = l - u
    a = u / l

    return (1.0 - (l*(1.0 - exp(-r * T))) / \
                  (l - u * exp(-r * T))) ** n



def num_topology_histories(node, leaves=None):
    """
    Returns the number of labeled histories for a topology

    The topology is specified by a root 'node' and a set of leaves 'leaves'.
    If leaves are not specified, the leaves of 'node' will be used.
    """

    # TODO: can simplify
    
    if leaves is None:
        leaves = node.leaves()
    leaves = set(leaves)

    prod = [1.0]

    def walk(node):        
        if node in leaves:
            return 0
        else:
            internals = map(walk, node.children)
            prod[0] *= stats.choose(sum(internals), internals[0])
            return 1 + sum(internals)
    walk(node)

    return prod[0]



#=============================================================================
# sampling


def sample_birth_death_count(n, t, birth, death):
    """
    Sample the gene count at time 't' with birth death rates
    'birth' and 'death' with starting count 'n'
    """

    t2 = 0.0
    bd = float(birth + death)

    while t2 < t and n > 0:
        rate = bd * n
        next_t = random.expovariate(rate)

        if t2 + next_t < t:
            if random.random() < birth / bd:
                # birth
                n += 1
            else:
                # death
                n -= 1

        t2 += next_t

    return n


def sample_birth_wait_time(n, T, birth, death):
    """
    Sample the next birth event from a reconstructed birthdeath process.
    Let there be 'n' lineages at time 0 that evolve until time 'T' with
    'birth' and 'death' rates.

    Conditioned that a birth will occur
    """
    
    # TODO: could make this much more efficient
    
    # uses rejection sampling
    start_y = birth_wait_time(0, n, T, birth, death)
    end_y = birth_wait_time(T, n, T, birth, death)
    M = max(start_y, end_y)
    
    while True:
        t = random.uniform(0, T)
        f = birth_wait_time(t, n, T, birth, death)

        if random.uniform(0, 1) <= f / M:
            return t


def sample_birth_death_tree(T, birth, death, tree=None, node=None,
                            keepdoom=False):
    """Simulate a reconstructed birth death tree"""
    
    # create tree if one is not given
    if tree is None:
        tree = treelib.Tree()
    
    # create starting node if one is not given
    if node is None:
        tree.make_root() 
        node = tree.root
    else:
        node = tree.add_child(node, tree.new_node())

    bd_rate = float(birth + death)
    doom = set()
    
    def walk(T, node):
        if bd_rate == 0.0:
            next_t = util.INF
        else:
            next_t = random.expovariate(bd_rate)
        
        if next_t > T:
            # finish branch
            node.dist = T
            
        elif random.random() < birth / bd_rate:
            # birth
            node.dist = next_t

            node2 = tree.add_child(node, tree.new_node())
            walk(T - next_t, node2)

            node2 = tree.add_child(node, tree.new_node())
            walk(T - next_t, node2)
            
        else:
            # death
            node.dist = next_t
            doom.add(node)
    walk(T, node)

    if not keepdoom:
        leaves = set(tree.leaves()) - doom
        treelib.subtree_by_leaves(tree, leaves)

        if len(leaves) == 0:
            doom.add(tree.root)
    
    return tree, doom


def sample_birth_death_gene_tree(stree, birth, death, 
                                 genename=lambda sp, x: sp + "_" + str(x),
                                 removeloss=True):
    """Simulate a gene tree within a species tree with birth and death rates"""
    
    # initialize gene tree
    tree = treelib.Tree()
    tree.make_root()
    recon = {tree.root: stree.root}
    events = {tree.root: "spec"}
    losses = set()
    
    def walk(snode, node):
        if snode.is_leaf():
            tree.rename(node.name, genename(snode.name, node.name))
            events[node] = "gene"
        else:
            for child in snode:
                # determine if loss will occur
                tree2, doom = sample_birth_death_tree(
                    child.dist, birth, death, 
                    tree=tree, node=node, keepdoom=True)
                
                # record reconciliation
                next_nodes = []
                def walk2(node):
                    node.recurse(walk2)
                    recon[node] = child
                    if node in doom:
                        losses.add(node)
                        events[node] = "gene"
                    elif node.is_leaf():
                        events[node] = "spec"
                        next_nodes.append(node)
                    else:
                        events[node] = "dup"
                walk2(node.children[-1])
                
                # recurse
                for leaf in next_nodes:
                    walk(child, leaf)
            
            # if no child for node then it is a loss
            if node.is_leaf():
                losses.add(node)
    walk(stree.root, tree.root)
    
    
    # remove lost nodes
    if removeloss:
        treelib.remove_exposed_internal_nodes(tree,
                                              set(tree.leaves()) - losses)
        treelib.remove_single_children(tree, simplify_root=False)
        
        delnodes = set()
        for node in recon:
            if node.name not in tree.nodes:
                delnodes.add(node)
        for node in delnodes:
            del recon[node]
            del events[node]

    if len(tree.nodes) <= 1:
        tree.nodes = {tree.root.name : tree.root}
        recon = {tree.root: stree.root}
        events = {tree.root: "spec"}
    
    return tree, recon, events    






#=============================================================================
# testing functions
# These are functions that are not efficient but implement the distributions
# in a more literal way, thus they are handy to test against.



def sample_birth1_literal(T, birth, death):
    """
    Sample the next birth from a reconstructed birth death process
    starting with only 1 lineage.

    This function does not condition on the survival of the lineage.    

    T     -- stopping time
    birth -- rate of birth
    death -- rate of death

    Returns (t, alive)
    t is a float of the first birth or None if none occurs.
    alive is True if the lineage is still alive, False if extinct.

    NOTE: This function uses a very literal way of performing the sampling.
    It is only good for testing purposes.
    """
    
    # sample events
    t1 = random.expovariate(birth)
    t2 = random.expovariate(death)

    # both events occur after stopping time T, simulation is done
    if t1 >= T and t2 >= T:
        return None, True

    if t2 < t1:
        # death occurs
        return None, False
    else:
        # birth occurs

        # recurse
        t3, alive3 = sample_birth1_literal(T - t1, birth, death)
        t4, alive4 = sample_birth1_literal(T - t1, birth, death)

        if alive3 and alive4:
            # if both lineages are alive then our birth is in the
            # reconstructed process
            return t1, True
        elif alive3:
            # lineage 3 contains the first birth of the recon proc
            if t3 is not None:
                return t1 + t3, True
            else:
                return None, True
        elif alive4:
            # lineage 4 contains the first birth of the recon proc
            if t4 is not None:
                return t1 + t4, True
            else:
                return None, True
        else:
            # both lineages died, so we are dead
            return None, False




def sample_birth_literal(n, T, birth, death):
    """
    Sample the next birth from a reconstructed birth death process
    
    n     -- number of current lineages
    T     -- stopping time
    birth -- rate of birth
    death -- rate of death
    """
    
    tmin = util.INF
        
    for i in xrange(n):
        # require each lineage to be alive
        while True:
            t, alive = sample_birth1_literal(T, birth, death)
            if alive:
                break

        if t is not None:
            tmin = min(tmin, t)

    if tmin < T:
        return tmin
    else:
        return None


