#
# HMM emission related functions
#

from math import exp, log

#=============================================================================


def parsimony_ancestral_seq(tree, seqs, pos):
    """Calculates ancestral sequence for a local tree using parsimony"""

    ancestral = {}
    sets = {}

    # do unweight parsimony
    for node in tree.postorder():
        if node.is_leaf():
            sets[node] = set([seqs[node.name][pos]])
        else:
            lset = sets[node.children[0]]
            rset = sets[node.children[1]]
            intersect = lset & rset
            if len(intersect) > 0:
                sets[node] = intersect
            else:
                sets[node] = lset | rset

    # traceback
    for node in tree.preorder():
        s = sets[node]
        if len(s) == 1 or not node.parents:
            # NOTE: this technique is used to make assignment deterministic
            ancestral[node.name] = ("A" if "A" in s else
                                    "C" if "C" in s else
                                    "G" if "G" in s else
                                    "T")
        else:
            pchar = ancestral[node.parents[0].name]
            if pchar in s:
                ancestral[node.name] = pchar
            else:
                ancestral[node.name] = ("A" if "A" in s else
                                        "C" if "C" in s else
                                        "G" if "G" in s else
                                        "T")

    return ancestral



def calc_emission(tree, model, pos, new_name):
    """
    Calculates emissions for all states at positions 'pos'
    """

    mu = model.mu
    seqs = model.seqs
    mintime = model.time_steps[0]
    emit = []

    for node_name, timei in model.states[pos]:
        node = tree[node_name]
        time = model.times[timei]

        local_site = parsimony_ancestral_seq(tree, seqs, pos)

        # v = new chromosome
        # x = current branch
        # p = parent of current branch

        if node.parents:
            parent = node.parents[0]
            parent_age = parent.age

            if not parent.parents:
                # unwrap top branch
                c = parent.children
                sib = (c[1] if node == c[0] else c[0])

                v = seqs[new_name][pos]
                x = local_site[node.name]
                p = local_site[sib.name]

                # modify (x,p) length to (x,p) + (sib,p)
                parent_age = 2 * parent_age - sib.age

            else:
                v = seqs[new_name][pos]
                x = local_site[node.name]
                p = local_site[parent.name]

        else:
            parent = None
            parent_age = None

            # adjust time by unwrapping branch
            time = 2 * time - node.age

            v = seqs[new_name][pos]
            x = local_site[node.name]
            p = x

        time = max(time, mintime)

        if v == x == p:
            # no mutation
            emit.append(- mu * time)

        elif v != p == x:
            # mutation on v
            emit.append(log(.3333 - .3333 * exp(-mu * time)))

        elif v == p != x:
            # mutation on x
            t1 = max(parent_age - node.age, mintime)
            t2 = max(time - node.age, mintime)

            emit.append(log((1 - exp(-mu *t2)) / (1 - exp(-mu * t1))
                            * exp(-mu * (time + t2 - t1))))

        elif v == x != p:
            # mutation on (y,p)
            t1 = max(parent_age - node.age, mintime)
            t2 = max(parent_age - time, mintime)

            emit.append(log((1 - exp(-mu * t2)) / (1 - exp(-mu * t1))
                            * exp(-mu * (time + t2 - t1))))

        else:
            # two mutations (v,x)
            # mutation on x
            if parent:
                t1 = max(parent_age - node.age, mintime)
                t2a = max(parent_age - time, mintime)
            else:
                t1 = max(model.times[-1] - node.age, mintime)
                t2a = max(model.times[-1] - time, mintime)
            t2b = max(time - node.age, mintime)
            t2 = max(t2a, t2b)
            t3 = time

            emit.append(log((1 - exp(-mu *t2)) * (1 - exp(-mu *t3))
                            / (1 - exp(-mu * t1))
                            * exp(-mu * (time + t2 + t3 - t1))))

    return emit
