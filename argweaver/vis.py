
import os
import subprocess

import argweaver
import argweaver.smc

from compbio import arglib


def inorder_tree(tree):
    queue = [("queue", tree.root)]

    while queue:
        cmd, node = queue.pop()

        if cmd == "visit":
            yield node
        elif cmd == "queue":
            if node.is_leaf():
                yield node
            else:
                queue.extend(
                    [("queue", node.children[1]),
                     ("visit", node),
                     ("queue", node.children[0])])


def layout_tree_leaves(tree, age_func=lambda node: node.age):
    layout = {}
    y = 0

    for node in inorder_tree(tree):
        if node.is_leaf():
            layout[node.name] = y
        else:
            y += (age_func(node) / 1e3) + 1

    vals = layout.values()
    mid = (max(vals) + min(vals)) / 2.0
    for k, v in layout.items():
        layout[k] = (v - mid)

    return layout


def orient_tree(tree, last_tree, recomb_node):
    """
    Flip nodes in tree to match last_tree.
    """
    order = dict((n, i) for i, n in enumerate(last_tree.leaf_names()))

    # make all leaves of recomb node have max order
    for leaf in recomb_node.leaf_names():
        order[leaf] = len(order)

    # flip nodes in tree
    for node in tree.postorder():
        if node.is_leaf():
            continue

        assert len(node.children) == 2
        if order[node.children[0].name] > order[node.children[1].name]:
            node.children = [node.children[1], node.children[0]]

        order[node.name] = min(order[node.children[0].name],
                               order[node.children[1].name])

    # flip node above subtree if needed
    inorder = dict((n.name, i) for i, n in enumerate(inorder_tree(last_tree)))
    node = tree[recomb_node.name].parent
    if inorder[node.children[0].name] > inorder[node.children[1].name]:
        node.children = [node.children[1], node.children[0]]


def iter_layout_smc(smc, names=None):
    age_func = lambda node: node.data["age"]
    last_tree = None
    spr = None

    for item in smc:
        if item["tag"] == "NAMES":
            names = item["names"]

        if item["tag"] == "SPR":
            spr = item

        elif item["tag"] == "TREE":
            assert names is not None
            tree = item["tree"]
            if isinstance(tree, basestring):
                raise Exception("Trees need to be parsed")

            block = [item["start"]-1, item["end"]]

            if last_tree:
                orient_tree(tree, last_tree, last_tree[spr["recomb_node"]])
            layout = layout_tree_leaves(tree, age_func=age_func)
            layout2 = dict((name, layout[i]) for i, name in enumerate(names))

            yield block, layout2

            last_tree = tree


class ArgLayout(object):
    """
    """
    def __init__(self):
        self.chrom = "chr"
        self.blocks = []
        self.leaf_layout = []

    def layout_smc(self, smc):
        names = smc.header["names"]
        self.chrom = smc.header["chrom"]
        self.blocks = []
        self.leaf_layout = []

        for block, leaf_layout in iter_layout_smc(smc, names=names):
            self.blocks.append(block)
            self.leaf_layout.append(leaf_layout)

    def read(self, filename):
        self.blocks = []
        self.leaf_layout = []

        for block, leaf_layout in iter_arg_layout(filename):
            self.chrom = block[0]
            self.blocks.append([block[1] - 1, block[2]])
            self.leaf_layout.append(leaf_layout)


def iter_arg_layout(filename):
    """
    Iterate through an ARG layout file.
    """
    with argweaver.open_stream(filename, compress='bgzip') as infile:
        for line in infile:
            tokens = line.rstrip().split("\t")
            block = [tokens[0], int(tokens[1]), int(tokens[2])]
            leaf_layout = {}
            for i in range(3, len(tokens), 2):
                leaf_layout[tokens[i]] = float(tokens[i+1])
            yield block, leaf_layout


def index_arg_layout(filename):
    subprocess.call(["tabix", "-s", "1", "-b", "2", "-e", "3", "-f", filename])


def query_arg_layout(filename, chrom, start, end):
    cmd = ["tabix", filename, "%s:%d-%d" % (chrom, start, end)]
    null = open(os.devnull, 'w')
    infile = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                              stderr=null).stdout
    return iter_arg_layout(infile)


def layout_arg(arg, start=None, end=None):

    if start is None:
        start = arg.start
    if end is None:
        end = arg.end

    tree = arg.get_marginal_tree(start)
    arglib.remove_single_lineages(tree)
    last_pos = start
    blocks = []
    leaf_layout = []

    layout_func = layout_tree_leaves

    for spr in arglib.iter_arg_sprs(
            arg, start=start, end=end, use_leaves=True):
        blocks.append([last_pos, spr[0]])
        leaf_layout.append(layout_func(tree))
        inorder = dict((n, i) for i, n in enumerate(inorder_tree(tree)))

        # determine SPR nodes
        rnode = arglib.arg_lca(tree, spr[1][0], spr[0])
        cnode = arglib.arg_lca(tree, spr[2][0], spr[0])

        # determine best side for adding new sister
        left = (inorder[rnode] < inorder[cnode])

        # apply spr
        arglib.apply_spr(tree, rnode, spr[1][1], cnode, spr[2][1], spr[0])

        # adjust sister
        rindex = rnode.parents[0].children.index(rnode)
        if left and rindex != 0:
            rnode.parents[0].children.reverse()

        last_pos = spr[0]

    blocks.append([last_pos, end])
    leaf_layout.append(layout_func(tree))

    return blocks, leaf_layout
