"""
  SMC input/output (sequentially Markov coalescent)

  The SMC data structure is used to describe an ARG using a series of
  local trees and SPR (subtree pruning and regrafting) operations.
"""

from contextlib import closing
import StringIO

import argweaver

from compbio import arglib
from rasmus import treelib
from rasmus import util


class SMCReader (object):
    """
    Reads an SMC file.
    """

    def __init__(self, filename, parse_trees=False, apply_spr=False):

        if isinstance(filename, basestring):
            # read SMC from file
            self._filename = filename
            self._infile = iter_smc_file(filename,
                                         parse_trees=parse_trees,
                                         apply_spr=apply_spr)
            self._stream = util.PushIter(self._infile)
        else:
            # read SMC from list of items
            self._filename = None
            self._infile = None
            self._stream = util.PushIter(filename)

        # read header
        header_tags = ("NAMES", "REGION")
        self.header = {}
        while self._stream.peek({"tag": ""})["tag"] in header_tags:
            self.header.update(self._stream.next())
            del self.header["tag"]

    def __iter__(self):
        return self

    def next(self):
        return self._stream.next()

    def peek(self, default=None):
        return self._stream.peek(default)

    def parse_tree(self, text):
        tree = parse_tree(text, self.header["names"])
        return tree

    def close(self):
        if self._infile:
            self._infile.close()


def read_smc(filename, parse_trees=False, apply_spr=False):
    """
    Read an SMC file.

    Read an entire SMC file into a list.
    """
    return list(iter_smc_file(filename, parse_trees=parse_trees,
                              apply_spr=apply_spr))


def iter_smc_file(filename, parse_trees=False, apply_spr=False,
                  region=None):
    """
    Iterates through a SMC file.

    parse_trees: If True, parses local trees.
    apply_spr: If True, avoids reading each tree by applying the SPR
        operation to the current tree.
    region: If given, returns only trees and SPRs within region=(start, end).

    Yields item, where item can be one of the following:
        {'tag': 'NAMES',
         'names': names_of_sequences}

        {'tag': 'REGION',
         'chrom': name_of_chromosome,
         'start': start_coordinate_of_region,
         'end': end_coordinate_of_region}

        {'tag': 'TREE',
         'start': start_coordinate_of_local_region,
         'end': end_coordinate_of_local_region,
         'tree': local_tree}

        {'tag': 'SPR',
         'pos': coordinate of recombination point,
         'recomb_node': name_of_recombination_node,
         'recomb_time': time_of_recombination,
         'coal_node': name_of_branch_with_recoalescence,
         'coal_time': time_of_recoalescence}
    """

    if region:
        tree = None
        spr = None

        for item in iter_subsmc(iter_smc_file(filename), region):

            if item["tag"] == "SPR":
                spr = item
            elif item["tag"] == "TREE":
                if parse_trees:
                    if apply_spr and tree is not None and spr is not None:
                        smc_apply_spr(tree, spr)
                    else:
                        tree = parse_tree(item["tree"])
                    item["tree"] = tree
            yield item
        return

    with closing(argweaver.open_stream(filename)) as infile:
        spr = None
        tree = None

        for line in infile:
            line = line.rstrip()
            tokens = line.split("\t")

            if tokens[0] == "NAMES":
                yield {"tag": "NAMES", "names": tokens[1:]}

            elif tokens[0] == "REGION":
                yield {"tag": "REGION",
                       "chrom": tokens[1],
                       "start": int(tokens[2]),
                       "end": int(tokens[3])}

            elif tokens[0] == "RANGE":
                raise Exception("deprecated RANGE line, use REGION instead")

            elif tokens[0] == "TREE":
                tree_text = tokens[3]
                if parse_trees:
                    if apply_spr and tree is not None and spr is not None:
                        smc_apply_spr(tree, spr)
                    else:
                        tree = parse_tree(tree_text)
                else:
                    tree = tree_text

                yield {"tag": "TREE",
                       "start": int(tokens[1]),
                       "end": int(tokens[2]),
                       "tree": tree}

            elif tokens[0] == "SPR":
                spr = {"tag": "SPR",
                       "pos": int(tokens[1]),
                       "recomb_node": int(tokens[2]),
                       "recomb_time": float(tokens[3]),
                       "coal_node": int(tokens[4]),
                       "coal_time": float(tokens[5])}
                yield spr


def iter_subsmc(smc, region):
    """
    Iterate through a region of an SMC stream.
    """
    for item in smc:
        if item["tag"] == "NAMES":
            yield item

        elif item["tag"] == "REGION":
            item["start"] = max(region[0], item["start"])
            item["end"] = min(region[1], item["end"])
            yield item

        elif item["tag"] == "SPR":
            if region[0] <= item["pos"] < region[1]:
                yield item

        elif item["tag"] == "TREE":
            if item["start"] <= region[1] and item["end"] >= region[0]:
                item["start"] = max(region[0], item["start"])
                item["end"] = min(region[1], item["end"])

                yield item


def smc_apply_spr(tree, spr):
    """
    Apply an SPR operation to a local tree.
    """

    recomb = tree[spr["recomb_node"]]
    coal = tree[spr["coal_node"]]
    broken = recomb.parent
    broken_dist = broken.dist
    assert broken is not None

    # remove recomb subtree from its parent
    broken.children.remove(recomb)

    # adjust recoal if coal branch and broken are the same
    if coal == broken:
        coal = broken.children[0]

    # remove broken from tree
    broken_child = broken.children[0]
    broken_child.parent = broken.parent
    if broken.parent:
        util.replace(broken.parent.children, broken, broken_child)
    broken_child.dist += broken_dist

    # reuse broken node as new coal node
    new_node = broken
    new_node.data["age"] = spr["coal_time"]
    new_node.children = [recomb, coal]
    new_node.parent = coal.parent
    if new_node.parent:
        new_node.dist = new_node.parent.data["age"] - new_node.data["age"]
    else:
        new_node.dist = 0.0
    recomb.parent = new_node
    recomb.dist = new_node.data["age"] - recomb.data["age"]
    coal.parent = new_node
    coal.dist = new_node.data["age"] - coal.data["age"]
    if new_node.parent:
        if coal in new_node.parent.children:
            util.replace(new_node.parent.children, coal, new_node)
        else:
            assert new_node in new_node.parent.children

    # change root
    while tree.root.parent is not None:
        tree.root = tree.root.parent


def write_smc(filename, smc):
    """Writes a SMC file"""

    out = argweaver.open_stream(filename, "w")

    for item in smc:
        if item["tag"] == "NAMES":
            util.print_row("NAMES", *item["names"], out=out)

        elif item["tag"] == "REGION":
            util.print_row("REGION",
                           item["chrom"], item["start"], item["end"], out=out)

        elif item["tag"] == "TREE":
            if not isinstance(item["tree"], basestring):
                tree = format_tree(item["tree"])
            else:
                tree = item["tree"]

            util.print_row("TREE", item["start"], item["end"], tree, out=out)

        elif item["tag"] == "SPR":
            util.print_row("SPR", item["pos"],
                           item["recomb_node"], item["recomb_time"],
                           item["coal_node"], item["coal_time"], out=out)
    out.close()


def parse_tree_data(node, data, *args, **kwargs):
    """Default data reader: reads optional bootstrap and branch length"""

    # detect comments
    data = treelib.read_nhx_data(node, data)

    # parse age
    if "age" in node.data:
        node.data["age"] = float(node.data["age"])

    # parse name and dist
    if ":" in data:
        name, dist = data.split(":")
        node.dist = float(dist)

        if len(name) > 0:
            node.name = int(name)
    else:
        data = data.strip()

        # treat as name
        if data:
            node.name = int(data)


def rename_tree(tree, names=None):
    """Rename the leaves of a tree from ints to names"""

    # rename leaves to integers
    for node in list(tree):
        if node.is_leaf():
            tree.rename(node.name, int(node.name))

    if names:
        # rename leaves according to given names
        for i, name in enumerate(names):
            tree.rename(i, name)

    return tree


def parse_tree(text, names=None):
    """Parse a newick string into a tree"""
    tree = treelib.parse_newick(text, read_data=parse_tree_data)
    return rename_tree(tree, names)


def format_tree(tree):
    """Format a tree into a newick string"""
    def write_data(node):
        if node.is_leaf():
            return ":%f%s" % (node.dist, treelib.format_nhx_comment(node.data))
        else:
            return "%d:%f%s" % (node.name, node.dist,
                                treelib.format_nhx_comment(node.data))

    stream = StringIO.StringIO()
    tree.write(stream, writeData=write_data, oneline=True, rootData=True)
    return stream.getvalue()


def smc2sprs(smc):
    """
    Convert SMC iterable into a series of SPRs

    NOTE: SMC is 1-index and SPRs are 0-index
    """

    names = None
    region = None
    init_tree = None
    tree = None

    for item in smc:
        if item["tag"] == "NAMES":
            names = item["names"]

        elif item["tag"] == "REGION":
            chrom = item["chrom"]
            # convert to 0-index
            region = [item["start"]-1, item["end"]]

        elif item["tag"] == "TREE":
            tree = item["tree"]
            if isinstance(tree, basestring):
                tree = parse_tree(tree)

            if not init_tree:
                init_tree = tree.copy()

                # rename leaves
                for i, name in enumerate(names):
                    init_tree.rename(i, name)

                if "age" in init_tree.root.data:
                    times = dict((node, float(node.data["age"]))
                                 for node in init_tree)
                else:
                    times = None

                arg = arglib.make_arg_from_tree(init_tree, times)

                arg.chrom = chrom
                arg.start = region[0]
                arg.end = region[1]
                yield arg

        elif item["tag"] == "SPR":
            rleaves = [names[i] for i in
                       tree.leaf_names(tree[item["recomb_node"]])]
            cleaves = [names[i] for i in
                       tree.leaf_names(tree[item["coal_node"]])]

            yield (item["pos"], (rleaves, item["recomb_time"]),
                   (cleaves, item["coal_time"]))


def smc2arg(smc):
    """
    Convert SMC to ARG.

    NOTE: SMC is 1-index and ARG is 0-index
    """

    it = smc2sprs(smc)
    tree = it.next()
    arg = arglib.make_arg_from_sprs(tree, it)

    return arg


def arg2smc(arg, names=None, chrom="chr", start=None, end=None,
            format_trees=True):
    """
    Convert ARG to SMC.

    NOTE: SMC is 1-index and ARG is 0-index
    """

    if names:
        leaf_names = names
    else:
        leaf_names = list(arg.leaf_names())

    # convert from 0-index to 1-index
    if start:
        region_start = start + 1
    else:
        region_start = arg.start + 1
    if end:
        region_end = end
    else:
        region_end = arg.end

    # yield SMC header
    yield {"tag": "NAMES",
           "names": leaf_names}
    yield {"tag": "REGION",
           "chrom": chrom,
           "start": region_start,
           "end": region_end}

    last_node_index = None
    for block, tree, last_tree, spr in argweaver.iter_arg_sprs(
            arg, start=start, end=end):
        mapping = argweaver.get_local_node_mapping(tree, last_tree, spr)

        # get node numbering
        node_index = {}
        if last_tree is None:
            # initialize node numbering
            index = 0
            for leaf_name in leaf_names:
                node_index[leaf_name] = index
                index += 1
            for node in tree:
                if not node.is_leaf():
                    node_index[node.name] = index
                    index += 1
        else:
            # determine node numbering using mapping.
            broken_index = None
            for last_node in last_tree:
                node_name = mapping[last_node.name]
                if node_name:
                    node_index[node_name] = last_node_index[last_node.name]
                else:
                    broken_index = last_node_index[last_node.name]

            # find out which node does not have an index and give it the
            # broken index
            for node in tree:
                if node.name not in node_index:
                    node_index[node.name] = broken_index
                    break

            # yield SPR
            yield {"tag": "SPR",
                   "pos": block[0],
                   "recomb_node": last_node_index[spr[0][0]],
                   "recomb_time": spr[0][1],
                   "coal_node": last_node_index[spr[1][0]],
                   "coal_time": spr[1][1]}

        # add age data and rename nodes using index
        tree2 = tree.get_tree()
        ages = treelib.get_tree_ages(tree2)
        for node in tree2:
            node.data["age"] = ages[node]
        rename_tree_nodes(tree2, [(node.name, node_index[node.name])
                                  for node in tree2])
        if format_trees:
            tree2 = format_tree(tree2)

        yield {"tag": "TREE",
               "start": block[0] + 1,
               "end": block[1],
               "tree": tree2}

        last_node_index = node_index


def rename_tree_nodes(tree, renames):
    """Rename all the nodes in a tree."""

    nodes = {}
    for old_name, new_name in renames:
        nodes[old_name] = tree[old_name]
        del tree.nodes[old_name]

    for old_name, new_name in renames:
        tree[new_name] = nodes[old_name]

    return tree


def read_arg(smc_filename, region=None):
    """Read an ARG from an SMC file"""
    return smc2arg(iter_smc_file(smc_filename, parse_trees=True,
                                 apply_spr=True, region=region))


def iter_smc_trees(smc, pos):
    """
    Iterate through local trees at positions 'pos' in filename 'smc_file'
    """
    need_close = False
    if isinstance(smc, basestring):
        need_close = True
        smc = SMCReader(smc)
    try:
        piter = iter(pos)
        item = smc.next()
        for p in piter:
            while (item["tag"] != "TREE" or
                   item["end"] < p):
                item = smc.next()
            if item["start"] <= p:
                yield smc.parse_tree(item["tree"])
    except StopIteration:
        pass
    if need_close:
        smc.close()


#=============================================================================
# multiple SMC files


def get_smc_sample_iter(filename):
    """Returns the iteration number of an SMC filename"""
    if filename.endswith(".smc.gz"):
        filename = filename[:-len(".smc.gz")]
    elif filename.endswith(".gz"):
        filename = filename[:-len(".gz")]
    i = filename.rindex(".")
    return int(filename[i+1:])


def list_smc_files(path):
    """Lists all SMC files in a directory"""
    files = util.list_files(path, ".smc.gz")
    files.sort(key=get_smc_sample_iter)
    return files
