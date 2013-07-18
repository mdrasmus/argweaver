
from itertools import izip
from collections import defaultdict

from compbio import arglib
from compbio import phylo
from rasmus import util


def get_conflicts(splits):
    n = len(splits)
    conflicts = []

    for i in range(n):
        for j in range(i+1, n):
            if not arglib.is_split_compatible(splits[i], splits[j]):
                conflicts.append((i, j))

    return conflicts


def get_unique_conflicts(splits):
    """Ignore redundant conflicts"""

    n = len(splits)
    conflicts = []
    right = [n] * n # nearest conflict to right
    left = [-1] * n  # nearest conflict to left

    # visit conflict edges in increase itervals
    for k in range(1, n):
        for i in range(n - k):
            j = i + k
            if right[i] < left[j]:
                # redundant, skip
                continue
            if not arglib.is_split_compatible(splits[i], splits[j]):
                # conflict
                conflicts.append((i, j))
                right[i] = min(right[i], j)
                left[j] = max(left[j], i)

    return conflicts, left, right


def find_break_points(conflicts):
    # build endpoints list
    points = []
    for a,b in conflicts:
        points.append((a, "start", (a,b)))
        points.append((b, "end", (a,b)))
    points.sort()

    # process endpoints
    done = set()
    group = []
    breaks = []
    for i, kind, edge in points:
        if kind == "start":
            group.append(edge)
        elif kind == "end" and edge not in done:
            breaks.append(i-.5)
            done.update(group)
            group = []

    return breaks


def place_break_points(breaks, start, end, mut_pos):
    """
    Place break points in the middle of mutations
    """

    for b in breaks:
        yield mutindex2pos(b, start, end, mut_pos)


def mutindex2pos(x, start, end, mut_pos):
    """
    Convert mutation indices back to genome positions
    """

    if x == -.5:
        return start
    elif x == len(mut_pos) - .5:
        return end
    else:
        i = int(x - .5)
        j = int(x + .5)
        return (mut_pos[i] + mut_pos[j]) / 2.0



def iter_blocks(start, end, breaks):

    ibreaks = iter(breaks)
    prev = start

    # yield blocks
    for b in ibreaks:
        yield (prev, b)
        prev = b
    yield (prev, end)




def get_split_tracks(conflicts, left, right, breaks):

    n = len(left)
    assert n == len(right)

    all_lefts = defaultdict(lambda: [])
    for i, j in conflicts:
        all_lefts[j].append(i)

    # init tracks
    tracks = [[-.5, n-.5] for i in xrange(n)]

    # set right endpoints
    for i in xrange(n):
        if right[i] == n:
            # no one is to the right
            tracks[i][1] = n-.5
        else:
            tracks[i][1] = max(b for b in breaks if b < right[i])


    # set left endpoints
    for i in xrange(n):
        if left[i] == -1:
            # no one is to the left
            tracks[i][0] = -.5
        else:
            tracks[i][0] = max(tracks[k][1] for k in all_lefts[i])

    return tracks



def get_split_sets(splits, tracks, breaks):
    """
    Group compatiable split into non-recombining blocks
    """

    blocks = []
    for b in breaks:
        blocks.append([])
    blocks.append([])

    seen = set()
    for split, track in izip(splits, tracks):
        if (split, track[0], track[1]) in seen:
            continue
        seen.add((split, track[0], track[1]))

        # find start and end blocks
        start_break = util.binsearch(breaks, track[0])[0]
        if start_break is None:
            start_block = 0
        else:
            start_block = start_break + 1

        end_break = util.binsearch(breaks, track[1])[1]
        if end_break is None:
            end_block = len(blocks) - 1
        else:
            end_block = end_break

        for i in xrange(start_block, end_block+1):
            blocks[i].append(split)

    return blocks



def assert_splits(blocks, split_sets, mut_splits, mut_pos):

    i = 0

    # are splits in split_sets mutually compatiable?
    for block, splits in izip(blocks, split_sets):
        for i in xrange(len(splits)):
            for j in xrange(i+1, len(splits)):
                assert arglib.is_split_compatible(splits[i], splits[j]), \
                    (splits[i], splits[j])

    # are splits in split_sets compatible with mutation splits?
    for block, splits in izip(blocks, split_sets):
        while i < len(mut_pos) and mut_pos[i] < block[1]:
            for split in splits:
                assert arglib.is_split_compatible(mut_splits[i], split), \
                    (mut_splits[i], split)
            i += 1



def infer_parsimonious_splits(splits, mut_pos, start, end):

    # get conflict edges
    conflicts, left, right = get_unique_conflicts(splits)

    # find max parsimonious break points
    breaks = find_break_points(conflicts)

    # determine split tracks
    tracks = get_split_tracks(conflicts, left, right, breaks)

    # build split sets and block boundaries
    split_sets = get_split_sets(splits, tracks, breaks)
    breaks2 = place_break_points(breaks, start, end, mut_pos)
    blocks = list(iter_blocks(start, end, breaks2))

    return blocks, split_sets



def iter_parsimonous_tree(split_sets, leaves, rooted=True):

    # make inferred trees
    for split_set in split_sets:
        splits2 = [(map(str, s), tuple(str(i) for i in leaves
                                       if i not in s))
                   for s in split_set]
        tree = phylo.splits2tree(splits2, rooted=rooted)

        # remove boot information
        for node in tree:
            if "boot" in node.data:
                del node.data["boot"]

        yield tree






def get_mutation_split_tracks(arg, mut_splits, mut_pos):

    mut_splits_set = set(mut_splits)
    mut_split_tracks = defaultdict(lambda: [])

    # find regions for splits
    i = 0
    for block, tree in izip(arglib.iter_recomb_blocks(arg),
                            arglib.iter_marginal_trees(arg)):
        for node in tree:
            if len(node.children) != 2 or node.children[0] == node.children[1]:
                continue
            split = tuple(sorted(tree.leaf_names(node)))
            if split in mut_splits_set:
                regions = mut_split_tracks[split]
                if len(regions) > 0 and regions[-1][1] == block[0]:
                    # extend region
                    regions[-1] = (regions[-1][0], block[1])
                else:
                    # add new region
                    regions.append(block)

    # keep only tracks who have a mutation in their interval
    mut_tracks = []
    for i in xrange(len(mut_pos)):
        for region in mut_split_tracks[mut_splits[i]]:
            if region[0] < mut_pos[i] < region[1]:
                a = util.binsearch(mut_pos, region[0])[0]
                a = a if a is not None else -.5
                b = util.binsearch(mut_pos, region[1])[1]
                b = b if b is not None else len(mut_pos)-.5
                mut_tracks.append((a,b))
                break
        else:
            assert False, i

    return mut_tracks

