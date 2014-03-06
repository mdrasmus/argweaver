"""

    Visualization of ARGs (Ancestral Recombination Graphs)

"""

# python imports
from itertools import chain, izip
import random
from math import *

# rasmus imports
from rasmus import util, treelib, stats, sets

# compbio imports
from compbio import arglib

# summon imports
import summon
from summon import sumtree
from summon.shapes import box
from summon.core import *




#=============================================================================
# visualization


def minlog(x, default=10):
    return log(max(x, default))



def layout_arg_leaves(arg):
    """Layout the leaves of an ARG"""

    basetree = treelib.Tree()
    nodes = list(arg.postorder())
    nodes.sort(key=lambda x: x.age)
    lookup = {}

    for node in nodes:
        if node.is_leaf():
            lookup[node] = basetree.new_node(node.name)
        else:
            basechildren = []
            for child in node.children:
                basechild = lookup[child]
                while basechild.parent:
                    basechild = basechild.parent
                basechildren.append(basechild)
            basechildren = util.unique(basechildren)
            if len(basechildren) > 1:
                lookup[node] = basenode = basetree.new_node(node.name)
                for basechild in basechildren:
                    basetree.add_child(basenode, basechild)
            else:
                lookup[node] = basechildren[0]
    basetree.root = lookup[nodes[-1]]

    # assign layout based on basetree layout
    # layout leaves
    return dict((arg[name], i) for i, name in enumerate(basetree.leaf_names()))


def layout_arg(arg, leaves=None, yfunc=lambda x: x):
    """Layout the nodes of an ARG"""

    layout = {}

    # layout leaves
    if leaves is None:
        leafx = layout_arg_leaves(arg)
    else:
        leafx = util.list2lookup(leaves)

    for node in arg.postorder():
        if node.is_leaf():
            layout[node] = [leafx[node], yfunc(node.age)]
        else:
            layout[node] = [
                stats.mean(layout[child][0] for child in node.children),
                yfunc(node.age)]

    return layout


def map_layout(layout, xfunc=lambda x: x, yfunc=lambda x: x):

    for node, (x, y) in layout.items():
        layout[node] = [xfunc(x), yfunc(y)]

    return layout





def get_branch_layout(layout, node, parent, side=0, recomb_width=.4):
    """Layout the branches of an ARG"""

    nx, ny = layout[node]
    px, py = layout[parent]

    if node.event == "recomb":
        if len(node.parents) == 2 and node.parents[0] == node.parents[1]:
            step = recomb_width * [-1, 1][side]
        else:
            step = recomb_width * [-1, 1][node.parents.index(parent)]
        return [nx+step, ny, nx+step, py]
    else:
        return [nx, ny, nx, py]



def show_arg(arg, layout=None, leaves=None, mut=None, recomb_width=.4,
             win=None):
    """Visualize an ARG"""

    if win is None:
        win = summon.Window()
    else:
        win.clear_groups()

    # ensure layout
    if layout is None:
        layout = layout_arg(arg, leaves)

    # callbacks
    def branch_click(node, parent):
        print node.name, parent.name

    # draw ARG
    win.add_group(draw_arg(arg, layout, recomb_width=recomb_width,
                           branch_click=branch_click))

    # draw mutations
    if mut:
        g = group()
        for node, parent, pos, t in mut:
            x1, y1, x2, y2 = get_branch_layout(layout, node, parent)
            g.append(group(draw_mark(x1, t, col=(0,0,1)), color(1,1,1)))
        win.add_group(g)
    return win


def draw_arg(arg, layout, recomb_width=.4, branch_click=None):

    def branch_hotspot(node, parent, x, y, y2):
        def func():
            branch_click(node, parent)
        return hotspot("click", x-.5, y, x+.5, y2, func)

    # draw branches
    g = group(color(1,1,1))
    for node in layout:
        if not node.is_leaf():
            x, y = layout[node]
            for i, child in enumerate(node.children):
                cx, cy = layout[child]
                x1, y1, x2, y2 = get_branch_layout(
                    layout, child, node, i, recomb_width=recomb_width)
                g.append(line_strip(x, y, x2, y2, x1, y1, cx, cy))
                if branch_click:
                    g.append(branch_hotspot(child, node, x1, y1, y2))

    # draw recomb
    for node in layout:
        if node.event == "recomb":
            x, y = layout[node]
            g.append(draw_mark(x, y, col=(1, 0, 0)))

    return g



def show_marginal_trees(arg, mut=None):

    win = summon.Window()
    x = 0
    step = 2
    treewidth = len(list(arg.leaves())) + step

    def trans_camera(win, x, y):
        v = win.get_visible()
        win.set_visible(v[0]+x, v[1]+y, v[2]+x, v[3]+y, "exact")

    win.set_binding(input_key("]"), lambda : trans_camera(win, treewidth, 0))
    win.set_binding(input_key("["), lambda : trans_camera(win, -treewidth, 0))

    blocks = arglib.iter_recomb_blocks(arg)

    for tree, block in izip(arglib.iter_marginal_trees(arg), blocks):
        pos = block[0]
        print pos

        leaves = sorted((x for x in tree.leaves()), key=lambda x: x.name)
        layout = layout_arg(tree, leaves)
        win.add_group(
            translate(x, 0, color(1,1,1),
                      draw_tree(tree, layout),
                      text_clip(
                    "%d-%d" % (block[0], block[1]),
                    treewidth*.05, 0,
                    treewidth*.95, -max(l[1] for l in layout.values()),
                    4, 20,
                    "center", "top")))

        # mark responsible recomb node
        for node in tree:
            if pos != 0.0 and node.pos == pos:
                nx, ny = layout[node]
                win.add_group(draw_mark(x + nx, ny))

        # draw mut
        if mut:
            for node, parent, mpos, t in mut:
                if (node.name in tree and node.name != tree.root.name and
                    block[0] < mpos < block[1]):
                    nx, ny = layout[tree[node.name]]
                    win.add_group(draw_mark(x + nx, t, col=(0,0,1)))
                if node.name in tree and tree[node.name].parents:
                    nx, ny = layout[tree[node.name]]
                    py = layout[tree[node.name].parents[0]][1]
                    start = arg[node.name].data["ancestral"][0][0]
                    win.add_group(lines(color(0,1,0),
                                        x+nx, ny, x+nx, py,
                                        color(1,1,1)))


        x += treewidth

    win.set_visible(* win.get_root().get_bounding() + ("exact",))

    return win


def show_tree_track(tree_track, mut=None, show_labels=False,
                    use_blocks=False, branch_click=None):
    """
    tree_track = [((start, end), tree), ...]
    """

    def draw_labels(tree, layout):
        return group(*
                [text_clip(leaf.name, layout[leaf][0], layout[leaf][1],
                          1, layout[leaf][1] + 1e4, 4, 20, "middle", "left")
                 for leaf in tree.leaves()])

    def branch_hotspot(node, parent, x, y, y2):
        def func():
            branch_click(node, parent)
        return hotspot("click", x-.5, y, x+.5, y2, func)

    def print_branch(node, parent):
        print "node", node.name


    tree_track = iter(tree_track)
    if mut:
        mut = util.PushIter(mut)
    block, tree = tree_track.next()
    if branch_click is True:
        branch_click = print_branch

    win = summon.Window()
    treex = 0
    step = 2
    treewidth = len(list(tree.leaves())) + step

    def trans_camera(win, x, y):
        v = win.get_visible()
        win.set_visible(v[0]+x, v[1]+y, v[2]+x, v[3]+y, "exact")

    win.set_binding(input_key("]"), lambda : trans_camera(win, treewidth, 0))
    win.set_binding(input_key("["), lambda : trans_camera(win, -treewidth, 0))

    for block, tree in chain([(block, tree)], tree_track):
        pos = block[0]
        print pos

        layout = treelib.layout_tree(tree, xscale=1, yscale=1)
        treelib.layout_tree_vertical(layout, leaves=0)
        g = win.add_group(
            translate(treex, 0, color(1,1,1),
                      sumtree.draw_tree(tree, layout,
                                        vertical=True),
                      (draw_labels(tree, layout) if show_labels else group()),
                      text_clip(
                    "%d-%d" % (block[0], block[1]),
                    treewidth*.05, 0,
                    treewidth*.95, -max(l[1] for l in layout.values()),
                    4, 20,
                    "center", "top")))


        clicking = group()
        g.append(clicking)

        # hotspots
        if branch_click:
            for node in tree:
                if node.parent:
                    x, y = layout[node]
                    x2, y2 = layout[node.parent]
                    clicking.append(branch_hotspot(node, node.parent, x, y, y2))
        #win.add_group(clicking)


        # draw mut
        if mut:
            for mpos, age, chroms in mut:
                if block[0] < mpos < block[1]:
                    node = arglib.split_to_tree_branch(tree, chroms)
                    parent = node.parent
                    if node and parent:
                        t = random.uniform(layout[node][1], layout[parent][1])
                        nx, ny = layout[node]
                        win.add_group(draw_mark(treex + nx, t, col=(0,0,1)))
                elif mpos > block[1]:
                    mut.push((mpos, age, chroms))
                    break


        treex += treewidth

    #win.set_visible(* win.get_root().get_bounding() + ("exact",))
    win.home("exact")

    return win



def show_coal_track(tree_track):

    win = summon.Window()

    bgcolor = (1, 1, 1, .1)
    cmap = util.rainbow_color_map(low=0.0, high=1.0)

    maxage = 0
    for (start, end), tree in tree_track:
        print start
        l = []
        times = treelib.get_tree_timestamps(tree)
        nleaves = len(tree.leaves())
        maxage2 = 0
        for node in tree:
            if len(node.children) > 1:
                age = times[node]
                freq = len(node.leaves()) / float(nleaves)
                l.extend([color(*cmap.get(freq)), start, age, end, age])
                if age > maxage2:
                    maxage2 = age
        win.add_group(group(lines(*l), color(*bgcolor),
                      box(start, 0, end, maxage2, fill=True)))
        if maxage2 > maxage:
            maxage = maxage2

    # hotspot
    def func():
        x, y = win.get_mouse_pos()
        print "pos=%s age=%f" % (util.int2pretty(int(x)), y)
    win.add_group(hotspot("click", 0, 0, end, maxage,
                          func))

    win.home("exact")


    return win



def show_smc(smc, mut=None, show_labels=False, branch_click=None,
             use_names=False):
    """
    """

    def draw_labels(tree, layout):
        return group(*
                [text_clip(names[leaf.name],
                           layout[leaf][0] - .4, layout[leaf][1],
                           layout[leaf][0] + .4, layout[leaf][1] - 1e4,
                           4, 20, "top", "center")
                 for leaf in tree.leaves()])

    def branch_hotspot(node, parent, x, y, y2):
        def func():
            branch_click(node, parent)
        return hotspot("click", x-.5, y, x+.5, y2, func)

    def print_branch(node, parent):
        print "node", node.name

    def trans_camera(win, x, y):
        v = win.get_visible()
        win.set_visible(v[0]+x, v[1]+y, v[2]+x, v[3]+y, "exact")

    def on_scroll_window(win):
        region = win.get_visible()
        print region

    def on_resize_window(win):
        region = win.get_visible()
        print region



    branch_color = (1, 1, 1)
    spr_color = (1, 0, 0, .5)
    recomb_color = (1, 0, 0)


    # create window
    win = summon.Window()
    win.set_binding(input_key("]"), lambda : trans_camera(win, treewidth, 0))
    win.set_binding(input_key("["), lambda : trans_camera(win, -treewidth, 0))
    win.add_view_change_listener(lambda : on_scroll_window(win))
    #win.remove_resize_listener(lambda : on_resize_window(win))

    treex = 0
    step = 2

    names = []
    seq_range = [0, 0]
    treewidth = 10
    tree = None
    layout = None

    for item in smc:
        if item["tag"] == "NAMES":
            names = item["names"]
            if not use_names:
                names = map(str, range(len(names)))

            treewidth = len(names)

        elif item["tag"] == "RANGE":
            seq_range = [item["start"], item["end"]]

        elif item["tag"] == "TREE":
            tree = item["tree"]

            layout = treelib.layout_tree(tree, xscale=1, yscale=1)
            treelib.layout_tree_vertical(layout, leaves=0)
            #map_layout(layout, yfunc=minlog)

            region_text = text_clip("%d-%d" % (item["start"], item["end"]),
                        treewidth*.05, 0,
                        treewidth*.95, -max(l[1] for l in layout.values()),
                        4, 20,
                        "center", "top")

            g = win.add_group(
                translate(treex, 0, color(1,1,1),
                          sumtree.draw_tree(tree, layout,
                                            vertical=True),
                          (draw_labels(tree, layout)
                           if show_labels else group()),
                          zoom_clamp(translate(0, -20, region_text),
                                     axis=(treewidth, 0),
                                     miny=1.0, maxy=1.0)
                          ))

            clicking = group()
            g.append(clicking)

        elif item["tag"] == "SPR":

            rx, ry = layout[tree[item["recomb_node"]]]
            ry = item["recomb_time"]
            cx, cy = layout[tree[item["coal_node"]]]
            cy = item["coal_time"]

            g.append(
                group(
                lines(color(*spr_color), rx, ry, cx, cy),
                mark_tree(tree, layout,
                          item["recomb_node"], time=item["recomb_time"],
                          col=recomb_color)))

            treex += treewidth + step


    '''
    tree_track = iter(tree_track)
    if mut:
        mut = util.PushIter(mut)
    block, tree = tree_track.next()
    if branch_click is True:
        branch_click = print_branch

    win = summon.Window()
    treex = 0
    step = 2
    treewidth = len(list(tree.leaves())) + step

    def trans_camera(win, x, y):
        v = win.get_visible()
        win.set_visible(v[0]+x, v[1]+y, v[2]+x, v[3]+y, "exact")

    win.set_binding(input_key("]"), lambda : trans_camera(win, treewidth, 0))
    win.set_binding(input_key("["), lambda : trans_camera(win, -treewidth, 0))

    for block, tree in chain([(block, tree)], tree_track):
        pos = block[0]
        print pos

        layout = treelib.layout_tree(tree, xscale=1, yscale=1)
        treelib.layout_tree_vertical(layout, leaves=0)
        g = win.add_group(
            translate(treex, 0, color(1,1,1),
                      sumtree.draw_tree(tree, layout,
                                        vertical=True),
                      (draw_labels(tree, layout) if show_labels else group()),
                      text_clip(
                    "%d-%d" % (block[0], block[1]),
                    treewidth*.05, 0,
                    treewidth*.95, -max(l[1] for l in layout.values()),
                    4, 20,
                    "center", "top")))


        clicking = group()
        g.append(clicking)

        # hotspots
        if branch_click:
            for node in tree:
                if node.parent:
                    x, y = layout[node]
                    x2, y2 = layout[node.parent]
                    clicking.append(branch_hotspot(node, node.parent, x, y, y2))
        #win.add_group(clicking)


        # draw mut
        if mut:
            for mpos, age, chroms in mut:
                if block[0] < mpos < block[1]:
                    node = arglib.split_to_tree_branch(tree, chroms)
                    parent = node.parent
                    if node and parent:
                        t = random.uniform(layout[node][1], layout[parent][1])
                        nx, ny = layout[node]
                        win.add_group(draw_mark(treex + nx, t, col=(0,0,1)))
                elif mpos > block[1]:
                    mut.push((mpos, age, chroms))
                    break


        treex += treewidth
    '''

    win.home("exact")


    return win



def show_coal_track3(tree_track):

    win = summon.Window()


    bgcolor = (1, 1, 1, .1)
    cmap = util.rainbow_color_map(low=0.5, high=1.0)

    maxage = 0
    for (start, end), tree in tree_track:
        print start
        l = []
        times = treelib.get_tree_timestamps(tree)
        nleaves = len(tree.leaves())
        maxage2 = 0
        for node in tree:
            if len(node.children) > 1:
                age = times[node]
                sizes = [len(x.leaves()) for x in node.children]
                bias = max(sizes) / float(sum(sizes))
                l.extend([color(*cmap.get(bias)), start, age, end, age])
                if age > maxage2:
                    maxage2 = age
        win.add_group(group(lines(*l), color(*bgcolor),
                      box(start, 0, end, maxage2, fill=True)))
        if maxage2 > maxage:
            maxage = maxage2

    def func():
        x, y = win.get_mouse_pos()
        print "pos=%s age=%f" % (util.int2pretty(int(x)), y)
    win.add_group(hotspot("click", 0, 0, end, maxage,
                          func))

    win.home("exact")


    return win


def show_coal_track2(tree_track):

    win = summon.Window()


    bgcolor = (1, 1, 1, .1)
    cmap = util.rainbow_color_map(low=0.0, high=1.0)
    tracks = {}

    maxage = 0
    for (start, end), tree in tree_track:
        print start
        l = []
        times = treelib.get_tree_timestamps(tree)
        nleaves = len(tree.leaves())
        maxage2 = 0
        for node in tree:
            if len(node.children) > 1:
                age = times[node]
                freq = len(node.leaves()) / float(nleaves)
                #sizes = [len(x.leaves()) for x in node.children]
                #m = max(sizes)
                #n = sum(sizes)
                #pval = 2 * (n - m) / float(n - 1)
                l.extend([color(*cmap.get(freq)), start, age, end, age])
                if age > maxage2:
                    maxage2 = age
        win.add_group(group(lines(*l), color(*bgcolor),
                      box(start, 0, end, maxage2, fill=True)))
        if maxage2 > maxage:
            maxage = maxage2

    def func():
        x, y = win.get_mouse_pos()
        print "pos=%s age=%f" % (util.int2pretty(int(x)), y)
    win.add_group(hotspot("click", 0, 0, end, maxage,
                          func))

    win.home("exact")


    return win


def show_coal_track2(tree_track):

    win = summon.Window()


    bgcolor = (1, 1, 1, .1)
    cmap = util.rainbow_color_map(low=0.0, high=1.0)

    maxage = 0
    for (start, end), tree in tree_track:
        print start
        l = []
        times = treelib.get_tree_timestamps(tree)
        nleaves = len(tree.leaves())
        maxage2 = 0
        for node in tree:
            if len(node.children) > 1:
                age = times[node]
                sizes = [len(x.leaves()) for x in node.children]
                m = max(sizes)
                n = sum(sizes)
                pval = 2 * (n - m) / float(n - 1)
                freq = len(node.leaves()) / float(nleaves)
                l.extend([color(*cmap.get(freq)), start, age, end, age])
                if age > maxage2:
                    maxage2 = age
        win.add_group(group(lines(*l), color(*bgcolor),
                      box(start, 0, end, maxage2, fill=True)))
        if maxage2 > maxage:
            maxage = maxage2

    def func():
        x, y = win.get_mouse_pos()
        print "pos=%s age=%f" % (util.int2pretty(int(x)), y)
    win.add_group(hotspot("click", 0, 0, end, maxage,
                          func))

    win.home("exact")


    return win



def draw_tree(tree, layout, orient="vertical"):

    vis = group()
    bends = {}

    for node in tree.postorder():
        # get node coordinates
        nx, ny = layout[node]
        px, py = layout[node.parents[0]] if node.parents else (nx, ny)

        # determine bend point
        if orient == "vertical":
            bends[node] = (nx, py)
        else:
            bends[node] = (px, ny)

        # draw branch
        vis.append(lines(nx, ny, bends[node][0], bends[node][1]))

        # draw cross bar
        if len(node.children) > 0:
            a = bends[node.children[-1]]
            b = bends[node.children[0]]
            vis.append(lines(a[0], a[1], b[0], b[1]))

    return vis


def draw_mark(x, y, col=(1,0,0), size=.5, func=None):
    """Draw a mark at (x, y)"""

    if func:
        h = hotspot("click", x-size, y-size, x+size, y+size, func)
    else:
        h = group()

    return zoom_clamp(
        color(*col),
        box(x-size, y-size, x+size, y+size, fill=True),
        h,
        color(1,1,1),
        origin=(x, y),
        minx=10.0, miny=10.0, maxx=20.0, maxy=20.0,
        link=True)



def mark_tree(tree, layout, name, y=None, time=None,
              col=(1, 0, 0), yfunc=lambda y: y, size=.5):
    nx, ny = layout[tree[name]]
    if y is not None:
        y += ny
    else:
        y = time
    return draw_mark(nx, yfunc(y), col=col, size=size)


def draw_branch_mark(arg, layout, node=None, parent=None, pos=None,
                     chroms=None, age=None, col=(0,0,1)):
    """Draw a mark on a branch of an ARG"""

    if node is None:
        node = arglib.split_to_arg_branch(arg, pos, chroms)
    if parent is None:
        assert pos is not None
        parent = arg.get_local_parent(node, pos)

    if node and parent:
        if age is None:
            t = random.uniform(layout[node][1], layout[parent][1])
        else:
            t = layout[node][1] + (age - node.age)
        nx, ny = layout[node]
        return draw_mark(nx, t, col=col)
    else:
        return group()



def draw_branch(arg, layout, node=None, parent=None, chroms=None,
                pos=None, col=None):
    """Draw a mark on a branch of an ARG"""


    if node is None:
        node = arglib.split_to_arg_branch(arg, pos, chroms)
    if parent is None:
        assert pos is not None
        parent = arg.get_local_parent(node, pos)

    if node and parent:
        x1, y1, x2, y2 = get_branch_layout(layout, node, parent)
        if col is None:
            return lines(x1, y1, x2, y2)
        else:
            return lines(color(*col), x1, y1, x2, y2)
    else:
        return group()




#=============================================================================
# haplotype visualization


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


def layout_tree_leaves_even(tree):
    layout = {}
    y = 0

    for node in inorder_tree(tree):
        if node.is_leaf():
            layout[node.name] = y
        else:
            y += 1

    return layout


def layout_tree_leaves(tree):
    layout = {}
    y = 0

    for node in inorder_tree(tree):
        if node.is_leaf():
            layout[node.name] = y
        else:
            #y += 1
            y += (node.age / 1e3) + 1
            #y += exp(node.age / 5e2) + 1
            #y += log(node.age + 1) ** 3

    vals = layout.values()
    mid = (max(vals) + min(vals)) / 2.0
    for k, v in layout.items():
        layout[k] = (v - mid)

    return layout


def layout_chroms(arg, start=None, end=None):

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
    #layout_func = layout_tree_leaves_even

    for spr in arglib.iter_arg_sprs(arg, start=start, end=end, use_leaves=True):
        print "layout", spr[0]
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


def layout_tree_block(tree, names):
    layout = {}

    def walk(node):
        if node.is_leaf():
            x = names.index(node.name)
            layout[node.name] = (x, x-.25, x+.25, node.age)
            return x-.25, x+.25
        else:
            assert len(node.children) == 2
            low1, high1 = walk(node.children[0])
            low2, high2 = walk(node.children[1])
            x = (min(high1, high2) + max(low1, low2)) / 2.0
            low = min(low1, low2)
            high = max(high1, high2)
            layout[node.name] = (x, low, high, node.age)
            return low, high
    walk(tree.root)
    return layout


def mouse_click(win):
    print win.get_mouse_pos("world")


def chrom_click(win, chrom, block):
    def func():
        if win:
            print chrom, block, win.get_mouse_pos("world")[0]
    return func


def draw_arg_threads(arg, blocks, layout, sites=None,
                     chrom_colors=None, chrom_color=[.2,.2,.8,.8],
                     snp_colors={"compat": [1, 0, 0],
                                 "noncompat": [0, 1, 0]},
                     spr_alpha=1,
                     spr_trim=10,
                     compat=False,
                     draw_group=None,
                     win=None):

    leaf_names = set(arg.leaf_names())

    # TEST:
    rnodes = dict((r.pos, r) for r in arg if r.event == "recomb")

    if draw_group is None:
        draw_group = group()

    # set chromosome color
    if chrom_colors is None:
        chrom_colors = {}
        for name in leaf_names:
            chrom_colors[name] = chrom_color

    spr_colors = {}
    for name in leaf_names:
        spr_colors[name] = list(chrom_colors[name])
        if len(spr_colors[name]) < 4:
            spr_colors[name].append(1.0)
        spr_colors[name][3] *= spr_alpha


    trims = []

    for k, (x1, x2) in enumerate(blocks):
        # calc trims
        length = x2 - x1
        minlen =  0
        spr_trim2 = min(spr_trim, (length - minlen) / 2.0)
        trims.append((x1 + spr_trim2, x2 - spr_trim2))
        trim = trims[-1]

        # horizontal lines
        l = []
        for name in leaf_names:
            c = chrom_colors[name]
            y = layout[k][name]
            l.extend([color(*c), trim[0], y, trim[1], y])
        draw_group.append(lines(*l))

        # SPRs
        if k > 0:
            l = []

            # TEST:
            #rnode = rnodes.get(x1, None)
            #young = (rnode is not None and rnode.age < 500)

            for name in leaf_names:
                #c = [1,0,0] if young else spr_colors[name]
                c = spr_colors[name]
                y1 = layout[k-1][name]
                y2 = layout[k][name]
                l.extend([color(*c), trims[k-1][1], y1, trims[k][0], y2])

            draw_group.append(lines(*l))

        # hotspots
        g = group()
        for name in leaf_names:
            y = layout[k][name]
            g.append(hotspot("click", x1+spr_trim, y+.4, x2-spr_trim, y-.4,
                             chrom_click(win, name, (x1, x2))))
        draw_group.append(g)

        # SNPs
        tree = None
        if sites:
            l = []
            for pos, col in sites.iter_region(x1, x2):
                split = set(sites.get_minor(pos)) & leaf_names
                if len(split) == 0:
                    continue
                if compat:
                    if tree is None:
                        tree = arg.get_marginal_tree((x1+x2)/2.0)
                        arglib.remove_single_lineages(tree)
                    node = arglib.split_to_arg_branch(tree, pos-.5, split)
                    if node is not None:
                        derived = list(tree.leaf_names(node))
                        c = color(*snp_colors["compat"])
                    else:
                        c = color(*snp_colors["noncompat"])
                        derived = split
                else:
                    c = color(*snp_colors["compat"])
                    derived = split

                for d in derived:
                    if d in layout[k]:
                        y = layout[k][d]
                        l.extend([c, pos, y+.4, pos, y-.4])
            draw_group.append(lines(*l))

    return draw_group
