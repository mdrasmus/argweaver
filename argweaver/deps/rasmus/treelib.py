#
# Tree data structures 
#
# Contains special features for representing phylogeny.  
# See compbio.phylo for more.
#
#


# python libs
import copy
import math
import random
import sys
import os
import StringIO

# rasmus libs
try:
    from rasmus import util
except ImportError:
    import util
try:
    from rasmus import textdraw
except ImportError:
    pass


# ply parsing support
try:
    from rasmus import treelib_parser
except ImportError:
    try:
        import treelib_parser
    except ImportError:
        treelib_parser = None


#============================================================================
# Tree data structures
#


class TreeNode (object):
    """A class for nodes in a rooted Tree
    
    Contains fields for branch length 'dist' and custom data 'data'
    """

    def __init__(self, name=None):
        self.name = name
        self.children = []
        self.parent = None
        self.dist = 0
        self.data = {}

    
    def __iter__(self):
        """Iterate through child nodes"""
        return iter(self.children)
    
    
    def copy(self, parent=None, copyChildren=True):
        """Returns a copy of a TreeNode and all of its children"""
        
        node = TreeNode(self.name)
        node.name = self.name
        node.dist = self.dist
        node.parent = parent
        node.data = copy.copy(self.data)
        
        if copyChildren:
            for child in self.children:
                node.children.append(child.copy(node))
        
        return node
        
    def is_leaf(self):
        """Returns True if the node is a leaf (no children)"""
        return len(self.children) == 0
    
    def recurse(self, func, *args):
        """Applies a function 'func' to the children of a node"""
        for child in self.children:
            func(child, *args)
    
    def leaves(self):
        """Returns the leaves beneath the node in traversal order"""
        leaves = []

        def walk(node):
            if node.is_leaf():
                leaves.append(node)
            for child in node.children:
                walk(child)
        walk(self)
          
        return leaves
    
    def leaf_names(self):
        """Returns the leaf names beneath the node in traversal order"""
        return [x.name for x in self.leaves()]
    
    def write_data(self, out):
        """Writes the data of the node to the file stream 'out'"""
        out.write(str(self.dist))

    def __repr__(self):
        """Returns a representation of the node"""

        return "<node %s>" % self.name



class BranchData (object):
    """A class for managing branch specific data for a Tree
    
    By default, this class implements bootstrap data for TreeNode's.  
    
    To incorporate new kinds of branch data, do the following.  Subclass this
    class (say, MyBranchData).  Create Tree's with
    Tree(branch_data=MyBranchData()).  This will ensure your new branch data
    manager is used for manipulations to the tree.  Any tree's copied from the
    tree (via tree.copy()) will also use the same branch manager.
    """

    def __init__(self):
        pass

    def get_branch_data(self, node):
        """Returns branch specific data from a node"""
        if "boot" in node.data:
            return {"boot": node.data["boot"]}
        else:
            return {}
    
    def set_branch_data(self, node, data):
        """Set the branch specific data from 'data' to node.data"""
        if "boot" in data:
            node.data["boot"] = data["boot"]
    
    def split_branch_data(self, node):
        """Split a branch's data into two copies"""
        if "boot" in node.data:
            return {"boot": node.data["boot"]}, {"boot": node.data["boot"]}
        else:
            return {}, {}
    
    def merge_branch_data(self, data1, data2):
        """Merges the branch data from two neighboring branches into one"""
        if "boot" in data1 and "boot" in data2:
            assert data1["boot"] == data2["boot"], (data1, data2)
            return {"boot": data1["boot"]}
        else:
            return {}
    
    

class Tree (object):
    """
    Basic rooted tree
    
    Well suited for phylogenetic trees
    """

    def __init__(self, nextname=1, branch_data=BranchData()):
        self.nodes = {}
        self.root = None
        self.nextname = nextname
        self.default_data = {}
        self.data = {}
        self.branch_data = branch_data


    def copy(self):
        """Returns a copy of the tree"""
        tree = Tree(nextname = self.nextname)
        
        # copy structure
        if self.root != None:
            # copy all nodes
            tree.root = self.root.copy()
            
            # set all names
            def walk(node):
                tree.nodes[node.name] = node
                for child in node.children:
                    walk(child)
            walk(tree.root)
        
        # copy extra data
        tree.copy_data(self)
        tree.copy_node_data(self)
        
        return tree


    #=========================================
    # iterators
    
    def __iter__(self):
        """Iterate through nodes of tree"""
        return self.nodes.itervalues()


    def __len__(self):
        """Returns number of nodes in tree"""
        return len(self.nodes)


    def __getitem__(self, key):
        """Returns node by name"""
        return self.nodes[key]


    def __setitem__(self, key, node):
        """Adds a node to the tree"""
        node.name = key
        self.add(node)


    def __contains__(self, name):
        """Returns True if tree has node with name 'name'"""
        return name in self.nodes
    

    def preorder(self, node=None, is_leaf=lambda x: x.is_leaf()):
        """Iterate through nodes in pre-order traversal"""
        
        if node is None:
            node = self.root

        queue = [node]

        while len(queue) > 0:
            node = queue.pop()
            yield node

            if not is_leaf(node):
                for child in reversed(node.children):
                    queue.append(child)

                
    def postorder(self, node=None, is_leaf=lambda x: x.is_leaf()):
        """Iterate through nodes in post-order traversal"""
        
        if node is None:
            node = self.root

        stack = [[node, 0]]

        while len(stack) > 0:
            node, i = stack[-1]

            if i < len(node.children) and not is_leaf(node):
                stack.append([node.children[i], 0])
                stack[-2][1] += 1
            else:
                yield node
                stack.pop()


    def inorder(self, node=None, is_leaf=lambda x: x.is_leaf()):
        """Iterate through nodes with in-order traversal"""

        if node is None:
            node = self.root

        stack = [[node, 0]]
        
        while len(stack) > 0:
            node, i = stack[-1]

            if node.is_leaf():
                yield node
                stack.pop()

            elif i < len(node.children) and not is_leaf(node):
                assert len(node.children) == 2

                if i == 1:
                    # left has been visited
                    # yield current node then visit right
                    yield node
                
                # recurse into children
                stack.append([node.children[i], 0])
                stack[-2][1] += 1
            else:
                stack.pop()
    
  
    #=============================
    # structure functions

    def make_root(self, name = None):
        """Create a new root node"""
        if name is None:
            name = self.new_name()
        self.root = TreeNode(name)
        return self.add(self.root)


    def add(self, node):
        """Add a node to the tree
           Does not add node to any specific location (use add_child instead).
        """
        self.nodes[node.name] = node
        return node


    def add_child(self, parent, child):
        """Add a child node to an existing node 'parent' in the tree"""
        assert parent != child
        self.nodes[child.name] = child
        self.nodes[parent.name] = parent
        child.parent = parent
        parent.children.append(child)
        return child


    def new_node(self, name=None):
        """Add a new node with name 'name' to the tree"""
        if name is None:
            name = self.new_name()
        return self.add(TreeNode(name))


    def remove(self, node):
        """
        Removes a node from a tree.
        Notifies parent (if it exists) that node has been removed.
        """
        
        if node.parent:
            node.parent.children.remove(node)
        del self.nodes[node.name]
    
    
    def remove_tree(self, node):
        """
        Removes subtree rooted at 'node' from tree.
        Notifies parent (if it exists) that node has been removed.
        """
        
        def walk(node):
            if node.name in self.nodes:
                del self.nodes[node.name]
            for child in node.children:
                walk(child)
        walk(node)
        
        if node.parent:
            node.parent.children.remove(node)
    
    
    def rename(self, oldname, newname):
        """Rename a node in the tree"""
        node = self.nodes[oldname]
        del self.nodes[oldname]
        self.nodes[newname] = node
        node.name = newname
    
    
    def new_name(self):
        """Returns a new node name that should be unique in the tree"""
        name = self.nextname
        self.nextname += 1
        return name

    
    def unique_name(self, name, names, sep="_"):
        """Create a new unique name not already in names"""
        i = 1
        name2 = name
        while name2 in names:
            name2 = name + sep + str(i)
            i += 1
        names.add(name2)
        return name2

    
    def add_tree(self, parent, childTree):
        """Add a subtree to the tree."""
        
        # Merge nodes and change the names of childTree names if they conflict
        # with existing names
        self.merge_names(childTree)        
        self.add_child(parent, childTree.root)
    
    
    def replace_tree(self, node, childTree):
        """Remove node and replace it with the root of childTree"""
    
        # merge nodes and change the names of childTree names if they conflict
        # with existing names
        self.merge_names(childTree)
        parent = node.parent
        if parent:
            index = parent.children.index(node)
            parent.children[index] = childTree.root
            childTree.root.parent = parent
            del self.nodes[node.name]
    
    
    def merge_names(self, tree2):
        """Merge the node names from tree2 into this tree.
           Change any names that conflict"""
    
        for name in tree2.nodes:
            if name in self.nodes:
                name2 = self.new_name()
                self.nodes[name2] = tree2.nodes[name]
                self.nodes[name2].name = name2
            else:
                # make sure I do not issue a name that matches this one
                if isinstance(name, int):
                    if name >= self.nextname:
                        self.nextname = name + 1
                self.nodes[name] = tree2.nodes[name]
        
    
    def clear(self):
        """Clear all nodes from tree"""
        self.nodes = {}
        self.root = None
    
    
    def leaves(self, node=None):
        """Return the leaves of the tree in order"""
        if node is None:
            node = self.root
            if node is None:
                return []
        return node.leaves()
    
    
    def leaf_names(self, node = None):
        """Returns the leaf names of the tree in order"""
        return map(lambda x: x.name, self.leaves(node))
    
    
    #===============================
    # data functions

    def has_data(self, dataname):
        """Does the tree contain 'dataname' in its extra data"""
        return dataname in self.default_data
    
    
    def copy_data(self, tree):
        """Copy tree data to another"""
        self.branch_data = tree.branch_data
        self.default_data = copy.copy(tree.default_data)
        self.data = copy.copy(tree.data)
    
    
    def copy_node_data(self, tree):
        """Copy node data to another tree"""
        for name, node in self.nodes.iteritems():
            if name in tree.nodes:
                node.data = copy.copy(tree.nodes[name].data)
        self.set_default_data()

    
    def set_default_data(self):
        """Set default values in each node's data"""
        for node in self.nodes.itervalues():
            for key, val in self.default_data.iteritems():
                node.data.setdefault(key, val)
    
    
    def clear_data(self, *keys):
        """Clear tree data"""
        for node in self.nodes.itervalues():
            if len(keys) == 0:
                node.data = {}
            else:
                for key in keys:
                    if key in node.data:
                        del node.data[key]
    
    
    #======================================================================
    # branch data functions
    # forward branch data calles to branch data manager        
    
    def get_branch_data(self, node):
        """Returns branch specific data from a node"""
        return self.branch_data.get_branch_data(node)

    def set_branch_data(self, node, data):
        """Set the branch specific data from 'data' to node.data"""
        return self.branch_data.set_branch_data(node, data)
    
    def split_branch_data(self, node):
        """Split a branch's data into two copies"""
        return self.branch_data.split_branch_data(node)
    
    def merge_branch_data(self, data1, data2):
        """Merges the branch data from two neighboring branches into one"""    
        return self.branch_data.merge_branch_data(data1, data2)
    
    
    #=======================================================================
    # input and output
    #
    
    def read_data(self, node, data):
        """Default data reader: reads optional bootstrap and branch length"""

        # also parse nhx comments
        data = read_nhx_data(node, data)
        
        if ":" in data:
            boot, dist = data.split(":")
            node.dist = float(dist)
            
            if len(boot) > 0:
                if boot.isdigit():
                    node.data["boot"] = int(boot)
                else:
                    try:
                        node.data["boot"] = float(boot)
                    except ValueError:
                        # treat as node name
                        name = boot.strip()
                        if name and node.name is None:
                            node.name = name
        else:
            data = data.strip()

            # treat as name
            if data:
                node.name = data
    
    
    def write_data(self, node):
        """Default data writer: writes optional bootstrap and branch length"""
        
        string = ""
        if "boot" in node.data and \
           not node.is_leaf() and \
           self.root != node:
            if isinstance(node.data["boot"], int):
                string += "%d" % node.data["boot"]
            else:
                string += "%f" % node.data["boot"]
        else:
            # see if internal node names exist
            if not node.is_leaf() and isinstance(node.name, str):
                string += node.name

        string += ":%f" % node.dist
        return string
    
    
    def read_newick(self, filename, readData=None):
        """
        Reads newick tree format from a file stream
        
        You can specify a specialized node data reader with 'readData'
        """

        return read_tree(filename, read_data=readData, tree=self)


    def write(self, out=sys.stdout, writeData=None, oneline=False,
              rootData=False):
        """Write the tree in newick notation"""
        self.write_newick(out, writeData=writeData, 
                          oneline=oneline, rootData=rootData)
    
    
    def write_newick(self, out=sys.stdout, writeData=None, oneline=False,
                    rootData=False):
        """Write the tree in newick notation"""
        write_newick(self, util.open_stream(out, "w"), 
                     writeData=writeData, oneline=oneline,
                     rootData=rootData)
            
    
    def get_one_line_newick(self, root_data=False, writeData=None):
        """Get a presentation of the tree in a oneline string newick format"""
        stream = StringIO.StringIO()
        self.write(stream, oneline=True,
                   writeData=writeData, rootData=root_data)
        return stream.getvalue()



#============================================================================
# Input/Output functions

def read_tree(infile, read_data=None, tree=None):
    """Read a tree from a file stream"""
    infile = util.open_stream(infile)
    return parse_newick(infile, read_data=read_data, tree=tree)


def read_newick(infile, read_data=None, tree=None):
    """Read a tree from a file stream"""
    infile = util.open_stream(infile)
    return parse_newick(infile, read_data=read_data, tree=tree)


def iter_trees(treefile):
    """read multiple trees from a tree file"""
    
    infile = util.open_stream(treefile)

    yield read_tree(infile)
    try:
        while True:
            yield read_tree(infile)
    except Exception, e:
        pass


def tokenize_newick(infile):
    """
    Iterates through the tokens in a stream in newick format

    infile -- a string or file stream
    """

    def iter_stream(infile):
        while True:
            yield infile.read(1)

    if not isinstance(infile, basestring):
        infile = iter_stream(infile)
    else:
        infile = iter(infile)
    
    running = True
    word = []
    for c in infile:
        if c == "":
            # EOF encountered
            break

        elif c in " \t\n":
            # skip white space
            if word:
                yield "".join(word)
                word[:] = []
        
        elif c in ";(),:[]":
            # special tokens
            if word:
                yield "".join(word)
                word[:] = []

            if c == "[":
                # parse comment
                word.append(c)
                for c in infile:
                    word.append(c)
                    if c == "]":
                        break
                yield "".join(word)
                word[:] = []
            else:
                yield c
        else:
            # word token
            word.append(c)
            
    if word:
        yield "".join(word)
        word[:] = []


def parse_newick(infile, read_data=None, tree=None):
    """
    Parse a newick string or stream

    infile    -- a string or file stream
    read_data -- an optional function for reading node data fields
    tree      -- an optional tree to populate
    """
    
    # node stack
    ancestors = []

    # create tree
    if tree is None:
        tree = Tree()
    if read_data is None:
        read_data = tree.read_data

    # create root
    node = TreeNode()
    tree.root = node
    nodes = [node]

    # process token stream
    tokens = tokenize_newick(infile)
    token = None
    data = []
    empty = True
    try:
        while True:
            prev_token = token
            token = tokens.next()
            empty = False

            if token == '(': # new branchset
                if data:
                    read_data(node, "".join(data))
                    data = []
                child = TreeNode()
                nodes.append(child)
                child.parent = node
                node.children.append(child)
                ancestors.append(node)
                node = child

            elif token == ',': # another branch
                if data:
                    read_data(node, "".join(data))
                    data = []
                parent = ancestors[-1]
                child = TreeNode()
                nodes.append(child)

                child.parent = parent
                parent.children.append(child)
                node = child

            elif token == ')': # optional name next
                if data:
                    read_data(node, "".join(data))
                    data = []
                node = ancestors.pop()

            elif token == ':': # optional length next
                data.append(token)

            elif token == ';': # end of tree
                if data:
                    read_data(node, "".join(data))
                    data = []
                break

            else:
                if prev_token in '(,':
                    node.name = token
                    
                elif prev_token in '):':
                    data.append(token)

                else:
                    data.append(token)
                    
    except StopIteration:
        if empty:
            raise Exception("Empty tree")

    except Exception, e:
        raise # Exception("Malformed newick: " + repr(e))

    # setup node names
    names = set()
    for node in nodes:
        if node.name is None:
            node.name = tree.new_name()
        node.name = tree.unique_name(node.name, names)
        tree.nodes[node.name] = node

    # test for bootstrap presence
    for node in nodes:
        if "boot" in node.data:
            tree.default_data["boot"] = 0
            break
    tree.set_default_data()

    
    return tree


def write_newick(tree, out=sys.stdout, writeData=None, oneline=False,
                 rootData=False):
    """Write the tree in newick notation"""
    write_newick_node(tree, tree.root, util.open_stream(out, "w"), 
                      writeData=writeData, oneline=oneline,
                      rootData=rootData)

    
def write_newick_node(tree, node, out=sys.stdout, 
                      depth=0, writeData=None, oneline=False,
                      rootData=False):
    """Write the node in newick format to the out file stream"""

    # default data writer
    if writeData is None:
        writeData = tree.write_data

    if not oneline:
        out.write(" " * depth)

    if len(node.children) == 0:
        # leaf
        out.write(str(node.name))
    else:
        # internal node
        if oneline:
            out.write("(")
        else:
            out.write("(\n")
        for child in node.children[:-1]:
            write_newick_node(tree, child, out, depth+1, 
                              writeData=writeData, oneline=oneline)
            if oneline:
                out.write(",")
            else:
                out.write(",\n")
        write_newick_node(tree, node.children[-1], out, depth+1,
                          writeData=writeData, oneline=oneline)
        if oneline:
            out.write(")")
        else:            
            out.write("\n" + (" " * depth) + ")")

    # don't print data for root node
    if depth == 0:
        if rootData:
            out.write(writeData(node))
        if oneline:
            out.write(";")
        else:
            out.write(";\n")
    else:
        out.write(writeData(node))


'''
def read_tree(filename):
    """Read a tree from a newick file"""
    tree = Tree()
    tree.read_newick(filename)
    return tree


def parse_newick(newick):
    """Read a tree from newick notation stored in a string"""
    tree = Tree()
    stream = StringIO.StringIO(newick)
    tree.read_newick(stream)
    return tree
'''


#=============================================================================
# alternate reading functions

def read_newick_ply(filename, readData=None, tree=None):
    """read with PLY"""

    if tree is None:
        tree = Tree()
    else:
        tree.clear()

    # default data reader
    if readData is None:
        readData = tree.read_data

    # get parse tree
    text = util.read_until(util.open_stream(filename), ";")[0] + ";"
    expr = treelib_parser.yacc.parse(text)
    
    # walk the parse tree and build the tree
    names = set()

    def walk(expr):
        children, name, data = expr
        assert ":" not in name, "bad name '%s'" % name

        # parse name
        if name == "":
            name = None
        node = TreeNode(name)

        # parse data
        readData(node, data)

        if node.name is None:
            node.name = tree.new_name()
            
        # ensure unique name
        node.name = tree.unique_name(node.name, names)

        # recurse
        for child in children:
            ret = walk(child)
            if ret:
                tree.add_child(node, ret)
        return node
    tree.root = walk(expr)
    tree.nodes[tree.root.name] = tree.root

    # test for bootstrap presence
    for node in tree.nodes.itervalues():
        if "boot" in node.data:
            tree.default_data["boot"] = 0
            break
    tree.set_default_data()

    return tree

    
def read_newick_recursive(filename, tree=None):
    """
    Reads a big newick file with a custom parser

    DEPRECATED
    """
    
    infile = util.open_stream(filename) #file(filename)    
    opens = [0]
    names = set()

    if tree is None:
        tree = Tree()

    def readchar():
        while True:
            char = infile.read(1)
            if not char or char not in " \t\n": break
        if char == "(": opens[0] += 1
        if char == ")": opens[0] -= 1
        return char

    def read_until(chars):
        token = ""
        while True:
            #char = readchar()
            while True:
                char = infile.read(1)
                if not char or char not in " \t\n": break
            if char == "(": opens[0] += 1
            if char == ")": opens[0] -= 1

            if char in chars or char == "":
                return token, char
            token += char

    def read_dist():
        word = ""
        while True:
            #char = readchar()
            while True:
                char = infile.read(1)
                if not char or char not in " \t\n": break
            if char == "(": opens[0] += 1
            if char == ")": opens[0] -= 1

            if not char in "-0123456789.e":
                return float(word)
            else:
                word += char

    def read_name():
        token = ""
        while True:
            #char = readchar()
            while True:
                char = infile.read(1)
                if not char or char not in " \t\n": break
            if char == "(": opens[0] += 1
            if char == ")": opens[0] -= 1

            if char in ":)," or char == "":
                return token, char
            token += char

    def read_item():
        char1 = readchar()

        if char1 == "(":
            node = TreeNode(tree.new_name())
            depth = opens[0]
            while opens[0] == depth:
                tree.add_child(node, read_item())

            token, char = read_until("):,")
            if char == ":":
                node.dist = read_dist()
            return node
        else:                   
            #word, char = read_until(":),")
            word, char = read_name()
            word = char1 + word.rstrip()

            name = tree.unique_name(word, names)

            node = TreeNode(name)
            if char == ":":
                node.dist = read_dist()
            return node


    def read_root():
        word, char = read_until("(")

        assert char == "("

        node = TreeNode(tree.new_name())
        depth = opens[0]
        while opens[0] == depth:
            tree.add_child(node, read_item())
        return node

    tree.root = read_root()
    tree.add(tree.root)

    return tree


def read_parent_tree(treefile, labelfile=None, labels=None, tree=None):
    """Reads a parent array from a file"""

    if tree is None:
        tree = Tree()

    lines = util.open_stream(treefile).readlines()

    if labelfile:
        labels = util.read_strings(labelfile)

    elif labels is None:
        nitems = (len(lines) + 1)/ 2
        labels = map(str, range(nitems))

    tree.make_root()

    for i, line in enumerate(lines):
        parentid = int(line.split(" ")[0])

        # determine current child
        if i < len(labels):
            child = TreeNode(labels[i])
        else:
            if i in tree.nodes:
                child = tree.nodes[i]
            else:
                child = TreeNode(i)

        if parentid == -1:
            # keep track of all roots
            tree.add_child(tree.root, child)
        else:                
            if not parentid in tree.nodes:
                parent = TreeNode(parentid)
                tree.add(parent)
            else:
                parent = tree.nodes[parentid]

            try:
                tree.add_child(parent, child)
            except:
                print i, parentid

    # remove unused internal nodes
    labelset = set(labels)
    for child in list(tree.root.children):
        if child.is_leaf() and child.name not in labelset:
            tree.remove(child)

    # remove redunant root
    if len(tree.root.children) == 1:
        tree.root = tree.root.children[0]
        tree.remove(tree.root.parent)
        tree.root.parent = None

    return tree


def write_parent_tree(treefile, tree, labels=None):
    """Writes tree to the parent array format"""
    
    ids = {}

    if labels is None:
        labels = tree.leaf_names()

    # assign ids to leaves
    for leafname in labels:
        ids[tree.nodes[leafname]] = len(ids)

    # assign ids to internal nodes
    def walk(node):
        node.recurse(walk)
        if not node.is_leaf():
            ids[node] = len(ids)
    walk(tree.root)

    # build ptree array
    ptree = [0] * len(ids)
    for node, idname in ids.iteritems():
        if node.parent != None:
            ptree[idname] = ids[node.parent]
        else:
            ptree[idname] = -1

    util.write_list(treefile, ptree)


#=============================================================================
# NHX format

def parse_nhx_comment(comment):
    """Parse a NHX comment"""
    for pair in comment.split(":"):
        if "=" in pair:
            yield pair.split("=")

def format_nhx_comment(data):
    """Format a NHX comment"""
    return "[&&NHX:" + ":".join("%s=%s" % (k, v)
                                for k, v in data.iteritems()) + "]"
        

def parse_nhx_data(text):
    """Parse the data field of an NHX file"""
    data = None
    
    if "[" in text:
        data = {}
        i = text.find("[")
        j = text.find("]")
        comment = text[i+1:j]
        text = text[:i]

        if comment.startswith("&&NHX:"):
            for k, v in parse_nhx_comment(comment[6:]):
                data[k] = v

    return text, data


def read_nhx_data(node, text):
    """Read data function for parsing the data field of an NHX file"""

    text, data = parse_nhx_data(text)
    if data:
        node.data.update(data)
    return text


def write_nhx_data(node):
    """Write data function for writing th data field of an NHX file"""
    
    text = Tree().write_data(node)
    if node.data:
        text += format_nhx_comment(node.data)
    return text


#============================================================================
# Misc. functions for manipulating trees

def assert_tree(tree):
    """Assert that the tree data structure is internally consistent"""
    
    visited = set()
    def walk(node):
        assert node.name in tree.nodes
        assert node.name not in visited
        visited.add(node.name)
        if node.parent:
            assert node in node.parent.children
        for child in node.children:
            assert child.parent == node
        node.recurse(walk)
    walk(tree.root)
    
    assert tree.root.parent is None
    assert len(tree.nodes) == len(visited), "%d %d" % (len(tree.nodes), len(visited))



def lca(nodes):
    """Returns the Least Common Ancestor (LCA) of a list of nodes"""
    
    if len(nodes) == 1:
        return nodes[0]
    elif len(nodes) > 2:
        return lca([lca(nodes[:2])] + nodes[2:])
    elif len(nodes) == 2:
        node1, node2 = nodes
        set1 = set([node1])
        set2 = set([node2])
        
        while True:
            if node1 in set2:
                return node1
            if node2 in set1:
                return node2
            if node1.parent != None:
                node1 = node1.parent
            if node2.parent != None:
                node2 = node2.parent
            
            set1.add(node1)
            set2.add(node2)
    else:
        raise Exception("No nodes given")


def find_dist(tree, name1, name2):
    """Returns the branch distance between two nodes in a tree"""

    if not name1 in tree.nodes or \
       not name2 in tree.nodes:
        raise Exception("nodes '%s' and '%s' are not in tree" %
                        (name1, name2))
    
    # find root path for node1
    node1 = tree.nodes[name1]
    path1 = [node1]    
    while node1 != tree.root:
        node1 = node1.parent
        path1.append(node1)
    
    # find root path for node2
    node2 = tree.nodes[name2]
    path2 = [node2]
    while node2 != tree.root:
        node2 = node2.parent
        path2.append(node2)
    
    # find when paths diverge
    i = 1
    while i <= len(path1) and i <= len(path2) and (path1[-i] == path2[-i]):
        i += 1
    
    dist = 0
    for j in range(i, len(path1)+1):
        dist += path1[-j].dist
    for j in range(i, len(path2)+1):
        dist += path2[-j].dist
    
    return dist
        

def descendants(node, lst=None):
    """Return a list of all the descendants beneath a node"""
    if lst is None:
        lst = []
    for child in node.children:
        lst.append(child)
        descendants(child, lst=lst)
    return lst


def count_descendants(node, sizes=None):
    """Returns a dict with number of leaves beneath each node"""
    if sizes is None:
        sizes = {}
    
    if len(node.children) > 0:
        sizes[node] = 0
        for child in node.children:
            count_descendants(child, sizes)
            sizes[node] += sizes[child]
    else:
        sizes[node] = 1
    
    return sizes


def subtree(tree, node):
    """Return a copy of a subtree of 'tree' rooted at 'node'"""
    
    # make new tree
    tree2 = Tree(nextname = tree.new_name())
    
    # copy nodes and data
    tree2.root = node.copy()
    tree2.copy_data(tree)
    
    # add nodes
    def walk(node):
        tree2.add(node)
        node.recurse(walk)
    walk(tree2.root)
    
    return tree2


def max_disjoint_subtrees(tree, subroots):
    """Returns a list of rooted subtrees with atmost one node from 
       the list 'subroots'
    """
    
    marks = {}
    
    # mark the path from each subroot to the root
    for subroot in subroots:
        ptr = subroot
        while ptr != None:
            lst = marks.setdefault(ptr, [])
            lst.append(subroot)
            ptr = ptr.parent

    # subtrees are those trees with nodes that have at most one mark
    subroots2 = []
    def walk(node):
        marks.setdefault(node, [])
        if len(marks[node]) < 2 and \
           (not node.parent or len(marks[node.parent]) >= 2):
            subroots2.append(node)
        node.recurse(walk)
    walk(tree.root)
    
    return subroots2


def tree2graph(tree):
    """Convert a tree to a graph data structure (sparse matrix)"""
    mat = {}
    
    # init all rows of adjacency matrix to 
    for name in tree.nodes:
        mat[name] = {}
    
    for name, node in tree.nodes.iteritems():
        for child in node.children:
            mat[name][child.name] = child.dist
        
        if node.parent:
            mat[name][node.parent.name] = node.dist
            
    return mat


def graph2tree(mat, root, closedset=None):
    """Convert a graph to a tree data structure"""
    
    if closedset is None:
        closedset = set()
    tree = Tree()

    def walk(name):
        node = TreeNode(name)
        node.dist = 0
        closedset.add(name)
        for child in mat[name]:
            if child not in closedset:
                child_node = walk(child)
                child_node.dist = mat[name][child]
                tree.add_child(node, child_node)
        return node            
    tree.root = walk(root)
    
    tree.nextname = max(name for name in tree.nodes if isinstance(name, int))
    
    return tree


def remove_single_children(tree, simplify_root=True):
    """
    Remove all nodes from the tree that have exactly one child
    
    Branch lengths are added together when node is removed.
    """
    
    # find single children
    removed = [node
               for node in tree
               if len(node.children) == 1 and node.parent]
    
    # actually remove children
    for node in removed:
        newnode = node.children[0]
        
        # add distance
        newnode.dist += node.dist
        
        # change parent and child pointers
        newnode.parent = node.parent
        index = node.parent.children.index(node)
        node.parent.children[index] = newnode
        
        # remove old node
        del tree.nodes[node.name]

    # remove singleton from root
    if simplify_root and tree.root and len(tree.root.children) == 1:
        oldroot = tree.root
        tree.root = tree.root.children[0]
        oldroot.children = []
        tree.remove(oldroot)
        tree.root.parent = None
        tree.root.dist += oldroot.dist
    
    return removed



def remove_exposed_internal_nodes(tree, leaves=None):
    """
    Remove all leaves that were originally internal nodes
    
    leaves -- a list of original leaves that should stay
    
    if leaves is not specified, only leaves with strings as names will be kept
    """
    
    if leaves != None:
        stay = set(leaves)
    else:
        # use the fact that the leaf name is a string to determine
        # wether to keep it
        stay = set()
        for leaf in tree.leaves():
            if isinstance(leaf.name, basestring):
                stay.add(leaf)
    
    # post order traverse tree
    def walk(node):
        # keep a list of children to visit, since they may remove themselves
        for child in list(node.children):
            walk(child)

        if node.is_leaf() and node not in stay:
            tree.remove(node)
    walk(tree.root)


def subtree_by_leaves(tree, leaves=None, keep_single=False,
                      simplify_root=True):
    """
    Remove any leaf not in leaves set
    
    leaves        -- a list of leaves that should stay
    keep_single   -- if False, remove all single child nodes
    simplify_root -- if True, basal branch is removed when removing single
                     children nodes
    """
    
    stay = set(leaves)    
    
    # post order traverse tree
    def walk(node):
        # keep a list of children to visit, since they may remove themselves
        for child in list(node.children):
            walk(child)

        if node.is_leaf() and node not in stay:
            tree.remove(node)
    if len(stay) == 0:
        tree.clear()
    else:
        walk(tree.root)

    if not keep_single:
        remove_single_children(tree, simplify_root=simplify_root)

    return tree


def subtree_by_leaf_names(tree, leaf_names, keep_single=False, newCopy=False):
    """Returns a subtree with only the leaves specified"""
    
    if newCopy:
        tree = tree.copy()
    return subtree_by_leaves(tree, [tree.nodes[x] for x in leaf_names],
                             keep_single=keep_single)


def reorder_tree(tree, tree2, root=True):
    """Reorders the branches of tree to match tree2"""

    if root:
        # reroot tree to match tree2
        root_branches = [set(n.leaf_names()) for n in tree2.root.children]

        def walk(node):
            if node.is_leaf():
                leaves = set([node.name])
            else:
                leaves = set()
                for child in node.children:
                    l = walk(child)
                    if l is None:
                        return None
                    leaves = leaves.union(l)

            if leaves in root_branches:
                # root found, terminate search
                reroot(tree, node.name, newCopy=False)
                return None

            return leaves
        walk(tree.root)
    

    # reorder tree to match tree2
    leaf_lookup = util.list2lookup(tree2.leaf_names())

    def mean(lst):
        return sum(lst) / float(len(lst))
    
    def walk(node):
        if node.is_leaf():
            return set([node.name])
        else:
            leaf_sets = []

            for child in node.children:
                leaf_sets.append(walk(child))

            scores = [mean(util.mget(leaf_lookup, l)) for l in leaf_sets]
            rank = util.sortindex(scores)
            node.children = util.mget(node.children, rank)

            # return union
            ret = leaf_sets[0]
            for l in leaf_sets[1:]:
                ret = ret.union(l)
            return ret

    walk(tree.root)


def set_tree_topology(tree, tree2):
    """
    Changes the topology of tree to match tree2

    trees must have nodes with the same names
    """

    nodes = tree.nodes
    nodes2 = tree2.nodes
    
    for node in tree:
        node2 = nodes2[node.name]

        # set parent
        if node2.parent:
            node.parent = nodes[node2.parent.name]
        else:
            node.parent = None
        
        # set children
        if node.is_leaf():
            assert node2.is_leaf()
        else:
            # copy child structure
            node.children[:] = [nodes[n.name] for n in node2.children]

    tree.root = nodes[tree2.root.name]



#=============================================================================
# Rerooting functions
#


def is_rooted(tree):
    """Returns True if tree is rooted"""
    return len(tree.root.children) <= 2



def unroot(tree, newCopy=True):
    """Return an unrooted copy of tree"""
    
    if newCopy:
        tree = tree.copy()

    nodes = tree.root.children
    if len(nodes) == 2 and not (nodes[0].is_leaf() and nodes[1].is_leaf()):
        dist = nodes[0].dist + nodes[1].dist
        data = tree.merge_branch_data(nodes[0].data, nodes[1].data)
        if len(nodes[0].children) < 2:
            nodes.reverse()
        tree.add_child(nodes[0], nodes[1])
        nodes[1].dist = dist
        tree.set_branch_data(nodes[1], data)
        nodes[0].dist = 0
        tree.set_branch_data(nodes[0], {})
        nodes[0].parent = None
        
        # replace root
        del tree.nodes[tree.root.name]
        tree.root = nodes[0]
    return tree


def reroot(tree, newroot, onBranch=True, newCopy=True):
    """
    Change the rooting of a tree
    """
    
    # TODO: remove newCopy (or assert newCopy=False)
    if newCopy:
        tree = tree.copy()
    

    # handle trivial case
    if (not onBranch and tree.root.name == newroot) or \
       (onBranch and newroot in [x.name for x in tree.root.children] and \
        len(tree.root.children) == 2):
        return tree

    assert not onBranch or newroot != tree.root.name, "No branch specified"

    unroot(tree, newCopy=False)

    # handle trivial case
    if not onBranch and tree.root.name == newroot:
        return tree
    
    if onBranch:
        # add new root in middle of branch
        newNode = TreeNode(tree.new_name())
        node1 = tree.nodes[newroot]
        rootdist = node1.dist
        rootdata1, rootdata2 = tree.split_branch_data(node1)
        node1.dist = rootdist / 2.0
        tree.set_branch_data(node1, rootdata1)
        newNode.dist = rootdist / 2.0
        tree.set_branch_data(newNode, rootdata2)
        
        node2 = node1.parent
        node2.children.remove(node1)
        tree.add_child(newNode, node1)
        tree.add_child(node2, newNode)
        
        ptr = node2
        ptr2 = newNode
        newRoot = newNode
    else:
        # root directly on node
        ptr2 = tree.nodes[newroot]
        ptr = ptr2.parent
        newRoot = ptr2
    
    newRoot.parent = None
    
    # reverse parent child relationship of all nodes on path node1 to root
    oldroot = tree.root    
    nextDist = ptr2.dist
    nextData = tree.get_branch_data(ptr2)
    ptr2.dist = 0
    while True:
        nextPtr = ptr.parent
        ptr.children.remove(ptr2)
        tree.add_child(ptr2, ptr)
        
        tmp = ptr.dist
        tmpData = tree.get_branch_data(ptr)
        ptr.dist = nextDist
        tree.set_branch_data(ptr, nextData)
        nextDist = tmp
        nextData = tmpData
        
        ptr2 = ptr
        ptr = nextPtr
        
        if nextPtr is None:
            break
    tree.root = newRoot
    
    return tree


def midpoint_root(tree):
    """
    Reroot a tree using midpoint rerooting
    """

    # get maximum distance from leaves to each node
    depths = {}
    for node in tree.postorder():
        if node.is_leaf():
            depths[node] = (0.0, node)
        else:
            depths[node] = max((c.dist + depths[c][0], depths[c][1])
                               for c in node.children)

    # find maximum path
    dists = []
    for node in tree:
        if node.is_leaf():
            continue
        assert len(node.children) != 1
        tmp = sorted([(c.dist + depths[c][0], depths[c][1], c)
                      for c in node.children])
        dists.append((tmp[-1][0] + tmp[-2][0], node,
                      tmp[-1][2], tmp[-1][1],
                      tmp[-2][2], tmp[-2][1]))
    
    maxdist, top, child1, leaf1, child2, leaf2 = max(dists)
    middist = maxdist / 2.0


    # find longer part of path
    if depths[child1][0] + child1.dist >= middist:
        ptr = leaf1
    else:
        ptr = leaf2

    # find branch that contains midpoint
    dist = 0.0
    while ptr != top:        
        if ptr.dist + dist >= middist:
            # reroot tree
            reroot(tree, ptr.name, onBranch=True, newCopy=False)

            # fixup branch lengths and return
            pdist = sum(c.dist for c in tree.root.children)
            other = filter(lambda x: x != ptr, tree.root.children)[0]
            ptr.dist = middist - dist
            other.dist = pdist - ptr.dist
            return tree

        dist += ptr.dist
        ptr = ptr.parent
    
    
    assert 0 # shouldn't get here


#=============================================================================
# ages (previous known as timestamps)
#
# Methods for calculating the ages (timestamps) of nodes in the tree given
# the branch lengths.
#

def get_tree_ages(tree, root=None, leaves=None, times=None):
    """
    Use the branch lengths of a tree to set timestamps for each node
    Assumes ultrametric tree.

    Leaves have time 0
    """

    if root is None:
        root = tree.root

    esp = .001
    if times is None:
        times = {}

    def walk(node):
        if node.is_leaf() or (leaves and node in leaves):
            t = times.get(node, 0.0)
        else:
            t2 = None
            for child in node.children:
                t = walk(child)

                # ensure branch lengths are ultrametrix
                if t2:                    
                    assert abs(t - t2)/t < esp, (node.name, t, t2)
                t2 = t

        times[node] = t
        return t + node.dist
    walk(root)
    
    return times
get_tree_timestamps = get_tree_ages # backwards compatiability


def set_dists_from_ages(tree, times):
    """
    Sets the branch lengths of a tree using a timestamp dict
    """

    for node in tree:
        if node.parent:
            node.dist = times[node.parent] - times[node]
        else:
            node.dist = 0.0
set_dists_from_timestamps = set_dists_from_ages # backwards compatiability


def check_ages(tree, times):
    """Asserts that timestamps are consistent with tree"""
    
    for node in tree:
        if node.parent:
            if times[node.parent] - times[node] < 0.0 or \
               abs(((times[node.parent] - times[node]) -
                    node.dist)/node.dist) > .001:
                draw_tree_names(tree, maxlen=7, minlen=7)
                util.printcols([(a.name, b) for a, b in times.items()])
                print
                print node.name, node.dist, times[node.parent] - times[node]
                raise Exception("negative time span")
check_timestamps = check_ages # backwards compatiability



#=============================================================================
# parent tables

def tree2parent_table(tree, data_cols=[]):
    """Converts tree to a parent table

    This parent table will have a special numbering for the internal nodes,
    such that their id is also their row in the table.
    
    parent table is a standard format of the Compbio Lab as of 02/01/2007.
    It is a list of triples (node_name, parent_name, dist, ...)
        
    * parent_name indicates the parent of the node.  If the node is a root
      (has no parent), then parent_name is -1
    
    * dist is the distance between the node and its parent.

    * additional columns can be added using the data_cols argument.  The
      values are looked up from node.data[col]
    """

    ptable = []

    for node in tree:
        if node.parent:
            pname = node.parent.name
        else:
            pname = -1
        row = [node.name, pname, node.dist]
        for col in data_cols:
            row.append(node.data[col])
        ptable.append(row)
    
    return ptable


def parent_table2tree(ptable, data_cols=[], convert_names=True):
    """Converts a parent table to a Tree

    if convert_names is True, names that are strings that look like integers
    are converted to ints.
    
    See tree2parent_table for details
    """
    
    tree = Tree()

    parents = {}

    # create nodes
    for row in ptable:
        name, parent = row[:2]
        if name.isdigit():
            name = int(name)
        if parent.isdigit() or parent == "-1":
            parent = int(parent)
        
        node = TreeNode(name)
        node.dist = row[2]
        tree.add(node)
        parents[node] = parent

        for col, val in zip(data_cols, row[3:]):
            node.data[col] = val
            
    # link up parents
    for node, parent_name in parents.iteritems():
        if parent_name == -1:
            tree.root = node
        else:
            parent = tree.nodes[parent_name]
            tree.add_child(parent, node)
    
    return tree



def tree2parent_table_ordered(tree, leaf_names=None):
    """Converts tree to a parent table

    This parent table will have a special numbering for the internal nodes,
    such that their id is also their row in the table.
    
    parent table is a standard format of the Compbio Lab as of 02/01/2007.
    It is a list of triples (node_name, parent_name, dist)
    
    * If the node is a leaf node_name is the leaf name (a string)
    * If the node is internal node_name is an int representing which row
      (0-based) the node is in the table.
    
    * parent_name indicates the parent of the node.  If the parent is root, a 
      -1 is used as the parent_name.
    
    * dist is the distance between the node and its parent.
    
    Arguments:
    leaf_names -- specifies that a tree with only a subset of the leaves 
                  should be used
    
    NOTE: root is not given a row, because root does not have a distance
    the nodeid of the root is -1
    """
    
    if leaf_names != None:
        tree = subtree_by_leaf_names(tree, leaf_names, newCopy=True)
    else:
        leaf_names = tree.leaf_names()

    # assign a numbering to the leaves as specified
    nodeid = 0
    nodeids = {}
    nodes = []
    for leaf in leaf_names:
        nodeids[tree.nodes[leaf]] = leaf
        nodes.append(tree.nodes[leaf])
        nodeid += 1
    
    # assign a numbering to the internal nodes
    for node in tree:
        if node.is_leaf():
            continue
        if node == tree.root:
            nodeids[node] = -1
        else:
            nodeids[node] = nodeid
            nodeid += 1
            nodes.append(node)

    # make parent table
    parentTable = []
    for node in nodes:
        parentTable.append([nodeids[node], nodeids[node.parent], node.dist])
    
    return parentTable
    

def parent_table2tree_ordered(ptable):
    """Converts a parent table to a Tree
    
    See tree2parentTable for details
    """

    # TODO: allow named internal nodes
    
    tree = Tree()
    
    # create nodes
    maxint = 0
    for name, parent_name, dist in parentTable:
        node = TreeNode(name)
        node.dist = dist
        tree.add(node)
        
        if isinstance(name, int):
            maxint = max(name, maxint)
    
    # make a root node
    tree.nextname = maxint + 1
    tree.make_root()

    # link up parents
    for name, parent_name, dist in parentTable:
        if parent_name == -1:
            parent = tree.root
        else:
            parent = tree.nodes[parent_name]
        tree.add_child(parent, tree.nodes[name])
    
    return tree


def write_parent_table(ptable, out=sys.stdout):
    """Writes a parent table to out
    
    out can be a filename or file stream
    """
    
    out = util.open_stream(out, "w")
    for row in ptable:
        out.write("\t".join(map(str, row)) + "\n")



def read_parent_table(filename):
    """Reads a parent table from the file 'filename'
    
    filename can also be an open file stream
    """
    
    infile = util.open_stream(filename)
    ptable = []
    
    for line in infile:
        row = line.rstrip("\n").split("\t")
        name, parent, dist = row[:3]
        
        if name.is_digit():
            name = int(name)
        if parent.is_digit() or parent == "-1":
            parent = int(parent)
            
        ptable.append([name, parent, float(dist)] + row[3:])
    
    return ptable



#=============================================================================
# conversion to other formats

def make_ptree(tree):
    """Make parent tree array from tree"""
    
    nodes = []
    nodelookup = {}
    ptree = []
    
    def walk(node):
        for child in node.children:
            walk(child)
        nodes.append(node)
    walk(tree.root)
    
    def leafsort(a, b):
        if a.is_leaf():
            if b.is_leaf():
                return 0
            else:
                return -1
        else:
            if b.is_leaf():
                return 1
            else:
                return 0
    
    # bring leaves to front
    nodes.sort(cmp=leafsort)
    nodelookup = util.list2lookup(nodes)
    
    for node in nodes:
        if node == tree.root:
            ptree.append(-1)
        else:
            ptree.append(nodelookup[node.parent])
    
    assert nodes[-1] == tree.root
    
    return ptree, nodes, nodelookup

    

#=============================================================================
# Tree visualization
   
def layout_tree(tree, xscale, yscale, minlen=-util.INF, maxlen=util.INF,
               rootx=0, rooty=0):
    """\
    Determines the x and y coordinates for every branch in the tree.
    
    Branch lengths are determined by node.dist
    """    
    
    """
       /-----   ] 
       |        ] nodept[node]
    ---+ node   ]
       |
       |
       \---------
    """
    
    # first determine sizes and nodepts
    coords = {}
    sizes = {}          # number of descendants (leaves have size 1)
    nodept = {}         # distance between node y-coord and top bracket y-coord
    def walk(node):
        # calculate new y-coordinate for node
                
        # compute node sizes
        sizes[node] = 0
        for child in node.children:
            sizes[node] += walk(child)
        
        if node.is_leaf():
            sizes[node] = 1
            nodept[node] = yscale - 1
        else:
            top = nodept[node.children[0]]
            bot = (sizes[node] - sizes[node.children[-1]])*yscale + \
                  nodept[node.children[-1]]
            nodept[node] = (top + bot) / 2.0
        
        return sizes[node]
    walk(tree.root)
    
    # determine x, y coordinates
    def walk(node, x, y):
        xchildren = x+min(max(node.dist*xscale, minlen), maxlen)        
        coords[node] = [xchildren, y + nodept[node]]
                
        if not node.is_leaf():
            ychild = y
            for child in node.children:
                walk(child, xchildren, ychild)
                ychild += sizes[child] * yscale
    walk(tree.root, rootx, rooty)
    
    return coords



def layout_tree_hierarchical(tree, xscale, yscale,
                             minlen=-util.INF, maxlen=util.INF,
                             rootx=0, rooty=0,
                             use_dists=True):
    """\
    Determines the x and y coordinates for every branch in the tree.
    
    Leaves are drawn to line up.  Best used for hierarchical clustering.
    """
    
    """
       /-----   ] 
       |        ] nodept[node]
    ---+ node   ]
       |
       |
       \---------
    """
    
    # first determine sizes and nodepts
    coords = {}
    sizes = {}          # number of descendants (leaves have size 1)
    depth = {}          # how deep in tree is node
    nodept = {}         # distance between node y-coord and top bracket y-coord
    def walk(node):
        # calculate new y-coordinate for node
        
        # recurse: compute node sizes
        sizes[node] = 0        
        for child in node.children:
            sizes[node] += walk(child)
        
        if node.is_leaf():
            sizes[node] = 1
            nodept[node] = yscale - 1
            depth[node] = 0
        else:
            top = nodept[node.children[0]]
            bot = (sizes[node] - sizes[node.children[-1]])*yscale + \
                  nodept[node.children[-1]]
            nodept[node] = (top + bot) / 2.0
            depth[node] = max(depth[child] + 1 for child in node.children)
        
        return sizes[node]
    walk(tree.root)
    
    # determine x, y coordinates
    maxdepth = depth[tree.root]
    def walk(node, x, y):
        xchildren = x + xscale * (maxdepth - depth[node])
        coords[node] = [xchildren, y + nodept[node]]
        
        if not node.is_leaf():
            ychild = y
            for child in node.children:
                walk(child, x, ychild)
                ychild += sizes[child] * yscale
    walk(tree.root, rootx, rooty)
    
    return coords


def layout_tree_vertical(layout, offset=None, root=0, leaves=None,
                         ydir=-1):
    """
    Make layout vertical
    """

    if offset is None:
        if leaves is not None:
            for node in layout:
                if node.is_leaf():
                    offset = leaves - ydir*layout[node][0]
                    break
        else:
            for node in layout:
                if node.parent is None:
                    offset = root - ydir*layout[node][0]
                    break

    for node, (x, y) in layout.iteritems():
        layout[node] = [y, offset + ydir*x]
    return layout



#=============================================================================
# Tree color map

def tree_color_map(leafmap=lambda x: (0, 0, 0)):
    """Returns a simple color mixing colormap"""

    def func(tree):
        def walk(node):
            if node.is_leaf():
                node.color = leafmap(node)
            else:
                colors = []
                for child in node.children:
                    walk(child)
                    colors.append(child.color)
                node.color = color_mix(colors)
        walk(tree.root)
    return func
    
    
def color_mix(colors):
    """Mixes together several color vectors into one"""
    
    sumcolor = [0, 0, 0]
    for c in colors:
        sumcolor[0] += c[0]
        sumcolor[1] += c[1]
        sumcolor[2] += c[2]                
    for i in range(3):
        sumcolor[i] /= float(len(colors))
    return sumcolor


def make_expr_mapping(maps, default_color=(0, 0, 0)):
    """Returns a function that maps strings matching an expression to a value
    
       maps -- a list of pairs (expr, value)
    """

    # find exact matches and expressions
    exacts = {}
    exps = []
    for key, val in maps:
        if "*" not in key:
            exacts[key] = val
        else:
            exps.append((key, val))
    
    # create mapping function
    def mapping(key):
        if key in exacts:
            return exacts[key]    
        
        # return default color
        if not isinstance(key, str):
            return default_color
        
        # eval expressions first in order of appearance
        for exp, val in exps:
            if exp[-1] == "*":
                if key.startswith(exp[:-1]):
                    return val
            elif exp[0] == "*":
                if key.endswith(exp[1:]):
                    return val
        
        raise Exception("Cannot map key '%s' to any value" % key)
    return mapping


def read_tree_color_map(filename):
    """Reads a tree colormap from a file"""
    
    infile = util.open_stream(filename)
    maps = []
    
    for line in infile:
        expr, red, green, blue = line.rstrip().split("\t")
        maps.append([expr, map(float, (red, green, blue))])
    
    name2color = make_expr_mapping(maps)
    
    def leafmap(node):
        return name2color(node.name)

    return tree_color_map(leafmap)


#=========================================================================
# Draw Tree ASCII art 

def draw_tree(tree, labels={}, scale=40, spacing=2, out=sys.stdout,
             canvas=None, x=0, y=0, display=True, labelOffset=-1,
             minlen=1,maxlen=10000):
    """
    Print a ASCII Art representation of the tree
    """
    if canvas is None:
        canvas = textdraw.TextCanvas()
    
    xscale = scale
    yscale = spacing

    
    # determine node sizes
    sizes = {}
    nodept = {}
    def walk(node):
        if node.is_leaf():
            sizes[node] = 1
            nodept[node] = yscale - 1 
        else:
            sizes[node] = 0
        for child in node.children:
            sizes[node] += walk(child)
        if not node.is_leaf():
            top = nodept[node.children[0]]
            bot = (sizes[node] - sizes[node.children[-1]])*yscale + \
                  nodept[node.children[-1]]
            nodept[node] = (top + bot) / 2
        return sizes[node]
    walk(tree.root)
    
    
    def walk(node, x, y):
        # calc coords
        xchildren = int(x+min(max(node.dist*xscale,minlen),maxlen))
        
        # draw branch
        canvas.line(x, y+nodept[node], xchildren, y+nodept[node], '-')
        if node.name in labels:
            branchlen = xchildren - x
            lines = str(labels[node.name]).split("\n")
            labelwidth = max(map(len, lines))
            
            labellen = min(labelwidth, 
                           max(int(branchlen-1),0))
            canvas.text(x + 1 + (branchlen - labellen)/2., 
                        y+nodept[node]+labelOffset, 
                        labels[node.name], width=labellen)
        
        if node.is_leaf():
            canvas.text(xchildren +1, y+yscale-1, str(node.name))
        else:
            top = y + nodept[node.children[0]]
            bot = y + (sizes[node]-sizes[node.children[-1]]) * yscale + \
                      nodept[node.children[-1]]
        
            # draw children
            canvas.line(xchildren, top, xchildren, bot, '|')
            
            ychild = y
            for child in node.children:
                walk(child, xchildren, ychild)
                ychild += sizes[child] * yscale

            
            canvas.set(xchildren, y+nodept[node], '+')
            canvas.set(xchildren, top, '/')
            canvas.set(xchildren, bot, '\\')
        canvas.set(x, y+nodept[node], '+')
    walk(tree.root, x+0, 0)
    
    if display:
        canvas.display(out)


def draw_tree_lens(tree, *args, **kargs):
    labels = {}
    for node in tree.nodes.values():
        labels[node.name] = "%f" % node.dist
    
    draw_tree(tree, labels, *args, **kargs)


def draw_tree_boot_lens(tree, *args, **kargs):
    if not tree.has_data("boot"):
        draw_tree_lens(tree, *args, **kargs)
        return

    labels = {}
    for node in tree.nodes.values():
        if node.is_leaf():
            labels[node.name] = "%f" % node.dist
        else:
            if isinstance(node.data["boot"], int):
                labels[node.name] = "(%d) %f" % (node.data["boot"], node.dist)
            else:
                labels[node.name] = "(%.2f) %f" % (node.data["boot"], node.dist)
    
    draw_tree(tree, labels, *args, **kargs)


def draw_tree_names(tree, *args, **kargs):
    labels = {}
    for node in tree.nodes.values():
        if not node.is_leaf():
            labels[node.name] = "%s" % node.name
    
    draw_tree(tree, labels, *args, **kargs)


def draw_tree_name_lens(tree, *args, **kargs):
    labels = {}
    for node in tree.nodes.values():
        if not node.is_leaf():
            labels[node.name] = "%s " % node.name
        else:
            labels[node.name] =""
        labels[node.name] += "%f" % node.dist
    
    draw_tree(tree, labels, *args, **kargs)



