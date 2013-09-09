"""

    algorithms for sets

"""



class UnionFind:
    """An implementation of the UNINON/FIND algorithm"""

    def __init__(self, items=[]):
        self._parent = None    
        self._items = set(items)

    def __contains__(self, item):
        """Returns True if item is in the set"""
        return item in self.root()._items

    def __len__(self):
        """Returns the size of the set"""
        return len(self.root()._items)
    
    def __iter__(self):
        """Returns an iterator through the items of the set"""
        return iter(self.root()._items)
    
    
    def add(self, item):
        """Adds an item to the set"""
        self.root()._items.add(item)
    
    def root(self):
        """Returns the root node of a set"""
        node = self
        while node._parent:
            node = node._parent
        if node != self:
            self._parent = node
        return node
    
    def same(self, other):
        """Returns True if two sets are the same"""
        return self.root() == other.root()
    
    def union(self, other):
        """Produces a union between this set and another set"""
        root1 = self.root()
        root2 = other.root()
        if root1 == root2:
            return
        
        root1._items.update(root2._items)
        root2._items = set()
        root2._parent = root1
    
    def members(self):
        """Returns members of the set"""
        return self.root()._items
    
    def has(self, item):
        """Returns True if item is in set"""
        return item in self.members()
    
    def size(self):
        """Returns size of the set"""
        return len(self.root()._items)



def connected_components(components):

    sets = {}

    for comp in components:
        comp_sets = []
        for item in comp:
            if item not in sets:
                s = UnionFind([item])
                sets[item] = s
                comp_sets.append(s)
            else:
                comp_sets.append(sets[item])

        if len(comp_sets) > 1:
            for s in comp_sets[1:]:
                comp_sets[0].union(s)


    # yield unique sets
    done = set()
    for s in sets.itervalues():
        if s.root() not in done:
            done.add(s.root())
            yield s.members()







if __name__ == "__main__":
    
    if True:
        set1 = UnionFind()
        set2 = UnionFind()
        set3 = UnionFind()

        set1.add(1)
        set1.add(2)
        print len(set1), 1 in set1, set1.has(-1)
        set2.add(3)
        set2.add(4)
        set2.add(5)    
        print len(set2)
        set3.add(5)
        set3.add(6)
        set3.add(7)
        print len(set3)
        print set1.same(set2)
        set1.union(set2)
        print set1.same(set2)
        
        set1.union(set3)
        print set1.members()
        print len(set1), len(set2)

