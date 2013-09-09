


class LinkedNode (object):
    """A node in a doubly linked list"""
    
    def __init__(self, item):
        self.next = None
        self.prev = None
        self.item = item

        

class LinkedList (object):
    """A doubly linked list"""
    
    def __init__(self, items=[]):
        self._head = None
        self._tail = None
        self._size = 0

        self.extend(items)
        

    def __len__(self):
        """Return size of list"""
        return self._size


    def __iter__(self):
        """Iterate over the items in a linked list"""
        
        ptr = self._head
        while ptr is not None:
            yield ptr.item
            ptr = ptr.next

    def __reversed__(self):
        """Iterate backwards over list"""
        
        ptr = self._tail
        while ptr is not None:
            yield ptr.item
            ptr = ptr.prev

    def get_head(self):
        return self._head

    def get_tail(self):
        return self._tail

    def get_first(self):
        if self._head is None:
            raise IndexError("No elements in list")
        self._head.item

    def get_last(self):
        if self._last is None:
            raise IndexError("No elements in list")
        self._tail.item
        

    def iter_nodes(self):
        """Iterate over the linked nodes in a list"""

        node = self._head
        while node is not None:
            next = node.next
            yield node
            node = next

    def iter_nodes_reversed(self):
        """Iterate over the linked nodes in a list in reverse"""
        
        node = self._tail
        while node is not None:
            prev = ndoe.prev
            yield node
            node = prev

    def remove_node(self, node):
        """Remove node from list"""

        if node.prev is not None:
            node.prev.next = node.next
        else:
            # first in list
            self._head = node.next
            if self._head:
                self._head.prev = None

        if node.next is not None:
            node.next.prev = node.prev
        else:
            # last in list
            self._tail = node.prev
            if self._tail:
                self._tail.next = None

        self._size -= 1
        

    def append(self, item):
        """Append item to end of list"""

        node = LinkedNode(item)
        
        if self._tail is None:
            # append first node
            self._head = node
            self._tail = self._head
        else:
            # append to end of list
            self._tail.next = node
            node.prev = self._tail
            self._tail = node

        self._size += 1
        return node


    def prepend(self, item):
        """Prepend item to front of list"""

        node = LinkedNode(item)

        if self._head is None:
            # append first node
            self._head = node
            self._tail = self._head
        else:
            # append to front of list
            self._head.prev = node
            node.next = self._head
            self._head = node

        self._size += 1
        return node

    def extend(self, items):
        """Append many items to end of list"""

        for item in items:
            self.append(item)        


    def extend_front(self, items):
        """Prepend many items to front of list"""

        for item in items:
            self.prepend(item)


    def pop(self):
        """Pop item from end of list"""

        if self._tail is None:
            raise IndexError("pop from empty list")
        
        item = self._tail.item
        self._tail = self._tail.prev

        if self._tail is None:
            # list is empty
            self._head = None
        else:
            self._tail.next = None

        self._size -= 1

        return item

    def pop_front(self):
        """Pop item from front of list"""

        if self._head is None:
            raise IndexError("pop from empty list")

        item = self._head.item
        self._head = self._head.next

        if self._head is None:
            # list is empty
            self._tail = None
        else:
            self._head.prev = None

        self._size -= 1

        return item


    def insert_after(self, node, item):
        """Insert a new item after a node in the list"""

        if node is None:
            self.prepend(item)

        # create new node
        node2 = LinkedNode(item)
        node2.prev = node
        node2.next = node.next

        # link surrounding nodes
        node.next = node2
        if node2.next:
            node2.next.prev = node2
        else:
            self._tail = node2

        self._size += 1
        
    
    def clear(self):
        """Clear the list of all items"""

        self._head = None
        self._tail = None
        self._size = 0
            
