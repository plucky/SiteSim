# Walter Fontana, 2022
"""
This module defines a 'heap' structure for more efficient reaction selection.
"""

# This is an adaptation of an implementation provided by Rolf Fagerberg (SDU, Odense, Denmark).

import math


class Heap:
    def __init__(self, data, ident=None):
        self.id = ident
        self.data = data  # 'data' should not be modified outside the class!
        # the number of data entries
        self.n_entries = len(self.data)
        # height of the tree (0, 1, ...) that accommodates at least n_entries (as leaves).
        if self.n_entries == 0:
            self.height = 0
        else:
            self.height = int(math.log(self.n_entries, 2)) + 1
        # define the variables of the class; explained in make_tree()
        self.n_leaves = 0
        self.n_internal_nodes = 0
        self.available = 0
        self.tree = []
        self.first_data_item = 0
        self.last_data_item = 0

        self.make_tree()

    def make_tree(self):
        # the number of leaves in a full and complete tree of height h
        if self.height == 0:
            self.n_leaves = 0
        else:
            self.n_leaves = 2 ** self.height
        # the number of internal nodes in a full and complete tree of height h
        self.n_internal_nodes = max(self.n_leaves - 1, 0)
        # n_entries of the leaves carry data. This, then, is the available space before we need to add another level.
        self.available = self.n_leaves - self.n_entries
        # the data "structure". (All the "structure" is virtual by how we address the elements of the list.)
        self.tree = [0] * self.n_internal_nodes + self.data + [0] * self.available
        # indices of first and last data items
        self.first_data_item = self.n_internal_nodes
        self.last_data_item = self.n_internal_nodes + self.n_entries - 1
        # make the tree structure
        self.initialize_tree()

    def initialize_tree(self):
        # initialize the internal node values, by 'descending' from the parent of the last data item towards the root
        i = (self.last_data_item - 1) // 2  # parent of the last data item
        while i >= 0:
            self.tree[i] = self.tree[2 * i + 1] + self.tree[2 * i + 2]  #self.update_node(i)
            i = i - 1

    def update_node(self, i):
        # the value at node i is the sum of the values of its left (2*i+1) and right (2*i+2) child.
        self.tree[i] = self.tree[2 * i + 1] + self.tree[2 * i + 2]

    def update_from_leaf(self, i):
        while i > 0:
            i = (i - 1) // 2  # parent of i
            self.tree[i] = self.tree[2 * i + 1] + self.tree[2 * i + 2]  # self.update_node(i)

    def add_layer(self):
        # Note: Python is 'pass-by-object-reference', so self.data points to our original data content.
        # If that content changes, as it does during simulation, it is reflected
        # in self.data of our heap object. We don't need to pass any items.

        # update the shape values of the heap
        self.height += 1
        self.n_leaves = 2 ** self.height
        self.n_internal_nodes = self.n_leaves - 1
        self.available = self.n_leaves - self.n_entries
        # make a new heap
        self.data = self.tree[self.first_data_item:self.last_data_item + 1]
        self.tree = [0] * self.n_internal_nodes + self.data + [0] * self.available
        self.first_data_item = self.n_internal_nodes
        self.last_data_item = self.n_internal_nodes + self.n_entries - 1
        self.initialize_tree()

    def insert(self, data_item):
        # The function assumes that the last item has been appended to 'data'.
        # Note: self.data points to the same contents as the data list of the caller,
        # so we don't need to pass an item.
        if self.n_entries == self.n_leaves:  # level is full, add another level
            self.add_layer()
        self.n_entries = self.n_entries + 1
        self.last_data_item = self.n_internal_nodes + self.n_entries - 1
        self.available = self.n_leaves - self.n_entries
        # we update the entry in the current heap
        self.tree[self.last_data_item] = data_item
        self.update_from_leaf(self.last_data_item)

    def delete(self, index):
        # 'index' is the position that the item to be deleted has in mixture.complexes[] _before_ removal.
        i = index + self.first_data_item
        j = self.last_data_item
        # we did that with mixture.complexes[], now we follow with tree[]
        self.tree[i] = self.tree[j]
        self.tree[j] = 0
        self.n_entries = self.n_entries - 1
        self.last_data_item -= 1
        self.update_from_leaf(i)
        self.update_from_leaf(j)

    def modify(self, data_item, index):
        # Call this _after_ modifying properties of a molecular species.
        i = index + self.first_data_item
        self.tree[i] = data_item
        # the below is self.update_from_leaf(i)
        while i > 0:
            i = (i - 1) // 2  # parent of i
            self.tree[i] = self.tree[2 * i + 1] + self.tree[2 * i + 2]  # self.update_node(i)

    def draw_node(self, rv):
        i = 0  # index of the root
        stop = (self.last_data_item - 1) // 2  # parent of the last data item
        while i <= stop:
            left = 2 * i + 1
            right = 2 * i + 2
            l_value = self.tree[left]
            if rv < l_value:  # go left
                i = left
            else:  # go right
                rv = rv - l_value
                i = right
        # This is the index of the entity associated with the selected data item
        return i - self.first_data_item

    def show(self, pp_width=40):
        """
        Print heap stats.
        """
        ident = f'HEAP {self.id}'
        info = f"\n{ident.ljust(70, '-')}\n\n"
        info += f'{"height":>25}: {self.height}\n'
        info += f'{"total nodes":>25}: {self.n_internal_nodes + self.n_leaves}\n'
        info += f'{"internal nodes":>25}: {self.n_internal_nodes}\n'
        info += f'{"occupied data leaves":>25}: {self.n_entries}\n'
        info += f'{"available data leaves":>25}: {self.available}\n'
        info += f'{"start index for data":>25}: {self.first_data_item}\n'
        info += f'{"end index for data":>25}: {self.last_data_item}\n'
        info += '\n'
        return info

    def __str__(self):
        """
        Print heap stats.
        """
        info = f"HEAP {self.id} -> height|{self.height} nodes|{self.n_internal_nodes + self.n_leaves} "
        info += f"occ|{self.n_entries / self.n_leaves: .2f} root value| {self.tree[0]:1.2E}"
        info += '\n'
        return info
