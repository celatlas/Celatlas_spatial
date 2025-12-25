from numba import jit, int64, types, float64
from numba.experimental import jitclass
from typing import Optional
import numpy as np


@jit(nopython=True)
def hamming_distance(x: int, y: int):
    """Calculate the hamming distance (number of bits different) between the
    two integers given.

    """
    z = x ^ y
    distance = 0
    while z:
        distance += 1
        z &= z - 1
    return distance


class AVLNode:
    def __init__(self, key: int):
        self.key = key
        self.height = 1
        self.left = None
        self.right = None


class AVLTree:
    """AVL tree data structure that allows fast querying of items within a given
    Hamming distance of a query item using a distance function (e.g., Hamming
    distance or Levenshtein distance).

    """
    def __init__(self, distance_func=hamming_distance, items=None):
        """Initialize AVL tree with optional items list"""
        self.root = None
        self.distance_func = distance_func
        if items is not None:
            for item in items:
                self.insert(item)

    def height(self, node):
        if not node:
            return 0
        return node.height

    def balance_factor(self, node):
        if not node:
            return 0
        return self.height(node.left) - self.height(node.right)

    def update_height(self, node):
        if not node:
            return
        node.height = max(self.height(node.left), self.height(node.right)) + 1

    def right_rotate(self, y):
        x = y.left
        T2 = x.right
        x.right = y
        y.left = T2
        self.update_height(y)
        self.update_height(x)
        return x

    def left_rotate(self, x):
        y = x.right
        T2 = y.left
        y.left = x
        x.right = T2
        self.update_height(x)
        self.update_height(y)
        return y

    def insert(self, key: int):
        """Insert a new key into the AVL tree"""

        def _insert(root: Optional[AVLNode], key: int) -> AVLNode:
            if not root:
                return AVLNode(key)

            if key < root.key:
                root.left = _insert(root.left, key)
            elif key > root.key:
                root.right = _insert(root.right, key)
            else:
                return root  # Duplicate keys not allowed

            self.update_height(root)
            balance = self.balance_factor(root)

            # Left-Left Case
            if balance > 1 and key < root.left.key:
                return self.right_rotate(root)

            # Right-Right Case
            if balance < -1 and key > root.right.key:
                return self.left_rotate(root)

            # Left-Right Case
            if balance > 1 and key > root.left.key:
                root.left = self.left_rotate(root.left)
                return self.right_rotate(root)

            # Right-Left Case
            if balance < -1 and key < root.right.key:
                root.right = self.right_rotate(root.right)
                return self.left_rotate(root)

            return root

        self.root = _insert(self.root, key)

    def find_exact(self, key: int) -> bool:
        """Find exact match in the tree"""

        def _find(root: Optional[AVLNode], key: int) -> bool:
            if not root:
                return False
            if root.key == key:
                return True
            if key < root.key:
                return _find(root.left, key)
            return _find(root.right, key)

        return _find(self.root, key)

    def find(self, key: int, max_distance: int) -> list:
        """Find all items within given Hamming distance"""
        results = []

        # Slight speed optimization -- avoid lookups inside the loop
        _results_append = results.append
        _distance_func = self.distance_func

        def _search(node: Optional[AVLNode]):
            if not node:
                return

            dist = _distance_func(key, node.key)
            if dist <= max_distance:
                _results_append((dist, node.key))

            # prune search
            key_diff = abs(key - node.key)

            # Search both subtrees since we need to find all matches within distance
            if node.left and key_diff - max_distance <= node.key:
                _search(node.left)
            if node.right and key_diff + max_distance >= node.key:
                _search(node.right)

        _search(self.root)
        return sorted(results, key=lambda x: x[0])


node_spec = [
    ('key', int64),
    ('height', int64),
    ('left_idx', int64),
    ('right_idx', int64),
]

@jit(nopython=True)
def get_height(nodes, idx):
    if idx == -1:
        return 0
    return nodes[idx].height

@jit(nopython=True)
def get_balance(nodes, idx):
    if idx == -1:
        return 0
    return get_height(nodes, nodes[idx].left_idx) - get_height(nodes, nodes[idx].right_idx)

@jit(nopython=True)
def update_height(nodes, idx):
    if idx != -1:
        nodes[idx].height = max(
            get_height(nodes, nodes[idx].left_idx),
            get_height(nodes, nodes[idx].right_idx)
        ) + 1

@jit(nopython=True)
def right_rotate(nodes, y_idx):
    x_idx = nodes[y_idx].left_idx
    T2_idx = nodes[x_idx].right_idx

    nodes[x_idx].right_idx = y_idx
    nodes[y_idx].left_idx = T2_idx

    update_height(nodes, y_idx)
    update_height(nodes, x_idx)
    return x_idx

@jit(nopython=True)
def left_rotate(nodes, x_idx):
    y_idx = nodes[x_idx].right_idx
    T2_idx = nodes[y_idx].left_idx

    nodes[y_idx].left_idx = x_idx
    nodes[x_idx].right_idx = T2_idx

    update_height(nodes, x_idx)
    update_height(nodes, y_idx)
    return y_idx

@jit(nopython=True)
def insert_key(nodes, used_nodes, root_idx, key):
    """Insert a key into the AVL tree and return new root index and next available index."""
    if root_idx == -1:
        node = NumbaNode(key)
        nodes[used_nodes]['key'] = node.key
        nodes[used_nodes]['height'] = node.height
        nodes[used_nodes]['left_idx'] = node.left_idx
        nodes[used_nodes]['right_idx'] = node.right_idx
        return used_nodes, used_nodes + 1

    if key < nodes[root_idx].key:
        new_left, used_nodes = insert_key(nodes, used_nodes, nodes[root_idx].left_idx, key)
        nodes[root_idx].left_idx = new_left
    elif key > nodes[root_idx].key:
        new_right, used_nodes = insert_key(nodes, used_nodes, nodes[root_idx].right_idx, key)
        nodes[root_idx].right_idx = new_right
    else:
        return root_idx, used_nodes

    update_height(nodes, root_idx)
    balance = get_balance(nodes, root_idx)

    # Left-Left Case
    if balance > 1 and key < nodes[nodes[root_idx].left_idx].key:
        return right_rotate(nodes, root_idx), used_nodes

    # Right-Right Case
    if balance < -1 and key > nodes[nodes[root_idx].right_idx].key:
        return left_rotate(nodes, root_idx), used_nodes

    # Left-Right Case
    if balance > 1 and key > nodes[nodes[root_idx].left_idx].key:
        nodes[root_idx].left_idx = left_rotate(nodes, nodes[root_idx].left_idx)
        return right_rotate(nodes, root_idx), used_nodes

    # Right-Left Case
    if balance < -1 and key < nodes[nodes[root_idx].right_idx].key:
        nodes[root_idx].right_idx = right_rotate(nodes, nodes[root_idx].right_idx)
        return left_rotate(nodes, root_idx), used_nodes

    return root_idx, used_nodes

@jit(nopython=True)
def find_exact(nodes, root_idx, key, max_distance=0):
    """Exact search in BST."""
    results = []
    _distance_func = hamming_distance

    current = root_idx
    while current != -1:
        current_key = nodes[current]['key']
        if key == current_key:
            distance = _distance_func(current_key, key)
            results.append((distance, key))  # find exact match, return immediately
            break
        elif key < current_key:
            current = nodes[current]['left_idx']
        else:
            current = nodes[current]['right_idx']
    return results


@jitclass(node_spec)
class NumbaNode:
    def __init__(self, key):
        self.key = key
        self.height = 1
        self.left_idx = -1  # use value -1 to indicate None
        self.right_idx = -1


class OptimizedAVLTree:
    """Optimized AVL tree implementation using numba."""
    def __init__(self, distance_func=hamming_distance, items=None):
        self.max_size = 500000 if items is None else len(items)
        self.nodes = np.zeros(self.max_size, dtype=[
            ('key', np.int64),
            ('height', np.int64),
            ('left_idx', np.int64),
            ('right_idx', np.int64)
        ])
        self.root_idx = -1
        self.used_nodes = 0
        self.distance_func = distance_func

        if items is not None:
            self.build_tree(items)

    def build_tree(self, items):
        """Batch build tree from items."""
        for item in items:
            self.root_idx, self.used_nodes = insert_key(
                self.nodes, self.used_nodes, self.root_idx, item
            )

    def find(self, key: int, max_distance: int) -> list:
        """Find all items within max_distance of key."""
        results = find_exact(self.nodes, self.root_idx, key, max_distance)
        if not results:
            return []
        return sorted(results, key=lambda x: x[0])

    def __repr__(self):
        return f'<{self.__class__.__name__} with {self.used_nodes} nodes>'


if __name__ == '__main__':
    import doctest
    doctest.testmod()
