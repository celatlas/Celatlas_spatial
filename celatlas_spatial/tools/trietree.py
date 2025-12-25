import numpy as np

from numba import jit, int64, boolean
from numba.experimental import jitclass
from typing import List, Optional, Dict
from celatlas_spatial.__init__ import BC_WIDTH


node_spec = [
    ('children', int64[:, :]),    # flattened child nodes array
    ('is_end', boolean[:]),       # terminal node flag
    ('next_idx', int64),          # next available index
    ('size', int64)               # current size of the trie
]

@jit(nopython=True)
def decimal_to_quaternary(decimal: int) -> np.ndarray:
    """Convert decimal number to quaternary digits (base-4)"""
    quaternary = np.zeros(BC_WIDTH, dtype=np.int64)
    for i in range(BC_WIDTH):
        pos = 2 * (BC_WIDTH - 1 - i)
        quaternary[i] = (decimal >> pos) & 0b11
    return quaternary

@jit(nopython=True)
def calculate_nodes_count(items) -> int:
    """calculate the number of nodes needed to store the given items in a trie"""
    temp = np.zeros((len(items) * BC_WIDTH, 4), dtype=np.int64) - 1
    next_idx = 1
    for item in items:
        quaternary = decimal_to_quaternary(item)
        current = 0

        for digit in quaternary:
            if temp[current, digit] == -1:
                temp[current, digit] = next_idx
                next_idx += 1
            current = temp[current, digit]
    return next_idx

@jit(nopython=True)
def insert_key(nodes, next_idx, quaternary):
    """
    Insert a key into the trie and return next available index

    Args:
        nodes: structured array containing both children arrays and is_end flags
               dtype=[('children', np.int64, 4), ('is_end', np.bool_)]   
        next_idx: next available node index
        quaternary: quaternary representation of the key to insert

    Returns:
        next_idx: updated next available node index
    """
    current = 0  # root node
    for digit in quaternary:
        if nodes['children'][current, digit] == -1:
            nodes['children'][current, digit] = next_idx
            next_idx += 1
        current = nodes['children'][current, digit]
    nodes['is_end'][current] = True
    return next_idx

@jit(nopython=True)
def find_exact_quaternary(nodes, n):
    """
    Find exact match in the trie

    Args:
        nodes: structured array containing both children arrays and is_end flags
        n: decimal value to search for in the trie
    """
    quaternary = decimal_to_quaternary(n)
    current = 0
    for digit in quaternary:
        if nodes['children'][current, digit] == -1:
            return False
        current = nodes['children'][current, digit]
    return nodes['is_end'][current]


class Trie:
    """
    High performance Trie implementation for ACGT sequences encoded as 0123.
    Sequences are fixed length (n bases) and stored/queried using decimal values.
    """
    def __init__(self, items: Optional[List[int]] = None):
        """Initialize trie with optional list of decimal values"""
        self.max_size = len(items)
        required_nodes = calculate_nodes_count(items)
        self.nodes = np.zeros(required_nodes, dtype=[
            ('children', np.int64, 4),
            ('is_end', np.bool_)
        ])

        self.nodes['children'].fill(-1)
        self.next_idx = 1  # 0 is reserved for root node

        if isinstance(items, np.ndarray) and items.size > 0:
            self.build_tree(items)

    def build_tree(self, items):
        for item in items:
            quaternary = decimal_to_quaternary(item)
            self.next_idx = insert_key(self.nodes, self.next_idx, quaternary)

    def find(self, n: int, max_distance: int) -> list:
        results = []
        if find_exact_quaternary(self.nodes, n):
            results.append((max_distance, n))
        return sorted(results, key=lambda x: x[0])

    def __repr__(self) -> str:
        """String representation of the trie"""
        return f'<{self.__class__.__name__} with {np.sum(self.nodes["is_end"])} sequences>'
