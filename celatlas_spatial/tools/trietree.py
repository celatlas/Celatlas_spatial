import numpy as np

from numba import jit, int64, boolean
from numba.experimental import jitclass
from typing import List, Optional, Dict
from celatlas_spatial.__init__ import BC_WIDTH


node_spec = [
    ('children', int64[:, :]),    # compatibility spec for the legacy jitclass
    ('is_end', boolean[:]),       # terminal node flag
    ('next_idx', int64),          # next available index
    ('size', int64)               # current size of the trie
]

@jit(nopython=True)
def decimal_to_quaternary(decimal: int, bc_width: int) -> np.ndarray:
    """Convert decimal number to quaternary digits (base-4)"""
    quaternary = np.zeros(bc_width, dtype=np.int64)
    for i in range(bc_width):
        pos = 2 * (bc_width - 1 - i)
        quaternary[i] = (decimal >> pos) & 0b11
    return quaternary

@jit(nopython=True)
def calculate_nodes_count(items, bc_width: int) -> int:
    """Calculate exact trie size without allocating a node matrix.

    Fixed-width base-4 barcodes have the same lexicographic and integer
    ordering.  After sorting, each barcode contributes only the suffix after
    its longest common prefix with the previous barcode.  This replaces the
    old ``len(items) * bc_width * 4`` temporary matrix, which was a large
    peak-memory cost for spatial whitelists.
    """
    if len(items) == 0:
        return 1

    ordered = np.sort(items)
    node_count = 1 + bc_width
    previous = ordered[0]
    for index in range(1, len(ordered)):
        current = ordered[index]
        common_prefix = 0
        for position in range(bc_width):
            shift = 2 * (bc_width - 1 - position)
            if ((previous >> shift) & 0b11) != ((current >> shift) & 0b11):
                break
            common_prefix += 1
        node_count += bc_width - common_prefix
        previous = current
    return node_count

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
def build_tree_from_items(nodes, items, bc_width: int) -> int:
    """Insert all encoded barcodes in one compiled loop."""
    next_idx = 1
    for item in items:
        current = 0
        for position in range(bc_width):
            shift = 2 * (bc_width - 1 - position)
            digit = (item >> shift) & 0b11
            child = nodes['children'][current, digit]
            if child == -1:
                child = next_idx
                nodes['children'][current, digit] = child
                next_idx += 1
            current = child
        nodes['is_end'][current] = True
    return next_idx

@jit(nopython=True)
def find_exact_quaternary(nodes, n, bc_width: int):
    """
    Find exact match in the trie

    Args:
        nodes: structured array containing both children arrays and is_end flags
        n: decimal value to search for in the trie
    """
    quaternary = decimal_to_quaternary(n, bc_width)
    current = 0
    for digit in quaternary:
        if nodes['children'][current, digit] == -1:
            return False
        current = nodes['children'][current, digit]
    return nodes['is_end'][current]


@jit(nopython=True)
def _find_unique_mismatch_trie(
    nodes,
    query_digits,
    current,
    position,
    errors,
    max_errors,
    encoded,
    bc_width,
):
    """Search the trie directly with a bounded Hamming distance.

    This prunes a candidate as soon as the corresponding trie branch does
    not exist.  It avoids materialising and testing every possible mismatch
    mask (1770 candidates for a 20bp barcode with two mismatches).
    """
    if position == bc_width:
        if nodes['is_end'][current]:
            return 1, encoded
        return 0, np.uint64(0)

    match_count = 0
    match_value = np.uint64(0)
    query_digit = query_digits[position]

    # Visit the matching base first.  Most valid corrections are close to
    # the query, and this also lets us terminate quickly on ambiguity.
    for digit in range(4):
        child = nodes['children'][current, digit]
        if child == -1:
            continue
        next_errors = errors + (1 if digit != query_digit else 0)
        if next_errors > max_errors:
            continue
        next_encoded = (encoded << 2) | np.uint64(digit)
        count, value = _find_unique_mismatch_trie(
            nodes,
            query_digits,
            child,
            position + 1,
            next_errors,
            max_errors,
            next_encoded,
            bc_width,
        )
        if count == 0:
            continue
        if count == 2:
            return 2, value
        match_count += count
        if match_count > 1:
            return 2, value
        match_value = value

    return match_count, match_value


@jit(nopython=True)
def find_unique_mismatch_quaternary(nodes, barcode_int, masks, start, end, bc_width: int):
    """Find whether a barcode has exactly one whitelist neighbour.

    The old barcode path called ``Trie.find`` once for every mismatch mask from
    Python.  For a 20 bp barcode with two mismatches that can mean 1770 Python
    calls for one read.  Keep the same ambiguity rule, but execute the mask
    loop and trie traversals inside one compiled Numba call.

    Returns ``(match_count, match_value)``.  ``match_count`` is capped at two,
    because callers only need to distinguish zero, one, and ambiguous matches.
    """
    match_count = 0
    match_value = np.uint64(0)
    for index in range(start, end):
        candidate = barcode_int ^ masks[index]
        if find_exact_quaternary(nodes, candidate, bc_width):
            match_count += 1
            if match_count > 1:
                return 2, match_value
            match_value = candidate
    return match_count, match_value


class Trie:
    """
    High performance Trie implementation for ACGT sequences encoded as 0123.
    Sequences are fixed length (n bases) and stored/queried using decimal values.
    """
    def __init__(self, items: Optional[List[int]] = None, bc_width: int = BC_WIDTH):
        """Initialize trie with optional list of decimal values"""
        if items is None:
            items = np.array([], dtype=np.uint64)
        elif not isinstance(items, np.ndarray):
            items = np.asarray(items, dtype=np.uint64)
        self.bc_width = bc_width
        self.max_size = len(items)
        required_nodes = max(1, calculate_nodes_count(items, self.bc_width))
        if required_nodes >= np.iinfo(np.int32).max:
            raise MemoryError(
                f'Trie requires {required_nodes:,} nodes, which exceeds the '
                'int32 node-index limit'
            )
        self.nodes = np.zeros(required_nodes, dtype=[
            ('children', np.int32, 4),
            ('is_end', np.bool_)
        ])

        self.nodes['children'].fill(-1)
        self.next_idx = 1  # 0 is reserved for root node

        if isinstance(items, np.ndarray) and items.size > 0:
            self.build_tree(items)

    def build_tree(self, items):
        self.next_idx = build_tree_from_items(self.nodes, items, self.bc_width)

    def find(self, n: int, max_distance: int) -> list:
        results = []
        if find_exact_quaternary(self.nodes, n, self.bc_width):
            results.append((max_distance, n))
        return sorted(results, key=lambda x: x[0])

    def contains(self, n: int) -> bool:
        """Return whether an encoded barcode is present without allocating a list."""
        return bool(find_exact_quaternary(self.nodes, n, self.bc_width))

    def find_unique_mismatch(self, barcode_int, masks, start, end, max_mismatch=None):
        """Return ``(count, barcode)`` for a mismatch-neighbour search."""
        if max_mismatch is not None:
            query_digits = decimal_to_quaternary(np.uint64(barcode_int), self.bc_width)
            return _find_unique_mismatch_trie(
                self.nodes,
                query_digits,
                0,
                0,
                0,
                int(max_mismatch),
                np.uint64(0),
                self.bc_width,
            )
        return find_unique_mismatch_quaternary(
            self.nodes,
            np.uint64(barcode_int),
            masks,
            int(start),
            int(end),
            self.bc_width,
        )

    def __repr__(self) -> str:
        """String representation of the trie"""
        return f'<{self.__class__.__name__} with {np.sum(self.nodes["is_end"])} sequences>'
