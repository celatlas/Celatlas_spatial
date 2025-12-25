# BK-tree implementation in Python
from numba import jit
from collections import deque
from operator import itemgetter


__all__ = ['hamming_distance', 'BKTree']
_getitem0 = itemgetter(0)


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


class BKTree(object):
    """BK-tree data structure that allows fast querying of matches that are
    "close" given a function to calculate a distance metric (e.g., Hamming
    distance or Levenshtein distance).

    Each node in the tree (including the root node) is a two-tuple of
    (item, children_dict), where children_dict is a dict whose keys are
    non-negative distances of the child to the current item and whose values
    are nodes.
    """
    def __init__(self, distance_func=hamming_distance, items=None):
        """Initialize a BKTree instance with given distance function
        (which takes two items as parameters and returns a non-negative
        distance integer). "items" is an optional list of items to add
        on initialization.

        """
        if items is None:
            items = []
        self.distance_func = distance_func
        self.tree = None
        self.nodes = 0

        _add = self.add
        for item in items:
            _add(item)

    def add(self, item):
        """Add given item to this tree.

        """
        node = self.tree
        if node is None:
            self.tree = (item, {})
            self.nodes += 1
            return

        # Slight speed optimization -- avoid lookups inside the loop
        _distance_func = self.distance_func

        while True:
            parent, children = node
            distance = _distance_func(item, parent)
            node = children.get(distance)
            if node is None:
                children[distance] = (item, {})
                self.nodes += 1
                break

    def find(self, item, n):
        """Find items in this tree whose distance is less than or equal to n
        from given item, and return list of (distance, item) tuples ordered by
        distance.

        """
        if self.tree is None:
            return []

        candidates = deque([self.tree])
        found = []

        # Slight speed optimization -- avoid lookups inside the loop
        _candidates_popleft = candidates.popleft
        _candidates_extend = candidates.extend
        _found_append = found.append
        _distance_func = self.distance_func

        while candidates:
            candidate, children = _candidates_popleft()
            distance = _distance_func(candidate, item)
            if distance <= n:
                _found_append((distance, candidate))

            if children:
                lower = distance - n
                upper = distance + n
                _candidates_extend(c for d, c in children.items() if lower <= d <= upper)

        found.sort(key=_getitem0)
        return found

    def __iter__(self):
        """Return iterator over all items in this tree; items are yielded in
        arbitrary order.

        """
        if self.tree is None:
            return

        candidates = deque([self.tree])

        # Slight speed optimization -- avoid lookups inside the loop
        _candidates_popleft = candidates.popleft
        _candidates_extend = candidates.extend

        while candidates:
            candidate, children = _candidates_popleft()
            yield candidate
            _candidates_extend(children.values())

    def __repr__(self):
        """Return a string representation of this BK-tree with a little bit of info.

        """
        return '<{} using {} with {} top-level nodes>'.format(
            self.__class__.__name__,
            self.distance_func.__name__,
            len(self.tree[1]) if self.tree is not None else 'no',
        )


if __name__ == '__main__':
    import doctest
    doctest.testmod()
