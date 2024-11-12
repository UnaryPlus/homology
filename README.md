# Homology
A Python program to compute the integral homology groups of a simplicial complex. The core of the algorithm is the function `smith_normal_form`, defined in `smith.py`, which computes the Smith normal form of an integer matrix.

Example:
```python
from homology import SimplicialSpace
triangle = SimplicialSpace([[0, 1], [0, 2], [1, 2]])
# alternatively: triangle = SimplicialSpace(["AB", "AC", "BC"])
H = triangle.homology()
print([ G.pretty() for G in H ])
# ['ℤ', 'ℤ']
```

I might add a class for constructing cell complexes in the future, once I figure out how to do that.