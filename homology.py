import numpy as np
from functools import reduce
from smith import smith_normal_form

class Group:
    def __init__(self):
        self.factors = {}
    
    def add_factor(self, k, times=1):
        assert k >= 0, "k must be nonnegative"
        assert times >= 0, "times must be nonnegative"
        if k in self.factors:
            self.factors[k] += times
        elif k != 1 and times > 0:
            self.factors[k] = times
        return self
    
    def get_factors(self):
        return self.factors

    def pretty(self, unicode=True):
        if len(self.factors) == 0:
            return "0"
        Z = "\u2124" if unicode else "Z"
        plus = "\u2295" if unicode else "+"
        parts = []
        for k in sorted(self.factors.keys()):
            pow = self.factors[k]
            paren1, paren2 = ("(", ")") if k > 0 and pow > 1 else ("", "")
            cgroup = Z if k == 0 else f"{Z}/{k}{Z}"
            end = f"^{pow}" if pow > 1 else ""
            parts.append(paren1 + cgroup + paren2 + end)
        return reduce(lambda x, y: f"{x} {plus} {y}", parts)

class ChainComplex:
    def __init__(self, matrices, top_dim=None):
        self.matrices = matrices
        if top_dim == None:
            self.top_dim = matrices[-1].shape[1]
        else:
            self.top_dim = top_dim

    def homology(self, reduced=False):
        groups = []
        next_rank = 1 if reduced else 0
        for A in self.matrices:
            v = smith_normal_form(A)
            r = len(v)
            group = Group()
            for alpha in v:
                group.add_factor(alpha)
            group.add_factor(0, A.shape[0] - r - next_rank)
            groups.append(group)
            next_rank = r
        group = Group().add_factor(0, self.top_dim - next_rank)
        groups.append(group)
        return groups

def binary_search(xs, x):
    lo = 0
    hi = len(xs) - 1
    while True:
        mid = (hi + lo) // 2
        if x < xs[mid]:
            hi = mid - 1
        elif x > xs[mid]:
            lo = mid + 1
        else: return mid

# Assumes xs is nonempty and sorted
def nub(xs):
    return [xs[0]] + [ xs[i] for i in range(1, len(xs)) if xs[i] != xs[i - 1] ]

def downward_closure(simplices):
    dim = max(map(len, simplices)) - 1
    stratified = [ [] for _ in range(dim + 1) ]
    for s in simplices:
        stratified[len(s) - 1].append(sorted(s))
    stratified[-1] = nub(sorted(stratified[-1]))
    for n in range(dim, 0, -1):
        for simplex in stratified[n]:
            for k in range(n + 1):
                face = simplex[:k] + simplex[k+1:]
                stratified[n - 1].append(face)
        stratified[n - 1] = nub(sorted(stratified[n - 1]))
    return stratified

# TODO: generalize to cellular homology?
class SimplicialSpace:
    def __init__(self, simplices):
        self.simplices = downward_closure(simplices)
        self.dim = len(self.simplices) - 1

    # Matrix mapping (n+1)-chains to n-chains
    def boundary_matrix(self, n):
        A = np.zeros((len(self.simplices[n]), len(self.simplices[n+1])), 'int64')
        for i, simplex in enumerate(self.simplices[n+1]):
            coef = 1
            for k in range(n + 2):
                face = simplex[:k] + simplex[k+1:]
                j = binary_search(self.simplices[n], face)
                A[j,i] = coef
                coef = -coef 
        return A

    def boundary_complex(self):
        return ChainComplex([ self.boundary_matrix(n) for n in range(self.dim) ], len(self.simplices[-1]))

    def homology(self, reduced=False):
        return self.boundary_complex().homology(reduced)

